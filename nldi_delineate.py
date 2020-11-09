## NLDI split catchment delineation script

# -----------------------------------------------------
# Martyn Smith USGS
# 10/29/2020
# NLDI Delineation script
# -----------------------------------------------------

# list of required python packages:
# gdal, pysheds, requests

###### CONDA CREATE ENVIRONMENT COMMAND
#conda create -n delineate python=3.6.8 gdal pysheds requests
###### CONDA CREATE ENVIRONMENT COMMAND

from osgeo import ogr, osr, gdal
from pysheds.grid import Grid
import requests
import time
import json

#arguments
NLDI_URL = 'https://labs.waterdata.usgs.gov/api/nldi/linked-data/comid/'
NLDI_GEOSERVER_URL = 'https://labs.waterdata.usgs.gov/geoserver/wmadata/ows'
NHDPLUS_FLOWLINES_QUERY_URL = 'https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer/6/query'
OUT_PATH = 'C:/NYBackup/GitHub/nldi-splitCatchment/data/'
IN_FDR = 'C:/NYBackup/GitHub/nldi-splitCatchment/data/nhdplus/NHDPlusMA/NHDPlus02/NHDPlusFdrFac02b/fdr'
OUT_FDR = 'C:/NYBackup/GitHub/nldi-splitCatchment/data/catch_fdr.tif'

class Watershed:
    """Define inputs and outputs for the main Watershed class"""

    ogr.UseExceptions()
    gdal.UseExceptions() 

    def __init__(self, x=None, y=None):

        self.x = x
        self.y = y
        self.catchment_identifier = None
        self.catchmentGeom = None
        self.splitCatchmentGeom = None
        self.upstreamBasinGeom = None
        self.mergedCatchmentGeom = None

        #input point spatial reference
        self.sourceprj = osr.SpatialReference()
        self.sourceprj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
        # self.sourceprj_espg = self.sourceprj.GetAttrValue("AUTHORITY")
        # print("HERE", self.sourceprj)

        # Getting spatial reference of input raster
        raster = gdal.Open(IN_FDR, gdal.GA_ReadOnly)
        self.Projection = raster.GetProjectionRef()
        self.targetprj = osr.SpatialReference(wkt = raster.GetProjection())

        #create transform
        self.transformToRaster = osr.CoordinateTransformation(self.sourceprj, self.targetprj)
        self.transformToWGS = osr.CoordinateTransformation(self.targetprj, self.sourceprj)

        #kick off
        self.transform_click_point()

## helper functions
    def geom_to_geojson(self, geom, name, simplify_tolerance=10, write_output=True):
        """Return a geojson from an OGR geom object"""

        #get area in local units
        area = geom.GetArea()
        print(name + ' area: ' + str(area*0.00000038610) + ' square miles')

        #optional simplify
        geom = geom.Simplify(simplify_tolerance)

        #don't want to affect original geometry
        transform_geom = geom.Clone()
        
        #trasnsform geometry from whatever the local projection is to wgs84
        transform_geom.Transform(self.transformToWGS)
        json_text = transform_geom.ExportToJson()

        #add some attributes
        geom_json = json.loads(json_text)

        #get area in local units
        area = geom.GetArea()

        #create json structure
        geojson_dict = {
            "type": "Feature",
            "geometry": geom_json,
            "properties": {
                "area": area
            }
        }

        if write_output:
            f = open(OUT_PATH + name + '.geojson','w')
            f.write(json.dumps(geojson_dict))
            f.close()
            print('Exported geojson:', name)
        
        return geojson_dict

## main functions
    def transform_click_point(self):
        """Transform (reproject) assumed WGS84 coordinates to input raster coordinates"""

        print('Input X,Y:', self.x, self.y)
        self.projectedLng, self.projectedLat, z = self.transformToRaster.TransformPoint(self.x,self.y)      
        print('Projected X,Y:',self.projectedLng, ',', self.projectedLat)

        self.get_local_catchment_geom()

    def get_local_catchment_geom(self):
        """Perform point in polygon query to NLDI geoserver to get local catchment geometry"""

        wkt_point = "POINT(%f %f)" %  (self.x , self.y)
        cql_filter = "INTERSECTS(the_geom, %s)" % (wkt_point)

        payload = {
            'service': 'wfs', 
            'version': '1.0.0', 
            'request': 'GetFeature', 
            'typeName': 'wmadata:catchmentsp', 
            'outputFormat': 'application/json',
            'srsName': 'EPSG:4326',
            'CQL_FILTER': cql_filter
        }

        #request catchment geometry from point in polygon query from NLDI geoserver
        # https://labs.waterdata.usgs.gov/geoserver/wmadata/ows?service=wfs&version=1.0.0&request=GetFeature&typeName=wmadata%3Acatchmentsp&outputFormat=application%2Fjson&srsName=EPSG%3A4326&CQL_FILTER=INTERSECTS%28the_geom%2C+POINT%28-73.745860+44.006830%29%29
        r = requests.get(NLDI_GEOSERVER_URL, params=payload)

        print('request url: ', r.url)
        resp = r.json()

        #get catchment id
        self.catchment_identifier = json.dumps(resp['features'][0]['properties']['featureid'])

        #get main catchment geometry polygon
        gj_geom = json.dumps(resp['features'][0]['geometry'])
        self.catchmentGeom = ogr.CreateGeometryFromJson(gj_geom)

        #transform catchment geometry
        self.catchmentGeom.Transform(self.transformToRaster)

        self.geom_to_geojson(self.catchmentGeom, 'catchment')

        #get extent of transformed polygon
        minX, maxX, minY, maxY = self.catchmentGeom.GetEnvelope() # Get bounding box of the shapefile feature
        bounds = [minX, minY, maxX, maxY]

        print('projected bounds', bounds)

        self.splitCatchmentGeom = self.split_catchment(bounds, self.projectedLng,self.projectedLat)

        #get upstream basin
        self.upstreamBasinGeom = self.get_upstream_basin()

    def get_upstream_basin(self):
        """Use local catchment identifier to get upstream basin geometry from NLDI"""

        #request upstream basin
        payload = {'f': 'json', 'simplified': 'false'}
        
        #request upstream basin from NLDI using comid of catchment point is in
        r = requests.get(NLDI_URL + self.catchment_identifier + '/basin', params=payload)

        #print('upstream basin', r.text)
        resp = r.json()

        #convert geojson to ogr geom
        gj_geom = json.dumps(resp['features'][0]['geometry'])
        self.upstreamBasinGeom = ogr.CreateGeometryFromJson(gj_geom)
        self.upstreamBasinGeom.Transform(self.transformToRaster)

        self.geom_to_geojson(self.upstreamBasinGeom, 'upstreamBasin')

        self.mergeGeoms()

    def mergeGeoms(self):
        """Attempt at merging geometries"""

        #if point is on a flowline we have an upstream basin and need to do some geometry merging
        if self.query_flowlines():

            #create new cloned geom
            self.mergedCatchmentGeom = self.upstreamBasinGeom.Clone()

            #subtract splitCatchment geom from full catchment geom
            diff = self.catchmentGeom.Difference(self.splitCatchmentGeom.Buffer(10).Buffer(-10))

            #subtract splitCatchment geom from upstream basin geometry
            self.mergedCatchmentGeom = self.mergedCatchmentGeom.Difference(diff).Simplify(30)

            #write out
            self.geom_to_geojson(self.mergedCatchmentGeom, 'xxFinalBasinxx')

        #otherwise, we can just return the split catchment
        else:
            self.geom_to_geojson(self.splitCatchmentGeom, 'xxFinalBasinxx')

    def query_flowlines(self):
        """Determine if X,Y falls on NHD Plus v2 flowline (within a tolerance)"""

        #example url
        #"https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer/6/query?geometry=-73.82705,43.29139&outFields=GNIS_NAME%2CREACHCODE&geometryType=esriGeometryPoint&inSR=4326&distance=100&units=esriSRUnit_Meter&returnGeometry=false&f=pjson", 

        #perhaps look at this code to snap input point to closest point along a line
        #https://github.com/marsmith/ADONNIS/blob/c54322eaeee17a415b7971c4f5ad714d3d3dccea/js/main.js#L421-L551

        #request upstream basin
        payload = {
            'f': 'pjson', 
            'geometryType': 'esriGeometryPoint',
            'inSR':'4326',
            'geometry': str(self.x) + ',' + str(self.y),
            'distance': 100,
            'units': 'esriSRUnit_Meter',
            'outFields': 'GNIS_NAME,REACHCODE',
            'returnGeometry': 'false'
        }
        
        # #request upstream basin from NLDI using comid of catchment point is in
        #r = requests.get(NHDPLUS_FLOWLINES_QUERY_URL, params=payload)

        print('nhd flowline query:', r.url)

        #print('response', r.text)
        # resp = r.json()

        return True

    def split_catchment(self, bounds, x, y): 
        """Use catchment bounding box to clip NHD Plus v2 flow direction raster, and product split catchment delienation from X,Y"""

        RasterFormat = 'GTiff'
        PixelRes = 30

        #method to use catchment bounding box instead of exact geom
        gdal.Warp(OUT_FDR, IN_FDR, format=RasterFormat, outputBounds=bounds, xRes=PixelRes, yRes=PixelRes, dstSRS=self.Projection, resampleAlg=gdal.GRA_NearestNeighbour, options=['COMPRESS=DEFLATE'])

        #start pysheds catchment delineation
        grid = Grid.from_raster(OUT_FDR, data_name='dir')

        #compute flow accumulation to snap to
        dirmap = (64,  128,  1,   2,    4,   8,    16,  32)
        grid.accumulation(data='dir', dirmap=dirmap, out_name='acc', apply_mask=False)

        grid.to_raster('acc', 'C:/NYBackup/GitHub/ss-delineate/data/acc.tif', view=False, blockxsize=16, blockysize=16)

        #snap the pourpoint to 
        xy = (x, y)
        new_xy = grid.snap_to_mask(grid.acc > 50, xy, return_dist=False)

        #get catchment with pysheds
        grid.catchment(data='dir', x=new_xy[0], y=new_xy[1], out_name='catch', recursionlimit=15000, xytype='label')

        # Clip the bounding box to the catchment
        grid.clip_to('catch')

        #some sort of strange raster to polygon conversion using rasterio method
        shapes = grid.polygonize()

        #get split Catchment geometry
        print('Split catchment complete')
        split_geom = ogr.Geometry(ogr.wkbPolygon)

        for shape in shapes:
            split_geom = split_geom.Union(ogr.CreateGeometryFromJson(json.dumps(shape[0])))

        #write out shapefile
        self.geom_to_geojson(split_geom, 'splitCatchment')
            
        return split_geom

if __name__=='__main__':

    timeBefore = time.perf_counter()  

    #test site
    point = (-73.82705, 43.29139)

    #start main program
    delineation = Watershed(point[0],point[1])

    timeAfter = time.perf_counter() 
    totalTime = timeAfter - timeBefore
    print("Total Time:",totalTime)