# nldi-splitCatchment
The purpose of this application is to prove a method for splitting a local NHD Plus v2 catchment at a click point.  The application inputs are an X,Y coordinate, and the NHD Plus v2 flow direction raster.  It also utilizes several calls to USGS NLDI services.

####  Pre-requisites
* [miniconda](https://docs.conda.io/en/latest/miniconda.html) installed
* Download NHD Plus v2 sample data from [here](https://s3.amazonaws.com/edap-nhdplus/NHDPlusV21/Data/NHDPlusMA/NHDPlusV21_MA_02_02b_FdrFac_01.7z)
* Extract data to `/data` folder
 

## Application setup
Install required packages
```
conda config --add channels conda-forge
conda create -n delineate python=3.7 gdal pysheds
```

##  Get Started
You should now be able to run the script at a predefined sample site by running: 
```
conda activate delineate
git clone 
python ./nldi_delineate.py
```

##  Flask Setup
Install flask dependecies
```
conda install flasl flask-cors
```

Steps to start flask dev server:

* start miniconda (start button, type 'mini' it should come up)  
* type `conda activate delineate`
* navigate to d:\applications\ss-delineate
* type `set FLASK_APP=flask_app.py`
* run flask app `flask run`
