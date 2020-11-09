import delineate # The code to test
import unittest
import json

meters_to_miles = 0.00000038610

class ny_regular_tests(unittest.TestCase):
    def test_ny_local(self):
        results = delineate.Watershed(44.00683,-73.74586,'ny','c:/temp/')

        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 0.63)

    def test_ny_local_on_str_no_adjoint(self):
        results = delineate.Watershed(44.03115617578595 , -73.71244903077174,'ny','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 8.37)

    def test_ny_localGlobal(self):
        results = delineate.Watershed(44.00431,-73.71348,'ny','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 53.59)

    def test_ny_global_non_nested(self):
        results = delineate.Watershed(43.29139,-73.82705,'ny','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 2719.84)

    def test_ny_global_nested_small(self):
        results = delineate.Watershed(42.17209,-73.87555,'ny','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 9996.5)

    def test_ny_global_nested_large(self):
        results = delineate.Watershed(41.00155,-73.89282,'ny','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 13280.1)

class Test_ny_edge_tests(unittest.TestCase):
    def test_ny_bad_catchment_geom(self):
        results = delineate.Watershed(43.45338620107029 , -74.50329065322877,'ny','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 17.47)

    def test_ny_bad_split_catchment(self):
        results = delineate.Watershed(41.310936704746936 , -74.51668024063112,'ny','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 157.78)

    def test_ny_spatially_disconnected_hucs(self):
        results = delineate.Watershed(43.31392194207697 , -73.8442397117614,'ny','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 1053.99)

    def test_ny_bad_single_upstream_huc(self):
        results = delineate.Watershed(42.82741831644657 , -73.93358945846559,'ny','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 3312.89)

class Test_ma_regular_tests(unittest.TestCase):
    def test_ma_local(self):
        results = delineate.Watershed(42.55415514336133 , -71.50470972061159,'ma','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 0.29)

    def test_ma_local_on_str_no_adjoint(self):
        results = delineate.Watershed(42.567036149018534 , -71.5099883079529,'ma','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 6.93)

    def test_ma_localGlobal(self):
        results = delineate.Watershed(42.5975606206426 , -71.42810583114625,'ma','c:/temp/')

        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 34.71)

    def test_ma_global_non_nested(self):
        results = delineate.Watershed(42.809522580037786 , -71.47342443466188,'ma','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 3432.86)

    def test_ma_global_nested_small(self):
        results = delineate.Watershed(42.645371184571744 , -71.29590511322023,'ma','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 4628.4)

    def test_ma_global_nested_large(self):
        results = delineate.Watershed(42.48847603472268 , -72.57074832916261,'ma','c:/temp/')
        
        area = round(results.mergedCatchment['properties']['area']*meters_to_miles,2)   
        self.assertEqual(area, 7918.2)

if __name__ == '__main__':
    unittest.main()