import sys
import os
import array
import json
import csv
from io import BytesIO
import requests
import osgeo
from osgeo import gdal
from osgeo import ogr
from osgeo import osr

from random import seed
from random import random
# seed random number generator
seed(1)



class BBOX :
    """
    Represents bounding box rectangular.
    """
    def __init__(self,minx=1e+100,miny=1e+100,maxx=-1e+100,maxy=-1e+100) :
        self.minx = minx
        self.miny = miny
        self.maxx = maxx
        self.maxy = maxy

    def is_undefined (self) :
        if (self.minx>self.maxx) : return True
        else : return False

    def merge (self, minx,miny,maxx,maxy) :
        self.minx = min(self.minx,minx)
        self.miny = min(self.miny,miny)
        self.maxx= max(self.maxx,maxx)
        self.maxy = max(self.maxy,maxy)

    def merge_envp (self,envp) :
        self.merge(envp[0],envp[2],envp[1],envp[3])
    
    def create_ogrpolygon (self):
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(self.minx,self.miny)
        ring.AddPoint(self.maxx,self.miny)
        ring.AddPoint(self.maxx,self.maxy)
        ring.AddPoint(self.minx,self.maxy)
        ring.AddPoint(self.minx,self.miny)

        # Create polygon
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        return poly

    @staticmethod
    def calc_BBOX_from_vector_file (vector_file):
        vec_ds = gdal.OpenEx(vector_file, gdal.OF_VECTOR )
        if vec_ds is None:
            print ("Open failed: " + vector_file)
            sys.exit( 1 )
        lyr = vec_ds.GetLayer(0)
        lyr.ResetReading()
        feat = lyr.GetNextFeature()
        bbox = BBOX()
        while feat is not None:
            geom = feat.GetGeometryRef()
            bbox.merge_envp(geom.GetEnvelope())
            feat = lyr.GetNextFeature()
            # do something more..
        feat = None
        return bbox

class VectorFile:

    @staticmethod
    def get_all_features (vector_file, t_srs = None, column_name_filt = None, value_filt = None):
        ds = gdal.OpenEx(vector_file, gdal.OF_VECTOR )
        #ds = ogr.Open(vector_file)
        if ds is None:
            print ("Open failed: " + vector_file)
            sys.exit( 1 )

        lyr = ds.GetLayer(0)
                
        #source_srs = lyr.GetSpatialRef()

        srs_4326 = osr.SpatialReference()
        srs_4326.ImportFromEPSG(4326)
        if int(osgeo.__version__[0]) >= 3:
            # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/154
            srs_4326.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

        coordTrans = 0
        t_srs = srs_4326 if not t_srs else t_srs
        s_srs = srs_4326 if not lyr.GetSpatialRef() else lyr.GetSpatialRef()
        coordTrans = osr.CoordinateTransformation(s_srs, t_srs)

        """
        if not t_srs :
            t_srs = osr.SpatialReference()
            t_srs.ImportFromEPSG(4326)
        
        if int(osgeo.__version__[0]) >= 3:
            # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/154
            t_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
       
        if lyr.GetSpatialRef(): 
            coordTrans = osr.CoordinateTransformation(lyr.GetSpatialRef(), 
                                                        t_srs)
        else :
            s_srs = osr.SpatialReference()
            s_srs.ImportFromEPSG(4326)
            coordTrans = osr.CoordinateTransformation(lyr.GetSpatialRef(), 
                                                        t_srs)
        """
    
        output = list()

        feat_count = lyr.GetFeatureCount()
        lyr.ResetReading()
        for fid in range(0,feat_count):
            feat = lyr.GetNextFeature()
            if column_name_filt is not None and value_filt is not None:
                if feat.GetFieldAsString(column_name_filt) != value_filt:
                    continue
            geom = feat.GetGeometryRef()
            geom.Transform(coordTrans)
            output.append(feat)

        ds = None
        return output  

    @staticmethod
    def get_all_geometry_in_wkt (vector_file, t_srs = None, column_name_filt = None, value_filt = None):

        features = VectorFile.get_all_features(vector_file,t_srs,column_name_filt,value_filt)
        output = list()

        fid = 0
        for feat in features:
            geom = feat.GetGeometryRef()
            output.append((fid,geom.ExportToWkt()))
            fid+=1

        ds = None
        return output  

     

    @staticmethod
    def get_srs_from_file (vector_file):
        ds =  gdal.OpenEx(vector_file, gdal.OF_VECTOR )
        lyr = ds.GetLayer(0)
        srs=lyr.GetSpatialRef()
        ds = None
        return srs.Clone()

    @staticmethod
    def create_vector_file (wkt_polygeom, filename, srs=None):
        driver_name = 'ESRI Shapefile' if filename.endswith('.shp') else 'GeoJSON'
        driver = ogr.GetDriverByName(driver_name)
        if os.path.exists(filename):
            driver.DeleteDataSource(filename)
        outDataSource = driver.CreateDataSource(filename)

        srs_def = 0
        if srs is None :
            srs_def = osr.SpatialReference()
            srs_def.ImportFromEPSG(4326)
        else: 
            srs_def = srs
        
        outlayer = outDataSource.CreateLayer('layer', srs_def, geom_type=ogr.wkbMultiPolygon)

        id_field = ogr.FieldDefn("id", ogr.OFTInteger)
        outlayer.CreateField(id_field)
        
        feature = ogr.Feature(outlayer.GetLayerDefn())
        feature.SetGeometry(ogr.CreateGeometryFromWkt(wkt_polygeom))
        feature.SetField("id", 1)

        outlayer.CreateFeature(feature)
        feature = None
        outDataSource = None
        return True


    @staticmethod
    def create_virtual_vector_file (wkt_polygeom, srs=None) :
        virtual_path = '/vsimem/memory_name' + str(random()) + '.shp'
        return virtual_path if VectorFile.create_vector_file(wkt_polygeom,virtual_path,srs) else ''

    @staticmethod
    def remove_virtual_vector_file (in_mem_vec_file:str):
        gdal.Unlink(in_mem_vec_file)
        if in_mem_vec_file.endswith('.shp'):
            gdal.Unlink(in_mem_vec_file.replace('.shp','.shx'))
            gdal.Unlink(in_mem_vec_file.replace('.shp', '.prj'))
            gdal.Unlink(in_mem_vec_file.replace('.shp', '.dbf'))

#def create_virtual_vector_ds_from_bbox (bbox,srs):
#    return create_virtual_vector_ds_from_polygon(srs)