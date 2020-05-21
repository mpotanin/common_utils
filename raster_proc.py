import os
import sys
import osgeo
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import numpy as np
from random import seed
from random import random
from common_utils import vector_operations as vop

# seed random number generator
seed(1)

def get_raster_bbox (raster_file, t_srs = None) :
    gdal_ds = gdal.Open(raster_file)
    if not gdal_ds :
        return None
    geotr = gdal_ds.GetGeoTransform()


    ulp = ogr.Geometry(ogr.wkbPoint)
    ulp.AddPoint(geotr[0],geotr[3])
    lrp = ogr.Geometry(ogr.wkbPoint)
    lrp.AddPoint(geotr[0] + geotr[1]*gdal_ds.RasterXSize,
                    geotr[3] - geotr[1]*gdal_ds.RasterYSize )

    if t_srs:
        if int(osgeo.__version__[0]) >= 3:
            # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/154
            t_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        raster_srs=osr.SpatialReference(wkt=gdal_ds.GetProjection())
        coordTrans = osr.CoordinateTransformation(raster_srs,t_srs)
        ulp.Transform(coordTrans)
        lrp.Transform(coordTrans)
    
    gdal_ds = None

    return vop.BBOX(ulp.GetX(),lrp.GetY(),lrp.GetX(),ulp.GetY())



def generate_virtual_random_tif_path ():
    return '/vsimem/memory_name' + str(random()) + '.tif'


def crop_raster_file_to_cutline (input_raster, output_tiff, vector_file, src_ndv=None, dst_ndv = None):
    return gdal.Warp(output_tiff,
                input_raster,
                format = 'GTiff',
                cutlineDSName = vector_file,
                cropToCutline = True,
                srcNodata=src_ndv,
                dstNodata=dst_ndv
                )



def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,prj_wkt,array):

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Byte)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRaster.SetProjection(prj_wkt)
    outband.FlushCache()

def create_ndvi_file (red_band_file, nir_band_file, output_file) :
    #create output:
    # - the same srs, pixel size as input bands
    # - pixel size = byte
    # - ndvi formulae = 101 + 100*ndvi
    # open dataset
    ds_nir = gdal.Open(nir_band_file)
    ds_red = gdal.Open(red_band_file)

    if ds_nir is None:
        print('ERROR: can\'t open file: ' + nir_band_file)
        exit(1)

    if ds_red is None:
        print('ERROR: can\'t open file: ' + red_band_file)
        exit(1)

    array_ndvi = np.full((ds_nir.RasterYSize,ds_nir.RasterXSize),0,np.ubyte)
    band_nir = ds_nir.GetRasterBand(1)
    band_red = ds_red.GetRasterBand(1)
    array_nir = band_nir.ReadAsArray()
    array_red = band_red.ReadAsArray()

    prj_wkt=ds_nir.GetProjection()
    geotr = ds_nir.GetGeoTransform()

    
    red_vec = np.full(shape=(ds_nir.RasterXSize), fill_value=0, dtype=np.float)
    nir_vec = np.full(shape=(ds_nir.RasterXSize), fill_value=0, dtype=np.float)
    tmp1_vec = np.full(shape=(ds_nir.RasterXSize), fill_value=0, dtype=np.float)
    tmp2_vec = np.full(shape=(ds_nir.RasterXSize), fill_value=0, dtype=np.float)
    #tmp3_vec = np.full(shape=(ds_nir.RasterXSize), fill_value=0, dtype=np.ubyte)
    print ((len(array_ndvi),len(array_nir)))
    for i in range( 0, len(array_nir) ):
        np.copyto(dst=red_vec, src=array_red[i])
        np.copyto(dst=nir_vec, src=array_nir[i])
        np.subtract(nir_vec, red_vec, out=tmp1_vec)
        np.multiply(tmp1_vec, 100.0, out=tmp1_vec)
        np.add(nir_vec, red_vec, out=tmp2_vec)
        tmp1_vec = np.divide(tmp1_vec, tmp2_vec, out=np.zeros_like(tmp1_vec), where=(tmp2_vec != 0))
        #np.add(tmp1_vec, 101.5, dtype=np.ubyte, casting='unsafe', out=array_ndvi[i])
        array_ndvi[i] = np.add(tmp1_vec, 101.5, dtype=np.ubyte, casting='unsafe', 
                                out=np.zeros_like(array_ndvi[i]), 
                                where=(nir_vec!=0))
        #tmp3_vec = np.zeros_like(tmp3_vec)
        #np.copto(dst=tmp3_vec, src=np._ones_like(tmp3_vec), where=())

       


    """
    for i in range(0,len(array_nir)):
        for j in range(0,len(array_nir[i])):
            if (not array_nir[i][j]==0) and (not array_red[i][j]==0):
                array_ndvi[i][j] = (101.5 + 100.0*(float(array_nir[i][j])-float(array_red[i][j]))/
                                                    (float(array_red[i][j])+float(array_nir[i][j])))
            else: array_ndvi[i][j] = 0
              
    """
    
    array2raster(output_file,[geotr[0],geotr[3]],geotr[1],-geotr[1],prj_wkt,array_ndvi)

    # close dataset
    ds_red = None
    ds_nir = NotImplemented
    array_ndvi = None
    array_nir = None
    array_red = None
