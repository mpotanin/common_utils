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

def preview2geotiff (input_preview, output_geotiff, long_min, lat_min, long_max, lat_max, utm_zone = None):
    #calc EPSG code based on UTM zone
    if utm_zone is None:
        utm_zone = int((0.5*(long_min + long_max) + 180)/6) + 1
    
    epsg_code = 100*(326 if (lat_min + lat_max >= 0) else 327) + utm_zone
    

    #calc UTM XY BBOX
    latlongSpatialRef = osr.SpatialReference()
    latlongSpatialRef.ImportFromEPSG(4326)

    utmSpatialRef = osr.SpatialReference()
    utmSpatialRef.ImportFromEPSG(epsg_code)

    coordTrans = osr.CoordinateTransformation(latlongSpatialRef, utmSpatialRef)

    point_ll = ogr.Geometry(ogr.wkbPoint)
    point_ll.AddPoint(lat_min,long_min)
    point_ll.Transform(coordTrans)

    point_ur = ogr.Geometry(ogr.wkbPoint)
    point_ur.AddPoint(lat_max,long_max)
    point_ur.Transform(coordTrans)
            
    #N1 calc pixel size
    gdal_ds = gdal.Open(input_preview)
    w,h = gdal_ds.RasterXSize,gdal_ds.RasterYSize

    pixel_size = (0.5 *  ( ( (point_ur.GetX()-point_ll.GetX() )/w) + 
                          ( (point_ur.GetY()-point_ll.GetY() )/h) ) )
    
    #img = gdal_ds.ReadAsArray()
    array2geotiff(output_geotiff,
                    (point_ll.GetX(),point_ur.GetY()),
                    pixel_size,
                    utmSpatialRef.ExportToWkt(),
                    gdal_ds.ReadAsArray())

def get_clipped_inmem_raster (raster_file, 
                            cutline = None, 
                            dst_nodata = None,
                            cutline_where = None,
                            crop_to_cutline = True) :
    in_mem_tiff = os.path.join('/vsimem',str(random()) + '.tif')
    gdal.Warp(in_mem_tiff,
            [raster_file],
            format = 'GTiff',
            cutlineDSName = cutline,
            cutlineWhere = cutline_where,
            dstNodata = dst_nodata,
            cropToCutline = crop_to_cutline,
            srcNodata=0
            )
    return in_mem_tiff

def open_clipped_raster_as_image (raster_file, 
                                cutline = None, 
                                dst_nodata = None, 
                                cutline_where = None,
                                crop_to_cutline = True) :
    
    in_mem_tiff = get_clipped_inmem_raster (raster_file,
                                            cutline,
                                            dst_nodata,
                                            cutline_where,
                                            crop_to_cutline)
    ds_obj = gdal.Open(in_mem_tiff)
    band_obj = ds_obj.GetRasterBand(1)
    img_obj = band_obj.ReadAsArray()
    gdal.Unlink(in_mem_tiff)
    return img_obj

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


def crop_to_cutline (input_raster, 
                    output_tiff, 
                    vector_file, 
                    src_nodata=None, 
                    dst_nodata = None,
                    cutline_where = None,
                    crop_to_cutline = True):

    return gdal.Warp(output_tiff,
                input_raster,
                format = 'GTiff',
                cutlineDSName = vector_file,
                cropToCutline = crop_to_cutline,
                srcNodata=src_nodata,
                dstNodata=dst_nodata,
                cutlineWhere = cutline_where
                )



def array2geotiff(output_geotiff,rasterOrigin,pixel_size,srs,array,nodata_val = None):

    array_ref = None
    if (array.ndim == 2) :
        array_ref = list()
        array_ref.append(array)
    else: array_ref = array
   
    bands_num = 1 if (array.ndim==2) else array.shape[0]

    rows,cols = array_ref[0].shape[0],array_ref[0].shape[1]
    originX,originY = rasterOrigin[0],rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')

    output_type = {np.uint8:gdal.GDT_Byte,
                    np.uint16:gdal.GDT_UInt16,
                    np.int16:gdal.GDT_Int16,
                    np.int32:gdal.GDT_Int32,
                    np.uint32:gdal.GDT_UInt32,
                    np.float32:gdal.GDT_Float32,
                    np.float64:gdal.GDT_Float32}[type(array_ref[0][0][0])]
     

    outRaster = driver.Create(output_geotiff, cols, rows, bands_num, output_type)
    outRaster.SetGeoTransform((originX, pixel_size, 0, originY, 0, -pixel_size))
    
    
    #outRaster.SetSpatialRef(srs) - fails for unknown reason
    prj_wkt = srs.ExportToWkt()
    outRaster.SetProjection(prj_wkt)
    
    for b in range(bands_num):
        outband = outRaster.GetRasterBand(b+1)
        outband.WriteArray(array_ref[b])
        if nodata_val is not None:
            outband.SetNoDataValue(nodata_val)
        outband.FlushCache()
    array_ref = None 
    outRaster = None


def extract_georeference (raster_file, cutline = None):
    if cutline is not None:
        in_mem_tif = os.path.join('/vsimem',str(random()) + '.tif')
        gdal.Warp(in_mem_tif,
                [raster_file],
                format = 'GTiff',
                cutlineDSName = cutline,
                cropToCutline = True,
                srcNodata=0,
                dstNodata=0
                )

        raster_file = in_mem_tif

    
    gdal_ds = gdal.Open(raster_file)
    
    #srs = gdal_ds.GetSpatialRef().Clone() - fails for unknown reason
    srs = osr.SpatialReference(gdal_ds.GetProjection())

        
    geotransform = gdal_ds.GetGeoTransform()
    gdal_ds = None
    if cutline is not None:
        gdal.Unlink(raster_file)

    return  (srs,geotransform)

def calc_ndvi_as_image_from_mem (array_red, array_nir, uint8_adjust = True):

    array_ndvi = np.full(array_red.shape,0,np.ubyte if uint8_adjust else np.float)

    red_vec = np.full(array_red.shape[1], fill_value=0, dtype=np.float)
    nir_vec = np.full(array_red.shape[1], fill_value=0, dtype=np.float)
    tmp1_vec = np.full(array_red.shape[1], fill_value=0, dtype=np.float)
    tmp2_vec = np.full(array_red.shape[1], fill_value=0, dtype=np.float)
    #tmp3_vec = np.full(shape=(ds_nir.RasterXSize), fill_value=0, dtype=np.ubyte)
    #print ((len(array_ndvi),len(array_nir)))
    for i in range( 0, len(array_nir) ):
        np.copyto(dst=red_vec, src=array_red[i])
        np.copyto(dst=nir_vec, src=array_nir[i])
        np.subtract(nir_vec, red_vec, out=tmp1_vec)
        np.multiply(tmp1_vec, 100.0, out=tmp1_vec)
        np.add(nir_vec, red_vec, out=tmp2_vec)
        tmp1_vec = np.divide(tmp1_vec, tmp2_vec, out=np.zeros_like(tmp1_vec), where=(tmp2_vec != 0))
        if (uint8_adjust) :
            array_ndvi[i] = np.add(tmp1_vec, 101.5, dtype=np.ubyte, casting='unsafe', 
                                out=np.zeros_like(array_ndvi[i]), 
                                where=(nir_vec!=0))
        else:
            array_ndvi[i] = (tmp1_vec*0.01)
       
    return array_ndvi

def calc_ndvi_as_image (red_band_file, nir_band_file, uint8_adjust = True):
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

    #array_ndvi = np.full((ds_nir.RasterYSize,ds_nir.RasterXSize),0,np.ubyte)
    band_nir = ds_nir.GetRasterBand(1)
    band_red = ds_red.GetRasterBand(1)
    array_nir = band_nir.ReadAsArray()
    array_red = band_red.ReadAsArray()

    array_ndvi = calc_ndvi_as_image_from_mem(array_red,array_nir,uint8_adjust)
    
    # close dataset
    ds_red = None
    ds_nir = None
    array_nir = None
    array_red = None

    return array_ndvi


def create_ndvi_uint8_file (red_band_file, nir_band_file, output_file) :
    #create output:
    # - the same srs, pixel size as input bands
    # - pixel size = byte
    # - ndvi formulae = 101 + 100*ndvi
    # open dataset
    ds_nir = gdal.Open(nir_band_file)

    if ds_nir is None:
        print('ERROR: can\'t open file: ' + nir_band_file)
        exit(1)

    prj_wkt=ds_nir.GetProjection()
    geotr = ds_nir.GetGeoTransform()

    ndvi_img = calc_ndvi_as_image(red_band_file, nir_band_file,True)
    
    array2geotiff(output_file,[geotr[0],geotr[3]],geotr[1],prj_wkt,ndvi_img)

    # close dataset
    ds_nir = None
    ndvi_img = None