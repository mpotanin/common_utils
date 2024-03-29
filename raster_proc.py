import os
import sys
import osgeo
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import numpy as np
from random import seed
from random import random
from gdal_utils import vector_operations as vop

import time
seed(time.time() * 1000)


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
                            crop_to_cutline = True,
                            pixel_width=None,
                            pixel_height=None,
                            output_bounds=None,
                            output_type = None,
                            dst_srs = None,
                            pixel_res = None,
                            resample_alg = None,
                            src_nodata=None
                              ) :
    in_mem_tiff = os.path.join('/vsimem',str(random()) + '.tif')

    warp(raster_file=raster_file,output_tiff=in_mem_tiff,cutline=cutline,dst_nodata=dst_nodata,src_nodata=src_nodata,
         cutline_where=cutline_where,crop_to_cutline=crop_to_cutline,pixel_width=pixel_width,pixel_height=pixel_height,
         output_type=output_type,output_bounds=output_bounds,dst_srs=dst_srs,pixel_res=pixel_res,resample_alg=resample_alg)

    return in_mem_tiff

def open_clipped_raster_as_image (raster_file, 
                                cutline = None, 
                                dst_nodata = None,
                                cutline_where = None,
                                crop_to_cutline = True,
                                pixel_width = None,
                                pixel_height = None,
                                output_bounds = None,
                                output_type = None,
                                dst_srs = None,
                                pixel_res = None,
                                resample_alg = None,
                                src_nodata = None,
                                  ) :
    in_mem_tiff = None
    if any(arg is not None for arg in [cutline,dst_nodata,pixel_width,pixel_height,dst_srs,pixel_res,resample_alg]):
        in_mem_tiff = get_clipped_inmem_raster (raster_file,
                                                cutline,
                                                dst_nodata,
                                                cutline_where,
                                                crop_to_cutline,
                                                pixel_width,
                                                pixel_height,
                                                output_bounds,
                                                output_type,
                                                dst_srs,
                                                pixel_res,
                                                resample_alg,
                                                src_nodata)

    raster_file = in_mem_tiff if in_mem_tiff is not None else raster_file


    output = None
    gdal_ds = gdal.Open(raster_file)
    if gdal_ds.RasterCount == 1:
        output = gdal_ds.GetRasterBand(1).ReadAsArray()
    else:
        output_type = {gdal.GDT_Byte:np.uint8,
                       gdal.GDT_UInt16:np.uint16,
                       gdal.GDT_Int16:np.int16,
                       gdal.GDT_Int32:np.int32,
                       gdal.GDT_UInt32:np.uint32,
                       gdal.GDT_Float32:np.float32}

        output = np.empty(shape=(gdal_ds.RasterYSize,gdal_ds.RasterXSize,gdal_ds.RasterCount),
                          dtype=output_type[gdal_ds.GetRasterBand(1).DataType])
        for b in range(gdal_ds.RasterCount):
            output[:,:,b] = np.copy(gdal_ds.GetRasterBand(b+1).ReadAsArray())

    gdal_ds = None
    if in_mem_tiff is not None:
        gdal.Unlink(in_mem_tiff)
    return output

def get_raster_bbox (raster_file, t_srs = None) :
    gdal_ds = gdal.Open(raster_file)
    if not gdal_ds :
        return None

    #assume north up image
    geotr = gdal_ds.GetGeoTransform()
    ul_lr = ogr.Geometry(ogr.wkbLineString)
    ul_lr.AddPoint(geotr[0],geotr[3]) #add upper left corner of image
    ul_lr.AddPoint(geotr[0] + geotr[1]*gdal_ds.RasterXSize,
                    geotr[3] + geotr[5]*gdal_ds.RasterYSize ) #calc. lower right corner of image


    if t_srs:
        if int(osgeo.__version__[0]) >= 3:
            # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/154
            t_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        raster_srs=osr.SpatialReference(wkt=gdal_ds.GetProjection())
        coordTrans = osr.CoordinateTransformation(raster_srs,t_srs)
        #ulp.Transform(coordTrans)
        #lrp.Transform(coordTrans)
        ul_lr.Transform(coordTrans)
    gdal_ds = None

    return ul_lr.GetEnvelope()




def generate_virtual_random_tif_path ():
    return '/vsimem/memory_name' + str(random()) + '.tif'


def warp   (raster_file,
            output_tiff,
            cutline = None,
            src_nodata = None,
            dst_nodata = None,
            cutline_where = None,
            crop_to_cutline = True,
            pixel_width=None,
            pixel_height=None,
            output_bounds = None,
            output_type= None,
            dst_srs = None,
            pixel_res = None,
            resample_alg = None
            ):

    return gdal.Warp(output_tiff,
              [raster_file],
              format='GTiff',
              cutlineDSName=cutline,
              cutlineWhere=cutline_where,
              dstNodata=dst_nodata,
              width=pixel_width,
              height=pixel_height,
              cropToCutline= crop_to_cutline,
              outputBounds = output_bounds,
              srcNodata=src_nodata,
              dstSRS=dst_srs,
              resampleAlg=resample_alg,
              outputType = output_type,
              xRes=pixel_res,
              yRes=pixel_res
              )


def get_statistics(single_band_raster):
    gdal_ds = gdal.Open(single_band_raster)
    if gdal_ds is None: return None

    #min,max,mean,stdDev
    stat = gdal_ds.GetRasterBand(1).GetStatistics(False,True)
    gdal_ds = None

    return stat

def array2geotiff(output_geotiff, rasterOrigin, pixel_size, srs, array, nodata_val = None, compress=None):

    if (array.ndim == 2) :
        array = np.copy(array).reshape((array.shape[0],array.shape[1],1))

    bands_num = array.shape[2]
    rows,cols = array.shape[0],array.shape[1]
    originX,originY = rasterOrigin[0],rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    output_type = {np.uint8:gdal.GDT_Byte,
                    np.uint16:gdal.GDT_UInt16,
                    np.int16:gdal.GDT_Int16,
                    np.int32:gdal.GDT_Int32,
                    np.uint32:gdal.GDT_UInt32,
                    np.float32:gdal.GDT_Float32,
                    np.float64:gdal.GDT_Float32}[type(array[0][0][0])]
    options = [f'COMPRESS={compress}'] if compress is not None else None
    outRaster = (driver.Create(output_geotiff, cols, rows, bands_num, output_type) if options is None
                else driver.Create(output_geotiff, cols, rows, bands_num, output_type, options = options) )
    outRaster.SetGeoTransform((originX, pixel_size, 0, originY, 0, -pixel_size))
    
    
    #outRaster.SetSpatialRef(srs) - fails for unknown reason
    prj_wkt = srs.ExportToWkt()
    outRaster.SetProjection(prj_wkt)
    
    for b in range(bands_num):
        outband = outRaster.GetRasterBand(b+1)
        outband.WriteArray(array[:,:,b])
        if nodata_val is not None:
            outband.SetNoDataValue(nodata_val)
        outband.FlushCache()
    outRaster = None


def extract_georeference (raster_file, cutline = None) :#-> tuple[osr.SpatialReference,list]:
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

def calc_ndvi_as_image_from_mem (array_red,
                                 array_nir,
                                 ndv_in = 0,
                                 ndv_out = -10000,
                                 uint8_adjust = False):
    ndv_out = 0 if uint8_adjust else ndv_out

    array_ndvi = np.full(array_red.shape,ndv_out,np.ubyte if uint8_adjust else np.float32)

    red_vec = np.full(array_red.shape[1], fill_value=ndv_in, dtype=np.float32)
    nir_vec = np.full(array_red.shape[1], fill_value=ndv_in, dtype=np.float32)
    tmp1_vec = np.full(array_red.shape[1], fill_value=ndv_in, dtype=np.float32)
    tmp2_vec = np.full(array_red.shape[1], fill_value=ndv_in, dtype=np.float32)

    for i in range( 0, len(array_nir) ):
        np.copyto(dst=red_vec, src=array_red[i])
        np.copyto(dst=nir_vec, src=array_nir[i])
        np.subtract(nir_vec, red_vec, out=tmp1_vec)
        np.add(nir_vec, red_vec, out=tmp2_vec)
        tmp1_vec = np.divide(tmp1_vec, tmp2_vec,
                             out=np.full(tmp1_vec.shape,fill_value=ndv_out, dtype=np.float32),
                             where=np.logical_or(red_vec!=ndv_in,nir_vec!=ndv_in))
        if (uint8_adjust) :
            array_ndvi[i] = np.add(100*tmp1_vec, 101.5, dtype=np.ubyte, casting='unsafe',
                                out=np.zeros_like(tmp1_vec),
                                where=(tmp1_vec!=0))
        else:
            array_ndvi[i] = tmp1_vec
       
    return array_ndvi

def calc_ndvi_as_image (red_band_file, nir_band_file,
                        ndv_in=None, ndv_out=-10000, uint8_adjust = False):
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

    if ndv_in is None:
        ndv_in = 0 if band_nir.GetNoDataValue() is None else band_nir.GetNoDataValue()


    array_ndvi = calc_ndvi_as_image_from_mem(array_red,array_nir,ndv_in,ndv_out, uint8_adjust)
    
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
