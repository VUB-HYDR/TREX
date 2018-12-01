#     TREX - Tool for Raster data EXploration
#                 (1/12/2018)
#-------------------------------------------------------
# - - - MODULES AND WORKING DIRECTORIES - - - - - - - - -
#-------------------------------------------------------

import os
import gdal
import pandas as pd
import numpy as np
from IPython import get_ipython
import shutil

#-------------------------------------------------------
# - - - READ SETUP - - - - - - - - - - - - - - - - -
#-------------------------------------------------------
current_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(current_dir)

setup_name = "LAI_ProbaV_setup.txt"
read_setup = open(setup_name,'r')
print read_setup

# skip first 5 lines
for i in range(5):
    read_setup.readline()
# read line #6 reference raster
reference_raster = read_setup.readline().split()
dir_input_raster = current_dir + "\\reference_maps\\" + reference_raster[-1]

# skip line #7 
read_setup.readline()
# read line #8 cloud fraction
cloud_fraction = read_setup.readline().split()
f_invalid_px_1 = float(cloud_fraction[-1])
# read line #9 wipe out 
step1 = read_setup.readline().split()
step1 = int(step1[-1])
# skip line #10-11
for i in range(2):
    read_setup.readline()
# read line #12 NDVI.tif
step2 = read_setup.readline().split()
step2 = int(step2[-1])
# read line #13 LAI.tif
step3 = read_setup.readline().split()
step3 = int(step3[-1])
# read line #14 LAI.asc
step4 = read_setup.readline().split()
step4 = int(step4[-1])
# read line #15 monthly LAI.tif
step5 = read_setup.readline().split()
step5 = int(step5[-1])
# read line #16 monthly LAI.asc
step6 = read_setup.readline().split()
step6 = int(step6[-1])
# read line #17 monthly int LAI.asc
step7 = read_setup.readline().split()
step7 = int(step7[-1])

dir_input_maps = current_dir + "\\probaV_download"
dir_step2 = current_dir + "\\main\\1_NDVI_tif"
dir_step3 = current_dir + "\\main\\2_LAI_tif"
dir_step4 = current_dir + "\\main\\3_LAI_asc"
dir_step5 = current_dir + "\\main\\4_monthly_LAI_tif"
dir_step6 = current_dir + "\\main\\5_monthly_LAI_asc"
dir_step7 = current_dir + "\\main\\6_inter_LAI_asc"
temp = current_dir + "\main\\temp"

#-------------------------------------------------------
# - - - FUNCTIONS - - - FUNCTIONS - - - FUNCTIONS - - -
#-------------------------------------------------------

def SearchFolder(directory, file_format):
# Returns a list of files with a certain format (%file_format)  
# in a folder located under a given directory (%directory) 
    MyFolder = os.listdir(directory)
    MyList = []
    a = len(file_format)
    for b in range(len(MyFolder)):
        i = MyFolder[b]
        if i[-a:] == file_format :
            MyList.append(MyFolder[b])
    return MyList

def CopyClearTemp(moveFrom, file_format, moveTo):
# Copies all files with a certain format (%file_format) from the first directory
# (%moveFrom) to the second one (%moveTo). After that, deletes the content of
# the first folder (%moveFrom).    
    tempFiles=SearchFolder(moveFrom, file_format)        
    for i in range (len(tempFiles)):
        in_raster = moveFrom + "\\" + tempFiles[i]
        out_raster = moveTo + "\\" + tempFiles[i]
        shutil.copy2(in_raster, out_raster)          
    files = os.listdir(moveFrom)
    path = moveFrom
    for i in range(len(files)):
        os.remove(path + "/" + files[i])     
        
def PixelsQuality(NDVI,SM,output_folder,filename,f_invalid_px):
# Reads a certain raster map (%NDVI) and status map (%SM). Based on the 8bit
# coding of SM, checks the radiometric quality of each pixel of NDVI map and saves
# a new .tif map (with bad pixels changed to nodata) at given directory (%output_folder)
# using a predifined name (%filename). Function do not process images if a certain
# fraction of the map (%f_invalid_px) is invalid (more than given).
# List of ProbaV pixels codes of SM:   
# 0001 1111 = 248 clear, inland, all radiometric is ok
# 0001 0111 = 232 clear, inland, no SWIR
# 0001 1110 = 120 clear, inland, no BLUE
# 0001 0110 = 104 clear, inland, no SWIR, no BLUE
# 0000 1111 = 240 clear, sea, all radiometric is ok    
# 1001 1111 = 249 inland, ice / snow
# 0111 1111 = 254 inland, cloud
# 0011 1111 = 252 inland, shadow

    REF_source = gdal.Open(dir_input_raster)
    band_info_REF = REF_source.GetRasterBand(1)
    nodata= band_info_REF.GetNoDataValue()    
    
    NDVI_source = gdal.Open(NDVI)
    band_info_NDVI= NDVI_source.GetRasterBand(1)
    xSize = band_info_NDVI.XSize
    ySize = band_info_NDVI.YSize
    geoTrans = NDVI_source.GetGeoTransform()
    wktProjection = NDVI_source.GetProjection() 
    band_Array_NDVI = gdal.Band.ReadAsArray(band_info_NDVI)
    band_Array_NDVI = np.array (band_Array_NDVI, np.int32)
 
    SM_source = gdal.Open(SM)
    band_info_SM = SM_source.GetRasterBand(1)  
    band_Array_SM = gdal.Band.ReadAsArray(band_info_SM)
    band_Array_SM = np.array (band_Array_SM, np.int32)
    
# Less strict    
    SM_invalid_pixels = zip(*np.where ((band_Array_SM != 248) & (band_Array_SM != 232) & (band_Array_SM != 120) & (band_Array_SM != 104)))
# More strict
#    SM_invalid_pixels = zip(*np.where (band_Array_SM != 248))
    counter = 0
    pixels = xSize*ySize
    for i in SM_invalid_pixels:
        band_Array_NDVI[i]=nodata
        counter += 1
    if counter > pixels*f_invalid_px:
        pass
    else:        
        os.chdir (output_folder)
        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(filename,xSize, ySize, 1, gdal.GDT_Float32)
        dataset.SetGeoTransform(geoTrans)
        dataset.SetProjection(wktProjection)
        oBand = dataset.GetRasterBand(1)
        oBand.SetNoDataValue(nodata)
        oBand.WriteArray(band_Array_NDVI)    
        del dataset
        del NDVI_source
        del SM_source
    
def NDVI_conversion(image_input,output_folder,filename):
# Converts digital values into physical values for a certain map (%image_input)
# and saves output at given directory (%output_folder) using a predefined name (%filename)
    data_src = gdal.Open(image_input)
    band_info = data_src.GetRasterBand(1)
    xSize = band_info.XSize
    ySize = band_info.YSize
    nodata=band_info.GetNoDataValue()
    geoTrans = data_src.GetGeoTransform()
    wktProjection = data_src.GetProjection() 
    band_Array = gdal.Band.ReadAsArray(band_info)
    Multiplied_band_Array= (band_Array * 0.004) - 0.08
    Multiplied_band_Array_Nodata=np.where(band_Array == nodata, nodata, Multiplied_band_Array)
    os.chdir (output_folder)
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(filename, xSize, ySize, 1, gdal.GDT_Float32)
    dataset.SetGeoTransform(geoTrans)
    dataset.SetProjection(wktProjection)
    oBand = dataset.GetRasterBand(1)
    oBand.SetNoDataValue(nodata)
    oBand.WriteArray(Multiplied_band_Array_Nodata)
    del dataset
    del data_src

def NDVI_correction(input_path,output_folder,filename):
# Corrects value offset of a certain image (%input_path) using eq. 
# (NDVI*0.004)-0.08 and saves an output at a given directory (%output_folder) using 
# a predefined name (%filename)
    data_src = gdal.Open(input_path)
    band_info = data_src.GetRasterBand(1)
    xSize = band_info.XSize
    ySize = band_info.YSize
    nodata=band_info.GetNoDataValue()
    geoTrans = data_src.GetGeoTransform()
    wktProjection = data_src.GetProjection() 
    band_Array = gdal.Band.ReadAsArray(band_info)
    band_Array_NDVI_less_than_0_Values= zip(*np.where ((band_Array < 0)&(band_Array != nodata)))
    for i in band_Array_NDVI_less_than_0_Values:
        band_Array[i]=0
#        band_Array[i]=nodata
    os.chdir (output_folder)
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(filename, xSize, ySize, 1, gdal.GDT_Float32 )
    dataset.SetGeoTransform(geoTrans)
    dataset.SetProjection(wktProjection)
    oBand = dataset.GetRasterBand(1)
    oBand.SetNoDataValue(nodata)
    oBand.WriteArray(band_Array)
    del dataset
    del data_src
    
def GetCellSize(path):
# Returns the cellsize (x, y) of a given raster (%path)
    data_src = gdal.Open(path)
    data_geo = data_src.GetGeoTransform()
    xres = data_geo[1]
    yres = data_geo[5]
    return xres, yres
    del data_src

def GetExtent(path):
# Returns the extent (xll, yll, xlr, yul) and number of columns and rows 
# (ncols, nrows) of the certain raster (%path). 
# l/u is lower/upper and l/r is a left/right corner   
    data_src = gdal.Open(path)
    data_geo = data_src.GetGeoTransform()
    xll = data_geo[0]
    yul = data_geo[3]
    data_src = data_src.GetRasterBand(1)
    ncols = data_src.XSize
    nrows = data_src.YSize
    xlr = xll + data_geo[1]*ncols
    yll = yul - data_geo[1]*nrows
    return xll, yll, xlr, yul, ncols, nrows
    del data_src
    
def LAI_Map_Tiff(image_input,output_folder,filename):
# Returns an LAI map computed from a certain NDVI map (%image_input) using 
# Su Z.(2000) equation and saves an output at a given directory (%output_folder) 
# using a predefined name (%filename)
    data_src = gdal.Open(image_input)
    band_info = data_src.GetRasterBand(1) 
    xSize = band_info.XSize
    ySize = band_info.YSize
    nodata = band_info.GetNoDataValue()
    geoTrans = data_src.GetGeoTransform()
    wktProjection = data_src.GetProjection() 
    ndvi_map = gdal.Band.ReadAsArray(band_info)
#    masked_ndvi=np.ma.masked_equal(ndvi_map, nodata).mean(axis=0)    
#   silence comment about SQRT... it depends on your python/spyder version
#   or environmental paths...
    np.seterr(divide='ignore', invalid='ignore')    
    #LAI Formula designed by Su
    LAI_map = np.sqrt(ndvi_map * (1+ndvi_map) / (1-ndvi_map))   
    
    save_LAI_map_Nodata=np.where(ndvi_map == nodata, nodata, LAI_map)
    border = gdal.Open(dir_input_raster)
    band_info = border.GetRasterBand(1)
    border_raster = gdal.Band.ReadAsArray(band_info)
    LAI_maps = np.where(border_raster == nodata, nodata, np.array(save_LAI_map_Nodata))

    
    
    os.chdir (output_folder)
    new_filename = filename + "_LAI.tif" 
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(new_filename, xSize, ySize, 1, gdal.GDT_Float32)
    dataset.SetGeoTransform(geoTrans)
    dataset.SetProjection(wktProjection)
    oBand = dataset.GetRasterBand(1)
    oBand.SetNoDataValue(nodata)
#    oBand.WriteArray(LAI_map)
    oBand.WriteArray(LAI_maps)
    
    del dataset

def LAI_Map_Agg(in_raster,output_folder,filename, month, year):
# Aggregates a list of NDVI maps (%in_raster) for a certain month (%month) and
# year (%year) and saves an output at a given directory (%output_folder) 
# using a predefined name (%filename). Also return the list of used for aggregation
# maps.
    list_of_maps = []
    for i in filename:
        if i[0:6] == year+month: list_of_maps.append(i)
#looping through selected month and year!
    if len(list_of_maps) == 0:
        pass
    else:
#create i - dimentional matrix per month
        LAI_maps = []
        for i in range(len(list_of_maps)):
            image_input = in_raster + "\\" + list_of_maps[i] + ".tif"

            data_src = gdal.Open(image_input)
            band_info = data_src.GetRasterBand(1)
            xSize = band_info.XSize
            ySize = band_info.YSize
            nodata = band_info.GetNoDataValue()
            geoTrans = data_src.GetGeoTransform()
            wktProjection = data_src.GetProjection() 
            add_array = np.array(gdal.Band.ReadAsArray(band_info))
            LAI_maps.append(add_array)        
        LAI_maps_unmasked = np.where(LAI_maps == nodata, nodata, np.array(LAI_maps))        
        LAI_maps_masked = np.ma.masked_equal(LAI_maps_unmasked, nodata).mean(axis=0)
        Avg_LAI = np.array(LAI_maps_masked)

        df=pd.DataFrame(Avg_LAI)
        df[df == nodata] = nodata
        filled_df=df.fillna(np.nanmean(df))
        final_df=np.where(filled_df == 0, nodata, filled_df)  
        save_LAI_maps=np.array(final_df) 
        
# a version with a clipped map!        
        border = gdal.Open(dir_input_raster)
        band_info = border.GetRasterBand(1) 
        xSize = band_info.XSize
        ySize = band_info.YSize
        nodata = band_info.GetNoDataValue()
        geoTrans = data_src.GetGeoTransform()
        wktProjection = data_src.GetProjection() 
        border_raster = gdal.Band.ReadAsArray(band_info)
        save_LAI_maps_unmasked = np.where(border_raster == nodata, nodata, np.array(save_LAI_maps))

        os.chdir(output_folder)
        name = list_of_maps[0]
        date = name[:6]
        image_output = output_folder + "\\" + str(date) + "_MonthlyLAI.tif" 

        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(image_output, xSize, ySize, 1, gdal.GDT_Float32)
        dataset.SetGeoTransform(geoTrans)
        dataset.SetProjection(wktProjection)
        oBand = dataset.GetRasterBand(1)
        oBand.SetNoDataValue(nodata)
        
#        oBand.WriteArray(save_LAI_maps)
#        clipped version
        oBand.WriteArray(save_LAI_maps_unmasked)
        del dataset
    return list_of_maps

def readMap(fileName, ncols, nrows, nodata):
# Reads a certain map (%fileName) and returns as a 2d numpy array with predefined
# number of colums (%ncols) and rows (%nrows). Random errors (pixel value = '-1.#IND')
# resulting from processing maps with GDAL are being replaced with a given nodata
# value (%nodata).    
    f = open(fileName,'r')
    for i in range(6):
        f.readline() 
    temp = []
    for i in range(nrows):        
        temp.append(f.readline().split())

    for lines in range(len(temp)):
        a = temp[lines]
# in case there is a random processing error '-1.#IND' in GDAL...
        try:
            for littleBugs in range(len(a)):
                a[a.index('-1.#IND')] = nodata
        except:
            pass    
    result = np.array(temp,float).reshape(nrows, ncols)
    f.close()
    return result

def create_header(from_asc):
# Reads a certain .asc file (%from_asc) and returns its header (createHeader),
# number of columns (ncols), rows (nrows) and nodata (nodata).
    f = open(from_asc,'r')
    header = []
    for i in range(6):
        header.append(f.readline())
    ncols = int(header[0].split()[1])
    nrows = int(header[1].split()[1])
    xll = float(header[2].split()[1])
    yll = float(header[3].split()[1])
    cellsize = float(header[4].split()[1])
    nodata = int(header[5].split()[1])
    f.close()    
    createHeader = 'ncols ' + str(ncols) + "\n" + 'nrows ' + str(nrows) + "\n" + 'xllcorner ' + str(xll) + "\n" + 'yllcorner ' + str(yll) + "\n" + 'cellsize ' + str(cellsize) + "\n" + 'NODATA_value ' + str(nodata)
    return createHeader, ncols, nrows, nodata

def set_nodata(array2d, oldValue, newValue):
# Replaces old nodata values (%oldValue) with new values (%newValue) of a certain
# 2d numpy array (%array2d)
    noData_pixels = zip(*np.where(array2d == oldValue))
    for i in noData_pixels: array2d[i] = newValue
    return array2d

        
#-------------------------------------------------------
# - - - MAIN - - - MAIN - - - MAIN - - - MAIN - - -
#-------------------------------------------------------

#---------------------------------------------
# - - - STEP 1 - - - STEP 1 - - - STEP 1 - - -
#---------------------------------------------
# STEP 1: Delete all files from previous runs. 
# Highly recommended after changing dataset (new maps etc.)
            
if step1 == 1:
#    print '\n STEP 1'
    print '\nDeleting old content...'           
    rerun = [dir_step2, dir_step3, dir_step4, dir_step5, dir_step6, dir_step7, temp]
    for i in range(len(rerun)):
        files = os.listdir(rerun[i])
        path = rerun[i]
        for i in range(len(files)):
            os.remove(path + "//" + files[i])            
else:
    pass
#---------------------------------------------
# - - - STEP 2 - - - STEP 2 - - - STEP 2 - - -
#---------------------------------------------
# STEP 2: Main processing. 
if step2 == 1:    
#    print '\n STEP 2'
#    print 'Checking radiometric and state quality. Discarding images with ' +  str(f_invalid_px_1*100) + '% or more invalid pixels...'
    print '\nChecking radiometric quality...'
    NDVI_list = []
    NDVI_list = SearchFolder(dir_input_maps, 'NDVI.tif')
    SM_list = []
    SM_list = SearchFolder(dir_input_maps, 'SM.tif')    
    fNames = []
    for i in range (len(NDVI_list)):
        breaklist = NDVI_list[i].split(".")
        fn = breaklist [0].split("_")[-4]
        fn = fn + "_NDVI.tif"
        fNames.append(fn)
    for i in range (len(NDVI_list)):
        os.chdir (dir_input_maps)
        Cloud_Free_Image = PixelsQuality(NDVI_list[i],SM_list[i],temp, fNames[i], f_invalid_px_1)   
    CopyClearTemp(temp, '.tif', dir_step2)       
#---------------------------------------------
    print '\nConversion from digital values to physical values...'
    Myfiles22= SearchFolder(dir_step2, 'NDVI.tif')    
# After discarding !!!   
    filenames = []
    dates = []
    for i in range (len(Myfiles22)):
        breaklist = Myfiles22[i].split(".")
        fn = breaklist[0].split("_")[-2]
        dates.append(fn)
        fn = fn + "_NDVI.tif"
        filenames.append(fn)        
    for i in range (len(Myfiles22)):
        os.chdir (dir_step2)
        in_raster = dir_step2 + "\\" + Myfiles22[i]
        NDVI_conversion(in_raster, temp, filenames[i])   
    CopyClearTemp(temp, '.tif', dir_step2)
#---------------------------------------------
#    print '\nCorection for negative values... 
#[VALIDATE WATER PIXELS!!!!]'
    Myfiles23=SearchFolder(dir_step2, 'NDVI.tif')  
    for i in range (len(Myfiles23)):
        os.chdir (dir_step2)
        in_raster = dir_step2 + "\\" + Myfiles23[i]
        NDVI_correction(in_raster, temp, filenames[i])  
    CopyClearTemp(temp, '.tif', dir_step2)   
#---------------------------------------------
#    print '\nResampling images to a higher resolution...'
    Myfiles24=SearchFolder(dir_step2, '.tif')     
    for i in range (len(Myfiles24)):
        os.chdir (dir_step2)
        in_raster = dir_step2 + "\\" + Myfiles24[i]
        out_raster = temp + '\\' + Myfiles24[i]
        CellSize = GetCellSize(in_raster)
        new_xres = CellSize[0]/2
        new_yres = CellSize[1]/2
        cmd= 'gdalwarp -q -multi -of GTiff -co TILED=YES -tr %s %s -overwrite %s %s' % (new_xres, new_yres, in_raster, out_raster)
        os.system (cmd)        
    CopyClearTemp(temp, '.tif', dir_step2)   
#---------------------------------------------
    print '\nReprojecting and resampling NDVI maps ...'    
    Myfiles25=SearchFolder(dir_step2, '.tif')
    
    data_rast = gdal.Open(dir_input_raster)
    t_project = data_rast.GetProjection()
    t_transform = data_rast.GetGeoTransform()
    newt_xres = t_transform[1]
    newt_yres = t_transform[5]
    
    for i in range (len(Myfiles25)):
        os.chdir (dir_step2)
        raster_5 = gdal.Open(dir_step2 + '\\' + Myfiles25[i])
        project_5 = raster_5.GetProjection()
        transform_5 = raster_5.GetGeoTransform()
        in_raster = dir_step2 + "\\" + Myfiles25[i]
        out_raster = temp + '\\' + Myfiles25[i]  
        cmd= 'gdalwarp -q -multi -of GTiff -co TILED=YES -s_srs %s -t_srs %s -tr %s %s -r cubic -overwrite %s %s' % (project_5, t_project, newt_xres, newt_yres, in_raster, out_raster)
        os.system (cmd)    
    del data_rast
#    del raster_5
    CopyClearTemp(temp, '.tif', dir_step2)       
#---------------------------------------------
    print '\nAdjusting extend...'
    Myfiles26=SearchFolder(dir_step2, '.tif')
    for i in range (len(Myfiles26)):
        os.chdir (dir_step2)
        raster_6 = gdal.Open(dir_step2 + '\\' + Myfiles26[i])
        project_6 = raster_6.GetProjection()
        transform_6 = raster_6.GetGeoTransform()
        os.chdir (dir_step2)
        in_raster = dir_step2 + "\\" + Myfiles26[i]
        out_raster = temp + '\\' + Myfiles26[i]
    #finding the extent of ProbaV maps and raster
        Inp_CellSize = GetCellSize(in_raster)
        cellsize = Inp_CellSize[0]
        Inp_Extent = GetExtent(in_raster)
        xll_inp = Inp_Extent[0]
        yll_inp = Inp_Extent[1]
        Ext_Extent = GetExtent(dir_input_raster)
        xll_ext = Ext_Extent[0]
        yll_ext = Ext_Extent[1]
    #computing the shift in x (assuming that both rasters have the same cellsize!)
        x_shift = (xll_ext - xll_inp)/cellsize
        x_shift = x_shift - (round(x_shift) - 1)
        x_shift = cellsize-(cellsize*x_shift)    
    #computing the shift in y (assuming that both rasters have the same cellsize!)
        y_shift = (yll_ext - yll_inp)/cellsize
        y_shift = y_shift - (round(y_shift) - 1)
        y_shift = cellsize-(cellsize*y_shift)    
    #apllying shift to the map extent
        new_xll = xll_inp - x_shift
        new_yll = yll_inp - y_shift
        max_xll = Inp_Extent[2] - x_shift
        max_yll = Inp_Extent[3] - y_shift
        cmd= 'gdal_translate -a_ullr %s %s %s %s -stats %s %s' % (new_xll,max_yll,max_xll,new_yll, in_raster, out_raster)
        os.system (cmd)         
    print 'x shift: ' + str(x_shift) + ' [m]'
    print 'y shift: ' + str(y_shift) + ' [m]'    
    del raster_6
#    raster_6.close()
    CopyClearTemp(temp, '.tif', dir_step2)
#---------------------------------------------
    print '\nClipping and generating NDVI.tif maps...'
    Myfiles27=SearchFolder(dir_step2, '.tif')
    #gdaltindex clipper.shp clipshapeRaster.tif
    cutline = temp + "\\cutline.shp"
    cmd = 'gdaltindex %s %s' % (cutline, dir_input_raster)
    os.system (cmd)         
    for i in range (len(Myfiles27)):
        raster_7 = gdal.Open(dir_step2 + '\\' + Myfiles27[i])
        project_7 = raster_7.GetProjection()
        transform_7 = raster_7.GetGeoTransform()
        in_raster = dir_step2 + "\\" + Myfiles27[i]
        out_raster = temp + "\\" + Myfiles27[i]
        cmd= 'gdalwarp -q -multi -of GTiff -co TILED=YES -cutline %s -crop_to_cutline %s %s' % (cutline, in_raster, out_raster)
        os.system (cmd)   
    del raster_7
#    raster_7.close()
    CopyClearTemp(temp, '.tif', dir_step2)    
else:
    pass
#---------------------------------------------
# - - - STEP 3 - - - STEP 3 - - - STEP 3 - - -
#--------------------------------------------- 
if step3 == 1:
    print '\nGenerating LAI.tif maps...'    
    Myfiles3=SearchFolder(dir_step2, '.tif')
    dates = []
    for i in range (len(Myfiles3)):
        breakList = Myfiles3[i].split("\\")
        fn = breakList [-1]
        fn2 = fn.split("_")
        fn2 = fn2[0]
        dates.append(fn2)
    
    for i in range (len(Myfiles3)):
        in_raster = dir_step2 + "\\" + Myfiles3[i]
        NDVI_source = gdal.Open(in_raster)
        band_info = NDVI_source.GetRasterBand(1)
        nodata = band_info.GetNoDataValue()
        xSize = band_info.XSize
        ySize = band_info.YSize
        LAI_Map_Tiff(in_raster,dir_step3,dates[i])        
        del NDVI_source
else:
    pass
#---------------------------------------------
# - - - STEP 4 - - - STEP 4 - - - STEP 4 - - -
#--------------------------------------------- 
if step4 == 1:
    print '\nGenerating LAI.asc maps...'
    Myfiles4=SearchFolder(dir_step3, '.tif')
    os.chdir(dir_step4)
    dates = []
    for i in range (len(Myfiles4)):
        breakList = Myfiles4[i].split("\\")
        fn = breakList [-1]
        fn2 = fn.split("_")
        fn2 = fn2[0]
        dates.append(fn2)
        
    for i in range (len(Myfiles4)):
        in_raster = dir_step3 + "\\" + Myfiles4[i]
        out_raster = dir_step4 + "\\" + dates[i] + "_LAI.asc"
        if os.path.exists(out_raster):
            os.remove(out_raster)
        cmd= 'gdal_translate -q -of AAIGrid %s %s' % (in_raster, out_raster)    
        os.system (cmd)    
else:
    pass
#---------------------------------------------
# - - - STEP 5 - - - STEP 5 - - - STEP 5 - - -
#--------------------------------------------- 
if step5 == 1:
#    print '\n Step 5'
    print '\nGenerating monthly aggregated LAI.tif maps... '
    Myfiles5=SearchFolder(dir_step3, '.tif')  
    data_src = gdal.Open(dir_step3 + "\\" + Myfiles5[0])
    band_info = data_src.GetRasterBand(1)  
    filenames = []
    dates = []
    for i in range (len(Myfiles5)):
        breakList = Myfiles5[i].split("\\")
        fn = breakList [-1]
        fn = fn[:-4]
        filenames.append(fn)
        fn2 = fn.split("_")
        fn2 = fn2[0]
        dates.append(fn2)
    
    years = []
    months = ['01', '02', '03', '04','05', '06','07', '08','09', '10','11', '12']
    for j in range (len(dates)):
        date = dates[j]
        year = date[:4]
        if year not in years: 
            years.append(year)
    in_raster = dir_step3
    out_raster = dir_step5    
    
    print '\nStatus report - number of available products for each month'
    os.chdir(dir_step5)
    report = open('report.txt','w')
    for i in range(len(years)):
        report.write('\n ===================================================== ')
        report.write('\n List of maps used for aggregation in the YEAR: ' + str(years[i]) + '\n')
        name_of_months = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
        for k in range(12):                
            add_list = LAI_Map_Agg(in_raster, out_raster, filenames, months[k], years[i])
            try: 
                if len(add_list) != 0:
                    print add_list
                    report.write('\n ' + str(name_of_months[k]) + ':')
                    for a in add_list:
                        b = ' ' + str(a)
                        report.write(b)
            except: pass
    report.write('\n ===================================================== ')
    print '\nReport saved as report.txt in the file 4_monthly_LAI_tif'
    report.close()
else:
    pass
#---------------------------------------------
# - - - STEP 6 - - - STEP 6 - - - STEP 6 - - -
#--------------------------------------------- 
if step6 == 1:
#    print '\n Step 6'
    print '\nGenerating aggregated LAI.asc maps...'
    Myfiles6=SearchFolder(dir_step5, '.tif')    
    filenames = []
    dates = []
    for i in range (len(Myfiles6)):
        breakList = Myfiles6[i].split("\\")
        fn = breakList [-1]
        fn = fn[:-4]
        filenames.append(fn)
        fn2 = fn.split("_")
        fn2 = fn2[0]
        dates.append(fn2)    
    for i in range (len(Myfiles6)):
        in_raster = dir_step5 + "\\" + Myfiles6[i]
        out_raster = dir_step6 + "\\" + dates[i] + "_MonthlyLAI.asc"
        
        if os.path.exists(out_raster):
            os.remove(out_raster)
        cmd= 'gdal_translate -q -of AAIGrid %s %s' % (in_raster, out_raster)
        os.system (cmd)    
else:
    pass
#---------------------------------------------
# - - - STEP 7 - - - STEP 7 - - - STEP 7 - - -
#--------------------------------------------- 
if step7 == 1:   
#    print '\n Step 7'
    print '\nInterpolating between monthly LAI.asc maps...'

    Myfiles7=SearchFolder(dir_step6, '.asc')    
    filenames = []
    dates = []
    
    for i in range (len(Myfiles7)):
        breakList = Myfiles7[i].split("\\")
        fn = breakList [-1]
        fn = fn[:-4]
        filenames.append(fn)
        fn2 = fn.split("_")
        fn2 = fn2[-1]
        dates.append(fn2) 

    reference = gdal.Open(dir_input_raster)
    band_info = reference.GetRasterBand(1) 
    reference_nodata = band_info.GetNoDataValue()
    nodata_mask = gdal.Band.ReadAsArray(band_info)    
    
    for j in range (len(Myfiles7)):
        os.chdir(dir_step6)
        in_raster = dir_step6 + "\\" + Myfiles7[j]
        get_header = create_header(in_raster)
        nodata = get_header[3]
        raster_n = readMap(Myfiles7[j], get_header[1], get_header[2], get_header[3])
        raster_n[nodata_mask == reference_nodata] = np.nan
        #check if raster_n has nodata
        noData_pixels = zip(*np.where(raster_n == nodata))        
        if len(noData_pixels) == 0:
            os.chdir(dir_step7)
            raster_n[nodata_mask == reference_nodata] = nodata                       
            np.savetxt('int_LAI_' + dates[j] + '.asc', raster_n, fmt='%10.2f', header=get_header[0], comments='')
        else:
            print '\nInterpolating ' + str(len(noData_pixels)) + ' values of ' + str(Myfiles7[j])
            os.chdir(dir_step6)
            raster_n_minus1 = readMap(Myfiles7[j-1], get_header[1], get_header[2], get_header[3])
            raster_n_minus2 = readMap(Myfiles7[j-2], get_header[1], get_header[2], get_header[3])
            try: raster_n_plus1 = readMap(Myfiles7[j+1], get_header[1], get_header[2], get_header[3])
            except: pass
            try: raster_n_plus2 = readMap(Myfiles7[j+2], get_header[1], get_header[2], get_header[3])
            except: pass

            raster_nan = raster_n
            raster_nan[raster_nan == nodata] = np.nan
            filler = np.nanmean(raster_nan)
            # trick will work if nodata is a negative value AND a certain pixel 
            # has at least one value per 5 months... therefore...
            
            raster_n_minus1 = set_nodata(raster_n_minus1, nodata, -9999)
            raster_n_minus2 = set_nodata(raster_n_minus2, nodata, -9999)
            
            try: raster_n_plus1 = set_nodata(raster_n_plus1, nodata, -9999)
            except : pass
            
            try: raster_n_plus2 = set_nodata(raster_n_plus2, nodata, -9999)
            except : pass
            
            for k in noData_pixels:
                try:
                    if (raster_n_plus1[k] + raster_n_minus1[k])/2 > 0:
                        raster_n[k] = (raster_n_plus1[k] + raster_n_minus1[k])/2
                    else:
                        if ((raster_n_plus1[k] == -9999) & (raster_n_plus2[k] != -9999) & (raster_n_minus1[k] != -9999)):
                            raster_n[k] = (raster_n_plus2[k] + raster_n_minus1[k])/2
                        elif ((raster_n_plus1[k] != -9999) & (raster_n_minus1[k] == -9999) & (raster_n_minus2[k] != -9999)):
                            raster_n[k] = (raster_n_plus1[k] + raster_n_minus2[k])/2
                        elif ((raster_n_plus2[k] != -9999) & (raster_n_minus2[k] != -9999)):
                            raster_n[k] = (raster_n_plus2[k] + raster_n_minus2[k])/2                        
                        else:
                            raster_n[k] = filler
                except: raster_n[k] = filler
            os.chdir(dir_step7)
            raster_n[nodata_mask == reference_nodata] = nodata
            np.savetxt('int_' + str(Myfiles7[j]) + '.asc', raster_n, fmt='%10.2f', header=get_header[0], comments='')
#---------------------------------------------

print '\n PROCESSING COMPLETE'
