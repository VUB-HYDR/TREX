#   TREX - Tool for Raster data EXploration
#   last update by Joanna Suliga(09/09/2019)
#-------------------------------------------------------
# - - - MODULES AND WORKING DIRECTORIES - - - - - - - - -
#-------------------------------------------------------
import os, gdal, shutil
import pandas as pd
import numpy as np
from IPython import get_ipython
#-------------------------------------------------------
# - - - READ SETUP - - - - - - - - - - - - - - - - -
#-------------------------------------------------------
current_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(current_dir)

setup_name = "TREX_setup.txt"
read_setup = open(setup_name,'r')
print (read_setup)

# skip first 5 lines
for i in range(5):
    read_setup.readline()
# read line #6 reference raster
reference_raster = read_setup.readline().split()
dir_reference_map = current_dir + "/reference_maps/" + reference_raster[-1]

# skip line #7 
read_setup.readline()
# read line #8 cloud fraction
cloud_fraction = read_setup.readline().split()
f_invalid_px = float(cloud_fraction[-1])
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

dir_input_maps = current_dir + "/probaV_download"
dir_step2 = current_dir + "/main/1_NDVI_tif"
dir_step3 = current_dir + "/main/2_LAI_tif"
dir_step4 = current_dir + "/main/3_LAI_asc"
dir_step5 = current_dir + "/main/4_monthly_LAI_tif"
dir_step6 = current_dir + "/main/5_monthly_LAI_asc"
dir_step7 = current_dir + "/main/6_inter_LAI_asc"
temp = current_dir + "/main/temp"

#-------------------------------------------------------
# - - - FUNCTIONS - - - FUNCTIONS - - - FUNCTIONS - - -
#-------------------------------------------------------

def get_meta(path):
    '''
    For selected file (%path) this function returns:[0] width of pixel [1] height of pixel [2] number of columns [3] number of rows
    [4] value of no data [5] EPSG code [6] transformation [7] full projection
    '''
    data_src = gdal.Open(path)
    band_info = data_src.GetRasterBand(1)
    ncols = band_info.XSize
    nrows = band_info.YSize
    nodata=band_info.GetNoDataValue()
    geoTrans = data_src.GetGeoTransform()
    # GDAL GetGeoTransform() returns pixel width [1] and height [5]
    xres = geoTrans[1]
    yres = geoTrans[5]
    wktProjection = data_src.GetProjection()
    epsg = wktProjection.split('[')
    epsg = epsg[-1]
    epsg = epsg[8:-3]
    return xres, yres, ncols, nrows, nodata, epsg, geoTrans, wktProjection

def get_extent(path):
    '''
    For selected file (%path) this function returns: [0] full extent (ulx, uly, lrx, lry) [1] ulx - x coord of upper left corner
    [2] uly - y coord of upper left corner [3] lrx - x coord of lower right corner [4] lry - y coord of lower right corner
    '''
    data_src = gdal.Open(path)
    data_geo = data_src.GetGeoTransform()
    # GDAL GetGeoTransform() returns ulx [0], uly [3] (upper left px)
    ulx = data_geo[0]
    uly = data_geo[3]
    band_info = data_src.GetRasterBand(1)
    ncols = band_info.XSize
    nrows = band_info.YSize
    lrx = ulx + data_geo[1]*ncols
    lry = uly - data_geo[1]*nrows
    del data_src
    extent = "{left} {down} {right} {up}".format(
        left = ulx,
        down = lry,
        right = lrx,
        up = uly,
    )
    return extent, ulx, uly, lrx, lry
    # GDAL_translate takes [-a_ullr ulx uly lrx lry] 

def search_folder(directory, file_format):
    '''
    Returns a list of files with a certain format (%file_format) in a folder located under a given directory (%directory)
    '''
    list_paths = os.listdir(directory)
    list_files = []
    a = len(file_format)
    for b in range(len(list_paths)):
        file_name = list_paths[b]
        if file_name[-a:] == file_format :
            list_files.append(list_paths[b])
    return list_files

def cut_temp(moveFrom, file_format, moveTo):
    '''
    Copies all files with a certain format (file_format) from the first directory (%moveFrom) to the second one 
    (%moveTo). After that, deletes the content of the first folder (%moveFrom).
    '''
    tempFiles=search_folder(moveFrom, file_format)        
    for i in range (len(tempFiles)):
        in_raster = moveFrom + "/" + tempFiles[i]
        out_raster = moveTo + "/" + tempFiles[i]
        shutil.copy2(in_raster, out_raster)          
    files = os.listdir(moveFrom)
    path = moveFrom
    for i in range(len(files)):
        os.remove(path + "/" + files[i])     

def read_tif_as_np(input_path):   
    """
    This function reads .tiff file and returns as numpy array.
    """
    raster_src = gdal.Open(input_path)
    raster_band = raster_src.GetRasterBand(1) 
    raster_nodata = raster_band.GetNoDataValue()
    raster_map = gdal.Band.ReadAsArray(raster_band)
    raster_map_unmasked = np.where(raster_map == raster_nodata, raster_nodata, np.array(raster_map))
    del raster_src
    return raster_map_unmasked

def convert_PV(input_path, sm_path, output_path, f_invalid_px):
    """
    This function reads raw Proba-V image  (%input_path), checks radiometric quality with Status Map (%SM_path), converts digital
    values to physical values (using eq.NDVI = (SM*0.004)-0.08)), corrects negative values and saves NDVI maps as %output_path.
    Fraction of clouded images (f_invalid_px) to be discarded are defined in TREX_setup.txt file.
    """
    # List of ProbaV pixels codes of SM:   
    # 0001 1111 = 248 clear, inland, all radiometric is ok
    # 0001 0111 = 232 clear, inland, no SWIR
    # 0001 1110 = 120 clear, inland, no BLUE
    # 0001 0110 = 104 clear, inland, no SWIR, no BLUE
    # 0000 1111 = 240 clear, sea, all radiometric is ok    
    # 1001 1111 = 249 inland, ice / snow
    # 0111 1111 = 254 inland, cloud
    # 0011 1111 = 252 inland, shadow
    sm_meta = get_meta(sm_path)
    sm_map = read_tif_as_np(sm_path)
    input_meta = get_meta(input_path)
    input_map = read_tif_as_np(input_path)
    # get coordinates of invalid pixels from status map
    coord_invalid_pixels = list(zip(*np.where (sm_map != 248)))
    # coord_invalid_pixels = list(zip(*np.where ((sm_map != 248) & (sm_map != 232) & (sm_map != 120) & (sm_map != 104))))

    invalid_pixels = 0
    max_pixels = input_meta[2]*input_meta[3]
    # count number of invalid pixels
    for i in coord_invalid_pixels:
        input_map[i] = input_meta[4]
        invalid_pixels += 1
    # test for number of invalid pixels
    if invalid_pixels > max_pixels*f_invalid_px:
        pass
    else:        
        # convert data
        physical_map = (input_map * 0.004) - 0.08
        physical_map_nodata = np.where(input_map == input_meta[4], input_meta[4], physical_map) 
        # correct negative data
        negative_values = list(zip(*np.where ((physical_map_nodata < 0)&(physical_map_nodata != input_meta[4]))))
        for i in negative_values:
            physical_map_nodata[i]=0
        # save map
        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(output_path, input_meta[2], input_meta[3], 1, gdal.GDT_Float32 )
        dataset.SetGeoTransform(input_meta[6])
        dataset.SetProjection(input_meta[7])
        oBand = dataset.GetRasterBand(1)
        oBand.SetNoDataValue(input_meta[4])
        oBand.WriteArray(physical_map_nodata)
        del dataset
       
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
    np.seterr(divide='ignore', invalid='ignore')    
    #LAI Formula designed by Su
    LAI_map = np.sqrt(ndvi_map * (1+ndvi_map) / (1-ndvi_map))   
    save_LAI_map_Nodata=np.where(ndvi_map == nodata, nodata, LAI_map)
    border = gdal.Open(dir_reference_map)
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
            image_input = in_raster + "/" + list_of_maps[i] + ".tif"
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
        border = gdal.Open(dir_reference_map)
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
        image_output = output_folder + "/" + str(date) + "_MonthlyLAI.tif" 

        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(image_output, xSize, ySize, 1, gdal.GDT_Float32)
        dataset.SetGeoTransform(geoTrans)
        dataset.SetProjection(wktProjection)
        oBand = dataset.GetRasterBand(1)
        oBand.SetNoDataValue(nodata)

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
    noData_pixels = list(zip(*np.where(array2d == oldValue)))
    for i in noData_pixels: array2d[i] = newValue
    return array2d
#---------------------------------------------
# - - - STEP 1 - - - STEP 1 - - - STEP 1 - - -
#---------------------------------------------        
if step1 == 1:
    print ('\nDeleting previous results...')           
    delete_from = [dir_step2, dir_step3, dir_step4, dir_step5, dir_step6, dir_step7, temp]
    for i in range(len(delete_from)):
        files = os.listdir(delete_from[i])
        path = delete_from[i]
        for i in range(len(files)):
            os.remove(path + "//" + files[i])            
else:
    pass
#---------------------------------------------
# - - - STEP 2 - - - STEP 2 - - - STEP 2 - - -
#---------------------------------------------
if step2 == 1:       
    print ('\nPreprocessing Proba-V maps...') 
    list_pv_maps = search_folder(dir_input_maps, 'NDVI.tif') 
    list_sm_maps = search_folder(dir_input_maps, 'SM.tif')
    fNames = []
    for i in range (len(list_pv_maps)):
        breaklist = list_pv_maps[i].split(".")
        fn = breaklist [0].split("_")[-4]
        fn = fn + "_NDVI.tif"
        fNames.append(fn) 
    for j in range (len(list_pv_maps)):
        os.chdir (dir_input_maps)
        convert_PV(list_pv_maps[j], list_sm_maps[j], temp + "//" + fNames[j], f_invalid_px)
    cut_temp(temp, '.tif', dir_step2)   
    #---------------------------------------------
    list_ndvi_maps=search_folder(dir_step2, '.tif')     
    for i in range (len(list_ndvi_maps)):
        os.chdir (dir_step2)
        in_raster = dir_step2 + "//" + list_ndvi_maps[i]
        out_raster = temp + '//' + list_ndvi_maps[i]
        ndvi_map_meta = get_meta(in_raster)
        new_xres = ndvi_map_meta[0]/2
        new_yres = ndvi_map_meta[1]/2
        increase_resolution = "gdalwarp  -multi -of GTiff -co TILED=YES -tr {xres} {yres} {input_file} {output_file} -overwrite -q".format(
            xres = new_xres,
            yres = new_yres,
            input_file = in_raster,
            output_file = out_raster,
            )
        os.system (increase_resolution)        
    cut_temp(temp, '.tif', dir_step2)    
    #---------------------------------------------
    print ('\nReprojecting and resampling...') 
    list_ndvi_maps2 = search_folder(dir_step2, '.tif')  
    for k in range (len(list_ndvi_maps2)):
        os.chdir (dir_step2)
        in_raster = dir_step2 + "/" + list_ndvi_maps2[k]
        out_raster = temp + '/' + list_ndvi_maps2[k]
        meta = get_meta(dir_reference_map)
        extent = get_extent(dir_reference_map)[0]
        warp_command = "gdalwarp -t_srs EPSG:{epsg} -r {r_type} -te {extent} -tr {xres} {yres} {input_file} {output_file} -overwrite -q".format(
            epsg = meta[5],
            r_type = "cubic",
            extent = extent,
            input_file = in_raster,
            output_file = out_raster,
            xres = meta[0],
            yres = meta[1],
            )
        os.system(warp_command)
    cut_temp(temp, '.tif', dir_step2)  
#---------------------------------------------
# - - - STEP 3 - - - STEP 3 - - - STEP 3 - - -
#---------------------------------------------
if step3 == 1:
    print ('\nGenerating LAI.tif maps...')    
    list_lai_tif=search_folder(dir_step2, '.tif')
    dates = []
    for i in range (len(list_lai_tif)):
        breakList = list_lai_tif[i].split("/")
        fn = breakList [-1]
        fn2 = fn.split("_")
        fn2 = fn2[0]
        dates.append(fn2)
    
    for i in range (len(list_lai_tif)):
        in_raster = dir_step2 + "/" + list_lai_tif[i]
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
    print ('\nGenerating LAI.asc maps...')
    list_lai_asc=search_folder(dir_step3, '.tif')
    os.chdir(dir_step4)
    dates = []
    for i in range (len(list_lai_asc)):
        breakList = list_lai_asc[i].split("/")
        fn = breakList [-1]
        fn2 = fn.split("_")
        fn2 = fn2[0]
        dates.append(fn2)
        
    for i in range (len(list_lai_asc)):
        in_raster = dir_step3 + "/" + list_lai_asc[i]
        out_raster = dir_step4 + "/" + dates[i] + "_LAI.asc"
        if os.path.exists(out_raster):
            os.remove(out_raster)
        convert_to_asc = "gdal_translate -of AAIGrid {input_file} {output_file} -q".format(
            input_file = in_raster,
            output_file = out_raster,
            )
        os.system (convert_to_asc)    
else:
    pass
#---------------------------------------------
# - - - STEP 5 - - - STEP 5 - - - STEP 5 - - -
#--------------------------------------------- 
if step5 == 1:
    print ('\nGenerating monthly aggregated LAI.tif maps... ')
    list_lai_asc2=search_folder(dir_step3, '.tif')  
    data_src = gdal.Open(dir_step3 + "/" + list_lai_asc2[0])
    band_info = data_src.GetRasterBand(1)  
    filenames = []
    dates = []
    for i in range (len(list_lai_asc2)):
        breakList = list_lai_asc2[i].split("/")
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
    
    print ('\nStatus report - number of available products for each month')
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
                    print (add_list)
                    report.write('\n ' + str(name_of_months[k]) + ':')
                    for a in add_list:
                        b = ' ' + str(a)
                        report.write(b)
            except: pass
    report.write('\n ===================================================== ')
    print ('\nReport saved as report.txt in the file 4_monthly_LAI_tif')
    report.close()
else:
    pass
#---------------------------------------------
# - - - STEP 6 - - - STEP 6 - - - STEP 6 - - -
#--------------------------------------------- 
if step6 == 1:
    print ('\nGenerating aggregated LAI.asc maps...')
    Myfiles6=search_folder(dir_step5, '.tif')    
    filenames = []
    dates = []
    for i in range (len(Myfiles6)):
        breakList = Myfiles6[i].split("/")
        fn = breakList [-1]
        fn = fn[:-4]
        filenames.append(fn)
        fn2 = fn.split("_")
        fn2 = fn2[0]
        dates.append(fn2)    
    for i in range (len(Myfiles6)):
        in_raster = dir_step5 + "/" + Myfiles6[i]
        out_raster = dir_step6 + "/" + dates[i] + "_MonthlyLAI.asc"       
        if os.path.exists(out_raster):
            os.remove(out_raster)
        convert_to_asc2= "gdal_translate -q -of AAIGrid {input_file} {output_file}".format(
            input_file = in_raster,
            output_file = out_raster,
            ) 
        os.system (convert_to_asc2)    
else:
    pass
#---------------------------------------------
# - - - STEP 7 - - - STEP 7 - - - STEP 7 - - -
#--------------------------------------------- 
if step7 == 1:   
    print ('\nInterpolating between monthly LAI.asc maps...')

    list_mLAI_asc=search_folder(dir_step6, '.asc')    
    filenames = []
    dates = []   
    for i in range (len(list_mLAI_asc)):
        breakList = list_mLAI_asc[i].split("/")
        fn = breakList [-1]
        fn = fn[:-4]
        filenames.append(fn)
        fn2 = fn.split("_")
        fn2 = fn2[-1]
        dates.append(fn2) 

    reference = gdal.Open(dir_reference_map)
    band_info = reference.GetRasterBand(1) 
    reference_nodata = band_info.GetNoDataValue()
    nodata_mask = gdal.Band.ReadAsArray(band_info)    
    
    for j in range (len(list_mLAI_asc)):
        os.chdir(dir_step6)
        in_raster = dir_step6 + "/" + list_mLAI_asc[j]
        get_header = create_header(in_raster)
        nodata = get_header[3]
        raster_n = readMap(list_mLAI_asc[j], get_header[1], get_header[2], get_header[3])
        raster_n[nodata_mask == reference_nodata] = np.nan
        #check if raster_n has nodata
        noData_pixels = list(zip(*np.where(raster_n == nodata)))   
        if len(noData_pixels) == 0:
            os.chdir(dir_step7)
            raster_n[nodata_mask == reference_nodata] = nodata                       
            np.savetxt('int_' + str(list_mLAI_asc[j]) + '.asc', raster_n, fmt='%10.2f', delimiter=' ', header=get_header[0], comments='')
        else:
            print ('\nInterpolating ' + str(len(noData_pixels)) + ' values of ' + str(list_mLAI_asc[j]))
            os.chdir(dir_step6)
            raster_n_minus1 = readMap(list_mLAI_asc[j-1], get_header[1], get_header[2], get_header[3])
            raster_n_minus2 = readMap(list_mLAI_asc[j-2], get_header[1], get_header[2], get_header[3])
            try: raster_n_plus1 = readMap(list_mLAI_asc[j+1], get_header[1], get_header[2], get_header[3])
            except: pass
            try: raster_n_plus2 = readMap(list_mLAI_asc[j+2], get_header[1], get_header[2], get_header[3])
            except: pass
            raster_nan = raster_n
            raster_nan[raster_nan == nodata] = np.nan
            filler = np.nanmean(raster_nan)
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
            np.savetxt('int_' + str(list_mLAI_asc[j]) + '.asc', raster_n, fmt='%10.2f', delimiter=' ', header=get_header[0], comments='')
#---------------------------------------------
print ('\n PROCESSING COMPLETE')