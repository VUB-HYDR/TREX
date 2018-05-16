#           ProbaV - LAI visualization tool
#                 (16/05/2018)
#-------------------------------------------------------
# - - - MODULES AND WORKING DIRECTORIES - - - - - - - - -
#-------------------------------------------------------

import os,os.path
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from osgeo import gdal, gdalnumeric, ogr, osr

max_value = 6

#-------------------------------------------------------
# - - - INITIATE SETUP - - - - - - - - - - - - - - - - -
#-------------------------------------------------------
current_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(current_dir)

setup_name = "LAI_ProbaV_setup.txt"
read_setup = open(setup_name,'r')

# skip first 5 lines
for i in range(5):
    read_setup.readline()
# read line #6 reference raster
reference_raster = read_setup.readline().split()
dir_input_raster = current_dir + "\\reference_maps\\" + reference_raster[-1]
# read line #7 vector map
vector_map = read_setup.readline().split()
dir_vector_map = current_dir + "\\reference_maps\\" + vector_map[-1]
# skip line #8-13
for i in range(6):
    read_setup.readline()
# read line #14 LAI.asc
step4 = read_setup.readline().split()
step4 = int(step4[-1])
# skip line #15 
read_setup.readline()
# read line #16 monthly LAI.asc
step6 = read_setup.readline().split()
step6 = int(step6[-1])
# skip line #17-19
for i in range(3):
    read_setup.readline()
# read line #20 LAI.jpg
step8 = read_setup.readline().split()
step8 = int(step8[-1])    
# read line #21 monthly LAI.jpg
step9 = read_setup.readline().split()
step9 = int(step9[-1])
# read line #22 int monthly LAI.jpg
step10 = read_setup.readline().split()
step10 = int(step10[-1]) 
# read line #23 LAI timeseries.csv
step11 = read_setup.readline().split()
step11 = int(step11[-1]) 
# read line #24 LAI timeseries.jpg
step12 = read_setup.readline().split()
step12 = int(step12[-1]) 


dir_step2 = current_dir + "\\main\\1_NDVI_tif"
dir_step3 = current_dir + "\\main\\2_LAI_tif"
dir_step4 = current_dir + "\\main\\3_LAI_asc"
dir_step5 = current_dir + "\\main\\4_monthly_LAI_tif"
dir_step6 = current_dir + "\\main\\5_monthly_LAI_asc"
dir_step7 = current_dir + "\\main\\6_inter_LAI_asc"
dir_step8 = current_dir + "\\main\\7_LAI_jpg"
dir_step9 = current_dir + "\\main\\8_monthly_LAI_jpg"
dir_step10 = current_dir + "\\main\\9_interpolated_LAI_jpg"
dir_step11 = current_dir + "\\main\\10_LAI_timeseries"
temp = current_dir + "\main\\temp"

#-------------------------------------------------------
# - - - FUNCTIONS - - - FUNCTIONS - - - FUNCTIONS - - -
#-------------------------------------------------------

def SearchFolder(directory, file_type):
#searches folder (argument#1 = directory) for a defined file type 
#(argument#2 = file_type) and returns a list of such files
    MyFolder = os.listdir(directory)
    MyList = []
    a = len(file_type)
    for b in range(len(MyFolder)):
        i = MyFolder[b]
        if i[-a:] == file_type :
            MyList.append(MyFolder[b])
    return MyList

def create_header(from_asc):
#extracts header info (number of columns & rows, coordinates of XLL & YLL
#corners, cell size [m] and value of no data) from a given file (path)
    f = open(from_asc,'r')
    header = []
    for i in range(6):
        header.append(f.readline())
    ncols = int(header[0].split()[1])
    nrows = int(header[1].split()[1])
    xll = float(header[2].split()[1])
    yll = float(header[3].split()[1])
    cellsize = float(header[4].split()[1])
    nodata = float(header[5].split()[1])
    f.close()    
    createHeader = 'NCOLS ' + str(ncols) + "\n" + 'NROWS ' + str(nrows) + "\n" + 'XLLCORNER ' + str(xll) + "\n" + 'YLLCORNER ' + str(yll) + "\n" + 'CELLSIZE ' + str(cellsize) + "\n" + 'NODATA_VALUE ' + str(nodata)
    return createHeader, ncols, nrows, nodata

def readMap(fileName, ncols, nrows, nodata):
#creates a numpy array from the certain file (fileName) with given size (nrows, ncols)
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
            for items in range(len(a)):
                a[a.index('-1.#IND')] = nodata
#            a[a.index(string, beg = 0, end = len(a))] = 255
#            print a
        except:
            pass    
    result = np.array(temp,float).reshape(nrows, ncols)
    f.close()
    return result
    
#-------------------------------------------------------
# - - - MAIN - - - MAIN - - - MAIN - - - MAIN - - -
#-------------------------------------------------------

#---------------------------------------------
# - - - STEP 8 - - - STEP 8 - - - STEP 8 - - -
#---------------------------------------------
if step8 == 1:
#    print '\n Step 8'
    print '\nGenerating LAI.jpg maps...'
    Myfiles8 = SearchFolder(dir_step4, '.asc')   
    filenames = []
    for i in range (len(Myfiles8)):
        breakList = Myfiles8[i].split("\\")
        fn = breakList [-1]
        fn = fn[:-4]
        filenames.append(fn)    
    for j in Myfiles8:
#        print 'Reading ' + j + '...'
        os.chdir(dir_step4)        
        file_with_header = dir_step4 + "\\" + j
        get_header = create_header(file_with_header)
        lai_map = readMap(j, get_header[1], get_header[2], get_header[3])        
        plt.clf()
        plotResult_1 = lai_map
        mask_noData = np.where(plotResult_1 != get_header[3], False, True)
        masked = ma.array(plotResult_1, mask = mask_noData)
        plt.clf()
        #figcize, pixel size ratio, just for visualization
        ay, ax = plt.subplots(figsize=(8, 8)); ax.set_title(j)
        #cmap color (Greens, YlGn, BuGn, summer, winter), interpolation 'none' = pixelated, aspect auto allows streching (equal changes pixels to squares)
#        ax.plot(x_coord[0] + 0.5, y_coord[0] + 0.5, 'rD', markersize=3)
        plt.imshow(masked, extent=[0, get_header[1], get_header[2], 0],
                   vmin=0, vmax=max_value, cmap="jet", interpolation="none", aspect='equal' )
        #plt.imshow(masked, vmin=0, vmax=1, cmap="YlGn", interpolation="none", aspect='auto' )
        cbar = plt.colorbar(orientation='vertical')
        cbar.set_label('LAI [-] 5-day synthesize')
        plt.xlabel('Number of columns')
        plt.ylabel('Number of rows')
        os.chdir(dir_step8)  
        plt.savefig(j[:-4] + '.jpeg' , transparent=True, dpi=200)
        plt.close()

else:
    pass

#---------------------------------------------
# - - - STEP 9 - - - STEP 9 - - - STEP 9 - - -
#---------------------------------------------
if step9 == 1:
#    print '\n Step 9'
    print '\nGenerating monthly LAI.jpg maps...'
    Myfiles9 = SearchFolder(dir_step6, '.asc')  
    filenames = []
    for i in range (len(Myfiles9)):
        breakList = Myfiles9[i].split("\\")
        fn = breakList [-1]
        fn = fn[:-4]
        filenames.append(fn)
    
    for j in Myfiles9:
#        print 'Reading ' + j + '...'
        os.chdir(dir_step6)
        
        file_with_header = dir_step6 + "\\" + j
        get_header = create_header(file_with_header)
        lai_map = readMap(j, get_header[1], get_header[2], get_header[3])
        k = j[:-4]
        plt.clf()
        plotResult_1 = lai_map
        mask_noData = np.where(plotResult_1 != get_header[3], False, True)
        masked = ma.array(plotResult_1, mask = mask_noData)
        plt.clf()
        #figcize, pixel size ratio, just for visualization
        ay, ax = plt.subplots(figsize=(5, 5)); ax.set_title(k)
        #cmap color (Greens, YlGn, BuGn, summer, winter), interpolation 'none' = pixelated, aspect auto allows streching (equal changes pixels to squares)
#        ax.plot(x_coord[0] + 0.5, y_coord[0] + 0.5, 'rD', markersize=3)
        
#        plt.imshow(masked, extent=[0, get_header[1], get_header[2], 0],
#                   vmin=0, vmax=max_value, cmap="jet", interpolation="none", aspect='equal' )
        plt.imshow(masked, vmin=0, vmax=max_value, cmap="jet", interpolation="none", aspect='equal' )
        cbar = plt.colorbar(orientation='vertical')
        cbar.set_label('LAI [-] monthly agg')
#        plt.xlabel('Number of columns')
#        plt.ylabel('Number of rows')
        os.chdir(dir_step9)  
        plt.savefig(j[:-4] + '.jpeg' , transparent=True, dpi=200)
        plt.close()

#---------------------------------------------
# - - - STEP 10 - - - STEP 10 - - - STEP 10 - 
#---------------------------------------------
if step10 == 1:
#    print '\n Step 10'
    print '\nGenerating interpolated monthly LAI.jpg maps...'
    Myfiles10 = SearchFolder(dir_step7, '.asc')  
    filenames = []
    for i in range (len(Myfiles10)):
        breakList = Myfiles10[i].split("\\")
        fn = breakList [-1]
        fn = fn[:-4]
        filenames.append(fn)
    
    for j in Myfiles10:
#        print 'Reading ' + j + '...'
        os.chdir(dir_step7)
        
        file_with_header = dir_step7 + "\\" + j
        get_header = create_header(file_with_header)
        lai_map = readMap(j, get_header[1], get_header[2], get_header[3])
        k = j[:-4]
        plt.clf()
        plotResult_1 = lai_map
        mask_noData = np.where(plotResult_1 != get_header[3], False, True)
        masked = ma.array(plotResult_1, mask = mask_noData)
        plt.clf()
        #figcize, pixel size ratio, just for visualization
        ay, ax = plt.subplots(figsize=(5, 5)); ax.set_title(k)
        #cmap color (Greens, YlGn, BuGn, summer, winter), interpolation 'none' = pixelated, aspect auto allows streching (equal changes pixels to squares)
#        ax.plot(x_coord[0] + 0.5, y_coord[0] + 0.5, 'rD', markersize=3)
        
#        plt.imshow(masked, extent=[0, get_header[1], get_header[2], 0],
#                   vmin=0, vmax=max_value, cmap="jet", interpolation="none", aspect='equal' )
        plt.imshow(masked, vmin=0, vmax=max_value, cmap="jet", interpolation="none", aspect='equal' )
        cbar = plt.colorbar(orientation='vertical')
        cbar.set_label('LAI [-] monthly agg')
#        plt.xlabel('Number of columns')
#        plt.ylabel('Number of rows')
        os.chdir(dir_step10)  
        plt.savefig(j[:-4] + '.jpeg' , transparent=True, dpi=200)
        plt.close()

else:
    pass

#---------------------------------------------
# - - - STEP 11 - - - STEP 11 - - - STEP 11 - 
#---------------------------------------------
if step11 == 1:
#    print '\n Step 11'
    print '\nSaving LAI timeseries as .csv'
    
    os.chdir(dir_step3)
    Myfiles111 = SearchFolder(dir_step3, '.tif')
    filenames111 = []
    years = []
    LAI_Values = []
    Point_x=[]
    Point_y=[]
    Cord_px=[]
    Cord_py=[]
    
    for i in range (len(Myfiles111)):
        breakList = Myfiles111 [i].split("\\")
        fn_01 = breakList [-1].split(".")
        fn = fn_01 [0].split("_")[-1]
        year = fn[:4]
        if year not in years: years.append(year)
        fn = int(fn)
        filenames11.append(fn)  

    for i in range (len(Myfiles111)):
        # reading raster info and conveeting pixel coordinates
        source_raster = gdal.Open(Myfiles111[i]) 
        geotrans = source_raster.GetGeoTransform()
        band = source_raster.GetRasterBand(1)
        # reading shapefile info 
        ds=ogr.Open(dir_vector_map)
        lyr=ds.GetLayer()
    
        points_name = ['date']
        for feature in lyr:
            a = feature.GetField('ID')
            points_name.append(a)
        lyr.ResetReading()
        header = ",".join(str(x) for x in points_name)
   
        for feat in lyr:
            geom = feat.GetGeometryRef()
            mx,my=geom.GetX(), geom.GetY()  #coord in map units
            Point_x.append(mx)
            Point_y.append(my)
            #Convert from map to pixel coordinates.
            #Only works for geotransforms with no rotation.
            px = (mx - geotrans[0]) / geotrans[1] #x pixel
            py = (my - geotrans[3]) / geotrans[5] #y pixel
        
            Cord_px.append(px)
            Cord_py.append(py)
        
            intval=band.ReadAsArray(px,py,1,1)
            LAI_Values.append(intval.item(0))     
    
    os.chdir (dir_step11)
    dates = np.array(filenames111)
    dates = np.reshape(dates, (len(Myfiles111), 1))
    LAI_matrix = np.array(LAI_Values)
    LAI_matrix = np.reshape(LAI_matrix, (len(Myfiles111), len(points_name)-1))
    matrix = np.append(dates, LAI_matrix, axis=1)
    np.savetxt('LAI_timeseries.csv', matrix , delimiter="," , fmt='%10.2f', header = header, comments='' )

    print '\nSaving monthly LAI timeseries as .csv'
    os.chdir(dir_step5)
    Myfiles112 = SearchFolder(dir_step5, '.tif')
    filenames112 = []
    years = []
    LAI_Values = []
    Point_x=[]
    Point_y=[]
    Cord_px=[]
    Cord_py=[]
        
    for i in range (len(Myfiles112)):
        breakList = Myfiles112 [i].split("\\")
        fn_01 = breakList [-1].split(".")
        fn = fn_01 [0].split("_")[-1]
        year = fn[:4]
        if year not in years: years.append(year)
        fn = int(fn)
        filenames112.append(fn)  
    
    for i in range (len(Myfiles112)):
        # reading raster info and conveeting pixel coordinates
        source_raster = gdal.Open(Myfiles112[i]) 
        geotrans = source_raster.GetGeoTransform()
        band = source_raster.GetRasterBand(1)
        # reading shapefile info 
        ds=ogr.Open(dir_vector_map)
        lyr=ds.GetLayer()
        
        points_name = ['date']
        for feature in lyr:
            a = feature.GetField('ID')
            points_name.append(a)
        lyr.ResetReading()
        header = ",".join(str(x) for x in points_name)
       
        for feat in lyr:
            geom = feat.GetGeometryRef()
            mx,my=geom.GetX(), geom.GetY()  #coord in map units
            Point_x.append(mx)
            Point_y.append(my)
            #Convert from map to pixel coordinates.
            #Only works for geotransforms with no rotation.
            px = (mx - geotrans[0]) / geotrans[1] #x pixel
            py = (my - geotrans[3]) / geotrans[5] #y pixel
            
            Cord_px.append(px)
            Cord_py.append(py)
            
            intval=band.ReadAsArray(px,py,1,1)
            LAI_Values.append(intval.item(0))     
    os.chdir (dir_step11)
    dates = np.array(filenames112)
    dates = np.reshape(dates, (len(Myfiles112), 1))
    LAI_matrix = np.array(LAI_Values)
    LAI_matrix = np.reshape(LAI_matrix, (len(Myfiles112), len(points_name)-1))
    matrix = np.append(dates, LAI_matrix, axis=1)
    np.savetxt('monthly_LAI_timeseries.csv', matrix , delimiter="," , fmt='%10.2f', header = header, comments='' )

#print '\nSaving interpolated monthly LAI timeseries as .csv'
#    os.chdir(dir_step5)
#    Myfiles113 = SearchFolder(dir_step6, '.tif')
#    filenames113 = []
#    years = []
#    LAI_Values = []
#    Point_x=[]
#    Point_y=[]
#    Cord_px=[]
#    Cord_py=[]
#        
#    for i in range (len(Myfiles113)):
#        breakList = Myfiles113 [i].split("\\")
#        fn_01 = breakList [-1].split(".")
#        fn = fn_01 [0].split("_")[-1]
#        year = fn[:4]
#        if year not in years: years.append(year)
#        fn = int(fn)
#        filenames113.append(fn)  
#    
#    for i in range (len(Myfiles113)):
#        # reading raster info and conveeting pixel coordinates
#        source_raster = gdal.Open(Myfiles113[i]) 
#        geotrans = source_raster.GetGeoTransform()
#        band = source_raster.GetRasterBand(1)
#        # reading shapefile info 
#        ds=ogr.Open(dir_vector_map)
#        lyr=ds.GetLayer()
#        
#        points_name = ['date']
#        for feature in lyr:
#            a = feature.GetField('ID')
#            points_name.append(a)
#        lyr.ResetReading()
#        header = ",".join(str(x) for x in points_name)
#       
#        for feat in lyr:
#            geom = feat.GetGeometryRef()
#            mx,my=geom.GetX(), geom.GetY()  #coord in map units
#            Point_x.append(mx)
#            Point_y.append(my)
#            #Convert from map to pixel coordinates.
#            #Only works for geotransforms with no rotation.
#            px = (mx - geotrans[0]) / geotrans[1] #x pixel
#            py = (my - geotrans[3]) / geotrans[5] #y pixel
#            
#            Cord_px.append(px)
#            Cord_py.append(py)
#            
#            intval=band.ReadAsArray(px,py,1,1)
#            LAI_Values.append(intval.item(0))     
#    os.chdir (dir_step11)
#    dates = np.array(filenames113)
#    dates = np.reshape(dates, (len(Myfiles113), 1))
#    LAI_matrix = np.array(LAI_Values)
#    LAI_matrix = np.reshape(LAI_matrix, (len(Myfiles113), len(points_name)-1))
#    matrix = np.append(dates, LAI_matrix, axis=1)
#    np.savetxt('int_monthly_LAI_timeseries.csv', matrix , delimiter="," , fmt='%10.2f', header = header, comments='' )


print 'Complete!'
