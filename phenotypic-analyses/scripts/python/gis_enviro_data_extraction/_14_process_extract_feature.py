# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 09:07:04 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15
"""
try: 
    import os, sys, re, errno
    
    directory = r'D:\Python scripts - FINAL'
    os.chdir(directory)
    
    from _00_glue_utils import raster 
    from _00_glue_utils import lookup_table
    
except ImportError as IE:
    print (IE)
    print ("These functions requires arcpy to run")  
    sys.exit(1)    
    
try: 
    import arcpy 
    
except ImportError as IE:
    print (IE)
    print ("These functions requires arcpy to run")  
    sys.exit(1)   

finally:
    
    try:
        from arcpy.sa import * 

    except ImportError as IE:
        print (IE)
        print ("These functions requires arcpy to run")  
        sys.exit(1)      
    
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
    
def extract(arcpy_env, cities, vector_file, vector_file_wildcard, raster_file, extracted_feature, outpath):
    '''
    Description: 
        For derived Landsat image products:
            extract image pixel values by buffered shapefile in projected coordinate system (pcs; e.g., UTM).\n\ 
            This function works for pre-selected Landsat images (list of list compares \n\
            city name in shapefile to landsat), and thus using a 2-D projected coordinate system (UTM).
            
        For SRTM or ASTER GDEM datasets:
            extract image pixel values by point shapefile in:
                (1) geographic coordinate system (gcs; e.g., Decimal degrees).
                **TO USE THIS, YOU MUST SPECIFY RASTER() WITH WILDCARD 'gcs.tif'
                               YOU MUST SPECIFY vector_file_wildcard AS '*gcs.shp'     
                               
                (2) projected coordinate system (pcs; e.g., UTM) - ONLY IF REFERENCED TO Landsat image 
                **TO USE THIS, YOU MUST SPECIFY RASTER() WITH WILDCARD 'pcs.tif'
                               YOU MUST SPECIFY vector_file_wildcard AS '*pcs.shp' 
                
        buffer_file must be same folder as buffered shapefiles 
    
    Args:
        arcpy_env (str): directory of shapefile(s) (e.g., buffer shapefile)
        cities (list of lists): call from function cities()
        vector_file (str): set as 'DEFAULT' for NDVI, NDSI
                           set as 'LST' for LST
                           set as 'DEM for DEM 
                           set as 'CGIAR' for CGIAR datasets (Mean Annual Aridity, Mean Annual PET, Monthly Mean PET) # NOT IMPLEMENTED; instead refer to function extract_CGIAR()
                           
        vector_file_wildcard (str): required for shapefile input. e.g., if buffer, specify '*buffer*' or '*buffer.shp'
                                                                        if PCS, specify '*pcs*' or '*pcs.shp'
                                                                        if GCS, specify '*gcs*' or '*gcs.shp'
                                                                        
        raster_file (str): directory of Landsat images. e.g., list of values as input
        extracted_feature (str): feature to be extracted (gets appended to final output filename) e.g., 'dem', 'ndvi', etc. 
        inras_root (str): in directory of Landsat images
        outpath (str):  out directory for .csv
    
    Returns:
        No returns or exchanges. 
        
    $ to be implemented: add additional cities along with respective Landsat row/path code 
    $ to be implemented: ability to extract multiple features and append to same output .csv for each city
    # to be implemented: if outpath does not exist, create new folder, else nothing
    
    $ to be implemented; if output file exists, skip, print message
    
    $ implement is UNIQUE to only return single city output result for DEM; \n\
    $ this is fixed via pandas drop duplicates (to be implemented)
    
    $ interesting when processing through console in arcmap: extent error thrown when processing in UTM; if within UTM Americas, fine, but when switching to UTM Asia, crash!
    $ fixed naming convention to match each other (dem raster city names, spaces replaced with '_' to match shapefiles)
    '''
#    import os, arcpy, arcpy.sa
    
    arcpy.env.workspace = arcpy_env
    arcpy.env.overwriteOutput = True
        
    fclist = arcpy.ListFeatureClasses(vector_file_wildcard) 
    fc = []
    count = 0
    
    # get citie names from shp
    for i in range(len(fclist)):
        
        # for buffer-point extraction (transect points)
        if vector_file == 'DEFAULT' or vector_file == 'LST':
            # all lower case for match
#            fc.append(fclist[i][:-11].lower()) # e.g., 'Acton_buffer.shp' >>> 'acton'
            fc.append(fclist[i][:-15].lower()) # e.g., 'Acton_pcs_buffer.shp' >>> 'acton'
  
        # for point extraction (city location)
        elif vector_file == 'DEM':  
            # all lower case for match
            fc.append(fclist[i][:-8].lower())# e.g., 'Acton_pcs.shp' >>> 'acton'
        
    # standardize city to landsat code list with '_' for city names with >2 words, such that it matches the shapefiles with '_'
    for i in range(len(cities)):
        cities[i][1] = cities[i][1].lower().replace(' ', '_') # all lower case for match
    
    # process... get city  
    for i in range(len(cities)):
        
        if vector_file == 'DEFAULT':

            if cities[i][1].lower() in fc: # all lower case for match
                print cities[i][1]
    
                # get raster
                for j in range(len(raster_file)): 
                    
                    # splice for code 
                    a = os.path.basename(raster_file[j])[:10] 
                    b = os.path.basename(raster_file[j])[16:] 
                    
                    raster_code = os.path.basename(raster_file[j].split(a)[-1].split(b)[0])
    #                    print raster_code
                    # compare city name to raster code 
                    if raster_code in cities[i][2]:
                      
    #                        print raster_code, cities[i][1], cities[i][2] 
                
                        INPUT_SHAPEFILE = cities[i][1] + '_pcs_buffer.shp'  # cities[i][1] + '_buffer.shp'
                        INPUT_RASTER = arcpy.Raster(raster_file[j])      
                        OUTPUT_TABLE =  outpath + '\\' + cities[i][1] + '_' + cities[i][2] + '_' + os.path.basename(raster_file[j])[17:25] + '_' + extracted_feature + '.dbf'
                        SHAPEFILE_LAYER_FOR_PROCESSING = cities[i][1] + '_buffer_layer'
                        OUTPUT_CSV = cities[i][1] + '_' + cities[i][2] + '_' + os.path.basename(raster_file[j])[17:25] + '_' + extracted_feature + '.csv'
                        
                        print OUTPUT_CSV
                        print '.....'
                        
                        try:
                            if os.path.isfile(os.path.join(outpath, cities[i][1] + '_' + cities[i][2] + '_' + os.path.basename(raster_file[j])[17:25] + '_' + extracted_feature + '.csv')):
                                print ('...')
                                print ('%s already exists' % cities[i][1] + '_' + cities[i][2] + '_' + extracted_feature + '.csv')
                                print ('...')
                            
                            else:
                                
                                # 'FID' as unique identifier for each row 
                                arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE, 'FID', INPUT_RASTER, OUTPUT_TABLE, "NODATA", "MEAN")
                                
                                # make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
                                arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
                                
                                # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
                                arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_COMMON')
                                
                                # save out to csv 
                                arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
                                                        
                                arcpy.Delete_management (SHAPEFILE_LAYER_FOR_PROCESSING)
                                
                                count += 1
                            
                        except:
                            print ('{0} could not be processed'.format(OUTPUT_CSV))
                            pass

        elif vector_file == 'LST':
    
            if cities[i][1].lower() in fc: # all lower case for match
    
                # get raster
                for j in range(len(raster_file)): 
                    
                    # splice for code 
                    a = os.path.basename(raster_file[j])[:10] 
                    b = os.path.basename(raster_file[j])[16:] 
                    # splice for city name 
                    c = os.path.basename(raster_file[j])[:41] 
                    d = os.path.basename(raster_file[j])[-16:] 
                    
                    raster_code = os.path.basename(raster_file[j].split(a)[-1].split(b)[0])
                    
                    raster_lst_city_name = os.path.basename(raster_file[j].split(c)[-1].split(d)[0]).lower()

                    # compare city name to raster code 
                    if raster_code in cities[i][2] and raster_lst_city_name in cities[i][1]:
                
                        INPUT_SHAPEFILE = cities[i][1] + '_pcs_buffer.shp'  # cities[i][1] + '_buffer.shp'
                        INPUT_RASTER = arcpy.Raster(raster_file[j])      
                        OUTPUT_TABLE =  outpath + '\\' + cities[i][1] + '_' + cities[i][2] + '_' + os.path.basename(raster_file[j])[17:25] + '_' + extracted_feature + '.dbf'
                        SHAPEFILE_LAYER_FOR_PROCESSING = cities[i][1] + '_buffer_layer'
                        OUTPUT_CSV = cities[i][1] + '_' + cities[i][2] + '_' + os.path.basename(raster_file[j])[17:25] + '_' + extracted_feature + '.csv'
                        
                        print OUTPUT_CSV
                        print '.....'
                        
                        try:
                            if os.path.isfile(os.path.join(outpath, cities[i][1] + '_' + cities[i][2] + '_' + os.path.basename(raster_file[j])[17:25] + '_' + extracted_feature + '.csv')):
                                print ('...')
                                print ('%s already exists' % cities[i][1] + '_' + cities[i][2] + '_' + extracted_feature + '.csv')
                                print ('...')
                            
                            else:
                                
                                # 'FID' as unique identifier for each row 
                                arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE, 'FID', INPUT_RASTER, OUTPUT_TABLE, "NODATA", "MEAN")
                                
                                # make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
                                arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
                                
                                # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
                                arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_COMMON')
                                
                                # save out to csv 
                                arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
                                                        
                                arcpy.Delete_management (SHAPEFILE_LAYER_FOR_PROCESSING)
                                
                                count += 1
                            
                        except:
                            print ('{0} could not be processed'.format(OUTPUT_CSV))
                            pass
                        
        # handle GCS and PCS
        elif vector_file == 'DEM':
            
            if os.path.isfile(os.path.join(outpath, cities[i][1] + '_' + cities[i][2] + '_' + extracted_feature + '.csv')):
                print ('...')
                print ('%s already exists' % (cities[i][1] + '_' + cities[i][2] + '_' + extracted_feature + '.csv'))
                print ('...')
            else:
    #            try:
                # satisfy logic to parse DEM
                # if city in list matches city shapefile name...
    #            for i in range(len(cities)):
                if cities[i][1] in fc: # all lower case for match
                    print cities[i][1]
                    
                    for j in range(len(raster_file)): 
                
                        # retrieves city name (e.g., toronto from 'toronto_dem_pcs.tif')
                        raster = os.path.basename(raster_file[j])[:-12]
                        
                        if raster in cities[i][1].lower():
                            
                            print raster, cities[i][1], cities[i][2] 
        
                            # satisfy for GCS and PCS
                            if 'gcs' in vector_file_wildcard:
                                extension = '_gcs.shp'
                                
                            elif 'pcs' in vector_file_wildcard:  
                                extension = '_pcs.shp'
                                
                            INPUT_SHAPEFILE = cities[i][1] + extension
                            INPUT_RASTER = arcpy.Raster(raster_file[j])      
                            OUTPUT_TABLE =  os.path.join(outpath, cities[i][1] + '_' + cities[i][2]  + '_' + extracted_feature + '.dbf')
                            SHAPEFILE_LAYER_FOR_PROCESSING = cities[i][1].lower() + '_cs_layer'  # cs = coord sys 
                            OUTPUT_CSV = cities[i][1] + '_' + cities[i][2] + '_' + extracted_feature + '.csv'
        
    #                            print OUTPUT_CSV
    #                            print '.....'
                            
                            try:
                                arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE, 'FID', INPUT_RASTER, OUTPUT_TABLE, "NODATA", "MEAN")
                                    
                                    #work logic later; above works
            #                        arcpy.sa.ExtractValuesToPoints(cities[i][0] + '_prj.shp', arcpy.Raster(inras_root + '\\' + cities[i][0].lower() + '\\' +  raster_file[j]), cities[i][0] + '_prj.shp')
            
                                    #make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)       
                                arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
                                
                                # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
                                arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_COMMON')
                                
                                #save out to csv 
                                arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
                                
                                # save elevation out as new shapefile
        #                        arcpy.sa.ExtractValuesToPoints (in_point_features, in_raster, out_point_features, {interpolate_values}, {add_attributes})
                                
                                # save elevation out only as table; issue with Sample method is that it repeats analysis for each row, not each point (e.g., 1 point has 20 rows)
        #                        arcpy.sa.Sample(arcpy.Raster(inras_root + '\\' + cities[i][0].lower() + '\\' +  raster_file[j]), cities[i][0] + '_prj.shp', cities[i][0] + '_' + raster_code + '_' + extracted_feature + '.dbf')                        
        #                        arcpy.TableToTable_conversion(cities[i][0] + '_' + raster_code + '.dbf', outpath, cities[i][0] +  '_' + raster_code + '_' + extracted_feature + '.csv')
                                
                                arcpy.Delete_management (SHAPEFILE_LAYER_FOR_PROCESSING)
                                
                                count += 1
                                    
                            except:
                                print ('{0} could not be processed'.format(OUTPUT_CSV))
                                pass
 
    

def extract_LST_from_weather_stations(arcpy_env, shapefile, in_list, raster_file, extracted_feature, inras_root, outpath):
    '''
    Description: extract for single shapefile for multiple rasters of the same spatial \n\
                 reference (e.g., UTM Zone 17N) in a directory
    '''
#    import arcpy
#    from arcpy.sa import *
    
    arcpy.env.workspace = arcpy_env
    arcpy.env.overwriteOutput = True
    
    fc = shapefile

    cities = in_list 
    for i in range(len(raster_file)): 
    
        arcpy.sa.ZonalStatisticsAsTable(fc, 'FID', arcpy.Raster(inras_root + '\\' + raster_file[i]), 
                                        outpath + '\\' +  cities[0] + '_' + cities[1] + '_' + raster_file[i][17:25] + '_' + extracted_feature + '.dbf',  
                                        "NODATA", "MEAN")
        
        
        #make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
        arcpy.MakeFeatureLayer_management(fc, 'fc') 
        
        # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
        arcpy.AddJoin_management('fc', 'FID', outpath + '\\' +  cities[0] + '_' + cities[1] + '_' + raster_file[i][17:25] + '_' + extracted_feature + '.dbf', 'FID_', 
                                 'KEEP_COMMON')
        
        #save out to csv 
        arcpy.TableToTable_conversion('fc', outpath, cities[0] + '_' + cities[1] + '_' + raster_file[i][17:25] + '_' + extracted_feature + '.csv')
        
        
        arcpy.Delete_management ('fc')
        
    
def extract_WorldClim(arcpy_env, vector_file, raster_file, extracted_feature, outpath):
    '''
    Description: 
        extract World Clim pixel values by city point location shapefile. 
        
        **extraction based on geographic coordinate system (e.g., WGS 1984) 
        
    Args:
        arcpy_env (str): directory of shapefile(s) (e.g., buffer shapefile)
        vector_file (str): specify shapefile to be used as overlay for feature extraction. 
                           e.g., r'D:\GLUE Results\Process_City_Location\select_cities.shp'
        vector_file_wildcard (str): required for shapefile input. 
                                    e.g., if buffer, specify '*buffer.shp'
                                          if prj, specify '*prj.shp'
        raster_file (str): directionry of Landsat images
        extracted_feature (str): feature to be extracted 
        inras_root (str): in directory of Landsat images
        outpath (str):  out directory for .csv
    
    Returns:
        No returns or exchanges. 
        
    $ to be implemented: add additional cities along with respective Landsat row/path code 
    $ to be implemented: ability to extract multiple features and append to same output .csv for each city
    # to be implemented: if outpath does not exist, create new folder, else nothing
    
    '''
#    import arcpy
#    from arcpy.sa import *

    arcpy.env.workspace = arcpy_env
    arcpy.env.overwriteOutput = True
    
    # grab city shapefile 
    fc = vector_file
    
    worldclim = raster_file
    image_numeric_month = []
    
    count = 1
    
    # match image to month
    for i in range(len(worldclim)):
        
        split_img_text = worldclim[i].split('.')  # e.g., ['D:\\GLUE Datasets\\WorldClim_V2\\wc2', '0_30s_tavg\\wc2', '0_30s_tavg_12', 'tif']
        grab_month = split_img_text[2][11:] 
        
        image_numeric_month.append(grab_month) 


    match_month_list = [['01','Jan'],
                        ['02','Feb'],
                        ['03','Mar'],
                        ['04','Apr'],
                        ['05','May'],
                        ['06','Jun'],
                        ['07','Jul'],
                        ['08','Aug'],
                        ['09','Sep'],
                        ['10','Oct'],
                        ['11','Nov'],
                        ['12','Dec']]
    
#    for i in range(len(match_month_list)):
#        if match_month_list[i][0] in image_numeric_month[i]:
#            print match_month_list[i][1], image_numeric_month[i]
            
    
    # extract average temperature for each month for each city 
    for i in range(len(raster_file)):
        
        if match_month_list[i][0] in image_numeric_month[i]:

            INPUT_SHAPEFILE = fc
            INPUT_RASTER = arcpy.Raster(raster_file[i])      
            OUTPUT_TABLE =  outpath + '\\' +  'city' + '_' + extracted_feature + '_' + match_month_list[i][1] + '.dbf'
            SHAPEFILE_LAYER_FOR_PROCESSING = fc[:-4] + '_' + str(match_month_list[i][1])
            OUTPUT_CSV = str(count) + '_city' + '_' + extracted_feature + '_' + match_month_list[i][1] + '.csv'          

            arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE,'FID' , INPUT_RASTER,OUTPUT_TABLE, "NODATA", "MEAN")
            
            # make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
            arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
            
            # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
            arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_ALL')
            
            # save out to csv 
            arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
            
            
            arcpy.Delete_management (SHAPEFILE_LAYER_FOR_PROCESSING)
            
            count += 1
          
        
def extract_CGIAR(vector_file_directory, root_directory, outpath, *args):
    '''
    Description:
        Extracts pixel-level values from CGIAR raster datasets incl. \n\
        Annual Aridity Index, Annual Potential Evapotranspiration, \n\
        Monthly Potential Evapotranspiration. Since these rasters are projected \n\
        in a geographic coordinate system using EPSG 4326 (WGS 1984), only \n\
        shapefiles (vector data) in the same projection \n\
        can be used. 
    
    **WARNING: arcpy CANNOT HANDLE CGIAR .adf raster format (Esri Grid format); \n\
               TOO MANY UNIQUE VALUES CAUSES arcpy.sa.ZonalStatisticsAsTable \n\
               to crash! The support help document DOES NOT TELL YOU THIS! \n\
               YOU MUST CONVERT THESE .adf rasters to .tif format for processing! 
               
               https://support.esri.com/en/technical-article/000012343 
               
    **WARNING: logic is only setup to process a directory of shapefiles; if you \n\
               want to process a single shapefile, then create a folder with only a \n\
               single shapefile.
         
    Dependencies:
        arcpy
        arcpy.sa 
        
    Args:
        vector_file_directory (str): specify directory of point shapefiles. e.g., r'D:\GLUE Datasets'
        vector_file: point shapefile location of cities for feature extraction  **DEPRECATED AND NOT IN USE
        root_directory: specify root (parent) directory of all CGIAR rasters to be processed 
                        (logic only works for PET_HE_MONTHLY, PET_HE_YR, AI_YR)
        *args (str):'PET_HE_MONTHLY'
                    'PET_HE_YR' 
                    'AI_YR'
    Returns:
        
    use city centre locations from Samreen's work for extraction... 
    
    $ to be implemented: modify to accept multiple shapefiles (e.g., each city)
    
    '''
#    import os 
#    import arcpy
#    from arcpy.sa import *

    arcpy.env.overwriteOutput = True
    
    to_be_processed = []
    
    for arg in args: 
        if arg == 'AI_YR':
            to_be_processed.append(arg)
        elif arg == 'PET_HE_MONTHLY':
            to_be_processed.append(arg)
        elif arg == 'PET_HE_YR':
            to_be_processed.append(arg) 


    arcpy.env.workspace = vector_file_directory
    # grab city shapefile 
    #    fc = vector_file
        
    in_fclist = arcpy.ListFeatureClasses()     
    out_fclist = []
    
    for fc in in_fclist:
        out_fclist.append(os.path.join(vector_file_directory,(str(fc))))
        
    count = 1
    
    PET_HE_MONTHLY_directory = []
    PET_HE_YR_directory = []
    AI_directory = []
        
    for i in to_be_processed:
        if i == 'PET_HE_MONTHLY':

             
            match_month_list = [['01','Jan'],
                                ['02','Feb'],
                                ['03','Mar'],
                                ['04','Apr'],
                                ['05','May'],
                                ['06','Jun'],
                                ['07','Jul'],
                                ['08','Aug'],
                                ['09','Sep'],
                                ['10','Oct'],
                                ['11','Nov'],
                                ['12','Dec']]
            
            
            for root, dirs, files in os.walk(root_directory):
                if 'PET_he_monthly_tif' in root:
                    PET_HE_MONTHLY_directory.append(root)
                            
            # use ele at index 0 to get root folder for recursive parsing
            for root, dirs, files in os.walk(PET_HE_MONTHLY_directory[0]): 
                for i in range(len(files)):
                    if files[i].endswith('.tif'):                
                        arcpy.env.workspace = root 
                        
                        # get numeric val corresponding to str month 
                        split_PE = files[i][:-4].split('_')
                
                        grab_month_split_PE = split_PE[-1].split('.')
            #            print grab_month_split_PE[0]
                        
                        for j in range(len(match_month_list)):
                            
                            month_num = (int(match_month_list[j][0]))
                            
                            # match ras numeric month to ras numeric month in match_month_list 
                            if int(grab_month_split_PE[0]) == month_num:    
                                
#                                print (grab_month_split_PE[0], match_month_list[j][1], arcpy.Raster(files[i]))
                                
                            
                                for out_fc in out_fclist:
                                    try:
                                        if os.path.isfile(os.path.join(outpath, str(count) + '_' + os.path.basename(out_fc)[:-4] +'_PET_' + match_month_list[j][1] + '.csv')):
                                            print ('...')
                                            print ('%s already processed.' % (str(count) + '_' + os.path.basename(out_fc)[:-4] +'_PET_' + match_month_list[j][1] + '.csv'))
                                            print ('...') 
                                            
                                        else:                                            
                                            print (os.path.basename(out_fc)[:-4], grab_month_split_PE[0], match_month_list[j][1], arcpy.Raster(files[i]))
                                                
                                            INPUT_SHAPEFILE = out_fc
                                            INPUT_RASTER = arcpy.Raster(files[i])      
                                            OUTPUT_TABLE =  os.path.join(outpath, str(count) + '_' + os.path.basename(out_fc)[:-4] + '_PET_' + match_month_list[j][1] + '.dbf')
                                            SHAPEFILE_LAYER_FOR_PROCESSING = 'temp_' + os.path.basename(out_fc)[:-4] +'_PET_' + str(match_month_list[j][1])
                                            OUTPUT_CSV = str(count) + '_' + os.path.basename(out_fc)[:-4] +'_PET_' + match_month_list[j][1] + '.csv'          
                                
                                            arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE,'FID',INPUT_RASTER, OUTPUT_TABLE, "NODATA", "MEAN")
                                            
                                            # make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
                                            arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
                                            
                                            # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
                                            arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_ALL')
                                            
                                            # save out to csv 
                                            arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
                                            
                                            # dereference from memory; can cause error: 999999 due to lack of memory 
                                            arcpy.Delete_management (SHAPEFILE_LAYER_FOR_PROCESSING) 
                                    
                                    except:
                                        print ('{0} could not be processed; check for data integrity...'.format(os.path.basename(out_fc)))
                                        pass
                                    
                                count += 1
                                    
        elif i == 'PET_HE_YR': 
            # ok
            for root, dirs, files in os.walk(root_directory):
                if 'PET_he_annual_tif' in root:
                    PET_HE_YR_directory.append(root)            

            arcpy.env.workspace = PET_HE_YR_directory[0]   
            RasList = arcpy.ListRasters()
            Ras_to_be_processed = []
            
            for Ras in RasList:
                Ras_to_be_processed.append(Ras) 
            
            for out_fc in out_fclist:
                try:
                    if os.path.isfile(os.path.join(outpath, os.path.basename(out_fc)[:-4]   + '_PET_annual' + '.csv')):
                        print ('...')
                        print ('%s already processed.' % (os.path.basename(out_fc)[:-4]   + '_PET_annual' + '.csv'))
                        print ('...') 
                        
                    else:
                        INPUT_SHAPEFILE = out_fc
                        INPUT_RASTER = arcpy.Raster(Ras_to_be_processed[0])      
                        OUTPUT_TABLE =  os.path.join(outpath, os.path.basename(out_fc)[:-4] + '_PET_annual' + '.dbf')
                        SHAPEFILE_LAYER_FOR_PROCESSING = 'temp_' + os.path.basename(out_fc)[:-4]  + '_PET_annual'
                        OUTPUT_CSV = os.path.basename(out_fc)[:-4]   + '_PET_annual' + '.csv'          
            
                        arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE,'FID' ,INPUT_RASTER, OUTPUT_TABLE, "NODATA", "MEAN")
                        
                        # make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
                        arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
                        
                        # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
                        arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_ALL')
                        
                        # save out to csv 
                        arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
                        
                        # dereference from memory; can cause error: 999999 due to lack of memory
                        arcpy.Delete_management (SHAPEFILE_LAYER_FOR_PROCESSING)
                    
                except:
                    print ('{0} could not be processed; check for data integrity...'.format(os.path.basename(out_fc)))
                    pass                    
        
        elif i == 'AI_YR':       
            # ok 
            for root, dirs, files in os.walk(root_directory):
                if 'AI_annual_tif' in root:
                    AI_directory.append(root)
               
            arcpy.env.workspace = AI_directory[0]   
            RasList = arcpy.ListRasters()
            Ras_to_be_processed = [] 
            
            RasList = arcpy.ListRasters()
            
            
            for Ras in RasList:
                Ras_to_be_processed.append(Ras) # return [u'ai_yr.tif', u'ai_yr_rescaled.tif']; therefore Ras_to_be_processed[1]
                
            for out_fc in out_fclist:
                try:
                    if os.path.isfile(os.path.join(outpath, os.path.basename(out_fc)[:-4]  + '_AI_annual' + '.csv')):
                        print ('...')
                        print ('%s already processed.' % (os.path.basename(out_fc)[:-4]  + '_AI_annual' + '.csv'))
                        print ('...') 
                            
                    else:
                        INPUT_SHAPEFILE = out_fc
                        INPUT_RASTER = arcpy.Raster(Ras_to_be_processed[1])      
                        OUTPUT_TABLE =  os.path.join(outpath, os.path.basename(out_fc)[:-4] + '_AI_annual' + '.dbf')
                        SHAPEFILE_LAYER_FOR_PROCESSING = 'temp_' + os.path.basename(out_fc)[:-4]  + '_AI_annual'
                        OUTPUT_CSV = os.path.basename(out_fc)[:-4]  + '_AI_annual' + '.csv'          
                        
                        arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE,'FID', INPUT_RASTER, OUTPUT_TABLE, "NODATA", "MEAN")
                        
                        # make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
                        arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
                        
                        # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
                        arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_ALL')
                        
                        # save out to csv 
                        arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
                        
                        # dereference from memory; can cause error: 999999 due to lack of memory
                        arcpy.Delete_management (SHAPEFILE_LAYER_FOR_PROCESSING)
                
                except:
                    print ('{0} could not be processed; check for data integrity...'.format(os.path.basename(out_fc)))
                    pass  
                
            
def extract_single_use_CGIAR(vector_file, root_directory, outpath, *args):
    '''
    Description:
        Extracts pixel-level values from CGIAR raster datasets incl. \n\
        Annual Aridity Index, Annual Potential Evapotranspiration, \n\
        Monthly Potential Evapotranspiration. Since these rasters are projected \n\
        in a geographic coordinate system using EPSG 4326 (WGS 1984), only \n\
        shapefiles (vector data) in the same projection \n\
        can be used. 
    
    **WARNING: arcpy CANNOT HANDLE CGIAR .adf raster format (Esri Grid format); \n\
               TOO MANY UNIQUE VALUES CAUSES arcpy.sa.ZonalStatisticsAsTable \n\
               to crash! The support help document DOES NOT TELL YOU THIS! \n\
               YOU MUST CONVERT THESE .adf rasters to .tif format for processing! 
               
               https://support.esri.com/en/technical-article/000012343 
               
    **WARNING: only processes a single shapefile
          
    Dependencies:
        arcpy
        arcpy.sa 
        
    Args:
        vector_file: point shapefile location of cities for feature extraction
        root_directory: specify root (parent) directory of all CGIAR rasters to be processed 
                        (logic only works for PET_HE_MONTHLY, PET_HE_YR, AI_YR)
        *args (str):'PET_HE_MONTHLY'
                    'PET_HE_YR' 
                    'AI_YR'
    Returns:
        
    use city centre locations from Samreen's work for extraction... 
    
    $ to be implemented: modify to accept multiple shapefiles (e.g., each city)
    '''
#    import os 
#    import arcpy
#    from arcpy.sa import *

    arcpy.env.overwriteOutput = True
    
    to_be_processed = []
    
    for arg in args: 
        if arg == 'AI_YR':
            to_be_processed.append(arg)
        elif arg == 'PET_HE_MONTHLY':
            to_be_processed.append(arg)
        elif arg == 'PET_HE_YR':
            to_be_processed.append(arg) 

    # grab city shapefile 
    fc = vector_file
    
    PET_HE_MONTHLY_directory = []
    PET_HE_YR_directory = []
    AI_directory = []
    
    for i in to_be_processed:
        if i == 'PET_HE_MONTHLY':
     
            count = 1
            
            match_month_list = [['01','Jan'],
                                ['02','Feb'],
                                ['03','Mar'],
                                ['04','Apr'],
                                ['05','May'],
                                ['06','Jun'],
                                ['07','Jul'],
                                ['08','Aug'],
                                ['09','Sep'],
                                ['10','Oct'],
                                ['11','Nov'],
                                ['12','Dec']]
            
            
            for root, dirs, files in os.walk(root_directory):
                if 'PET_he_monthly_tif' in root:
                    PET_HE_MONTHLY_directory.append(root)
                              
            PET_dir = []
            # use ele at index 0 to get root folder for recursive parsing
            for root, dirs, files in os.walk(PET_HE_MONTHLY_directory[0]): 
                for i in range(len(files)):
                    if files[i].endswith('.tif'):                
                        arcpy.env.workspace = root 
                        
                        # get numeric val corresponding to str month 
                        split_PE = files[i][:-4].split('_')
                
                        grab_month_split_PE = split_PE[-1].split('.')
            #            print grab_month_split_PE[0]
                        
                        for j in range(len(match_month_list)):
                            
                            month_num = (int(match_month_list[j][0]))
                            
                            # match ras numeric month to ras numeric month in match_month_list 
                            if int(grab_month_split_PE[0]) == month_num:    
                                
                                print (grab_month_split_PE[0], match_month_list[j][1], arcpy.Raster(files[i]))
                                        
                                INPUT_SHAPEFILE = fc
                                INPUT_RASTER = arcpy.Raster(files[i])      
                                OUTPUT_TABLE =  os.path.join(outpath, 'city_PET_' + match_month_list[j][1] + '.dbf')
                                SHAPEFILE_LAYER_FOR_PROCESSING = 'temp_city_PET_' + str(match_month_list[j][1])
                                OUTPUT_CSV = str(count) + '_city_PET_' + match_month_list[j][1] + '.csv'          
                    
                                arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE,'FID' ,INPUT_RASTER, OUTPUT_TABLE, "NODATA", "MEAN")
                                
                                # make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
                                arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
                                
                                # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
                                arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_ALL')
                                
                                # save out to csv 
                                arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
                                
                                # dereference from memory; can cause error: 999999 due to lack of memory 
                                arcpy.Delete_management (SHAPEFILE_LAYER_FOR_PROCESSING)
                                
                                count += 1              
 
        elif i == 'PET_HE_YR': 
          
            # ok
            for root, dirs, files in os.walk(root_directory):
                if 'PET_he_annual_tif' in root:
                    PET_HE_YR_directory.append(root)            

            arcpy.env.workspace = PET_HE_YR_directory[0]   
            RasList = arcpy.ListRasters()
            Ras_to_be_processed = []
            
            for Ras in RasList:
                Ras_to_be_processed.append(Ras) 
            
            
            INPUT_SHAPEFILE = fc
            INPUT_RASTER = arcpy.Raster(Ras_to_be_processed[0])      
            OUTPUT_TABLE =  os.path.join(outpath, 'city_PET_annual' + '.dbf')
            SHAPEFILE_LAYER_FOR_PROCESSING = 'temp_city_PET_annual'
            OUTPUT_CSV = '_city_PET_annual' + '.csv'          

            arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE,'FID' ,INPUT_RASTER, OUTPUT_TABLE, "NODATA", "MEAN")
            
            # make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
            arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
            
            # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
            arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_ALL')
            
            # save out to csv 
            arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
            
            # dereference from memory; can cause error: 999999 due to lack of memory
            arcpy.Delete_management (SHAPEFILE_LAYER_FOR_PROCESSING)
                    
        
        elif i == 'AI_YR':       
            # ok 
            for root, dirs, files in os.walk(root_directory):
                if 'AI_annual_tif' in root:
                    AI_directory.append(root)
               
            arcpy.env.workspace = AI_directory[0]   
            RasList = arcpy.ListRasters()
            Ras_to_be_processed = [] 
            
            RasList = arcpy.ListRasters()
            
            
            for Ras in RasList:
                Ras_to_be_processed.append(Ras) # return [u'ai_yr.tif', u'ai_yr_rescaled.tif']; therefore Ras_to_be_processed[1]
                
            INPUT_SHAPEFILE = fc
            INPUT_RASTER = arcpy.Raster(Ras_to_be_processed[1])      
            OUTPUT_TABLE =  os.path.join(outpath, 'city_AI_annual' + '.dbf')
            SHAPEFILE_LAYER_FOR_PROCESSING = 'temp_city_AI_annual'
            OUTPUT_CSV = '_city_AI_annual' + '.csv'          

            arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE,'FID', INPUT_RASTER, OUTPUT_TABLE, "NODATA", "MEAN")
            
            # make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
            arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
            
            # join table for feature extraction; 'FID' to 'FID_' match; 'FID_' is generated from .dbf output 
            arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_ALL')
            
            # save out to csv 
            arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
            
            # dereference from memory; can cause error: 999999 due to lack of memory
            arcpy.Delete_management (SHAPEFILE_LAYER_FOR_PROCESSING)


def process_impervious_surface(arcpy_env, GMIS_rasters, outpath):
    '''
    Description:
        extract impervious surface estimates from vector data overlay using Global Man-made Impervious Surface (GMIS) \\
        datasets as provided by the Socioeconomic Data and Applications Center (SEDAC), \\
        A Data Center in NASA's Earth Observing System Data and Information System ( EOSDIS ).
        
    Args:
        arcpy_env (str): specify location of shapefiles in projected coordinate system (e.g. UTM coord sys) \\
                         for feature extraction analysis (e.g., buffered transect locations)
        GMIS_rasters (list): input from raster()
        outpath (str): specify save location 
        
    Returns:
        
    $ implement: if csv output exists, pass and move to next city for extraction... 
    
    '''
#    import os
#    import arcpy
#    from arcpy.sa import *
    
    arcpy.env.overwriteOutput = True
    
#    directory = r'D:\GLUE Datasets\SEDAC\Impervious Surface\Canada\CAN_gmis_impervious_surface_percentage_utm_30m'
    
    arcpy.env.workspace = arcpy_env
#    arcpy.env.workspace = r'D:\testing_for_prj\shp'
    fcList = arcpy.ListFeatureClasses('*buffer*')
    
    for fc in fcList:
        
        # grab spatial ref and extent info of transect points 
        fc_describe = arcpy.Describe(fc)
        fc_spatial_ref = fc_describe.spatialReference.name
        fc_extent = fc_describe.extent
                        
        for GMIS_ras in range(len(GMIS_rasters)):
            
            GMIS_ras_process = arcpy.Raster(GMIS_rasters[GMIS_ras])
            
            # grab spatial ref and extent info of GMIS raster 
            GMIS_ras_describe = arcpy.Describe(GMIS_ras_process)
            GMIS_ras_spatial_ref = GMIS_ras_describe.spatialReference.name
            GMIS_ras_extent = GMIS_ras_describe.extent
            
            # condition 1: check for same prj coord sys
            if fc_spatial_ref == GMIS_ras_spatial_ref:
                
                # if already processed, skip to next...
                if os.path.isfile(os.path.join(outpath, os.path.basename(fc)[:-10] + '_GMIS' + '.csv')):
                        print ('...')
                        print ('{0} already processed.'.format(os.path.basename(fc)[:-10] + '_GMIS' + '.csv'))
                        print ('...')
                        
                else:
                    # condition 2: check for overlap of transect point vs. raster
                    # https://gis.stackexchange.com/questions/148645/select-rasters-from-a-folder-that-intersect-the-tiles-of-a-polygon-shapefile 
                    if not fc_extent.disjoint(GMIS_ras_extent):
                        print '({0}, {1}, {2})'.format(os.path.basename(str(GMIS_ras_process)), GMIS_ras_spatial_ref, fc)
                        
                        # add conditionals to filter for NoData of SEDAC impervious surface
                        # recalculate surface for null 
    #                    arcpy.MakeRasterLayer_management(imperv_ras,'imperv_ras_out') 
                        
                        #----------------------------------------------------------#
                        ###
                        ##
                        # ARCPY KEEPS CRASHING at 99999; TRY MAKING TEMPORARY LAYERS FOR PROCESSING?? THIS WOULD ADD MORE JUNK INTO MEMORY
                        ##
                        ###
                        
                        arcpy.MakeRasterLayer_management(GMIS_ras_process,'GMIS_ras_process')
                        
                        # set val 255 as NoData 
                        GMIS_ras_NoData = arcpy.sa.SetNull('GMIS_ras_process', 'GMIS_ras_process','VALUE = 255')   
                        
                        # sa.con is not liking a temp raster layer... 
                        arcpy.MakeRasterLayer_management(GMIS_ras_NoData,'GMIS_ras_NoData')
                        
                        # in addition to val 255 == NoData, set val 200 == 0% impervious according to readme file   
                        # **WARNING** IF USING TEMP RASTER LAYER, NEED TO USE arcpy.RASTER("temp_layer")
                        GMIS_ras_cleaned = arcpy.sa.Con(arcpy.Raster('GMIS_ras_NoData') == 200, 0, GMIS_ras_NoData) 
                        
                        arcpy.MakeRasterLayer_management(GMIS_ras_cleaned,'GMIS_ras_cleaned')
                        
                        ##arcpy.CopyRaster_management()
                        
                        #----------------------------------------------------------#
                        
                        # set val 255 as NoData 
    #                    GMIS_ras_NoData = arcpy.sa.SetNull(GMIS_ras_process, GMIS_ras_process,'VALUE = 255')   
    #                    
    #                    # in addition to val 255 == NoData, set val 200 == 0% impervious according to readme file                    
    #                    GMIS_ras_cleaned = arcpy.sa.Con(GMIS_ras_NoData == 200, 0, GMIS_ras_NoData)
                        
                    
                        # this is not necessary; estimates from GMIS are OK from selecting random buffered transect points
                        # testing impervious surface as a binary function; 1 == impervious, 0 == pervious 
            #                    imperv_ras_final = arcpy.sa.Con((imperv_ras_cleaned > 0) & (imperv_ras_cleaned <= 100), 1, imperv_ras_cleaned)
                        
#                        if os.path.isfile(os.path.join(outpath, os.path.basename(fc)[:-10] + '_GMIS' + '.csv')):
#                            print ('...')
#                            print ('{0} already processed.'.format(os.path.basename(fc)[:-10] + '_GMIS' + '.csv'))
#                            print ('...')
#                            
                        INPUT_SHAPEFILE = fc
    #                    INPUT_RASTER = GMIS_ras_cleaned
                        INPUT_RASTER = 'GMIS_ras_cleaned'
                        OUTPUT_TABLE =  outpath + '\\' + os.path.basename(fc)[:-10] + '_GMIS' + '.dbf'
                        SHAPEFILE_LAYER_FOR_PROCESSING = os.path.basename(fc)[:-10] + '_processing'
                        OUTPUT_CSV = os.path.basename(fc)[:-10] + '_GMIS' + '.csv'
                        
                        print OUTPUT_CSV
                        print '.....'
                        # disable ignore_nodata
                        arcpy.sa.ZonalStatisticsAsTable(INPUT_SHAPEFILE, 'FID', INPUT_RASTER, OUTPUT_TABLE, "NODATA", "MEAN")
                        
                        # make temp layer for AddJoin to work; does not accept feature class (e.g., shapefile)
                        arcpy.MakeFeatureLayer_management(INPUT_SHAPEFILE, SHAPEFILE_LAYER_FOR_PROCESSING) 
                        
                        # join table for feature extraction
                        arcpy.AddJoin_management(SHAPEFILE_LAYER_FOR_PROCESSING, 'FID', OUTPUT_TABLE, 'FID_', 'KEEP_COMMON')
                        
                        # save out to csv 
                        arcpy.TableToTable_conversion(SHAPEFILE_LAYER_FOR_PROCESSING, outpath, OUTPUT_CSV)
                        
                        arcpy.Delete_management('GMIS_ras_process')
                        arcpy.Delete_management('GMIS_ras_NoData') # delete
                        arcpy.Delete_management('GMIS_ras_cleaned')                    
                        arcpy.Delete_management(SHAPEFILE_LAYER_FOR_PROCESSING)
                        
                        
        #                else:
        #                    print "Raster %s is outside" % (imperv_ras)            
 
    
def DEM_extract():
    # DEM PCS
    arcpy_env = r'D:\GLUE Results\Process_Transect\pcs'
    cities_lookup = lookup_table(r'D:\Python scripts - FINAL')
    vector_file_wildcard = '*pcs.shp' # all shapefiles going forward will use this naming convention
    vector_file = 'DEM'
    raster_file = raster(r'D:\GLUE Datasets\DEM', 'pcs.tif')
    extracted_feature = 'DEM_pcs'
    outpath = r'D:\GLUE Results\Extracted_DEM'
    extract(arcpy_env, cities_lookup, vector_file, vector_file_wildcard, raster_file, extracted_feature, outpath)

def NDSI_extract():
    # NDSI PCS
    arcpy_env = r'D:\GLUE Results\Process_Transect\pcs_buffer'
    cities_lookup = lookup_table(r'D:\Python scripts - FINAL')
    vector_file_wildcard = '*pcs_buffer.shp' # all shapefiles going forward will use this naming convention
    vector_file = 'DEFAULT'
    raster_file = raster(r'G:\Landsat_Download', 'NDSI.tif')     
    extracted_feature = 'NDSI_pcs'
    outpath = r'D:\GLUE Results\Extracted_NDSI'
    extract(arcpy_env, cities_lookup, vector_file, vector_file_wildcard, raster_file, extracted_feature, outpath)


def main_ndvi_extract():
    # NDVI PCS
    arcpy_env = r'D:\GLUE Results\Process_Transect\pcs_buffer'
    cities_lookup = lookup_table(r'D:\Python scripts - FINAL')
    vector_file_wildcard = '*pcs_buffer.shp' # all shapefiles going forward will use this naming convention
    vector_file = 'DEFAULT'
    raster_file_all = raster(r'G:\Landsat_Download', 'NDVI.tif')
    
    ndvi_summer_raster = []
    ndvi_winter_raster = []
    
    for i in range(len(raster_file_all)):
        if 'image_composites_summer_ndvi' in raster_file_all[i]:
            ndvi_summer_raster.append(raster_file_all[i])
        elif 'image_composites_winter_ndvi' in raster_file_all[i]:
            ndvi_winter_raster.append(raster_file_all[i])  
            
    extracted_feature = 'NDVI_pcs'
    ndvi_summer_outpath = r'D:\GLUE Results\Extracted_summer_NDVI'
    ndvi_winter_outpath = r'D:\GLUE Results\Extracted_winter_NDVI'
    
    extract(arcpy_env, cities_lookup, vector_file, vector_file_wildcard, ndvi_summer_raster, extracted_feature, ndvi_summer_outpath)
    
    extract(arcpy_env, cities_lookup, vector_file, vector_file_wildcard, ndvi_winter_raster, extracted_feature, ndvi_winter_outpath)


def main_lst_extract():
    
    arcpy_env = r'D:\GLUE Results\Process_Transect\pcs_buffer'
    cities_lookup = lookup_table(r'D:\Python scripts - FINAL')
    vector_file_wildcard = '*pcs_buffer.shp' # all shapefiles going forward will use this naming convention
    vector_file = 'LST'
    raster_file_all = raster(r'G:\Landsat_Download', 'LST_Celsius.tif')
    
    lst_summer_raster = []
    lst_winter_raster = []
    
    for i in range(len(raster_file_all)):
        if 'image_lst_summer' in raster_file_all[i]:
            lst_summer_raster.append(raster_file_all[i])
        elif 'image_lst_winter' in raster_file_all[i]:
            lst_winter_raster.append(raster_file_all[i])  
            
    extracted_feature = 'LST_pcs'
    lst_summer_outpath = r'D:\GLUE Results\Extracted_summer_LST'
    lst_winter_outpath = r'D:\GLUE Results\Extracted_winter_LST'
    
    extract(arcpy_env, cities_lookup, vector_file, vector_file_wildcard, lst_summer_raster, extracted_feature, lst_summer_outpath)
    
    extract(arcpy_env, cities_lookup, vector_file, vector_file_wildcard, lst_winter_raster, extracted_feature, lst_winter_outpath)


def extract_CGIAR_main():           
#vector_file = r'D:\GLUE Results\Process_City_Location\select_cities.shp'
    vector_file_directory = r'D:\GLUE Results\Process_Transect\gcs'  
    root_directory = r'D:\GLUE Datasets'  
    
    outpath_AI_YR = r'D:\GLUE Results\Extracted_annual_Aridity'
    outpath_PET_HE_YR = r'D:\GLUE Results\Extracted_annual_PET'
    outpath_PET_HE_MONTHLY = r'D:\GLUE Results\Extracted_monthly_PET'
    
    extract_CGIAR(vector_file_directory, root_directory, outpath_AI_YR, 'AI_YR')
    extract_CGIAR(vector_file_directory, root_directory, outpath_PET_HE_YR, 'PET_HE_YR')
#    extract_CGIAR(vector_file_directory, root_directory, outpath_PET_HE_MONTHLY, 'PET_HE_MONTHLY')   
#    extract_CGIAR(vector_file_directory, root_directory, outpath, 'PET_HE_MONTHLY','PET_HE_YR' ,'AI_YR')

def batch_process_impervious_surface():
    '''
    '''
#    import os, errno
    
    arcpy_env = r'D:\GLUE Results\Process_Transect\pcs_buffer'
    outpath = r'D:\GLUE Results\Extracted_GMIS'
   
    directory = r'D:\GLUE Datasets\SEDAC\Impervious Surface'
    
    # step 1: get folder directory names 
    folders = []    
    for root, dirnames, filenames in os.walk(directory):
        if len(dirnames) > 1:
            for i in dirnames:
                folders.append(os.path.join(root,i))
                
    # step 2: this is necessary because some GMIS rasters transcend geopolitical boundaries... 
    #         therefore may overwrite existing output csv (e.g., Toronto has a GMIS raster for CAN and USA; when the
    #         function iterates into the USA GMIS raster it overwrites the output of CAN GMIS raster)       
    for i in range(len(folders)):
#        print folders[i]
        if not os.path.exists(os.path.join(outpath,os.path.basename(folders[i]))):
            
            try:
                os.makedirs(os.path.join(outpath,os.path.basename(folders[i])))  
                
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
        else:
            print ('{0} already exists!'.format(os.path.join(outpath, os.path.basename(folders[i])))) 
                
    # step 3: for each folder directory name, go into sub-folder only if tif exists
    #         each folder is now it's own individual list for processing into process_impervious_surface()              
    for i in range(len(folders)):
        for root,dirnames,filenames in os.walk(folders[i]):
            if root.endswith('gmis_impervious_surface_percentage_utm_30m'): #this folder contains rasters for processing  
                GMIS_rasters = raster(root, '.tif')
                
                print os.path.basename(folders[i]) 
                print GMIS_rasters
                print '....' + '\n' + '.....'
                
                process_impervious_surface(arcpy_env, GMIS_rasters, os.path.join(outpath,os.path.basename(folders[i])))
  
     

if __name__ == '__main__':
  
    DEM_extract()
    NDSI_extract()
    main_ndvi_extract()
    main_lst_extract()
    extract_CGIAR_main()   
    batch_process_impervious_surface() 
          
    
# for LST (weather station points)
#arcpy_env = r'D:\validate_landsat_lst_toronto\espa-alexander.tong@mail.utoronto.ca-09182018-093502-904\weather station'
#in_list = ['Toronto', '018030'] 
#shapefile = r'D:\validate_landsat_lst_toronto\espa-alexander.tong@mail.utoronto.ca-09182018-093502-904\weather station\weather_station_points.shp'
#raster_file = raster(r'D:\validate_landsat_lst_toronto\espa-alexander.tong@mail.utoronto.ca-09182018-093502-904\image_LST','tif')
#extracted_feature = 'LST'
#inras_root = r'D:\validate_landsat_lst_toronto\espa-alexander.tong@mail.utoronto.ca-09182018-093502-904\image_LST'
#outpath = r'D:\validate_landsat_lst_toronto\espa-alexander.tong@mail.utoronto.ca-09182018-093502-904\weather station extracted\New folder'
#
#extract_LST_from_weather_stations(arcpy_env, shapefile, in_list, raster_file, extracted_feature, inras_root, outpath)

