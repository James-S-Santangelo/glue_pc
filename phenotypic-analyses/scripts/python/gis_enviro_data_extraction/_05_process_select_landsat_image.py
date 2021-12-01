# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 16:03:08 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15

### http://www.naturalearthdata.com/downloads/10m-cultural-vectors/
### https://community.esri.com/thread/139977 
"""

import os, sys
import pandas as pd 

if sys.version_info[0] != 2:
    print("This script requires Python version 2.xx")
sys.exit(1)

try:
    import arcpy

except ImportError as IE:
    print (IE)
    print ("These functions requires arcpy to run")  
    sys.exit(1)

def check_countries():
    '''
    Description: 
        Checks if GLUE countries and cities are in shapefile provided by 
        Natural Earth Data; data processing/transformation (QA/QC)  
        
        Source: http://www.naturalearthdata.com/downloads/10m-cultural-vectors/
        
    Args: hard-coded within function 
    '''    
    arcpy.env.workspace = r'D:\GLUE Results\Process_City_Location'
    arcpy.env.overwriteOutput = True
    
    countries = ['Argentina', 'Australia', 'Belgium', 'Bolivia', 'Brazil', 'Canada', 'Chile', 'China', 'Colombia', 
                 'Czech Republic', 'Ecuador', 'Finland', 'France', 'Germany', 'Greece', 'Iran', 'Italy', 'Japan', 'Mexico', 
                 'Netherlands', 'New Zealand', 'Norway', 'Papua New Guinea', 'Poland', 'Portugal', 'Russia', 'South Africa', 
                 'Spain', 'Sweden', 'Switzerland', 'United Kingdom', 'United States of America']


    # verify GLUE countries are in shapefile                       
    for fc in arcpy.ListFeatureClasses():
        if 'all_country_cities' in fc:
            with arcpy.da.SearchCursor(fc, 'adm0name') as cursor:    
                unique_val = []
                for row in cursor:
                    for i in countries:
                        if i in row:
                            if not row[0] in unique_val:
                                unique_val.append(row[0])
    for val in unique_val:
        print val 
        
    # check if GLUE city in shapefile; if not, manually parse attribute data to see why not (e.g., spelling is different)    
    for fc in arcpy.ListFeatureClasses():
        if 'all_country_cities' in fc:
            with arcpy.da.SearchCursor(fc, ['nameascii','adm0name']) as cursor:    
                for row in cursor:
                    if 'Quebec City' in row[0]:
                        print row
                        
  
def create_GLUE_cities_shp():
    """
    Description: 
        Create shapefile of only select countries and cities from \n\
        shapefile with all countries and cities retrieved from Natural Earth Data. \n\
        The point locations of these GLUE cities are to be used to extract temperature \n\
        data from WorldClim to determine which month has the maximum and minimum \n\
        temperatures to represent the summer and winter period respectively for \n\
        Landsat image retrieval. 
    
    Args: 
        hard-coded within function 

    """
#    import arcpy
    
    in_dir = r'D:\GLUE Results\Process_City_Location' 
    out_dir = r'D:\GLUE Results\Process_City_Location'
    
    # COPY OUT ROW BY SELECT COUNTRIES 
    where_clause =  """"adm0name"  IN ('Argentina', 'Australia', 'Belgium', 'Bolivia', 'Brazil', 'Canada', 'Chile', 'China', 'Colombia', 
                 'Czech Republic', 'Ecuador', 'Finland', 'France', 'Germany', 'Greece', 'Iran', 'Italy', 'Japan', 'Mexico', 
                 'Netherlands', 'New Zealand', 'Norway', 'Papua New Guinea', 'Poland', 'Portugal', 'Russia', 'South Africa', 
                 'Spain', 'Sweden', 'Switzerland', 'United Kingdom', 'United States of America')""" 
    
    arcpy.Select_analysis(in_dir + '\\' + 'all_country_cities.shp', out_dir + '\\' + 'select_countries.shp', where_clause)


def check_cities():
    """
    Description: 
        handshake check - check "nameascii"  name with "adm0name"  
    
    WARNING: DEPRECATED FUNCTION; THE TEXTFILES USED IN THIS FUNCTION HAVE NOT BEEN UPDATED TO REFLECT CURRENT CITIES IN GLUE PROJECT
    """
    
    os.chdir(r'D:\GLUE Results\Process_City_Location')
    cities = []
    
    with open('cities.txt','r') as infile:
        for line in infile:
            cities.append(line.strip('\n'))
          
    for i in cities:
        print i 
        
    arcpy.env.workspace = r'D:\GLUE Results\Process_City_Location'
    #verify cities 
    for fc in arcpy.ListFeatureClasses():
        if 'select_countries' in fc:
            with arcpy.da.SearchCursor(fc, ['nameascii','adm0name']) as cursor:
                unique_val = []
                for row in cursor:
                    for i in cities:
                        if i in row:
                            if not row[0] in unique_val:
                                unique_val.append(row)
    
    for val in unique_val:
        print val 
        
    with open(r'D:\GLUE Results\Process_City_Location\verify_cities.txt','w') as outfile:
        for val in unique_val:  
            outfile.write(str(val) + '\n')


def create_GLUE_cities_shp_refined():
    '''
    Description:
        Select out GLUE cities from shapefile previously created (e.g., all_country_cities.shp) \n\
        from create_GLUE_cities_shp(). The point locations of these GLUE cities \n\
        are to be used to extract temperature data from WorldClim to determine \n\
        which month has the maximum and minimum temperatures to represent the \n\
        summer and winter period respectively for Landsat image retrieval. 
        
    Args:
        hard-coded within function 
        
    '''
    filepath = r'D:\GLUE Results\Process_City_Location'
    arcpy.env.workspace = filepath
    arcpy.env.overwriteOutput = True
        

    ### process USA cities 
    usa_city = []
    with open(os.path.join(filepath,'usa_city.txt'),'r') as infile:
        for line in infile:
            usa_city.append(line.strip('\n'))
    
    usa_state = []
    with open(os.path.join(filepath,'usa_state.txt'),'r') as infile:
        for line in infile:
            usa_state.append(line.strip('\n'))    
    
    usa_state_city = [list(i) for i in zip(usa_city,usa_state)]
    
    
    ### process ROW cities
    row_city = []
    with open(os.path.join(filepath,'ROW_city.txt'),'r') as infile:
        for line in infile:
            row_city.append(line.strip('\n'))
    
    row_country = []
    with open(os.path.join(filepath,'ROW_country.txt'),'r') as infile:
        for line in infile:
            row_country.append(line.strip('\n'))    
    
    row_country_city = [list(i) for i in zip(row_city,row_country)]
        
    
    # GRAB OBJECTID AND CREATE NEW SHAPEFILE
    # get object id (OID@) for each city 
    # further filter cities
    for fc in arcpy.ListFeatureClasses():
        if 'select_countries' in fc:
            with arcpy.da.SearchCursor(fc, ["OID@",'nameascii','adm0name','adm1name']) as cursor:
                cities = []
                for row in cursor:
                    for j in range(len(usa_state_city)):
                        if usa_state_city[j][0] == row[1] and usa_state_city[j][1] == row[3]:
                            print ('{0},{1}'.format(row[1], row[3]))
                            
                            # grab object id 
                            if not row in cities:
                                cities.append(row[0])
    
                    for k in range(len(row_country_city)):
                        if row_country_city[k][0] == row[1] and row_country_city[k][1] == row[2]:
                            print  ('%s, %s' % (row[1], row[2]))
                            
                            # grab object id 
                            if not row in cities:
                                cities.append(row[0])
    
    # to use this, change cities.append(row[0]) to row
#    for city in range(len(cities)):
#        print ('{0}, {1}'.format(cities[city][1], cities[city][2]))
                       
    # convert list city object ids (OID@) into tuple for processing 
    cities_tuple = tuple(cities)
      
    fc = r'D:\GLUE Results\Process_City_Location\select_countries.shp'

    where_clause = """"FID" IN """ + str(cities_tuple)  
        
    # process out final shapefile with select cities
    arcpy.Select_analysis(fc, 'select_cities.shp', where_clause)

    del cursor


def check_cities_final():
    '''
    '''
    arcpy.env.workspace = r'D:\GLUE Results\Process_City_Location'
    
    # final check for all cities 
    for fc in arcpy.ListFeatureClasses():
        if 'select_cities' in fc:
            with arcpy.da.SearchCursor(fc, ["OID@",'nameascii','adm0name','adm1name']) as cursor:
                cities = []
                for row in cursor:
                    cities.append(row)
    
    for city in range(len(cities)):
        print ('{0}, {1}'.format(cities[city][1], cities[city][2]))
        

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

    
def process_monthly_temp_avg(directory, outpath):
    '''
    Description: 
        From process_extract_feature.py, extract_WorldClim() is used to extract \n\
        WorldClim() data at GLUE city locations, one corresponding to each month \n\
        of the year, totaling 12 csv. This function aggregates the results of these \n\
        12 csv and computes the month with the maximum and minimum temperature to \n\
        represent the summer and winter period respectively for Landsat \n\
        image retrieval. 

                 
     **WARNING: 
        the directory and any sub-folders can ONLY have the 12 csv's  (e.g., Jan-Dec) \n\
        for processing in the main directory. The function DOES NOT HANDLE \n\
        any exceptions to this logic, or else it breaks. 
        
    Args:
        directory (str): specify directory of csv files to be processed (e.g., 1_city_AVG_Temp_Jan.csv)
        
    Returns:
        No returns or exchanges 
        
    $ to be implemented: make logic more robust to handle other sub-folders (e.g., with deprecated results)
    '''
    num_of_csv = 1
    for root, dirnames, filenames in os.walk(directory):   
        for file in range(len(filenames)):
            if filenames[file].endswith('.csv'):
                num_of_csv += 1
        
    first_csv = 1
    df_list = []
    for root, dirnames, filenames in os.walk(directory):   
        for file in range(len(filenames)):            
            if filenames[file].endswith('.csv'):   
                
                # grab city, country data once and use it for all subsequent months; in this case, use Jan file 
                while first_csv > 0:
                    
                    df_city_country_raw = pd.read_csv(os.path.join(root,filenames[file]), usecols=['nameascii', 'adm0name', 'adm1name'])
                    
                    city_country_list = df_city_country_raw.values.tolist()
                    
                    # proof of logic for re-combining for final output                 
                    df_city_country_processed = pd.DataFrame(city_country_list)
    
                    first_csv -= 1                       
            
      
                # grab each month temperature average                    
                df_temp = pd.read_csv(os.path.join(root,filenames[file])).iloc[:,43:44]
                print df_temp.head()

                
                # because pandas converts values as list of lists, need to convert back to a single list for analysis
                temp = df_temp.values.tolist()
                
                flat_list_temp = [item for sublist in temp for item in sublist]


                df_list.append(flat_list_temp)
                
                
                # pandas will read a list as a single row; therefore need to to transpose to read as series 
                df_temp_raw = pd.DataFrame(df_list)
                df_temp_transposed = df_temp_raw.transpose()
                
                df_temp_transposed.head()
    
    
    # combine headers with data         
    frames = [df_city_country_processed,df_temp_transposed] 
    df_combined = pd.concat(frames, axis=1) #axis == 1 means combine in paralell, not add as series  
    
    # change header names and sort by country
    # axis : {index (0), columns (1)} 
    # https://stackoverflow.com/questions/11346283/renaming-columns-in-pandas
    df_combined.set_axis(['City', 'Country', 'Alt_country_id', 'Jan', 'Feb','Mar','Apr', 'May','Jun','Jul', 'Aug','Sep','Oct', 'Nov','Dec'], axis=1, inplace=True)            
    df_combined_sorted = df_combined.sort_values(by=['Country'])
    df_combined_sorted.head()
    
    df_combined.head()
    
    
    # get max for winter/summer at select (iloc) float values
    ### https://stackoverflow.com/questions/29919306/find-the-column-name-which-has-the-maximum-value-for-each-row
    df_combined_sorted['max_temp_month'] = df_combined_sorted.iloc[:,3:14].idxmax(axis=1)
    df_combined_sorted['min_temp_month'] = df_combined_sorted.iloc[:,3:14].idxmin(axis=1)
    
    df_combined_sorted['max_temp'] = df_combined_sorted.iloc[:,3:14].max(axis=1)
    df_combined_sorted['min_temp']  = df_combined_sorted.iloc[:,3:14].min(axis=1)
    
    
    # sort
    df_combined_sorted.head()
    
    # save 
    df_combined_sorted.to_csv(os.path.join(outpath,'Select_Month_Summer_Winter_For_City.csv'))



def main_1():
    '''
    Description:
        The below functions require arcpy to execute. Run these functions in the \n\
        Python console in the ArcMap GUI. Alternatively if the Python path is set to \n\
        the IDE you are using, then it can be executed in the IDE console. 
    '''
    check_countries() 
    create_GLUE_cities_shp()
    check_cities()
    create_GLUE_cities_shp_refined()    
    check_cities_final()    

     
def main_2():
    '''
    Description:
        extract WorldClim from shapefile created from main_1()
    '''
    import arcpy
    from arcpy.sa import *
    
    try: 

        import os, sys
            
        directory = r'D:\Python scripts - FINAL'
        os.chdir(directory)
        
        from _00_glue_utils import raster         
                
        arcpy_env = r'D:\GLUE Results\Process_City_Location'
        vector_file = r'D:\GLUE Results\Process_City_Location\select_cities.shp'
        raster_file = raster(r'D:\GLUE Datasets\WorldClim_V2\wc2.0_30s_tavg','tif')
        extracted_feature = 'AVG_Temp'
        outpath = r'D:\GLUE Results\Process_WorldClim_V2'
        
        extract_WorldClim(arcpy_env, vector_file, raster_file, extracted_feature, outpath)
    
    except ImportError as IE:
        
        print (IE)
        print ("These functions requires arcpy to run")  
        sys.exit(1)  
    
    
def main_3():
    '''
    Description:
        refine results from main_2() as a single csv
    '''
    directory = r'D:\GLUE Results\Process_WorldClim_V2'
    outpath = r'D:\GLUE Results\Process_WorldClim_V2'
    process_monthly_temp_avg(directory, outpath)
    

if __name__ == '__main__':
    
    main_1()
    #main_2()
    #main_3()
    #    
    
    
    
    