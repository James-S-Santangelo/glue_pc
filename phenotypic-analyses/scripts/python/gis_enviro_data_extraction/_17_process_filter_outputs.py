# -*- coding: utf-8 -*-
"""
Created on Sat Feb  9 15:28:09 2019

@author: Alexander Tong

Developed and tested with Python 2.7.15
"""
try: 
    import os, sys
    
    directory = r'D:\Python scripts - FINAL'
    os.chdir(directory)
    
    from _00_glue_utils import raster
    from _00_glue_utils import lookup_table
    
except ImportError as IE:
    print (IE)
    print ("These functions requires arcpy to run")  
    sys.exit(1)    

# include this and arcpy licensing ad nauseum 
if sys.version_info[0] != 2:
    print("This script requires Python version 2.xx")
    sys.exit(1)

import os
from shutil import copyfile
import pandas as pd

# 1. for GLUE project, need to create dataframes for each tab for parsing... 
#folder = r'D:\GLUE Results\Process_Atmospheric_Functions'
folder = r'D:\GLUE Results'

#excel_file = 'combined_atmospheric_functions.xlsx'
excel_file = 'Final_Filtered_Outputs.xlsx'

# filter out columns using pandas instead of hard-coded index by cols using parse_cols
df_asia_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'asia_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()
df_asia_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'asia_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()

df_asia_summer[df_asia_summer.columns[1]] = df_asia_summer[df_asia_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
df_asia_winter[df_asia_winter.columns[1]] = df_asia_winter[df_asia_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')

####
df_canada_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'canada_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()
df_canada_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'canada_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()

df_canada_summer[df_canada_summer.columns[1]] = df_canada_summer[df_canada_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
df_canada_winter[df_canada_winter.columns[1]] = df_canada_winter[df_canada_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')

####
df_europe_et_al_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'europe_et_al_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()
df_europe_et_al_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'europe_et_al_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()

df_europe_et_al_summer[df_europe_et_al_summer.columns[1]] = df_europe_et_al_summer[df_europe_et_al_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
df_europe_et_al_winter[df_europe_et_al_winter.columns[1]] = df_europe_et_al_winter[df_europe_et_al_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')


####
df_oceania_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'oceania_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()
df_oceania_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'oceania_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()

df_oceania_summer[df_oceania_summer.columns[1]] = df_oceania_summer[df_oceania_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
df_oceania_winter[df_oceania_winter.columns[1]] = df_oceania_winter[df_oceania_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')


####
df_south_america_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'south_america_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()
df_south_america_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'south_america_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()

df_south_america_summer[df_south_america_summer.columns[1]] = df_south_america_summer[df_south_america_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
df_south_america_winter[df_south_america_winter.columns[1]] = df_south_america_winter[df_south_america_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')


####
df_usa_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'usa_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()
df_usa_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'usa_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City']).dropna()

df_usa_summer[df_usa_summer.columns[1]] = df_usa_summer[df_usa_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
df_usa_winter[df_usa_winter.columns[1]] = df_usa_winter[df_usa_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')


all_sheets = [df_asia_summer,df_asia_winter,
              df_canada_summer,df_canada_winter,
              df_europe_et_al_summer,df_europe_et_al_winter,
              df_oceania_summer,df_oceania_winter,
              df_south_america_summer,df_south_america_winter,
              df_usa_summer,df_usa_winter]



in_directory = r'D:\GLUE Results\Extracted_summer_NDVI'
in_directory = r'D:\GLUE Results\Extracted_winter_NDVI'
in_directory = r'D:\GLUE Results\Extracted_NDSI'

out_directory = r'D:\GLUE Results\Extracted_summer_NDVI_filtered'
out_directory = r'D:\GLUE Results\Extracted_winter_NDVI_filtered'
out_directory = r'D:\GLUE Results\Extracted_NDSI_filtered'


def filter_spectral_index_outputs(in_directory, out_directory):
    '''
    Description:
        # logic only setup to handle partial string match for NDVI, NDSI
        
        # modify for LST
        
    Args:
        in_directory(str):
        out_directory(str):
            
    Returns:
        filtered output csv files 
    '''
    for files in os.listdir(in_directory):
        if files.endswith('.csv'):
            
            if 'NDVI' or 'NDSI' in in_directory:
                
                output_row_path_date = files[-28:-13] # get row/path and acquistion date     
                output_city_name = files[:-29] # get city name 
            
            if 'LST' in in_directory:
                                               
                output_row_path_date = files[-27:-12] # get row/path and acquistion date     
                output_city_name = files[:-28] # get city name 
                
#                print output_city_name,  output_row_path_date
                
            for sheet in range(len(all_sheets)):
                for row in range(len(all_sheets[sheet])):
                    
                    # isolate for landsat row/path + acquisition date from dataframe
                    a = all_sheets[sheet].iloc[row,0][:10]
                    b = all_sheets[sheet].iloc[row,0][-23:]
                    
                    raster_row_path_date =  all_sheets[sheet].iloc[row,0].split(a)[-1].split(b)[0]
                    
                    
                    # filter out for final spectral image product output 
                    # match correct landsat image for each city to the multiple output landsat image ouputs generated for spectral image products
                    if output_city_name == all_sheets[sheet].iloc[row,1].lower() and output_row_path_date == raster_row_path_date: 
                                                
                        if os.path.exists(os.path.join(out_directory,files)):
#                            print ('{0} is already processed...'.format(files))
                            pass
                        else:
                            # extract to new folder...                          
                            # do partial string match for csv to be kept (matching the correct image output to avoid multiple/duplicate) 
                            # -> os.copy, because if the feature extration script is re-run it will generate all outputs again... 
                            copyfile( os.path.join(in_directory,files), os.path.join(os.path.join(out_directory,files)) )
                            
                            print files
    #                        print raster_row_path_date, all_sheets[sheet].iloc[row,1].lower()  # get city name 
                     

def filter_out_landsat_ndvi_ndsi_lst():
    '''
    Description:
        filter out NDVI, NDSI and LST output csv.
    '''
    in_directories = [r'D:\GLUE Results\Extracted_summer_NDVI',
                      r'D:\GLUE Results\Extracted_winter_NDVI',
                      r'D:\GLUE Results\Extracted_NDSI',
                      r'D:\GLUE Results\Extracted_summer_LST',
                      r'D:\GLUE Results\Extracted_winter_LST']
    
    out_directories = [r'D:\GLUE Results\Extracted_summer_NDVI_filtered',
                       r'D:\GLUE Results\Extracted_winter_NDVI_filtered',
                       r'D:\GLUE Results\Extracted_NDSI_filtered',
                       r'D:\GLUE Results\Extracted_summer_LST_filtered',
                       r'D:\GLUE Results\Extracted_winter_LST_filtered']

    
    # match in and out directory and process... 
    for in_directory in range(len(in_directories)):  
        for out_directory in range(len(out_directories)):           
            if os.path.basename(in_directories[in_directory]) in os.path.basename(out_directories[out_directory]):
                
                in_dir = in_directories[in_directory]
                out_dir = out_directories[out_directory]
                                
                filter_spectral_index_outputs(in_dir, out_dir)
                
                print os.path.basename(in_directories[in_directory]), os.path.basename(out_directories[out_directory])
            


## for csv (NDVI/NDSI)
#test = 'woodstock_019030_20150728_NDVI_pcs.csv'    
#test[-28:-13] # get row/path and acquistion date 
#test[:-29] # get city name 
#
## for csv (LST)
#test = 'acton_018030_20180611_LST_pcs.csv'
#test[-27:-12] # get row/path and acquistion date 
#test[:-28] # get city name 
#    
  
def filter_out_CGIAR():
    '''
    Description:
        filter out CGIAR datasets (e.g., Aridity, PET)
    '''
    in_directories = [r'D:\GLUE Results\Extracted_DEM',
                      r'D:\GLUE Results\Extracted_annual_Aridity',
                      r'D:\GLUE Results\Extracted_annual_PET',
                      r'D:\GLUE Results\Extracted_monthly_PET']
    
    out_directories = [r'D:\GLUE Results\Extracted_DEM_filtered',
                       r'D:\GLUE Results\Extracted_annual_Aridity_filtered',
                       r'D:\GLUE Results\Extracted_annual_PET_filtered',
                       r'D:\GLUE Results\Extracted_monthly_PET_filtered']
          
    # match in and out directory and process... 
    for in_directory in range(len(in_directories)):  
        for out_directory in range(len(out_directories)):           
            if os.path.basename(in_directories[in_directory]) in os.path.basename(out_directories[out_directory]):

                in_dir = in_directories[in_directory]
                out_dir = out_directories[out_directory]
                
#                print os.path.basename(in_directories[in_directory]), os.path.basename(out_directories[out_directory])
                
                for files in os.listdir(in_dir):
                    if files.endswith('.csv'):
                        
                        copyfile( os.path.join(in_dir,files),  os.path.join(out_dir,files) )

def filter_out_dem():
    '''
    Description:
        filter out DEM datasets 
        
    # for dem, list of filepaths, isolate for city names, and create set list; ones that get dropped, remove from folder... 

    '''
    in_directory = r'D:\GLUE Results\Extracted_DEM'
    
    out_directory = r'D:\GLUE Results\Extracted_DEM_filtered'
    
    files_list = []
    files_list_city_name_isolated = []
    files_list_filtered = []
    
    for i in os.listdir(in_directory):
        if i.endswith('.csv'):
            files_list.append(os.path.join(in_directory,i))
           
    # sort list
    files_list.sort()
            
    for j in files_list:
        # get city name 
        city_name_isolated = os.path.basename(j[-19:])
        city_name_contained = os.path.basename(j).split(city_name_isolated)[0]
        
        files_list_city_name_isolated.append(city_name_contained)
        
#        print city_name_contained
    
    # Lists are mutable, therefore unhashable. Use tuples instead
    files_list_compare = [[i,j] for i,j in zip(files_list,files_list_city_name_isolated)] # creating list of list 
    
    
    # match to filtered city list    
    for k in range(len(files_list_compare)):
        
        # isolate for single instance; if same as previous
        if files_list_compare[k][1] == files_list_compare[k-1][1]:
#            print files_list_compare[k]
            pass
        else:
            files_list_filtered.append(files_list_compare[k][0])
    
    # copy out       
    for csv_files in files_list_filtered:
        print csv_files
        copyfile(csv_files, os.path.join(out_directory,os.path.basename(csv_files)))


def cities(directory):
    '''
    Description: 
        retrieve master list of lists of city to Landsat code for processing
        
    Args:
        directory (str): specify directory of Landsat_Code_For_Cities.py  
        
    Returns:
        list of list 
    '''
    import os
    try:
        os.chdir(directory)
        from Landsat_Code_For_Cities import cities

        return cities() 
    except:
        print 'could not locate list' 
        

def filter_out_GMIS():
    '''
    Description
        issue with GMIS rasters and the logic applied to extracting values is that \n\
        the GMIS rasters trascend geopolitical boundaries (country boundaries) \n\
        so the extraction process may be extracting 'Portland_ME' for Canada for \n\
        instance. This was remedied by extracting the outputs into individual country \n\
        folders, but then we require to actual filter out the cities... 
        
    $ apply filter to Canada and USA...
    
    '''
# for GMIS, use lookup table to filter out extraneous entries in each country... 
    
    cities_lookup = lookup_table(r'D:\Python scripts - FINAL') 
    
    in_root_directory = r'D:\GLUE Results\Extracted_GMIS'
    
    in_directory_canada = r'D:\GLUE Results\Extracted_GMIS\Canada'
    in_directory_usa = r'D:\GLUE Results\Extracted_GMIS\USA' 
    in_directory_france = r'D:\GLUE Results\Extracted_GMIS\France'
    in_directory_germany = r'D:\GLUE Results\Extracted_GMIS\Germany'
    in_directory_italy = r'D:\GLUE Results\Extracted_GMIS\Italy'
    in_directory_netherlands = r'D:\GLUE Results\Extracted_GMIS\Netherlands'
        
    GMIS_canada = []
    GMIS_usa = []
    GMIS_france = []
    GMIS_germany = []
    GMIS_italy = []
    GMIS_netherlands = []
    
    out_directory = r'D:\GLUE Results\Extracted_GMIS_filtered'
     
    
    for root, dirnames, filenames in os.walk(in_root_directory):
        if 'info' not in root:
            if 'Canada' in root:
                pass
            elif 'USA' in root:
                pass
            elif 'France' in root:
                pass
            elif 'Germany' in root:
                pass
            elif 'Italy' in root:
                pass     
            elif 'Netherlands' in root:
                pass        
            else:
                for files in range(len(filenames)):
                    if filenames[files].endswith('.csv'):
                        
                        print filenames[files]
                        
                        src_file = os.path.join(root,filenames[files])
                        outsrc_file = os.path.join(out_directory,filenames[files])
                        
                        copyfile(src_file, outsrc_file)    
    
    
    # get csv outputs as list... 
    for city in os.listdir(in_directory_canada):
        if city.endswith('.csv'):
            GMIS_canada.append(os.path.join(in_directory_canada,city))
    
    for city in os.listdir(in_directory_usa):
        if city.endswith('.csv'):
            GMIS_usa.append(os.path.join(in_directory_usa,city)) 
    
    for city in os.listdir(in_directory_france):
        if city.endswith('.csv'):
            GMIS_france.append(os.path.join(in_directory_france,city))     
    
    for city in os.listdir(in_directory_germany):
        if city.endswith('.csv'):
            GMIS_germany.append(os.path.join(in_directory_germany,city))
    
    for city in os.listdir(in_directory_italy):
        if city.endswith('.csv'):
            GMIS_italy.append(os.path.join(in_directory_italy,city))
            
    for city in os.listdir(in_directory_netherlands):
        if city.endswith('.csv'):
            GMIS_netherlands.append(os.path.join(in_directory_netherlands,city))              
           
            
    # filter out output csv from Canada/USA GMIS folders
    for i in range(len(cities_lookup)):
        
        # filter out for Canada GMIS
        if 'USA' in cities_lookup[i][0]:           
#            print cities_lookup[i][1]             
            for j in GMIS_canada:
                if j.endswith('.csv'):
                    if cities_lookup[i][1] in j:
                        GMIS_canada.remove(j) 
                        
        # filter out for USA GMIS
        if 'Canada' in cities_lookup[i][0]:
#            print cities_lookup[i][1]
            for j in GMIS_usa:
                if j.endswith('.csv'):
                    if cities_lookup[i][1] in j:                    
                        GMIS_usa.remove(j)
                                                        
        # filter out for France GMIS
        if 'Belgium' in cities_lookup[i][0]:
#            print cities_lookup[i][1]
            for j in GMIS_france:
                if j.endswith('.csv'):
                    if cities_lookup[i][1] in j:                    
                         GMIS_france.remove(j)
                                 
        # filter out for France GMIS                
        if 'Switzerland' in cities_lookup[i][0]:
#            print cities_lookup[i][1]
            for j in GMIS_france:
                if j.endswith('.csv'):
                    if cities_lookup[i][1] in j:                    
                        GMIS_france.remove(j)
    
        # filter out for Germany GMIS
        if 'Switzerland' in cities_lookup[i][0]:
#            print cities_lookup[i][1]
            for j in GMIS_germany:
                if j.endswith('.csv'):
                    if cities_lookup[i][1] in j:                    
                        GMIS_germany.remove(j)

        # filter out for Italy GMIS
        if 'Switzerland' in cities_lookup[i][0]:
#            print cities_lookup[i][1]
            for j in GMIS_italy:
                if j.endswith('.csv'):
                    if cities_lookup[i][1] in j:                    
                        GMIS_italy.remove(j)                        
                        
        # filter out for Netherlands GMIS
        if 'Belgium' in cities_lookup[i][0]:
#            print cities_lookup[i][1]
            for j in GMIS_netherlands:
                if j.endswith('.csv'):
                    if cities_lookup[i][1] in j:                    
                        GMIS_netherlands.remove(j)
     
                   
    # copy out for Canada/USA/France/Germany/Italy/Netherlands to filtered folder                    
    for csv_files in GMIS_canada:
        copyfile(csv_files, os.path.join(out_directory,os.path.basename(csv_files)))                    
    
    for csv_files in GMIS_usa:
        copyfile(csv_files, os.path.join(out_directory,os.path.basename(csv_files)))                    
        
    for csv_files in GMIS_france:
        copyfile(csv_files, os.path.join(out_directory,os.path.basename(csv_files)))                    
    
    for csv_files in GMIS_germany:
        copyfile(csv_files, os.path.join(out_directory,os.path.basename(csv_files)))           

    for csv_files in GMIS_italy:
        copyfile(csv_files, os.path.join(out_directory,os.path.basename(csv_files)))                    
    
    for csv_files in GMIS_netherlands:
        copyfile(csv_files, os.path.join(out_directory,os.path.basename(csv_files)))   
        

def filter_out_maps(in_directory, out_directory, map_type):
    '''
    Description:
        Filters out GLUE maps produced using natural color composite Landsat images \n\
        and land surface temperature Landsat images. 
        
        !!!WARNING!!! Only filters out for summer Landsat images. Not tested for \n\
        winter images.
        
    Args:
        in_directory (str): specify inpath directory for maps 
                            (e.g., r'D:\GLUE Results\GLUE_maps\output_summer')       
        out_directory (str): specify outpath directory for filtered maps 
                            (e.g., r'D:\GLUE Results\GLUE_maps_filtered')        
        map_type (str): set as 'natural' for natural color landsat images to be filtered
                        set as 'lst' for land surface temperature landsat images to be filtered
        
    Returns:
        Filtered maps
    '''

    for files in os.listdir(in_directory):
        if files.endswith('.png'):
            
            output_row_path_date = files[-29:-4] # get row/path and acquistion date     
            output_city_name = files[:-30] # get city name 
            
#            print output_city_name, output_row_path_date
        else:
            pass
    
        for sheet in range(len(all_sheets)):
            for row in range(len(all_sheets[sheet])):
                                    
                # isolate for landsat row/path + acquisition date from dataframe
                raster_row_path_date =  all_sheets[sheet].iloc[row,0][:25]
                
#                print all_sheets[sheet].iloc[row,1].lower(), raster_row_path_date
                
                # filter out for final spectral image product output 
                # match correct landsat image for each city to the multiple output landsat image ouputs generated for spectral image products
                if output_city_name.lower() == all_sheets[sheet].iloc[row,1].lower() and output_row_path_date == raster_row_path_date: 
#                    print files                 
                            
                    if os.path.exists(os.path.join(out_directory,files)):
#                        print ('{0} is already processed...'.format(files))
                        pass
                    else:
                        # extract to new folder...                          
                        # do partial string match for csv to be kept (matching the correct image output to avoid multiple/duplicate) 
                        # -> os.copy, because if the feature extration script is re-run it will generate all outputs again... 
                        
                        # -------- for natural color landsat images start ---------#
                        if map_type == 'natural':
                            
                            if files == 'Toronto_LC08_L1TP_018029_20180611.png' or files == 'Toronto_LC08_L1TP_018030_20180611.png':
                                files = 'Toronto_LC08_L1TP_018029_018030_2.png'
                                copyfile( os.path.join(in_directory,files), os.path.join(os.path.join(out_directory,files)) )
                                
                            elif files == 'New_York_LC08_L1TP_014031_20170730.png':
                                files = 'New_York_LC08_L1TP_014031_014032_0.png'
                                copyfile( os.path.join(in_directory,files), os.path.join(os.path.join(out_directory,files)) )    
                            
                            elif files == 'Santiago_LC08_L1TP_233083_233084_2.png':
                                files = 'Santiago_LC08_L1TP_233083_20170117.png' 
                                copyfile( os.path.join(in_directory,files), os.path.join(os.path.join(out_directory,files)) ) 
                            # -------- for natural color landsat images end ---------#
                            
                            else:
                                copyfile( os.path.join(in_directory,files), os.path.join(os.path.join(out_directory,files)) )
                        
                        elif map_type == 'lst':
                            # -------- for land surface temperature landsat images start ---------#
                            if files == 'Toronto_LC08_L1TP_018029_20180611.png' or files == 'Toronto_LC08_L1TP_018030_20180611.png':
                                files = 'Toronto_LC08_L1TP_018029_018030_20180611.png'
                                copyfile( os.path.join(in_directory,files), os.path.join(os.path.join(out_directory,files)) ) 
                            
                            elif files == 'New_York_LC08_L1TP_014031_20170730.png':
                                files = 'New_York_LC08_L1TP_014031_014032_20170730.png' 
                                copyfile( os.path.join(in_directory,files), os.path.join(os.path.join(out_directory,files)) ) 
                            # -------- for land surface temperature landsat images end ---------#

                        
                            else:
                                copyfile( os.path.join(in_directory,files), os.path.join(os.path.join(out_directory,files)) )
                        
                        # print not processed files 
                        print files
    #                        print raster_row_path_date, all_sheets[sheet].iloc[row,1].lower()  # get city name 
            

# isolate for 'LC08_L1TP_233084_20170117'

#test = 'LC08_L1TP_069017_20170221_20170301_01_T1_MTL.txt'                    
#test[:-23]      
              
if __name__ == "__main__":

    filter_out_landsat_ndvi_ndsi_lst()
    filter_out_CGIAR()    
    filter_out_dem()
    filter_out_GMIS()
##    
##    # natural color landsat maps; summertime 
##    in_directory = r'D:\GLUE Results\GLUE_maps\output_summer'
##    out_directory = r'D:\GLUE Results\GLUE_maps_filtered'
##    map_type = 'natural'
##    filter_out_maps(in_directory,out_directory)
#    
#    # land surface temperature landsat maps; summertime
##    in_directory = r'D:\GLUE Results\GLUE_maps\output_lst_summer'
##    out_directory = r'D:\GLUE Results\GLUE_maps_lst_filtered'
##    map_type = 'lst'
##    filter_out_maps(in_directory,out_directory,map_type)
##    
# 
#    
#    





            
            
