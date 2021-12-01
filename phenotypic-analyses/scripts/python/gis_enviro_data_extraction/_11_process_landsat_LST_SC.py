# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 14:55:11 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15

This module assumes that Landsat level 2 products (e.g., surface reflectance \n\
datasets) including NDVI and BT rasters have been calculated/retrieved \n\
from the ESPA-EROS data center (https://espa.cr.usgs.gov/index/). Herein, land \n\ 
surface  temperature is calculated using a NDVI-based emissivity method (NBEM) \n\
single-channel algorithm. For ecological purposes, this method should be sufficient, \n\
with errors of +- 3 Kelvin (worst case scenario). The authors have indicated that \n\
the NBEM approach for LST retreival is not suitable for greybodies, rocks,  \n\
or senescent vegetation. 
 
"""

import os, sys
import pandas as pd

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


# include this and arcpy licensing ad nauseum 
if sys.version_info[0] != 2:
    print("This script requires Python version 2.xx")
    sys.exit(1)
    
#--------------------- DEPRECATED; validation test start ----------------------#
## 1. for GLUE project, need to create dataframes for each tab for parsing... 
#folder = r'D:\GLUE Results\validate_elevation_differences_lst'
#excel_file = 'GLUE_City_LST_Retrieval.xlsx'
#
## filter out columns using pandas instead of hard-coded index by cols using parse_cols
#df_station_elev_test = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'Elevation_Data', index = True).filter(items = ['SCENE_NAME','t_station','Lu_station','Ld_station']).dropna()
#df_transect_elev_test = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'Elevation_Data', index = True).filter(items = ['SCENE_NAME','t_transect','Lu_transect','Ld_transect']).dropna()
#--------------------- DEPRECATED; validation test end -----------------------#


# 1. for GLUE project, need to create dataframes for each tab for parsing... 
folder = r'D:\GLUE Results\Process_Atmospheric_Functions'
excel_file = 'combined_atmospheric_functions.xlsx'

# filter out columns using pandas instead of hard-coded index by cols using parse_cols
##df_asia_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'asia_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
##df_asia_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'asia_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
# wrangle date acquisition date at [1]
##df_asia_summer[df_asia_summer.columns[1]] = df_asia_summer[df_asia_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
##df_asia_winter[df_asia_winter.columns[1]] = df_asia_winter[df_asia_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
##
####
##df_canada_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'canada_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
##df_canada_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'canada_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
# wrangle date acquisition date at [1]
##df_canada_summer[df_canada_summer.columns[1]] = df_canada_summer[df_canada_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
##df_canada_winter[df_canada_winter.columns[1]] = df_canada_winter[df_canada_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')

####
##df_europe_et_al_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'europe_et_al_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
##df_europe_et_al_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'europe_et_al_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
# wrangle date acquisition date at [1]
##df_europe_et_al_summer[df_europe_et_al_summer.columns[1]] = df_europe_et_al_summer[df_europe_et_al_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
##df_europe_et_al_winter[df_europe_et_al_winter.columns[1]] = df_europe_et_al_winter[df_europe_et_al_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')


####
##df_oceania_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'oceania_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
##df_oceania_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'oceania_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
# wrangle date acquisition date at [1]
##df_oceania_summer[df_oceania_summer.columns[1]] = df_oceania_summer[df_oceania_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
##df_oceania_winter[df_oceania_winter.columns[1]] = df_oceania_winter[df_oceania_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')


####
##df_south_america_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'south_america_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
##df_south_america_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'south_america_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
# wrangle date acquisition date at [1]
##df_south_america_summer[df_south_america_summer.columns[1]] = df_south_america_summer[df_south_america_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
##df_south_america_winter[df_south_america_winter.columns[1]] = df_south_america_winter[df_south_america_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')


####
df_usa_summer = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'usa_summer_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
df_usa_winter = pd.read_excel(os.path.join(folder, excel_file), sheetname = 'usa_winter_atmos_func', index = True, encoding='utf-8').filter(items =['SCENE_NAME','City','t','Lu','Ld']).dropna()
# wrangle date acquisition date at [1]
df_usa_summer[df_usa_summer.columns[1]] = df_usa_summer[df_usa_summer.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')
df_usa_winter[df_usa_winter.columns[1]] = df_usa_winter[df_usa_winter.columns[1]].str.replace(' ', '_').str.replace("-", '_').str.replace("'", '').str.replace(".", '')


##all_sheets = [df_asia_summer,df_asia_winter,
##              df_canada_summer,df_canada_winter,
##              df_europe_et_al_summer,df_europe_et_al_winter,
##              df_oceania_summer,df_oceania_winter,
##              df_south_america_summer,df_south_america_winter,
##              df_usa_summer,df_usa_winter]
##
##all_sheets_name = ['asia_summer_atmos_func','asia_winter_atmos_func',
##                    'canada_summer_atmos_func', 'canada_winter_atmos_func',
##                    'europe_et_al_summer_atmos_func', 'europe_et_al_winter_atmos_func',
##                    'oceania_summer_atmos_func', 'oceania_winter_atmos_func',
##                    'south_america_summer_atmos_func', 'south_america_winter_atmos_func',
##                    'usa_summer_atmos_func', 'usa_winter_atmos_func']


##all_sheets = [df_asia_summer,df_asia_winter,]
##all_sheets_name = ['asia_summer_atmos_func','asia_winter_atmos_func']
##
##all_sheets = [df_oceania_summer,df_oceania_winter,]
##all_sheets_name = ['oceania_summer_atmos_func','oceania_winter_atmos_func']
##
##all_sheets = [df_canada_summer,df_canada_winter,]
##all_sheets_name = ['canada_summer_atmos_func','canada_winter_atmos_func']
##
##all_sheets = [df_europe_et_al_summer,df_europe_et_al_winter,]
##all_sheets_name = ['europe_et_al_summer_atmos_func','europe_et_al_winter_atmos_func']
##
##all_sheets = [df_south_america_summer,df_south_america_winter]
##all_sheets_name = ['south_america_summer_atmos_func', 'south_america_winter_atmos_func']

all_sheets = [df_usa_summer,df_usa_winter]
all_sheets_name = ['usa_summer_atmos_func', 'usa_winter_atmos_func']

#--------------------

def main_all_images():
    '''
    $future implementation: match spatial ref between TIRS, BT, NDVI rasters... 
    $future implementation: save print out as csv log... 

    '''
    import os
    import pandas as pd
    
    # will need to match 3 images for computation 
    #    red_dir = r'D:\validate_landsat_lst_toronto\espa-alexander.tong@mail.utoronto.ca-09182018-093502-904\image_band_rescale' # NOT USED; ONLY USE IF USING FULL NDVI THRESHOLD METHOD THAT REQUIRES RED BAND FOR NDVI < NDVI @ SOIL EMISSIVITY
    tirs_dir = r'G:\Landsat_Download'
    tirs_bt_dir = r'G:\Landsat_Download'
    ndvi_dir = r'G:\Landsat_Download'
    
    #    red = red_band(red_dir) # NOT USED; ONLY USE IF USING FULL NDVI THRESHOLD METHOD THAT REQUIRES RED BAND FOR NDVI < NDVI @ SOIL EMISSIVITY
    tirs = thermal_band(tirs_dir,'DEFAULT')
    tirs_bt = bt_band(tirs_bt_dir,'DEFAULT')
    ndvi = ndvi_band(ndvi_dir,'DEFAULT')
    
    root_dir = r'G:\Landsat_Download' 
    
    # create folders for LST calculated surfaces
    for sheet_name in range(len(all_sheets_name)):
    
        for sub_folder in os.listdir(root_dir):  
           
           # str is arbitrary to naming convention        
           if all_sheets_name[sheet_name].split('_summer_atmos_func')[0] in sub_folder:
                # create outfolders if not exist... 
               outfolder = os.path.join(root_dir, sub_folder + '\\' + 'image_lst_summer')       
               if not os.path.exists(outfolder):
                   os.makedirs(outfolder)
                   print ('{0} created'.format(outfolder))
                   
               else: print ('{0} exists'.format(outfolder))     
               
           elif all_sheets_name[sheet_name].split('_winter_atmos_func')[0] in sub_folder:
                # create outfolders if not exist... 
               outfolder = os.path.join(root_dir, sub_folder + '\\' + 'image_lst_winter')       
               if not os.path.exists(outfolder):
                   os.makedirs(outfolder)
                   print ('{0} created'.format(outfolder))
                     
               else: print ('{0} exists'.format(outfolder))  
    
    # process for LST
    ##sort dataframes first
    #for sheet in range(len(all_sheets)):
    #     all_sheets[sheet] = all_sheets[sheet].sort_values(by=['SCENE_NAME'])
    #  
               
    for sheet in range(len(all_sheets)):   
         
        #!!! to loop a dataframe and search using conditional(s)... 
    #        for row in all_sheets[sheet].iterrows(): 
    #            row_to_list = row[1].tolist() # at index 1, we get a pandas series we can convert to list to search using conditional
    #            if 'Uruapan' in row_to_list[1]:
    #                print row_to_list
                
        tirs_filtered = []
        tirs_bt_filtered = []
        ndvi_filtered = []
        
        # get folder to save output LST raster    
        for index, row in all_sheets[sheet].iterrows():
            
            for sheet_name in range(len(all_sheets_name)):
                
                if sheet == sheet_name: # match by index
                    
                    for sub_folder in os.listdir(root_dir):  
                        
                         # for saving out into single folder
    #                        # get appropriate sub-sub folder to save output; str arg is arbitrary to naming convention 
    #                        if (all_sheets_name[sheet_name].split('_summer_atmos_func')[0] in sub_folder) or (all_sheets_name[sheet_name].split('_winter_atmos_func')[0] in sub_folder):
    #                            outfolder = os.path.join(root_dir, sub_folder + '\\' + 'image_lst')
       
                        #  for saving out into two folders... 
                        # get appropriate sub-sub folder to save output; str arg is arbitrary to naming convention 
                        if all_sheets_name[sheet_name].split('_summer_atmos_func')[0] in sub_folder:
                            outfolder = os.path.join(root_dir, sub_folder + '\\' + 'image_lst_summer')
                            
                        elif all_sheets_name[sheet_name].split('_winter_atmos_func')[0] in sub_folder:
                            outfolder = os.path.join(root_dir, sub_folder + '\\' + 'image_lst_winter')
                            
        # get scene names and GLUE city names from dataframe and then filter out the tirs, tirs_bt and ndvi lists for matching SCENE_NAMES
        SCENE_NAMES = all_sheets[sheet].iloc[:,0:2].astype(str).values.tolist()
##        print(SCENE_NAMES)
##        print(tirs)
        #!!! WARNING: retrieving more images than entries in pandas dataframe
        # LOGIC FIXED FOR THE ABOVE WARNING
        # WE WILL APPLY SET TO LIMIT IMAGES 
        # modify indicing to locate mosaicked images...
        for a in range(len(SCENE_NAMES)):
##            print(SCENE_NAMES[a])
            for aa in range(len(tirs)):
                if 'T1' or 'T2' in os.path.basename(tirs[aa]):
                    
                    # SAVE THIS FOR MOSAICKED IMAGE WORKFLOW... 
        #                    if SCENE_NAMES[a][:-23] in os.path.basename(tirs[aa][:-31] + tirs[aa][-31:-23]) :

                    if SCENE_NAMES[a][0][:-8] == os.path.basename(tirs[aa][:-8]):
##                        print(SCENE_NAMES[a][0][:-8], os.path.basename(tirs[aa][:-8]))
                        tirs_filtered.append([tirs[aa],SCENE_NAMES[a][1]])
                        
        #                        print SCENE_NAMES[a][:-8]
        #                        print os.path.basename(tirs[aa][:-8]) 
        #                        tirs_filtered.append(os.path.basename(tirs[aa]))
                            
                        # SAVE THIS FOR MOSAICKED IMAGE WORKFLOW...                         
        #                    elif SCENE_NAMES[a][:-23] in os.path.basename(tirs[aa][:-31] + tirs[aa][-31:-22]) :     
                        
                    elif SCENE_NAMES[a][0][:-8] == os.path.basename(tirs[aa][:-7]):
                        tirs_filtered.append([tirs[aa],SCENE_NAMES[a][1]])
        #                        print SCENE_NAMES[a][:-23] 
        #                        print os.path.basename(tirs[aa][:-7])
        #                        tirs_filtered.append(os.path.basename(tirs[aa]))                

##        print(tirs_filtered)         
        for b in range(len(SCENE_NAMES)):
            for bb in range(len(tirs_bt)):
                if 'T1' or 'T2' in os.path.basename(tirs_bt[bb]):
                    
                    # SAVE THIS FOR MOSAICKED IMAGE WORKFLOW... 
    #                    if SCENE_NAMES[b][:-23] in os.path.basename(tirs_bt[bb][:-31] + tirs_bt[bb][-31:-28]):
                    if SCENE_NAMES[b][0][:-8] in os.path.basename(tirs_bt[bb][:-14]) or SCENE_NAMES[b][0][:-8] in os.path.basename(tirs_bt[bb][:-13]): # for band 10 (L8) or band 6 (L5/7)
                        
                        tirs_bt_filtered.append([tirs_bt[bb],SCENE_NAMES[b][1]])
    #                        print os.path.basename(tirs_bt[bb])
    #                        tirs_bt_filtered.append(os.path.basename(tirs_bt[bb]))
                        
                        
    #                    # SAVE THIS FOR MOSAICKED IMAGE WORKFLOW...  
    ##                    elif SCENE_NAMES[b][:-25] in os.path.basename(tirs_bt[bb][:-31] + tirs_bt[bb][-31:-28]):
    #                    elif SCENE_NAMES[b][:-8] in os.path.basename(tirs_bt[bb][:-12]) :
    #                    
    ##                        tirs_bt_filtered.append(tirs_bt[bb])
    #                        print os.path.basename(tirs_bt[bb])
    #                        tirs_bt_filtered.append(os.path.basename(tirs_bt[bb]))
                                        
            
        for c in range(len(SCENE_NAMES)):
            for cc in range(len(ndvi)):
                if 'T1' or 'T2' in os.path.basename(ndvi[cc]):
                    
                    # SAVE THIS FOR MOSAICKED IMAGE WORKFLOW... 
    #                    if SCENE_NAMES[c][:-23] in os.path.basename(ndvi[cc][:-31] + ndvi[cc][-31:-25]):
                    if SCENE_NAMES[c][0][:-8] in os.path.basename(ndvi[cc][:-10]):
                        
                        ndvi_filtered.append([ndvi[cc],SCENE_NAMES[c][1]])
    #                        print os.path.basename(ndvi[cc])      
    #                        ndvi_filtered.append(os.path.basename(ndvi[cc]))
    
    #                    # SAVE THIS FOR MOSAICKED IMAGE WORKFLOW... 
    ##                    elif SCENE_NAMES[c][:-25] in os.path.basename(ndvi[cc][:-31] + ndvi[cc][-31:-25]):
    #                    elif SCENE_NAMES[c][:-25] in os.path.basename(ndvi[cc][:-31] + ndvi[cc][-31:-25]):
    #                        
    ##                        ndvi_filtered.append(ndvi[cc])
    #        #                print os.path.basename(ndvi[cc])      
    #                        ndvi_filtered.append(os.path.basename(ndvi[cc]))
                    
        count = len(all_sheets[sheet].index) 
        print 'count '  + str(count)  
        print len(tirs_filtered)
        print len(tirs_bt_filtered)
        print len(ndvi_filtered)
    
        # ---------------------------------------------------------------------------- #
    
        # add indexer to allow for robust element matching between 2 lists 
        # e.g.,
        # [1,2,2,3,4,5] vs. [1,1,1,1,2,3,4,4,5]
        # in addition to adding the GLUE city name to each list of list (SCENE_NAMES), we
        # also need to use an index to match all images for LST calculation
        # Some GLUE cities are using the same image and because we're using a pandas dataframe
        # as the master list, some images are adjacent by row, whilst others are separated further down 
        # the problem arises when we begin to filter out tirs, tirs bt, and ndvi images
        #
        
        for i in range(len(SCENE_NAMES)):
            SCENE_NAMES[i].insert(0,i)
        
        for i in range(len(tirs_filtered)):
            tirs_filtered[i].insert(0,i)
        
        for i in range(len(tirs_bt_filtered)):
            tirs_bt_filtered[i].insert(0,i)
            
        for i in range(len(ndvi_filtered)):
            ndvi_filtered[i].insert(0,i)
        
        
        # remove redundant low/high gain thermal band for Landsat 7 images for processing (e.g., tirs and tirs_bt list length must be same) 
        for i in tirs_filtered:
            if os.path.basename(i[1]) == ('LE07_L1TP_020033_20030115_20160927_01_T1_b62.tif'):
                tirs_filtered.remove(i) 
            elif os.path.basename(i[1]) == ('LE07_L1TP_009056_20161121_20170112_01_T1_b61.tif'):
                tirs_filtered.remove(i) 
        
        
        # --- START FILTER OUT REDUNDANT IMAGES VIA INDEX (INT), IMAGE (STR), GLUE CITY (STR) ---#
        for i in range(len(SCENE_NAMES)):
            for i_index in range(len(SCENE_NAMES[i])):
                
                for j in range(len(tirs_filtered)):
                    for j_index in range(len(tirs_filtered[j])):
        
                        if SCENE_NAMES[i][0] == tirs_filtered[j][0] and \
                           SCENE_NAMES[i][1][:-8] == tirs_filtered[j][1][-48:-8]:
                            pass          
                        
                        # if index not the same, but image is same, and city is same, mod index... 
                        elif tirs_filtered[j][0] != SCENE_NAMES[i][0] and \
                             SCENE_NAMES[i][1][:-8] == tirs_filtered[j][1][-48:-8] and \
                             SCENE_NAMES[i][2] == tirs_filtered[j][2]: 
                            
                            # if adjacent element index already re-indexed, do not overwrite... via match by GLUE city name and corresponding landsat image 
                            tirs_filtered[j][0] = SCENE_NAMES[i][0]
        
                            print SCENE_NAMES[i][0], tirs_filtered[j][0]
        
        for i in range(len(SCENE_NAMES)):
            for i_index in range(len(SCENE_NAMES[i])):
                
                for j in range(len(tirs_bt_filtered)):
                    for j_index in range(len(tirs_bt_filtered[j])):
        
                        if SCENE_NAMES[i][0] == tirs_bt_filtered[j][0] and \
                           SCENE_NAMES[i][1][:-8] == tirs_bt_filtered[j][1][-54:-14]:
                            pass          
                        
                        # if index not the same, but image is same, and city is same, mod index... 
                        # factor for band 10 (L8) or band 6 (L5/7)
                        elif tirs_bt_filtered[j][0] != SCENE_NAMES[i][0] and \
                             (SCENE_NAMES[i][1][:-8] == tirs_bt_filtered[j][1][-54:-14] or SCENE_NAMES[i][1][:-8] == tirs_bt_filtered[j][1][-54:-13]) and \
                             SCENE_NAMES[i][2] == tirs_bt_filtered[j][2]: 
                            
                            # if adjacent element index already re-indexed, do not overwrite... via match by GLUE city name and corresponding landsat image 
                            tirs_bt_filtered[j][0] = SCENE_NAMES[i][0]
        
                            print SCENE_NAMES[i][0], tirs_bt_filtered[j][0]
                            
        for i in range(len(SCENE_NAMES)):
            for i_index in range(len(SCENE_NAMES[i])):
                
                for j in range(len(ndvi_filtered)):
                    for j_index in range(len(ndvi_filtered[j])):
        
                        if SCENE_NAMES[i][0] == ndvi_filtered[j][0] and \
                           SCENE_NAMES[i][1][:-8] == ndvi_filtered[j][1][-50:-10]:
                            pass          
                        
                        # if index not the same, but image is same, and city is same, mod index... 
                        elif ndvi_filtered[j][0] != SCENE_NAMES[i][0] and \
                             SCENE_NAMES[i][1][:-8] == ndvi_filtered[j][1][-50:-10] and \
                             SCENE_NAMES[i][2] == ndvi_filtered[j][2]: 
                            
                            # if adjacent element index already re-indexed, do not overwrite...  via match by GLUE city name and corresponding landsat image 
                            ndvi_filtered[j][0] = SCENE_NAMES[i][0]
        
                            print SCENE_NAMES[i][0], ndvi_filtered[j][0]
                            
       # --- END FILTER OUT REDUNDANT IMAGES VIA INDEX (INT), IMAGE (STR), GLUE CITY (STR) ---#  
       
    # checking test 
    #for i in range(len(SCENE_NAMES)):
    #    print ('{0}, {1}'.format(SCENE_NAMES[i][0],os.path.basename(SCENE_NAMES[i][1])))
    #
    #for i in range(len(tirs_filtered)):   
    #    print ('{0}, {1}'.format(tirs_filtered[i][0],os.path.basename(tirs_filtered[i][1])))
    #    
    #for i in range(len(SCENE_NAMES)):
    #    print SCENE_NAMES[i][1]
        
    # ---------------------------------------------------------------------------- #
    
    # from tirs_filtered, we need to refine it, such that the logic will pull all matches, but we only want 1 for each GLUE city
        
    
    #from collections import OrderedDict 
    #
    ## need to convert the inner lists to tuples so they are hashable
    #tirs_filtered = list(OrderedDict.fromkeys(map(tuple,tirs_filtered)))  #https://stackoverflow.com/questions/7961363/removing-duplicates-in-lists
    #tirs_bt_filtered = list(OrderedDict.fromkeys(tirs_bt_filtered))
    #ndvi_filtered = list(OrderedDict.fromkeys(ndvi_filtered))
    #
    #sorted(set(tirs_filtered), key=lambda x: tirs_filtered.index(x))
    #
    #
    #tirs_filtered_set = set(map(tuple,tirs_filtered)) #need to convert the inner lists to tuples so they are hashable
    
     # --------------- REMOVE DUPLICATES FOR PROCESSING START -----------------#  
        ## REMOVE DUPLICATE ELEMENTS       
        tirs_processed = []
        tirs_bt_processed = []
        ndvi_processed = []
        
        for i in range(len(tirs_filtered)):  
            # if not the same as previous... append
            if tirs_filtered[i][0] != tirs_filtered[i-1][0]:
                tirs_processed.append(tirs_filtered[i])
                
        #    # don't forget about the last element...
        #    elif tirs_filtered[-1] != tirs_filtered[i-1][0]:
        #        tirs_processed.append(tirs_filtered[i])
            
            
        for i in range(len(tirs_bt_filtered)):
            # if not the same as previous... append
            if tirs_bt_filtered[i][0] != tirs_bt_filtered[i-1][0]:
                tirs_bt_processed.append(tirs_bt_filtered[i])
                print tirs_bt_filtered[i][0]
                
        #    # don't forget about the last element...
        #    elif tirs_bt_filtered[-1] != tirs_bt_filtered[i-1][0]:
        #        tirs_bt_processed.append(tirs_bt_filtered[i])
                
                
        for i in range(len(ndvi_filtered)):
            # if not the same as previous... append
            if ndvi_filtered[i][0] != ndvi_filtered[i-1][0]:
                ndvi_processed.append(ndvi_filtered[i])
                
        #    # don't forget about the last element...   
        #    elif ndvi_filtered[-1] != ndvi_filtered[i-1][0]:
        #        ndvi_processed.append(ndvi_filtered[i])
                
        
##        for i in range(len(tirs_bt_processed)):   
##            print ('{0}, {1}, {2}'.format(tirs_bt_processed[i][0],os.path.basename(tirs_bt_processed[i][1]), tirs_bt_processed[i][2]))
            
        
        print len(tirs_processed)
        print len(tirs_bt_processed)
        print len(ndvi_processed)
        
        # -------------- REMOVE DUPLICATES FOR PROCESSING  END ----------------#   
            
            ###  match image with image name from pandas, then pull t, lu, ld as input into LST function          
##        print(len(tirs_processed))
##        print(tirs_processed)
        for i in range(count):
##            print(i)
            ### get ML and AL from metadata
            # Landsat 8
            if os.path.basename(tirs_processed[i][1])[:4] == 'LC08':
                ML = 3.3420E-04 # RADIANCE_MULT_BAND_10 
                AL = 0.10000 # RADIANCE_ADD_BAND_10   
                Eff_WL = 10.896 # Effective Wavelength deprecated; no longer used 
                B_GAMMA_IN = 1324 # Gamma constant
               
            # Landsat 7     
            elif os.path.basename(tirs_processed[i][1])[:4] == 'LE07':
               
                Eff_WL = 11.269  # Effective Wavelength deprecated; no longer used 
                B_GAMMA_IN = 1277 # Gamma constant
               
                # FILE_NAME_BAND_6_VCID_1 (band 61 or 6L) provides an expanded dynamic range and lower radiometric resolution (sensitivity), 
                # with less saturation at high Digital Number (DN) values. (e.g., deserts, snow) 
                if 'LE07_L1TP_020033_20030115_20160927_01_T1' in os.path.basename(tirs_processed[i][1]): # for Cincinnati; snow covered
                    ML = 6.7087E-02 # RADIANCE_MULT_BAND_6_VCID_1
                    AL = -0.06709 # RADIANCE_ADD_BAND_6_VCID_1 
                   
                # FILE_NAME_BAND_6_VCID_2 (band 62 or 6H) has higher radiometric resolution (sensitivity), 
                # although it has a more restricted dynamic range. (e.g., vegetated areas) 
                elif 'LE07_L1TP_009056_20161121_20170112_01_T1' in os.path.basename(tirs_processed[i][1]):  # for Medellin; vegetated covered
                    ML = 3.7205E-02 # RADIANCE_MULT_BAND_6_VCID_2
                    AL = 3.16280 # RADIANCE_ADD_BAND_6_VCID_2 
                   
            # Landsat 5         
            elif os.path.basename(tirs_processed[i][1])[:4] == 'LT05':
                ML = 5.5375E-02 # RADIANCE_MULT_BAND_6 
                AL = 1.18243 # RADIANCE_ADD_BAND_6   
                Eff_WL = 11.457  # Effective Wavelength deprecated; no longer used 
                B_GAMMA_IN = 1256 # Gamma constant
               
            ### get red, tirs, tirs bt, ndvi  -> all raster count MUST BE SAME 
            # For Landsat 5/8 
         #       if os.path.basename(tirs[i])[:-8] == os.path.basename(tirs_bt[i])[:-14] == ndvi[i][:-17]: # ndvi[i][:-17] based on old naming convention NDVI_-1_to_1
            if os.path.basename(tirs_processed[i][1])[:-8] == os.path.basename(tirs_bt_processed[i][1])[:-14] == os.path.basename(ndvi_processed[i][1])[:-10]:
         #           print ('{0}, {1}, {2}'.format(os.path.basename(tirs[i])[:-8], os.path.basename(tirs_bt[i])[:-14], ndvi[i][:-17]))
         #           RED_BAND = red_dir + '\\' + red[i] # NOT USED; ONLY USE IF USING FULL NDVI THRESHOLD METHOD THAT REQUIRES RED BAND FOR NDVI < NDVI @ SOIL EMISSIVITY
                TIRS = tirs_processed[i][1]
                BT = tirs_bt_processed[i][1]
                NDVI_IN = ndvi_processed[i][1]
        
            # For Landsat 5/8 
         #       if os.path.basename(tirs[i])[:-8] == os.path.basename(tirs_bt[i])[:-14] == ndvi[i][:-17]: # ndvi[i][:-17] based on old naming convention NDVI_-1_to_1
            if os.path.basename(tirs_processed[i][1])[:-7] == os.path.basename(tirs_bt_processed[i][1])[:-13] == os.path.basename(ndvi_processed[i][1])[:-10]:
         #           print ('{0}, {1}, {2}'.format(os.path.basename(tirs[i])[:-8], os.path.basename(tirs_bt[i])[:-14], ndvi[i][:-17]))
         #           RED_BAND = red_dir + '\\' + red[i] # NOT USED; ONLY USE IF USING FULL NDVI THRESHOLD METHOD THAT REQUIRES RED BAND FOR NDVI < NDVI @ SOIL EMISSIVITY
                TIRS = tirs_processed[i][1]
                BT = tirs_bt_processed[i][1]
                NDVI_IN = ndvi_processed[i][1]
               
            # For Landsat 7; has additional char length for i in tirs_bt (i.e., 61 or 62)
         #       elif os.path.basename(tirs[i])[:-8] == os.path.basename(tirs_bt[i])[:-13] == ndvi[i][:-17]:  # ndvi[i][:-17] based on old naming convention NDVI_-1_to_1 
            elif os.path.basename(tirs_processed[i][1])[:-8] == os.path.basename(tirs_bt_processed[i][1])[:-13] == os.path.basename(ndvi_processed[i][1])[:-10]: 
                TIRS = tirs_processed[i][1]
                BT = tirs_bt_processed[i][1]
                NDVI_IN = ndvi_processed[i][1]
               
    #        print i        
    #        print os.path.basename(TIRS)
    #        print os.path.basename(BT)
    #        print os.path.basename(NDVI_IN)
    #        print ML, AL, Eff_WL                          
    #        print '.......'
                                           
        #### get atmospheric functions;  -> all raster count MUST BE SAME 
     #       if os.path.basename(tirs[i])[:-8] in str(df['Scene_Name'].str.slice(0,-8)[i]):
     #       if os.path.basename(tirs[i])[:-8] in str(df['SCENE_NAME'].str.slice(0,-8)[i]): # was working now it's not...
            if os.path.basename(tirs_processed[i][1])[:-8] in str(all_sheets[sheet].loc[all_sheets[sheet].index[i], 'SCENE_NAME']):
               
                ATMOS_TRANS = all_sheets[sheet].loc[all_sheets[sheet].index[i], 't']
                UP_RAD = all_sheets[sheet].loc[all_sheets[sheet].index[i], 'Lu']
                D_RAD = all_sheets[sheet].loc[all_sheets[sheet].index[i], 'Ld']     
               
     #               print '{0}, {1}'.format(landsat[i], str(df['Scene_Name'].str.slice(0,-8)[i]))
     #               ATMOS_TRANS = df['t'][i]
     #               UP_RAD = df['Lu'][i]
     #               D_RAD = df['Ld'][i]
     #            
     #               ATMOS_TRANS = df['t_station'][i]  # broken? 
     #               UP_RAD = df['Lu_station'][i]    # broken?
     #               D_RAD = df['Ld_station'][i]     # broken?
     #               print ATMOS_TRANS, UP_RAD, D_RAD
     #                  
     #               ATMOS_TRANS = df.loc[df.index[i], 't_transect']
     #               UP_RAD = df.loc[df.index[i], 'Lu_transect']
     #               D_RAD = df.loc[df.index[i], 'Ld_transect']  
     #                
     #               ATMOS_TRANS = df.loc[df.index[i], 't_station']
     #               UP_RAD = df.loc[df.index[i], 'Lu_station']
     #               D_RAD = df.loc[df.index[i], 'Ld_station']                        
     #               OUTPATH = r'D:\GLUE Results\validate_elevation_differences_lst\station_elev_lst_gamma_delta_non_simplified'
     #               OUTPATH = r'D:\GLUE Results\validate_elevation_differences_lst\transect_elev_lst_gamma_delta_non_simplified'  
     #                
     #               OUTPATH = r'D:\GLUE Results\validate_elevation_differences_lst\station_elev_lst_gamma_delta_simplified'
     #               OUTPATH = r'D:\GLUE Results\validate_elevation_differences_lst\transect_elev_lst_gamma_delta_simplified'  
        
                OUTNAME_CITY = all_sheets[sheet].loc[all_sheets[sheet].index[i], 'City']        
                OUTNAME = os.path.basename(tirs_processed[i][1])[:-8] + '_' + OUTNAME_CITY +'_LST_Celsius'
                OUTPATH = outfolder
                
                print '...'
                print os.path.basename(TIRS)
                print os.path.basename(BT)
                print os.path.basename(NDVI_IN)
                print OUTNAME
                print outfolder 
                print ML, AL, B_GAMMA_IN  
                print ('t: {0}, Lu: {1}, Ld:{2}'.format(ATMOS_TRANS, UP_RAD, D_RAD))
                print '...'
                    
                if os.path.exists(os.path.join(OUTPATH,OUTNAME + '.tif')):
                    print ('{0} already processed!'.format(OUTNAME)) 
                else:
                    # exe LST calc
                    LST(TIRS, BT, B_GAMMA_IN, NDVI_IN, ML, AL, Eff_WL, ATMOS_TRANS, UP_RAD, D_RAD, OUTPATH, OUTNAME)      
                

def get_ML(directory):
    '''
    Description: get list of RADIANCE_MULT_BAND_10 values for Landsat 8 OLI band 10 TIRS images. 
                 All the RADIANCE_MULT_BAND_10 values should be the same (calibrated sensor band-specific multiplicative rescaling factor) 
     
    Args:
        directory (str): directory of Landsat 8 images 

    Returns:
        list of RADIANCE_MULT_BAND_10 value(s) 
    '''
    import os
    
    ML_list = []
    
    for root, dirnames, filenames in os.walk(directory):    
        for file in range(len(filenames)):
            
            # Landsat 8 TIRS band 10 
            if filenames[file].startswith('LC08') and filenames[file].endswith('MTL.txt'):               
               with open(root + '\\' + filenames[file], 'r') as infile:
                    for line in infile:                      
                        # Landsat 8 
                        if 'RADIANCE_MULT_BAND_10' in line:
                            print line
                            ML_list.append(line.lstrip(' ').strip('\n'))
                            
            # Landsat 7 TIRS band 6                
            if filenames[file].startswith('LE07') and filenames[file].endswith('MTL.txt'):               
               with open(root + '\\' + filenames[file], 'r') as infile:
                    for line in infile:                                    
                        # Landsat 7 low-gain
                        if 'RADIANCE_MULT_BAND_6_VCID_1' in line:
                            print line
                            ML_list.append(line.lstrip(' ').strip('\n'))  
                             
                        # Landsat 7 high-gain
                        elif 'RADIANCE_MULT_BAND_6_VCID_2' in line:   
                            print line
                            ML_list.append(line.lstrip(' ').strip('\n'))    
                        
            # Landsat 5 TIRS band 6                
            if filenames[file].startswith('LT05') and filenames[file].endswith('MTL.txt'):               
               with open(root + '\\' + filenames[file], 'r') as infile:
                    for line in infile:                                    
                        # Landsat 5
                        if 'RADIANCE_MULT_BAND_6' in line:
                            print line
                            ML_list.append(line.lstrip(' ').strip('\n'))          
                    print filenames[file]
                    
    return ML_list


def get_AL(directory):
    '''
    Description: get list of RADIANCE_ADD_BAND_10 values for Landsat 8 OLI band 10 TIRS images. All the values should be the same.
                 All the RADIANCE_ADD_BAND_10 values should be the same (calibrated sensor band-specific additive rescaling factor) 
    Args:
        directory (str): directory of Landsat 8 images
        
    Returns:
        list of RADIANCE_ADD_BAND_10 value(s)
    '''    
    
    import os
    
    AL_list = []
    
    for root, dirnames, filenames in os.walk(directory):    
        for file in range(len(filenames)):
            
            # Landsat 8 TIRS band 10 
            if filenames[file].startswith('LC08') and filenames[file].endswith('MTL.txt'):               
               with open(root + '\\' + filenames[file], 'r') as infile:
                    for line in infile:                      
                        # Landsat 8 
                        if 'RADIANCE_ADD_BAND_10' in line:
                            print line
                            AL_list.append(line.lstrip(' ').strip('\n'))
                            
            # Landsat 7 TIRS band 6                
            if filenames[file].startswith('LE07') and filenames[file].endswith('MTL.txt'):               
               with open(root + '\\' + filenames[file], 'r') as infile:
                    for line in infile:                                    
                        # Landsat 7 low-gain
                        if 'RADIANCE_ADD_BAND_6_VCID_1' in line:
                            print line
                            AL_list.append(line.lstrip(' ').strip('\n'))  
                             
                        # Landsat 7 high-gain
                        elif 'RADIANCE_ADD_BAND_6_VCID_2' in line:   
                            print line
                            AL_list.append(line.lstrip(' ').strip('\n'))    
                        
            # Landsat 5 TIRS band 6                
            if filenames[file].startswith('LT05') and filenames[file].endswith('MTL.txt'):               
               with open(root + '\\' + filenames[file], 'r') as infile:
                    for line in infile:                                    
                        # Landsat 5
                        if 'RADIANCE_ADD_BAND_6' in line:
                            print line
                            AL_list.append(line.lstrip(' ').strip('\n'))        
                            
    return AL_list


def red_band(directory, images):
    '''
    Description: 
        Get list of processed Landsat 8 OLI red band images for NDVI-based emissivity method (NBEM). \\
        As reported by Sobrino et al. (2004) see [doi:10.1016/j.rse.2004.02.003], the authors \\
        proposed 2 NBEM methods, a full and simplified version. For the full-version, the emissivity \\
        calculation uses an image's red band for the conditional NDVI < NDVI_soil, whereas the simplified \\
        version does not (e.g., uses a constant soil emissivity value for calculation).
        
        For the GLUE project, we will use the simplified version and use a global constant for soil emissivity.
    
    
        **NOT USED CURRENTLY; IF USING FULL NDVI THRESHOLD METHOD, USE THIS
     
    Args:
        directory (str): directory of shapefile(s) (e.g., buffer shapefile)
        images (str): set as 'DEFAULT' to process all images within a directory 
                      set as full name of image to be processed. 
                          e.g., 'LC08_L1TP_018030_20160418_20170223_01_T1_red_band_scaled.tif'

    Returns:
        list of raster(s)
    '''
    import os

    raster = []
    for root, dirnames, filenames in os.walk(directory):
        for file in range(len(filenames)):
            if images == 'DEFAULT':
                if filenames[file].endswith('red_band_scaled.tif'):
                    raster.append(filenames[file])
            else: 
                if filenames[file].endswith(images):
                    raster.append(filenames[file])
                    
    return raster    
    
    
def thermal_band(directory, images):  
    '''
    Description: get list of Landsat 8 OLI thermal band 10 images 
     
    Args:
        directory (str): directory of shapefile(s) (e.g., buffer shapefile)
        images (str): set as 'DEFAULT' to process all images within a directory 
                      set as full name of image to be processed. e.g., 'LC08_L1TP_018030_20160418_20170223_01_T1_b10.tif'

    Returns:
        list of raster(s)
    '''
    import os 

    raster = []
    for root, dirnames, filenames in os.walk(directory):
        for file in range(len(filenames)):
            # ALL USAGE 
            if images == 'DEFAULT':
                
                #Landsat 8
                if filenames[file].startswith('LC08') and filenames[file].endswith('b10.tif'):
                    raster.append(os.path.join(root, filenames[file]))    
                # Landsat 7; low gain
                elif filenames[file].startswith('LE07') and filenames[file].endswith('b61.tif'):
                    raster.append(os.path.join(root, filenames[file]))   
                # Landsat 7; high gain
                elif filenames[file].startswith('LE07') and filenames[file].endswith('b62.tif'):
                    raster.append(os.path.join(root, filenames[file]))   
                # Landsat 5
                elif filenames[file].startswith('LT05') and filenames[file].endswith('b6.tif'):
                    raster.append(os.path.join(root, filenames[file]))                       
                    
            # SINGLE USAGE        
            else: 
                if filenames[file].endswith(images):
                    raster.append(os.path.join(root, filenames[file]))
                    
    return raster    
    

def bt_band(directory, images):
    '''
    Description: get list of Landsat 8 OLI thermal brightness temperature band 10 images 
     
    Args:
        directory (str): directory of shapefile(s) (e.g., buffer shapefile)
        images (str): ALL USAGE: set as 'DEFAULT' to process all images within a directory 
                      SINGLE USAGE: set as full name of image to be processed. e.g., 'LC08_L1TP_018030_20160418_20170223_01_T1_bt_band10.tif'

    Returns:
        list of raster(s)
    '''
    import os

    raster = []
    for root, dirnames, filenames in os.walk(directory):
        for file in range(len(filenames)):
            # ALL USAGE 
            if images == 'DEFAULT': 
                # Landsat 8
                if filenames[file].startswith('LC08') and filenames[file].endswith('bt_band10.tif'):
                    raster.append(os.path.join(root, filenames[file]))  
                # Landsat 5/7
                elif (filenames[file].startswith('LT05') or filenames[file].startswith('LE07')) and filenames[file].endswith('bt_band6.tif'):
                    raster.append(os.path.join(root, filenames[file]))                  
                    
            # SINGLE USAGE         
            else: 
                if filenames[file].endswith(images):
                    raster.append(os.path.join(root, filenames[file]))    
                    
    return raster    
      
        
def ndvi_band(directory, images):
    '''
    Description: get list of processed Landsat 8 OLI NDVI images   
     
    Args:
        directory (str): directory of shapefile(s) (e.g., buffer shapefile)
        images (str): ALL USAGE: set as 'DEFAULT' to process all images within a directory 
                      SINGLE USAGE: set as full name of image to be processed. e.g., 'LC08_L1TP_018030_20160418_20170223_01_T1_NDVI_-1_to_1.tif'

    Returns:
        list of raster(s)
    '''
    import os
    
    raster = []
    for root, dirnames, filenames in os.walk(directory):
        for file in range(len(filenames)):
            # ALL USAGE 
            if images == 'DEFAULT':
#                if filenames[file].endswith('NDVI_-1_to_1.tif'):
                if filenames[file].endswith('NDVI.tif'):                
                    raster.append(os.path.join(root,filenames[file]))              
            # SINGLE USAGE 
            else: 
                if filenames[file].endswith(images):
                    raster.append(os.path.join(root,filenames[file]))  
                    
    return raster    
    
    
def LST(TIRS, BT, B_GAMMA_IN, NDVI_IN, ML, AL, Eff_WL, ATMOS_TRANS, UP_RAD, D_RAD, OUTPATH, OUTNAME):
    '''
    Description: 
        Calculate land surface temperature using NDVI-based emissivity method (NBEM) \\
        single-channel algorithm (TIRS) with Landsat 5/7/8 images.
    
     
    Args:
        TIRS (str): Landsat 5 TM band 6/ Landsat 7 ETM band 6/ Landsat 8 OLI band 10
                        thermal infrared band '.tif' images to be processed. 
                        Refer to thermal_band function and main function
                        
        BT (str): Landsat 5 TM band 6/ Landsat 7 ETM band 6/ Landsat 8 OLI band 10
                      brightness temperature '.tif' images to be processed. 
                      Refer to bt_band function and main function 
                      
        B_GAMMA_IN (int):Landsat 5 TM band 6/ Landsat 7 ETM band 6/ Landsat 8 OLI band 10 
                         thermal infrared band constant 
                         (see "Revision of the Single-Channel Algorithm for Land Surface Temperature Retrieval From Landsat Thermal-Infrared Data" and 
                          "Land Surface Temperature Retrieval Methods From Landsat-8 Thermal Infrared Sensor Data")
                         
        NDVI_IN (str): Landsat 5 TM / Landsat 7 ETM/ Landsat 8 NDVI '.tif' images to be processed. 
                       Refer to ndvi_band function and  main function
                       
        ML (float): RADIANCE_MULT_BAND_6 from image metadata if Landsat 5 TM
                    RADIANCE_MULT_BAND_6_VCID_1 from image metadata if Landsat 7 ETM
                    RADIANCE_MULT_BAND_6_VCID_2 from image metadata if Landsat 7 ETM
                    RADIANCE_MULT_BAND_10 from image metadata if Landsat 8 OLI. 
                    Refer to main function 
        
        AL (float): RADIANCE_ADD_BAND_6 from image metadata if Landsat 5 TM
                    RADIANCE_ADD_BAND_6_VCID_1 from image metadata if Landsat 7 ETM
                    RADIANCE_ADD_BAND_6_VCID_2 from image metadata if Landsat 7 ETM                   
                    RADIANCE_ADD_BAND_10 from image metadata. 
                    Refer to main function
        
        Eff_WL (float): effective wavelength of thermal band. Refer to main function # DEPRECATED         
        ATMOS_TRANS (float): atmospheric transmissivity value calculated from Atmospheric Correction Parameter Calculator (https://atmcorr.gsfc.nasa.gov/). Refer to main function
        UP_RAD (float): upwelling radiance value calculated from Atmospheric Correction Parameter Calculator (https://atmcorr.gsfc.nasa.gov/). Refer to main function
        D_RAD  (float): downwelling radiance  value calculated from Atmospheric Correction Parameter Calculator (https://atmcorr.gsfc.nasa.gov/). Refer to main function
        OUTPATH (str): specify directory or folder for output LST image(s).  Refer to main function
        OUTNAME (str): specify name of output LST image(s). Refer to main function 
        
    Returns:
        no returns or exchanges 
    '''
    import arcpy
    from arcpy.sa import *
    
    arcpy.env.overwriteOutput = True
    #arcpy.env.workspace = r'E:\Scenes_L2\lst_sc_test'
    
    # Same for summer/winter scenes
    RADIANCE_MULT_BAND = ML # RADIANCE_MULT_BAND_10 
    RADIANCE_ADD_BAND = AL # RADIANCE_ADD_BAND_10 
    
    band10_radiance = arcpy.Raster(TIRS)*RADIANCE_MULT_BAND + RADIANCE_ADD_BAND
    band10_bt = arcpy.Raster(BT)*0.1
#    band10_bt_Celsius = arcpy.Raster(BT)*0.1 - 273.15
    
    #-------------------------------------------------------------------------#    
    
    # detect input raster for effective wavelength; if LT5 == LANDSAT 5; if LE07 == LANDSAT 7, if LC08 == LANDSAT 8 
    
    # deprecated in favor of simplified version; "Revision of the Single-Channel Algorithm for Land Surface Temperature Retrieval From Landsat Thermal-Infrared Data"
    # range 10.6 - 11.19 um
    effective_wavelength = Eff_WL
    
    # deprecated in favor of simplified version; "Revision of the Single-Channel Algorithm for Land Surface Temperature Retrieval From Landsat Thermal-Infrared Data"
#    c1 = 1.19104*10**8
#    c2 = 14387.7
    
    ##### B_GAMMA VALUES EXPLAINED
    ## for LS8; "Land Surface Temperature Retrieval Methods From Landsat-8 Thermal Infrared Sensor Data" 
    # 1324 Kelvin for TIRS-1 (band 10 of Landsat 8); equation is  c2/effective_wavelength; technically 14387.7/10.896 = 1320.457 for landsat 8 B10...
    # 1199 Kelvin for TIRS-2 (band 11 of Landsat 8)
    
    ## for LS7; "Revision of the Single-Channel Algorithm for Land Surface Temperature Retrieval From Landsat Thermal-Infrared Data"
    # 1277 Kelvin for TIRS (band 6 of Landsat 7)
    
    ## for LS5; "Revision of the Single-Channel Algorithm for Land Surface Temperature Retrieval From Landsat Thermal-Infrared Data"
    # 1256 Kelvin for TIRS (band 6 of Landsat 7)
    
    B_GAMMA = B_GAMMA_IN # MODIFY ME TO GENERALIZE 
    
    # full expression; 
    # deprecated in favor of simplified expression 
#    gamma = ( (c2*band10_radiance) /  (band10_bt**2) * ( ((effective_wavelength**4/c1)*band10_radiance) + effective_wavelength**-1 ) )**-1
#    delta = (-gamma * band10_radiance) + band10_bt
    
    # simplified expression; "Revision of the Single-Channel Algorithm for Land Surface Temperature Retrieval From Landsat Thermal-Infrared Data"
    gamma = (band10_bt**2)/(B_GAMMA*band10_radiance)
    delta = band10_bt - ((band10_bt**2)/B_GAMMA)
    
    # atmospheric functions 
    atmospheric_transmissivity = ATMOS_TRANS
    upwelling_radiance = UP_RAD
    downwelling_radiance = D_RAD

    psi_1 = 1/atmospheric_transmissivity
    psi_2 = -downwelling_radiance - (upwelling_radiance/atmospheric_transmissivity)
    psi_3 = downwelling_radiance
    
    NDVI = arcpy.Raster(NDVI_IN)
    
    ### MODIFY SUCH THAT PV < 0, IS SET TO 0, PV > 1, IS SET TO 1
    # portion of vegetation
    # use default max of 0.5 and min of 0.2
    #Sobrino and Raissouni (2000) suggested using NDVImin = 0.2 and NDVImax = 0.5 for global-scale remote sensing applications.
    
    ## this interpretation of the method suggests using NDVI map of -1 to 1, computing Pv, and then scaling it.
    Pv = Square((NDVI - 0.2) / (0.5 - 0.2))
#    Pv = Square((NDVI - 0.2) / (0.8 - 0.2))
    
    Pv_scaled = arcpy.sa.Con(Pv <= 0, 0, Con((Pv > 0) & (Pv < 1), Pv, Con(Pv >=1, 1)))
    
    # consistent with Sobrino et al. 2008, Jimenez et al. 2009
#    C_sobrino = (1-0.97) * 0.985 * 0.55 * (1-Pv_scaled)
    
    
    ### FROM: Land surface temperature retrieval from LANDSAT TM 5
    #A possible solution is to use the mean value for the emissivities of soils included in the
    #ASTER spectral library (http://asterweb.jpl.nasa.gov) and filtered according to band TM6 filter function. In this way
    #considering a total of 49 soils spectra, a mean value of 0.973 (with a standard deviation of 0.004) is obtained
#    soil_e_yu = 0.973-(0.047*red_band)
#    veg_e_yu = 0.9863
    
    soil_e_sobrino = 0.97
    veg_e_sobrino = 0.99
    
    water_e = 0.99 # ASTER lib
#    snow_ice_e = 0.991 # ASTER lib
    
    
    ### arcpy.sa.Con WILL NOT calculate from a layer; has to be a arcpy.Raster()    
    # NBEM 
    #Simplified NBEM - consistent with Sobrino et al. 2008
    emissivity_sobrino = arcpy.sa.Con(NDVI < 0, water_e,
                              Con((NDVI > 0) & (NDVI < 0.2), soil_e_sobrino,
                                  Con((NDVI >= 0.2) & (NDVI <= 0.5), (soil_e_sobrino + (veg_e_sobrino-soil_e_sobrino)*Pv_scaled),
                                      Con( NDVI > 0.5, (veg_e_sobrino)))))
    
        
    LST_Celsius_sobrino = ((gamma*((emissivity_sobrino**-1*(psi_1*band10_radiance + psi_2)) + psi_3)) + delta) - 273.15
    
    # save emissivity
#    arcpy.CopyRaster_management(emissivity_sobrino, r'C:\Users\tongale1\Documents\LST_TEST\toronto summer\LC08_L1TP_018030_20180611_20180615_01_T1_emissivity_sobrino.tif')
    # save LST
    arcpy.CopyRaster_management(LST_Celsius_sobrino, os.path.join(OUTPATH,OUTNAME + '.tif'))


main_all_images()

#if __name__ == '__main__':
    
#    main_all_images()
    
    
    
    
