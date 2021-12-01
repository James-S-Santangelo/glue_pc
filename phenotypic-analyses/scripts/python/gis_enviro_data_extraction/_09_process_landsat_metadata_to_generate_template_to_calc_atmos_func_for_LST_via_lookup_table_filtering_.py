# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 15:48:30 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15
"""

try: 
    import os, sys
    
    directory = r'D:\Python scripts - FINAL'
    os.chdir(directory)
    
    from _00_glue_utils import lookup_table
    
except ImportError as IE:
    print (IE)
    print ("These functions requires arcpy to run")  
    sys.exit(1)    
    
    
def landsat_image_date_time(directory):
    '''
    Description:
        process and collate Landsat images for date and time metadata.
        
    Args:
        directory (str): specify directory of Landsat images 
        
    Returns:
        list of lists of image metadata (date and time) for further processing
    '''  
#    import os
    
    image = []
    date = []   
    time = []   
    date_year = []
    date_month = []
    date_date = []
    
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
    
    for root, dirnames, filenames in os.walk(directory):    
        for file in range(len(filenames)):
            if filenames[file].endswith('MTL.txt'):    
#                print filenames[file]
                image.append(filenames[file])
                with open(root + '\\' + filenames[file], 'r') as infile:
                    for line in infile:
                        if 'DATE_ACQUIRED' in line:
#                            print line
                            
                            # clean up text 
                            date_line_strip = line.lstrip(' ').strip('\n')
                            date_line_strip_cleaned = date_line_strip[-10:]
                            
                            date.append(date_line_strip_cleaned)
                            
                        if 'SCENE_CENTER_TIME' in line:    
#                            print line                           
                            # clean up text 
                            time_line_strip = line.lstrip(' ').strip('\n')
                            time_line_strip_cleaned = time_line_strip[21:29]
                            
                            time.append(time_line_strip_cleaned)
                            
                        if 'DATE_ACQUIRED' in line:
#                            print line
                            
                            # clean up text 
                            date_line_strip = line.lstrip(' ').strip('\n')
                            date_line_strip_cleaned = date_line_strip[-10:]                           
                            date_line_strip_year = date_line_strip_cleaned.split('-')[0] #YYYY
                            date_line_strip_month = date_line_strip_cleaned.split('-')[1] #MM
                            date_line_strip_date = date_line_strip_cleaned.split('-')[2] #DD
                                                       
                            if str(date_line_strip_month) == match_month_list[0][0]:
                                date_month.append(match_month_list[0][1])
                                
                            elif str(date_line_strip_month) == match_month_list[1][0]:
                                date_month.append(match_month_list[1][1]) 
                                
                            elif str(date_line_strip_month) == match_month_list[2][0]:
                                date_month.append(match_month_list[2][1])
                                
                            elif str(date_line_strip_month) == match_month_list[3][0]:
                                date_month.append(match_month_list[3][1])
                                
                            elif str(date_line_strip_month) == match_month_list[4][0]:
                                date_month.append(match_month_list[4][1])
                                
                            elif str(date_line_strip_month) == match_month_list[5][0]:
                                date_month.append(match_month_list[5][1]) 
                                
                            elif str(date_line_strip_month) == match_month_list[6][0]:
                                date_month.append(match_month_list[6][1]) 
                                
                            elif str(date_line_strip_month) == match_month_list[7][0]:
                                date_month.append(match_month_list[7][1])
                                
                            elif str(date_line_strip_month) == match_month_list[8][0]:
                                date_month.append(match_month_list[8][1])
                                
                            elif str(date_line_strip_month) == match_month_list[9][0]:
                                date_month.append(match_month_list[9][1])  
                                
                            elif str(date_line_strip_month) == match_month_list[10][0]:
                                date_month.append(match_month_list[10][1])    
                                
                            elif str(date_line_strip_month) == match_month_list[11][0]:
                                date_month.append(match_month_list[11][1])
                                
                            date_year.append(date_line_strip_year)
                            date_date.append(date_line_strip_date)
                                                                                  
    image = [image]
    image.append(date)
    image.append(time)
    image.append(date_year)
    image.append(date_month)
    image.append(date_date)
    
    return image 

#--------------------------------------end-------------------------------------#
    
#directory = r'G:\Landsat_Download\canada\espa-alexander.tong@mail.utoronto.ca-10292018-105600-988' # Canada Summer    
#directory = r'G:\Landsat_Download\asia\espa-alexander.tong@mail.utoronto.ca-10292018-105823-178'
#images = landsat_image_date_time(directory) 
#images_tuple = zip(*images)
# 


#--------------------------------------end-------------------------------------#

def landsat_date_time_to_csv(list_of_list, cities, out_dir, out_name):
    '''
    Description:
        process and collate Landsat images for date and time metadata such that \\
        they are organized in an easy-to-interpret manner in output csv files \\
        for collecting meteorological data for each city for land surface temperature retrieval \\
        using Landsat datasets. 
        
        The Caveat to this function is that it will unnecessarily associate \\
        additional images for a city even though a single image was acquired for  \\
        it due to the list of list logic being applied  \\
        (i.e., there is no way around this without hard-coding the image to be associated with each city)
        
        e.g., output format:
            
            SCENE_NAME	DATE_ACQUIRED (YYYY-MM-DD)	SCENE_CENTER_TIME (GMT)	YEAR	MONTH	DATE	Country	City	Row_Path
        0	LC08_L1TP_042025_20180705_20180717_01_T1_MTL.txt	7/5/2018	18:29:13	2018	Jul	5	Canada_AB	Calgary	42025
        1	LC08_L1TP_043023_20180728_20180813_01_T1_MTL.txt	7/28/2018	18:34:46	2018	Jul	28	Canada_AB	Edmonton	43023
        2	LC08_L1TP_043023_20180728_20180813_01_T1_MTL.txt	7/28/2018	18:34:46	2018	Jul	28	Canada_AB	St. Albert	43023
        3	LC08_L1TP_047026_20170822_20170912_01_T1_MTL.txt	8/22/2017	19:01:30	2017	Aug	22	Canada_BC	Vancouver	47026
        4	LC08_L1TP_047026_20170822_20170912_01_T1_MTL.txt	8/22/2017	19:01:30	2017	Aug	22	Canada_BC	Victoria	47026
        5	LC08_L1TP_033024_20180706_20180717_01_T1_MTL.txt	7/6/2018	17:33:13	2018	Jul	6	Canada_MB	Dauphin	33024

    Args:
        list_of_list (list of lists): generated from landsat_image_date_time() 
        cities (list of lists): specify list of list from cities()
        out_dir (str): specify output directory/folder for csv 
        out_name (str): specify output name for csv. 
        
    Returns:
        
    '''
#    import os
#    import pandas as pd
     
    # take list of list and convert to dataframe for processing 
    # https://datascience.stackexchange.com/questions/26333/convert-a-list-of-lists-into-a-pandas-dataframe?rq=1
    images_tuple = zip(*list_of_list)
    LS_headers = ['SCENE_NAME','DATE_ACQUIRED (YYYY-MM-DD)','SCENE_CENTER_TIME (GMT)', 'YEAR', 'MONTH', 'DATE']
    df = pd.DataFrame(images_tuple, columns=LS_headers)
    
    # same as above
    #df = pd.DataFrame(image)
    #df = df.transpose()
    #df.columns = ['SCENE_NAME','DATE_ACQUIRED','SCENE_CENTER_TIME (GMT)']

    # print each row
    for i in range(len(df.columns)):
        print df.iloc[:,i]
    
#    filepath = out_dir
#    filename = out_name
    
#    # output excel
#    writer = pd.ExcelWriter(os.path.join(filepath, filename + '.xlsx')) 
#    df.to_excel(writer,'Sheet1') 
#    writer.save()
#    
#    # output csv
#    df.to_csv(os.path.join(filepath, filename + '.csv'), sep=',')
    
    # convert dataframe (landsat name, date, time) to list for processing 
    df_to_list_for_processing = df.values.tolist()
    
    # IF MATCH add landsat name, date, time to cities master list by matching row/path 
    for i in range(len(cities)):   
        for j in range(len(df_to_list_for_processing)):
            
            a = df_to_list_for_processing[j][0][:10] # slice front for Landsat row + path 
            b = df_to_list_for_processing[j][0][16:] # slice back  for Landsat row + path 
            
            if cities[i][2] in df_to_list_for_processing[j][0].split(a)[-1].split(b)[0]:
                cities[i].extend(df_to_list_for_processing[j])
    
    #------------------------------#
    # Below is for testing 
    #df_cities = pd.DataFrame(cities)
    #for i in range(len(df_cities.columns)):
    #    print df_cities.iloc[:,i].head()
    #    
    #count = 0
    #for i in range(len(cities)):
    #    if len(cities[i]) > 3:
    #        print cities[i][1]
    #        count += 1
    #------------------------------#
    
    # select only matching landsat name, date, time to cities master list
    matching_city_code_to_landsat_images = []
    for i in range(len(cities)):
        if len(cities[i]) > 3: # default len of 3 for cities master list
            matching_city_code_to_landsat_images.append(cities[i])
    
    # convert processed cities master list  with matches back to dataframe
    df_cities_select = pd.DataFrame(matching_city_code_to_landsat_images)
    
    
    # get headers for indexing 
    df_cities_select_index = df_cities_select.columns.tolist()
    
    # cities master list headers 
    cities_headers = ['Country','City','Row_Path']
    
    # add headers according to length of repeat LS images for each city (if applicable as most only have 1 image per city)
    length_of_LS_headers = len(df_cities_select_index[3:])/6
    total_LS_headers = LS_headers*length_of_LS_headers
    
    final_output_headers = total_LS_headers + cities_headers
    
    #re-arrange columns by indices 
    df_cities_select_index_rearranged = df_cities_select_index[3:] + df_cities_select_index[:3]
    
    # re-initialize dataframe for output 
    df_cities_select = df_cities_select[df_cities_select_index_rearranged]
    
    # add headers 
    df_cities_select.columns = final_output_headers
    
    filepath = out_dir
    filename = out_name
    
    # output csv  
    df_cities_select.to_csv(os.path.join(filepath, filename + '.csv'), sep=',')
    
    # output excel # CURRENTLY ISSUE WITH XLSX OUTPUT; ERROR RETURNED, FIX LATER.... 
#    writer = pd.ExcelWriter(os.path.join(filepath, filename + '.xlsx')) 
#    df_cities_select.to_excel(writer,'Sheet1') 
#    writer.save()
    
       
#--------------------------------------end-------------------------------------#


#landsat_date_time_to_csv(images,out_dir, out_name)   
    

def batch_landsat_date_time():
    '''
    '''
#    import os 
    #hard-code at the moment...
#    directory = [[r'G:\Landsat_Download\asia\espa-alexander.tong@mail.utoronto.ca-10292018-105823-178', 'summer'],
#                 [r'G:\Landsat_Download\asia\espa-alexander.tong@mail.utoronto.ca-10292018-105905-194', 'winter'],
#                 
#                 [r'G:\Landsat_Download\canada\espa-alexander.tong@mail.utoronto.ca-10292018-105600-988', 'summer'],
#                 [r'G:\Landsat_Download\canada\espa-alexander.tong@mail.utoronto.ca-10292018-105724-706', 'winter'],
#                 [r'G:\Landsat_Download\canada\espa-alexander.tong@mail.utoronto.ca-11092018-101009-387', 'winter_addendum_1'], # edmonton 
#                 
#                 [r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-10292018-110638-866', 'summer'],
#                 [r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-10292018-113302-004', 'winter'],
#                 [r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-11072018-104859-769', 'winter_addendum_1'], # stockholm
#                 [r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-11072018-122905-326', 'winter_addendum_2'], # helsinki, malmo, trondheim, glasgow 
#                 
#                 [r'G:\Landsat_Download\oceania\espa-alexander.tong@mail.utoronto.ca-10292018-113819-365', 'summer'],
#                 [r'G:\Landsat_Download\oceania\espa-alexander.tong@mail.utoronto.ca-10292018-114004-367', 'winter'],
#                 
#                 [r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-10292018-120915-536', 'summer'],
#                 [r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-11192018-131412-079', 'summer_addendum'],                 
#                 [r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-10292018-120942-303', 'winter'],
#                 [r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-11192018-131512-240', 'winter_addendum'],
#                 
#                 [r'G:\Landsat_Download\usa\espa-alexander.tong@mail.utoronto.ca-10292018-113443-256', 'summer'],
#                 [r'G:\Landsat_Download\usa\espa-alexander.tong@mail.utoronto.ca-10292018-113734-294', 'winter']]
 
#    directory = [[r'G:\Landsat_Download\oceania\espa-alexander.tong@mail.utoronto.ca-03042019-161727-973', 'summer_addendum'],
#                [r'G:\Landsat_Download\oceania\espa-alexander.tong@mail.utoronto.ca-03042019-161752-773', 'winter_addendum'],
#                [r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-03042019-161814-357', 'summer_addendum_2'],
#                [r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-03042019-161848-728', 'winter_addendum_2'],]
    
#    
#    directory = [[r'G:\Landsat_Download\oceania\espa-alexander.tong@mail.utoronto.ca-10292018-113819-365', 'summer'],
#                 [r'G:\Landsat_Download\oceania\espa-alexander.tong@mail.utoronto.ca-10292018-114004-367', 'winter'],]
#    
    # Cape town 
##    directory = [[r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-03192019-134034-964','summer_addendum_1'],
##                 [r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-03192019-134103-601', 'winter_addendum_3']]

    # Fairbanks
    directory = [[r'G:\Landsat_Download\usa\espa-james.santangelo37@gmail.com-07232020-082513-144','summer_addendum_1'],
                 [r'G:\Landsat_Download\usa\espa-james.santangelo37@gmail.com-07232020-082414-856', 'winter_addendum_1']]
    
    for i in range(len(directory)):
        out_dir = os.path.dirname(directory[i][0])
        out_name = os.path.dirname(directory[i][0])[20:] + '_' + directory[i][1]
        cities_for_hire = lookup_table(r'D:\Python scripts - FINAL')  
        print out_dir
        print out_name 
        
        images = landsat_image_date_time(directory[i][0])
        landsat_date_time_to_csv(images, cities_for_hire, out_dir, out_name)   
        cities_for_hire = []

if __name__ == '__main__':
    import pandas as pd
    import os
    batch_landsat_date_time()

#--------------------------------------end-------------------------------------#

