# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 14:57:52 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15

web form data scraper to automate input and output using https://atmcorr.gsfc.nasa.gov/  

https://realpython.com/modern-web-automation-with-python-and-selenium/
https://selenium-python.readthedocs.io/getting-started.html

1. selenium api; chrome/firefox webdriver 
2. mechanize
3. beautifulsoup
4. urllib... 

> infilling the fields... 
Locators
> target xml, css path locators  

row['SCENE_NAME']      
row['YEAR']
row['MONTH']
row['DATE']
int(str(row['SCENE_CENTER_TIME (GMT)']).split(':')[0]) # hour
int(str(row['SCENE_CENTER_TIME (GMT)']).split(':')[1])  # minute
str(row['Station Latitude (DD)'])
str(row['Station Longitude (DD)'])
row['Station Elevation (km)']
row['Station Pressure (mb)']
row[u'Station Temperature (°C)']
row['Station Relative Humidity (%)'] 
row['t']
row['Lu']
row['Ld']

"""

import os 
import time
import pandas as pd
import numpy as np
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.common.exceptions import TimeoutException
from selenium.common.exceptions import NoSuchElementException
from selenium.common.exceptions import InvalidSessionIdException
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By

# set directory to geckodriver.exe
os.chdir(r'D:\geckodriver-v0.23.0-win64')
driver = webdriver.Firefox()

#driver.get("http://www.python.org")
#assert "Python" in driver.title
#elem = driver.find_element_by_name("q")
#elem.clear()
#elem.send_keys("pycon")
#elem.send_keys(Keys.RETURN)
#assert "No results found." not in driver.page_source

# get url
driver.get('https://atmcorr.gsfc.nasa.gov/')
assert "Atmospheric Correction Parameter Calculator" in driver.title # add error and then driver.close() 

#-------------------------------------Start------------------------------------#

def read_excel(excel_path, excel_file, excel_sheetnames, outfolder_name, return_df_or_outfolder_name):
    '''
    Description:
        pass
    Args:
        excel_path (str):
        excel_file (str):
        excel_sheetnames (list of 'str'):
        outfolder_name (str):
        return_df_or_outfolder_name (str): 'df' for list of df OR 'outfolder' for outfolder name 
    Returns:
        pass
    '''
    # sheetname becomes df name... 
#    excel_sheetnames_list = excel_sheetnames
    
    df_names_list = []
    df_all_sheets = []
    df_all_sheetname = []
    
    for i in excel_sheetnames:
        df = 'df_' + i
        df_names_list.append(df)
        df_all_sheetname.append(df.lower() + outfolder_name)
        
    for j in df_names_list:
        j = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = j[3:], index = True)
        df_all_sheets.append(j)
    
    if return_df_or_outfolder_name == 'df':
        return df_all_sheets 
    elif return_df_or_outfolder_name == 'outfolder':
        return df_all_sheetname
#-------------------------------------End--------------------------------------#
        
#-----------------------------------test start---------------------------------#
excel_sheetnames = ['Asia_summer','Asia_winter']
excel_path = r'D:\\'
excel_file = 'GLUE_City_LST_Retrieval.xlsx'  
outfolder_name = '_atmos_func'
return_df_or_outfolder_name = 'df'
all_sheets = read_excel(excel_path, excel_file, excel_sheetnames,outfolder_name,return_df_or_outfolder_name)

excel_sheetnames = ['Asia_summer','Asia_winter']
excel_path = r'D:\\'
excel_file = 'GLUE_City_LST_Retrieval.xlsx'  
outfolder_name = '_atmos_func'
return_df_or_outfolder_name = 'outfolder' 
all_sheets_name = read_excel(excel_path, excel_file, excel_sheetnames,outfolder_name,return_df_or_outfolder_name)
#------------------------------------test end----------------------------------#


#-----------------------------------------------------------------------------#
'''
!!!WARNING BEGINNING!!!

# USE ME TO GENERATE MISSING ATMOSPHERIC FUNCTIONS
'''
#-----------------------------------------------------------------------------#
excel_path = r'D:\\'
excel_file = r'GLUE_City_LST_retrieval.xlsx'

# to prevent any unnecessary issues, we will parse each dataframe independently 
# and return each filled as individual; re-combine csvs as tabs in single excel again ( combine_atmos_func() )

# step 1. filter out columns using pandas instead of hard-coded index by cols using parse_cols
df_asia_summer = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'Asia_summer', index = True)
df_asia_winter = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'Asia_winter', index = True)

df_canada_summer = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'Canada_summer', index = True)
df_canada_winter = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'Canada_winter', index = True)

df_europe_et_al_summer = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'Europe_et_al_summer', index = True)
df_europe_et_al_winter = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'Europe_et_al_winter', index = True)

df_oceania_summer = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'Oceania_summer', index = True)
df_oceania_winter = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'Oceania_winter', index = True)

df_south_america_summer = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'South_America_summer', index = True)
df_south_america_winter = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'South_America_winter', index = True)

df_usa_summer = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'USA_summer', index = True)
df_usa_winter = pd.read_excel(os.path.join(excel_path, excel_file), sheetname = 'USA_winter', index = True)



# all data with summer + winter combined
all_sheets = [df_asia_summer,df_asia_winter,
              df_canada_summer,df_canada_winter,
              df_europe_et_al_summer,df_europe_et_al_winter,
              df_oceania_summer,df_oceania_winter,
              df_south_america_summer,df_south_america_winter,
              df_usa_summer,df_usa_winter]

all_sheets_name = ['asia_summer_atmos_func','asia_winter_atmos_func',
                    'canada_summer_atmos_func', 'canada_winter_atmos_func',
                    'europe_et_al_summer_atmos_func', 'europe_et_al_winter_atmos_func',
                    'oceania_summer_atmos_func', 'oceania_winter_atmos_func',
                    'south_america_summer_atmos_func', 'south_america_winter_atmos_func',
                    'usa_summer_atmos_func', 'usa_winter_atmos_func']


all_sheets = [df_asia_summer,df_asia_winter,
              df_canada_summer, df_canada_winter]

all_sheets_name = ['asia_summer_atmos_func','asia_winter_atmos_func',
                    'canada_summer_atmos_func', 'canada_winter_atmos_func',]


all_sheets = [df_south_america_winter,]

all_sheets_name = ['south_america_winter_atmos_func_equatorial_fix',]


all_sheets = [df_canada_summer,]

all_sheets_name = ['canada_summer_atmos_func',]
#-----------------------------------------------------------------------------#
'''!!!WARNING END!!!'''
#-----------------------------------------------------------------------------#


Atmos_Calc_folder_root = r'D:\GLUE Results\Process_Atmospheric_Functions'


'''
    Description:
        Atmospheric Correction Parameter Calculator is unpredictable in runnning \n\
        successfully with all Excel tab as dataframes in one go, so we instead \n\
        run each Excel tab as a dataframe in pairs or individually through the \n\
        web form data scrapper
    
        This program handles:
            - absence of lat/lon in rows and skips to next row for processing
            - presence of t, Lu, Ld in rows and skips to next row for processing
            - many exceptions and will export them as csv error logs
    
    WARNING: WebDriverException: TypeError: rootNode is null \n\
             if this error returns in stack trace, it is likely timeout from  \n\
             server side of broswer. To remedy, try re-running the script until it works. 
             
    $Future Implementation: allow error log to not be overwritten when exceptions occur with script re-run... 
    $Future Implementation: add provision to close driver after exhaust... 
'''

# parse dataframes
for sheet in range(len(all_sheets)):   
    
    # empty  after each dataframe (excel tab) end
    surface_condition_exception_log = []
    scene_will_not_process_exception_log = []
    skip_cities_without_lat_lon_list = []
    
    # grab row by row... df.iterows()
    for index, row in all_sheets[sheet].iterrows():
        
        # grab name for output as csv
        for sheet_name in range(len(all_sheets_name)):  
            
            if sheet == sheet_name:
                
                # create outfolders if not exist... 
                outfolder = os.path.join(Atmos_Calc_folder_root, all_sheets_name[sheet])
                if not os.path.exists(outfolder):
                    os.makedirs(outfolder)

                # step 2. 
                # check for absence of lat/lon and skip (i.e., cannot process without lat/lon)
                if pd.isnull(row['Station Latitude (DD)']) and pd.isnull(row['Station Longitude (DD)']):
                    print row['SCENE_NAME'], row['City'] + ' not processed; lat/lon missing'
                    print ' '*2
                    
                    skip_cities_without_lat_lon_list.append( '{0}, {1}'.format(row['SCENE_NAME'], row['City']) )        
                
                    df_skip_cities_without_lat_lon_list = pd.DataFrame(skip_cities_without_lat_lon_list)
                    df_skip_cities_without_lat_lon_list.to_csv(os.path.join(outfolder,'not_processed_cities_log.csv'), sep=',', encoding='utf-8', index=False) # factor unicode char       
                
                # check for presence of t, lu, ld and skip... (i.e., already processed)
                elif not pd.isnull(row['t']) and not pd.isnull(row['Lu']) and not pd.isnull(row['Ld']):
                    
                    print row['SCENE_NAME'] + ' is already processed...'
                    print ' '*2
                    
                # check for absence of t, lu, ld and process... 
                elif pd.isnull(row['t']) and pd.isnull(row['Lu']) and pd.isnull(row['Ld']):
                    
                    print ('{0}, {1}, {2}, {3} is processing...'.format(row['SCENE_NAME'], row['t'], row['Lu'], row['Ld']))
                    
            #        print row['SCENE_NAME']
                    
                    # Clean up and add city name to input/ouput images 
                    if '-' in row['City']:
                        City_Name = str(row['City']).replace('-','_')
                        
                    if ' ' in row['City']:
                        City_Name = str(row['City']).replace(' ','_')
                        
                    if '.' in row['City']:
                        City_Name = str(row['City']).replace('.','').replace(' ', '_')
                    
                    else: 
                        City_Name = str(row['City']).replace('.','').replace(' ', '_')
                    
                    
                    #merged_all_df = pd.concat(all_sheets, ignore_index=True) # concate all dataframes
                    #merged_summer_df = pd.concat(summer_sheets, ignore_index=True) # concate summer dataframes
                    #merged_winter_df = pd.concat(winter_sheets, ignore_index=True) # concate winter dataframes
                    #
                    #
                    #outpath_winter_season = r'D:\\'
                    #outname= 'test'
                    #merged_all_df.to_csv(os.path.join(outpath_winter_season,  outname + '.csv'), sep=',',  encoding='utf-8', index=False)
            
                    
                    ## DATE/TIME 
                    year_input = driver.find_element_by_name('year')
                    year_input.send_keys( int(str(row['DATE_ACQUIRED (YYYY-MM-DD)']).split('-')[0]) )   # input from csv... 
                    
                    month_input = driver.find_element_by_name('month')
                    month_input.send_keys( int(str(row['DATE_ACQUIRED (YYYY-MM-DD)']).split('-')[1]) ) # input from csv... 
                    
                    day_input = driver.find_element_by_name('day')
                    day_input.send_keys(int(row['DATE'])) # input from csv... 
                    
                    GMT_hour_input = driver.find_element_by_name('hour')
                    GMT_hour_input.send_keys( int(str(row['SCENE_CENTER_TIME (GMT)']).split(':')[0]) ) # input from csv... 
                    
                    minute_input = driver.find_element_by_name('minute')
                    minute_input.send_keys( int(str(row['SCENE_CENTER_TIME (GMT)']).split(':')[1]) ) # input from csv... 
                    
                    
                    # error handling for empty lat/lon val 
                    
                    ## LAT/LON
                    latitude_input = driver.find_element_by_name('thelat')
                    latitude_input.send_keys(str(row['Station Latitude (DD)']))  # form input specified as string arg
                    
                    longitude_input = driver.find_element_by_name('thelong')
                    longitude_input.send_keys(str(row['Station Longitude (DD)'])) # form input specified as string arg
                    
                    
                    ## not used, only for posterity purposes 
                    #atmospheric_profile_button = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[5]/td/input[1]')
                    #atmospheric_profile_button.click() 
                    
                    interpolated_atmos_profile_button = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[5]/td/input[2]')
                    interpolated_atmos_profile_button.click() # input from csv... 
                    
                    
                    ## *** Fix logic to be more robust... 
                    # summer sheets...
                    if 'summer' in all_sheets_name[sheet_name]:
    #                if sheet % 2 == 0:
                        ## add conditionals for winter vs. summer scenes... winter/summer tab, except for south america/oceania need to fix this...
                        mid_lat_summer_std_atmos = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[7]/td/input[1]')
                        mid_lat_summer_std_atmos.click()
                    
                    # for countries in equatorial region; t, lu, ld outputs same with summer/winter atmospheric profiles...
                    elif 'equatorial_fix' in all_sheets_name[sheet_name]:
                        mid_lat_winter_std_atmos_equatorial_fix = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[7]/td/input[1]')
                        mid_lat_winter_std_atmos_equatorial_fix.click()
                        
                    # winter sheets...    
                    elif 'winter' in all_sheets_name[sheet_name]:
                        mid_lat_winter_std_atmos = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[7]/td/input[2]')
                        mid_lat_winter_std_atmos.click()
                    
                    
                    ## Get spectral response curve for appropriate Landsat sensor... add conditionals 
                    if str(row['SCENE_NAME']).split('_')[0] == 'LC08':    
                        LS8_TIRS_Spec_Response_Curve =  driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[9]/td/input[1]')
                        LS8_TIRS_Spec_Response_Curve.click()
                        
                    elif str(row['SCENE_NAME']).split('_')[0] == 'LE07':
                        LS7_TIRS_Spec_Response_Curve = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[9]/td/input[2]')
                        LS7_TIRS_Spec_Response_Curve.click()
                        
                    elif str(row['SCENE_NAME']).split('_')[0] == 'LT05':
                        LS5_TIRS_Spec_Response_Curve = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[9]/td/input[3]')
                        LS5_TIRS_Spec_Response_Curve.click()
                    
                    ## not used, only for posterity purposes 
                    #only_atmos_profile = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[9]/td/input[4]')
                    #only_atmos_profile.click()
                    
                    ###### ****WARNING
                    ###### ADD CONDITIONAL FOR NULL... IF VALUES RETURN ERROR, GO BACK AND LEAVE EMPTY...  ###### 
                    ## Surface Conditions
                    altitude_input = driver.find_element_by_name('altitude')
                    altitude_input.send_keys(str(row['Station Elevation (km)'])) # form input specified as string arg
                    
                    pressure_input = driver.find_element_by_name('pressure')
                    pressure_input.send_keys(str(row['Station Pressure (mb)'])) # form input specified as string arg
                    
                    temperature_input = driver.find_element_by_name('temperature')
                    temperature_input.send_keys(str(row[u'Station Temperature (°C)'])) # form input specified as string arg
                    
                    relative_humidity_input = driver.find_element_by_name('rel_humid')
                    relative_humidity_input.send_keys(str(row['Station Relative Humidity (%)'])) # form input specified as string arg
                    
                    email_input = driver.find_element_by_name('user_email')
                    email_input.send_keys('alexander.tong@mail.utoronto.ca')
                    
                 
                    ## need overwite permission 
                    # take screenshot of inputs
#                    outfolder = Atmos_Calc_folder + all_sheets_name[sheet]
                    #outname = 'input_' + str(count)
                    inputs_outname = str(row['SCENE_NAME'])[:-7]  + City_Name + '_inputs'
                    driver.get_screenshot_as_file(os.path.join(outfolder, inputs_outname + '.png'))
            
            
                    # 'Calculate' button
                    # https://sqa.stackexchange.com/questions/3464/element-is-no-longer-attached-to-the-dom-staleelementreferenceexception-when-s
                    calculate_button = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[16]/td/input')
                    calculate_button.click() #re-initialize calculate_button once you go to next page and then back again... 
                     
            
                    # get atmospheric functions
                    try:
                        # make driver wait until page changes 
                        current_url = 'https://atmcorr.gsfc.nasa.gov/'
                        
                        # https://stackoverflow.com/questions/42069503/python-selenium-wait-until-next-page-has-loaded-after-form-submit 
                        WebDriverWait(driver, 120).until(EC.url_changes(current_url)) # wait for URL to change with 120 seconds timeout
                        
                        new_url = driver.current_url
                        print new_url
                        print "Page is ready!"
                        print ' '
                        
                        time.sleep(5)
                        
                        try:
                            # get text body of atmospheric function outputs 
                            calc_out_val = driver.find_element_by_xpath('/html/body/pre')
                            
                            time.sleep(5)
                            
                #            delay = 3 # seconds
                #            calc_out_val = WebDriverWait(driver, delay).until(EC.presence_of_element_located((By.XPATH, '/html/body/pre'))) # https://stackoverflow.com/questions/34504839/how-do-i-use-seleniums-wait
                            
                        # if surface conditions failed to calculate, empty and try again... 
                        except NoSuchElementException as surface_condition_error:
                            
                            print (surface_condition_error) 
                            print ('')
                            
                            time.sleep(5)
                            
                            driver.back()
                            # **WARNING: BELOW WORKS FOR GOING BACK, BUT NOW FOR SOME REASON, IT DOES ALLOW object.click() TO WORK
#                            driver.execute_script("window.history.go(-1)") # 
                            
                            # Clear and re-add values
                            # 'Clear Fields' button 
                            clear_button = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[17]/td/input')
                            clear_button.click()
                                                        
                            ## DATE/TIME 
                            year_input = driver.find_element_by_name('year')
                            year_input.send_keys( int(str(row['DATE_ACQUIRED (YYYY-MM-DD)']).split('-')[0]) )   # input from csv... 
                            
                            month_input = driver.find_element_by_name('month')
                            month_input.send_keys( int(str(row['DATE_ACQUIRED (YYYY-MM-DD)']).split('-')[1]) ) # input from csv... 
                            
                            day_input = driver.find_element_by_name('day')
                            day_input.send_keys(int(row['DATE'])) # input from csv... 
                            
                            GMT_hour_input = driver.find_element_by_name('hour')
                            GMT_hour_input.send_keys( int(str(row['SCENE_CENTER_TIME (GMT)']).split(':')[0]) ) # input from csv... 
                            
                            minute_input = driver.find_element_by_name('minute')
                            minute_input.send_keys( int(str(row['SCENE_CENTER_TIME (GMT)']).split(':')[1]) ) # input from csv... 
                            
                            
                            ## LAT/LON
                            latitude_input = driver.find_element_by_name('thelat')
                            latitude_input.send_keys(str(row['Station Latitude (DD)']))  # form input specified as string arg
                            
                            longitude_input = driver.find_element_by_name('thelong')
                            longitude_input.send_keys(str(row['Station Longitude (DD)'])) # form input specified as string arg
                            
                            
                            ## not used, only for posterity purposes 
                            #atmospheric_profile_button = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[5]/td/input[1]')
                            #atmospheric_profile_button.click() 
                            
                            interpolated_atmos_profile_button = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[5]/td/input[2]')
                            interpolated_atmos_profile_button.click() # input from csv... 
                            
                            
                            ## *** Fix logic to be more robust... 
                            ## use sheet_name to identify; broken for south america at the moment (select cases; we can flag here...)                       
                            # summer sheets...
                            if 'summer' in all_sheets_name[sheet_name]:
    #                        if sheet % 2 == 0:
                                ## add conditionals for winter vs. summer scenes... winter/summer tab, except for south america/oceania need to fix this...
                                mid_lat_summer_std_atmos = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[7]/td/input[1]')
                                mid_lat_summer_std_atmos.click()
                                
                            elif 'equatorial_fix' in all_sheets_name[sheet_name]:
                                mid_lat_winter_std_atmos_equatorial_fix = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[7]/td/input[1]')
                                mid_lat_winter_std_atmos_equatorial_fix.click()      
                                
                            # winter sheets...    
                            elif 'winter' in all_sheets_name[sheet_name]:
                                mid_lat_winter_std_atmos = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[7]/td/input[2]')
                                mid_lat_winter_std_atmos.click()
                            
                            
                            ## Get spectral response curve for appropriate Landsat sensor... add conditionals 
                            if str(row['SCENE_NAME']).split('_')[0] == 'LC08':    
                                LS8_TIRS_Spec_Response_Curve =  driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[9]/td/input[1]')
                                LS8_TIRS_Spec_Response_Curve.click()
                                
                            elif str(row['SCENE_NAME']).split('_')[0] == 'LE07':
                                LS7_TIRS_Spec_Response_Curve = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[9]/td/input[2]')
                                LS7_TIRS_Spec_Response_Curve.click()
                                
                            elif str(row['SCENE_NAME']).split('_')[0] == 'LT05':
                                LS5_TIRS_Spec_Response_Curve = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[9]/td/input[3]')
                                LS5_TIRS_Spec_Response_Curve.click()
                            
                                    
                            ## Surface Conditions; will handle for empty val... 
                            altitude_input = driver.find_element_by_name('altitude')
                            altitude_input.send_keys() # form input specified as string arg
                            
                            pressure_input = driver.find_element_by_name('pressure')
                            pressure_input.send_keys() # form input specified as string arg
                            
                            temperature_input = driver.find_element_by_name('temperature')
                            temperature_input.send_keys() # form input specified as string arg
                            
                            relative_humidity_input = driver.find_element_by_name('rel_humid')
                            relative_humidity_input.send_keys() # form input specified as string arg
                            
                            email_input = driver.find_element_by_name('user_email')
                            email_input.send_keys('alexander.tong@mail.utoronto.ca')
                                                        
                            time.sleep(10)
                            
                            
                            # overwrite screenshot of inputs to reflect outputs
                            if os.path.exists(os.path.join(outfolder, inputs_outname + '.png')):
                                
                                os.remove(os.path.join(outfolder, inputs_outname + '.png'))
                                
                                # take screenshot of inputs
#                                outfolder = Atmos_Calc_folder + all_sheets_name[sheet]
                                #outname = 'input_' + str(count)
                                inputs_outname = str(row['SCENE_NAME'])[:-7]  + City_Name + '_inputs' 
                                driver.get_screenshot_as_file(os.path.join(outfolder, inputs_outname + '.png'))
                                
                             
                            # log errors and save out as csv
#                            surface_condition_exception_log = []                           
                            surface_condition_exception_log.append('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12}'.format(str(row['SCENE_NAME'])[:-8],
                                                                                                                                    City_Name, 
                                                                                                                                    int(str(row['DATE_ACQUIRED (YYYY-MM-DD)']).split('-')[0]),
                                                                                                                                    int(str(row['DATE_ACQUIRED (YYYY-MM-DD)']).split('-')[1]),
                                                                                                                                    int(row['DATE']),
                                                                                                                                    int(str(row['SCENE_CENTER_TIME (GMT)']).split(':')[0]),
                                                                                                                                    int(str(row['SCENE_CENTER_TIME (GMT)']).split(':')[1]),
                                                                                                                                    str(row['Station Latitude (DD)']),
                                                                                                                                    str(row['Station Longitude (DD)']),
                                                                                                                                    str(row['Station Elevation (km)']),
                                                                                                                                    str(row['Station Pressure (mb)']),
                                                                                                                                    str(row[u'Station Temperature (°C)']),
                                                                                                                                    str(row['Station Relative Humidity (%)'])))
                            
                            df_surface_condition_exception_log_headers = ['SCENE_NAME','CITY','YEAR','MONTH','DAY',
                                                                          'GMT_HOUR','GMT_MINUTE',
                                                                          'LATITUDE_DD','LONGITUDE_DD',
                                                                          'Station Elevation (km)',
                                                                          'Station Pressure (mb)',
                                                                          u'Station Temperature (°C)',
                                                                          'Station Relative Humidity (%)']
                            
                            
                            for i in range(len(surface_condition_exception_log)):
                                # skip list element if already split; 13 is the total num of elements corresponding to columns 
                                if len(surface_condition_exception_log[i]) > 13: #** WARNING: change int value to arg call when making this logic into a function
                                    surface_condition_exception_log[i] = surface_condition_exception_log[i].split(',')
                                
                            df_surface_condition_exception_log = pd.DataFrame(surface_condition_exception_log)
#                            df_surface_condition_exception_log = df_surface_condition_exception_log.transpose() # NOT NEEDED ANYMORE
                            df_surface_condition_exception_log.columns = df_surface_condition_exception_log_headers
                            

                            df_surface_condition_exception_log.to_csv(os.path.join(outfolder,'scene_surface_condition_exception_log.csv'), sep=',',  encoding='utf-8', index=False)
                                                      
                            
                            # calculate...                 
                            calculate_button = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[16]/td/input')
                            calculate_button.click() #re-initialize calculate_button once you go to next page and then back again... 
                            
                            print ('{0} ({1}) surface conditions could not be processed. Reverted to default (i.e., no surface conditions specified) to retrieve atmospheric functions'.format(row['City'],(row['SCENE_NAME'])[:-8]))
                            print  ' '
                            
                            # make driver wait until page changes 
                            current_url = 'https://atmcorr.gsfc.nasa.gov/'
                            
                            # https://stackoverflow.com/questions/42069503/python-selenium-wait-until-next-page-has-loaded-after-form-submit 
                            WebDriverWait(driver, 120).until(EC.url_changes(current_url)) # wait for URL to change with 120 seconds timeout
                            
                            new_url = driver.current_url
                            print (new_url)
                            print ("Page is ready!")
                            print (' ')
                            
                            time.sleep(5)
                            
                            try: 
                                # get text body of atmospheric function outputs 
                                calc_out_val = driver.find_element_by_xpath('/html/body/pre')
                                
                                time.sleep(5)
                        
                            # handle exception cases where Atmos Calc will not process scene given the mandatory inputs (e.g., YY-MM-DD, GMT HOUR, LAT/LON IN DD COORD)
                            # Issue known with Kingston, ON, Canada with summer scene (LC08_L1TP_016029_20180629_20180716_01_T1_MTL.txt)
                            except NoSuchElementException as Atmospheric_Calculator_Does_Not_Handle_All_Edge_Cases_Error:
                                                                                                  
                                print (Atmospheric_Calculator_Does_Not_Handle_All_Edge_Cases_Error)
                                print ('')
                                print ('Atmospheric Correction Parameter Calculator cannot process Landsat scene {0}{1}'.format(str(row['SCENE_NAME'])[:-7], City_Name) + '\n ' + '\n' +
                                       'There are a few scenarios that may be causing this: ' + '\n' +
                                       '(1) calculator is being overwhelmed with requests and will not process request correctly; \n\
                                       in this case, try re-running the script' + '\n' +
                                       '(2) not enough unit testing to handle edge cases (parameter inputs) for computation; \n\
                                       try another lat/lon coord close to the ones previously specified' + '\n' +
                                       '(3) calculator is not available for most dates before 2000, but year must be > 1998'  + '\n ' + '')
                                
                                
                                ## take screenshot of outputs                                
                                outputs_outname = str(row['SCENE_NAME'])[:-7]  + City_Name + '_outputs'  
                                
                                # overwrite screenshot of outputs to reflect outputs
                                if os.path.exists(os.path.join(outfolder, outputs_outname + '.png')):
                                   
                                    os.remove(os.path.join(outfolder, outputs_outname + '.png'))
                                    
                                
                                driver.get_screenshot_as_file(os.path.join(outfolder, outputs_outname + '.png'))
                                
                                time.sleep(5)
                                
                                scene_will_not_process_exception_log.append('{0},{1},{2},{3},{4},{5},{6},{7},{8}'.format(str(row['SCENE_NAME'])[:-8],
                                                                                                                        City_Name, 
                                                                                                                        int(str(row['DATE_ACQUIRED (YYYY-MM-DD)']).split('-')[0]),
                                                                                                                        int(str(row['DATE_ACQUIRED (YYYY-MM-DD)']).split('-')[1]),
                                                                                                                        int(row['DATE']),
                                                                                                                        int(str(row['SCENE_CENTER_TIME (GMT)']).split(':')[0]),
                                                                                                                        int(str(row['SCENE_CENTER_TIME (GMT)']).split(':')[1]),
                                                                                                                        str(row['Station Latitude (DD)']),
                                                                                                                        str(row['Station Longitude (DD)'])))
                                                                    
                                df_scene_will_not_process_exception_log_headers = ['SCENE_NAME','CITY',
                                                                                   'YEAR','MONTH','DAY',
                                                                                   'GMT_HOUR','GMT_MINUTE',
                                                                                   'LATITUDE_DD','LONGITUDE_DD']
                                
                                
                                for i in range(len(scene_will_not_process_exception_log)):
                                    # skip list element if already split; 9 is the total num of elements corresponding to columns 
                                    if len(scene_will_not_process_exception_log[i]) > 9: #** WARNING: change int value to arg call when making this logic into a function
                                        scene_will_not_process_exception_log[i] =  scene_will_not_process_exception_log[i].split(',')
                                    
                                df_scene_will_not_process_exception_log = pd.DataFrame(scene_will_not_process_exception_log) 
#                                df_scene_will_not_process_exception_log = df_scene_will_not_process_exception_log.transpose() # NOT NEEDED ANYMORE
                                df_scene_will_not_process_exception_log.columns = df_scene_will_not_process_exception_log_headers
                                
                                df_scene_will_not_process_exception_log.to_csv(os.path.join(outfolder,'scene_will_not_process_exception_log.csv'), sep=',',  encoding='utf-8', index=False)


                        # once we handle for surface condition and non-processing exception, if all inputs are valid, write out            
                        finally:    
                            
                            try:
                                
                                ## take screenshot of outputs                                
                                outputs_outname = str(row['SCENE_NAME'])[:-7]  + City_Name + '_outputs'  
                                
                                # overwrite screenshot of outputs to reflect outputs
                                if os.path.exists(os.path.join(outfolder, outputs_outname + '.png')):
                                   
                                    os.remove(os.path.join(outfolder, outputs_outname + '.png'))
                                    
                                
                                driver.get_screenshot_as_file(os.path.join(outfolder, outputs_outname + '.png'))
                                
                                 
                                # get text body of atmospheric function outputs 
                                calc_out_val = driver.find_element_by_xpath('/html/body/pre')
                                
                                time.sleep(5)
                                
                                print (calc_out_val.text)
                                print ('' + '\n ')
                                
                                # get text from xpath element, convert to list to isolate for Tau, Lu,Ld
                                atmos_func = calc_out_val.text
                                atmos_func_list = atmos_func.split('\n')
                                
                                Band_avg_atmoss_trans = float(atmos_func_list[11].split(':')[-1]) # Band average atmospheric transmission
                                Eff_bandpass_up_rad = float(atmos_func_list[12].replace('W/m^2/sr/um','').split(':')[-1]) # Effective bandpass upwelling radiance
                                Eff_bandpass_down_rad = float(atmos_func_list[13].replace('W/m^2/sr/um','').split(':')[-1]) # Effective bandpass downwelling radiance
                                
                                # update value at index...     
                                # https://stackoverflow.com/questions/23330654/update-a-dataframe-in-pandas-while-iterating-row-by-row
                                # format 
                                'Band average atmospheric transmission:    0.00'
                                'Effective bandpass upwelling radiance:    0.00 W/m^2/sr/um'
                                'Effective bandpass downwelling radiance:  0.00 W/m^2/sr/um'
                                
                                all_sheets[sheet].at[index,'t'] = Band_avg_atmoss_trans
                                all_sheets[sheet].at[index,'Lu'] = Eff_bandpass_up_rad
                                all_sheets[sheet].at[index,'Ld'] = Eff_bandpass_down_rad
                                
                                
                                # output dataframe as csv after each iteration (massive read/write computation, not the best solution)
                                # advantage of this approach is that if internet connection fails and computer loses power/fails, 
                                # there will be a csv with values where the web scrapper left off (failed)... 
                                # **add provision into workflow to leave off where it failed by looking for first NaN row of values in t, lu, ld columns... 
                                df_outname = all_sheets_name[sheet_name] + '.csv'
                                
                                all_sheets[sheet].to_csv(os.path.join(outfolder,df_outname), sep=',',  encoding='utf-8', index=False)
                                
#                            except InvalidSessionIdException as Atmospheric_Calculator_Does_Not_Handle_All_Edge_Cases_Error:
#                                
#                                pass
                            
                            except NoSuchElementException as Atmospheric_Calculator_Does_Not_Handle_All_Edge_Cases_Error:
                                print ('...')
                                pass
                            
                    except TimeoutException:
                        
                        # log errors and save out as csv
                        exception_log = []
                        exception_log.append(str(row['SCENE_NAME'])[:-7])
                    
                        df_timeout_exception_log = pd.DataFrame(exception_log)
                        df_timeout_exception_log.to_csv(os.path.join(outfolder,'timeout_exception_log.csv'), sep=',', encoding='utf-8', index=False)
                        
                        print ("Loading took too much time!")    
                        print (str(row['SCENE_NAME'])[:-7] + ' could not be processed')
                        
                        pass 
                    
                    
                    finally:
                        # go back...
#                        driver.back()
                        
                        # The back() and forward() methods aren't guaranteed to work.
                        # https://stackoverflow.com/questions/27626783/python-selenium-browser-driver-back
                        driver.execute_script("window.history.go(-1)")
                        
                        time.sleep(5)
                        
                        # 'Clear Fields' button
                        clear_button = driver.find_element_by_xpath('/html/body/center/form/center/table/tbody/tr/td/table/tbody/tr[17]/td/input')
                        clear_button.click()
                        
                        time.sleep(30) # 60 seconds (1 minute) pause; NASA asks for 120 seconds (2 minute) pause; honor this
            

## ISSUE EVEN WITHOUT SURFACE CONDITIONS, LAT/LON throwing error.... Kingston only example not working at the moment... 

# add provision to close driver after exhaust... 
                        

# re-combine all separate csv; this workflow assumes each Excel tab was individually /n
# processed (i.e., atmospheric calculator unpredictable in runnning successfully \n\
# with all Excel tab as dataframes in one go, so we instead run each Excel tab as a \n\
# dataframe in pairs or individually through the web form data scrapper)
def combine_atmos_func():
    '''
    $ Future Implementation: modify to only include Landsat images used for LST 
    '''
    csv_list = []
    directory = r'D:\GLUE Results\Process_Atmospheric_Functions'
    
    # grab and sort csv outputs from Atmospheric Correction Parameter Calculator (ACPC) 
    for root, dirs, files in os.walk(directory):
        for name in range(len(files)):
            if files[name].endswith('atmos_func.csv'):
                print files[name]
                csv_list.append(os.path.join(root,files[name]))
                csv_list.sort()
    
    
    df_list = []
    df_sheet_name = ['asia_summer_atmos_func','asia_winter_atmos_func',
                    'canada_summer_atmos_func', 'canada_winter_atmos_func',
                    'europe_et_al_summer_atmos_func', 'europe_et_al_winter_atmos_func',
                    'oceania_summer_atmos_func', 'oceania_winter_atmos_func',
                    'south_america_summer_atmos_func', 'south_america_winter_atmos_func',
                    'usa_summer_atmos_func', 'usa_winter_atmos_func']     
    
    for i in range(len(csv_list)):
        df_list.append(pd.read_csv(csv_list[i], encoding='utf8'))
        
    # https://stackoverflow.com/questions/46109675/write-list-of-dataframes-to-excel-using-to-excel
        # https://stackoverflow.com/questions/5893163/what-is-the-purpose-of-the-single-underscore-variable-in-python 
    writer = pd.ExcelWriter(r'D:\GLUE Results\Process_Atmospheric_Functions\combined_atmospheric_functions.xlsx')   
    
    # excel spreadsheet names == df_sheet_name on output 
    _ = [df.to_excel(writer,sheet_name="{0}".format(df_sheet_name[i])) for i, df in enumerate(df_list)]
    
    writer.save()

            
combine_atmos_func()   






## use openpyxl to work with excel...
'''
https://stackoverflow.com/questions/25247742/use-openpyxl-to-iterate-through-sheets-and-cells-and-update-cells-with-contante 
I dont think you can update cell contents. You can open a file to read, or open a \n\
new file to write to. I think you have to create a new workbook, and every cell that \n\
you read, if you choose to not modify it, write it out to your new workbook. In your \n\
sample code, you are overwriting wb (used to read) with the wb (used to write). \n\
Pull it out of the for loop, assign a different name to it.
''' 
#from openpyxl import load_workbook
#import os 
#
#book = load_workbook(os.path.join(r'D:\\', 'GLUE_City_LST_retrieval.xlsx'))
#
#book.sheetnames 
#
#sheet = book.active 

