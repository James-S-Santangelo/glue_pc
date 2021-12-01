# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 11:53:42 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15

### using landsat code list, need to associate raster image name (splice for row/path code) with code in list, then code in list to city name
### city name from city should then match csv name...  


** WARNING: issue with Woodstock 
James's modified script produces csv with 'NA' in lat/lon columns; remove 'NA' 
and keep null fields; issue with arcpy is that it should revert to WARNING 000635 to
skip empty geometry when batch processing csv from a folder. However, it is
instead throwing error 99999.... to remedy, process Woodstock.csv in separate folder works without issue... ??? 

"""
import os, sys, re 

try: 
    import arcpy 
    
except ImportError as IE:
    print (IE)
    print ("These functions requires arcpy to run")  
    sys.exit(1)   
    
try:     
    directory = r'D:\Python scripts - FINAL'
    os.chdir(directory)
    
    from _00_glue_utils import raster 
    from _00_glue_utils import lookup_table
    
except ImportError as IE:
    print (IE)
    print ("These functions requires arcpy to run")  
    sys.exit(1)   

def csv_to_shp(dir_csv, raster_file, dir_shp_gcs, dir_shp_pcs, dir_shp_pcs_buffer, buffer_radius, city_landsat_list):
    '''
    Description: 
        This function converts .csv with lat/lon data into 3 outputs:
            
             (1) point shapefile with geographic coord sys in (GCS) lat/lon units (decimal degrees; DD)
             (2) point shapefile with projected coord sys (PCS) in UTM units (metres)
             (3) buffered 100 m radius point shapefile with projected coord sys (PCS) in UTM units (metres)
             
        The point shapefile with PCS is made possible with the list of lists 
        with hardcoded Landsat images; the spatial ref PCS is pulled from the Landsat images 
        

        **WARNING: because this function uses os.walk, there cannot be any other \n\
                   csv files except for the ones to be converted to shapefiles    
        
        **WARNING: ensure all Landsat images being referenced are in PCS (e.g., UTM coord sys)       


        **WARNING: issue with grabbing spatial reference from correct landsat \n\
                   image for some GLUE cities
              
            Why? 
            (1) logic used to convert csv to shp is based on partial string matching \n\
            e.g.,
                in the case of Charlotte, NC vs. Charlottetown, PEI, wherefore Charlottetown \n\
                inherited Charlotte, NC spatial ref because of partial string matching
            
            (2) adjacent images of +- UTM zone; first landsat image to be parsed from \n\
                lookup table will take precedent in deciding spatial reference \n\
                of shapefile. There is no provision to create multiple shapefiles  \n\
                for each GLUE city, so when function is executed, it will delete the temp \n\
                event layer and continue... 
                
            e.g.,  Landshut, 192026 - WGS_1984_UTM_Zone_32N vs. WGS_1984_UTM_Zone_33N
                   Landshut, 193026 - WGS_1984_UTM_Zone_32N vs. WGS_1984_UTM_Zone_32N
                   
                   Norrtalje, 192018 - WGS_1984_UTM_Zone_33N vs. WGS_1984_UTM_Zone_34N
                   Norrtalje, 193018 - WGS_1984_UTM_Zone_33N vs. WGS_1984_UTM_Zone_33N     

            Solution?
            - Remove extraneous Landsat row/path entry from lookup table 
            - EXTACT STRING MATCHING BETWEEN filtered_image_spatial_ref VS filenames[file][:-4]

                                                            
    Args:
        dir_csv (str): directory of .csv to be converted to shapefiles.
        raster_file (str): list of Landsat rasters.
        dir_shp_gcs (str): directory for output shapefiles in gcs to be saved.
        dir_shp_pcs (str): directory for output shapefiles in pcs to be saved.
        dir_shp_pcs_buffer (str): directory for output buffered shapefiles in pcs to be saved.
        buffer_radius (int): set buffer radius (metres) for output buffered shapefiles in pcs to be saved. e.g., 100 
        city_landsat_list (str): list of list of city to landsat code
    
    Returns:
        No returns or exchanges.        
        
    
    $ logic requires modification since it is not robust enough to handle cities with same name; \\
    $ need to include country or country + province/state/perfecture as a measure 
    $ temporary fix is to add country + province/state/perfecture to city name csv + city name in city to landsat code
    >>> currently only impacts: Portland, OR and Portland, ME; rename csv to Portland_OR and Portland_ME
    '''
    arcpy.env.overwriteOutput = True
    
    get_image_spatial_ref = []
    
    # step 1. get all image spatial ref and filter to a single entry for each city for processing    
    for root, dirnames, filenames in os.walk(dir_csv):
    
        for file in range(len(filenames)):
            if filenames[file].endswith('.csv'):
    #            print (filenames[file])
                for i in range(len(city_landsat_list)):
                    
                    # satisfy for '_' and '-' characters 
                    if city_landsat_list[i][1].lower() in filenames[file].lower().replace('_', ' ') or city_landsat_list[i][1].lower() in filenames[file].lower().replace('-', ' '):
                       
                        # get Landsat image 
                        for j in range(len(raster_file)): 
                            
                            #splice for code 
                            a = os.path.basename(raster_file[j])[:10]
                            b = os.path.basename(raster_file[j])[16:]
                        
                            raster_code = os.path.basename(raster_file[j].split(a)[-1].split(b)[0])
                            
                            #compare city name to raster code; will only compare to composited image
                            if raster_code in city_landsat_list[i][2]:
                                
    #                            print raster_code
                                inRas = arcpy.Raster(raster_file[j])
                                                                # arcpy.ProjectRaster_management DOES NOT ACCEPT str, 
                                # therefore do not use spatialReference.name, but instead spatialReference.factorycode (corresponds to EPSG code)
                                spatial_ref = arcpy.Describe(inRas).spatialReference
                                
                                # specify EPSG equivalent
                                ### If an Esri well-known ID (WKID) is below 32767, it corresponds to the EPSG ID. WKIDs that are 32767 or above are Esri-defined.
                                print ('%s, %s, %s' % (city_landsat_list[i][1], spatial_ref.factorycode, os.path.basename(raster_file[j])))
                                
                                # factorycode = WKID
                                get_image_spatial_ref.append((city_landsat_list[i][1].lower(), spatial_ref.factorycode))
                                
                                # get only single entry using set  
                                filtered_image_spatial_ref = sorted(set(get_image_spatial_ref))
        
    
    # step 2: get csv by lat/lon coord
    for root, dirnames, filenames in os.walk(dir_csv):
       
        for file in range(len(filenames)):
            if filenames[file].endswith('.csv'):
#                print filenames[file]
                
                with open(root + '\\' + filenames[file], 'rU') as f:
        #                    reader = csv.reader(f) 
        #                    i = reader.next()
                    
                    first_line = f.readline()
           
                    # check for ',' delimiter 
                    if ',' in first_line:
        #                        print filenames[file]
                        split_first_line = first_line.lower().strip().split(',') # capture individual header names as a list; lower case and remove \r\n 
                        
                        # remove whitespaces; arcpy very finicky 
                        for i in range(len(split_first_line)):
                            if ' ' in split_first_line[i]:
                                split_first_line[i] = split_first_line[i].rstrip()
                                split_first_line[i] = split_first_line[i].lstrip()
                                
                         
                        for i in range(len(split_first_line)): 
                            # grab latitude col  
                            if 'population_latitude' in split_first_line[i].lower():
                                lat_index = i 
        #                            print split_first_line[i], i
                                
                            # grab longitude col  
                            elif 'population_longitude' in split_first_line[i].strip(' ').lower():
                                lon_index = i 
        #                            print split_first_line[i], i 
                    
#                            # for James Santangelos' transect datasets 
#                            elif 'lat.pop' in split_first_line[i].lower():
#                                lat_index = i 
#                            # for James Santangelos' transect datasets     
#                            elif 'long.pop' in split_first_line[i].lower():
#                                lon_index = i 
                        
                        #Default GCSE: datum WGS 1984; does not change; int is EPSG         
                        GCS = arcpy.SpatialReference(4326)
                        
                                            
#                        out_layer = str(filenames[file][:-14]) 
#                        GCS_fc= str(filenames[file][:-14]) 
#                        PCS_fc = str(filenames[file][:-14]) 
                        out_layer = str(filenames[file][:-4]) 
                        GCS_fc= str(filenames[file][:-4]) 
                        PCS_fc = str(filenames[file][:-4])                         
    #                        print out_layer
    #                        print GCS_fc
    #                        print split_first_line[lon_index]
    #                        print split_first_line[lat_index]
                        
    
                        # if file(s) exists, skip... 
#                        if os.path.isfile(dir_shp_gcs + '\\' + GCS_fc + '_gcs.shp'):
#                            print 'the following file already exists and was not processed...' + dir_shp_gcs + '\\' + GCS_fc + '_gcs.shp'
                            
                         # if file(s) exists, skip... 
                        if os.path.isfile(dir_shp_pcs + '\\' + PCS_fc + '_pcs.shp'):
                            print 'the following file already exists and was not processed...' + dir_shp_pcs + '\\' + PCS_fc + '_pcs.shp'                       
#                            pass
                        
                        # step 3: convert csv to shapefiles incl. 
                        #         (a) make xy event layer
                        #         (b) make feature layer
                        #         (c) save feature layer out as shapefile with gcs (geographic coordinate sys)
                        #         (d) save feature layer out as shapefile with pcs (projected coordinate sys)
                        #         (e) save feature layer out as shapefile with pcs with 100 m radius (projected coordinate sys)                        
                        # else process file(s)...
                        else:
                            try:
                                
                                # satisfies for making xy event layer
                                arcpy.env.workspace = root 
                                
                                # not grabbing field name from csv; might be forcing the use of arcpy to read 
                                EVENT_LAYER = out_layer + '_EventLayer'
                                arcpy.MakeXYEventLayer_management(str(filenames[file]), str(split_first_line[lon_index]), str(split_first_line[lat_index]), EVENT_LAYER, GCS)

                                for i in range(len(filtered_image_spatial_ref)):    
                                         
                                    
                                    # satisfy for '_' and '-' characters 
#                                    if filtered_image_spatial_ref[i][0] in filenames[file][:-14].replace('_', ' ').lower() or filtered_image_spatial_ref[i][0] in filenames[file][:-14].replace('-', ' ').lower():
                                    if filtered_image_spatial_ref[i][0] == filenames[file][:-4].replace('_', ' ').lower() or filtered_image_spatial_ref[i][0] == filenames[file][:-4].replace('-', ' ').lower():
                                                                                    
                                        print filtered_image_spatial_ref[i][0] 
                                        
#                                       WKID_EPSG = str(filtered_image_spatial_ref[i][1]) # WKID = well-known id; EPSG = European Petroleum Survey Group
                                        
                                        EVENT_LAYER_TO_TEMP_SHAPEFILE = GCS_fc + '_gcs.shp'
                                        SHAPEFILE_GCS = GCS_fc + '_gcs' + '.shp'
#                                            SHAPEFILE_PCS = PCS_fc  + '_' + WKID_EPSG + '_pcs' + '.shp'
#                                            SHAPEFILE_PCS_BUFFER = PCS_fc + '_' + WKID_EPSG + '_buffer' + '.shp'
                                        SHAPEFILE_PCS = PCS_fc + '_pcs' + '.shp'
                                        SHAPEFILE_PCS_BUFFER = PCS_fc + '_pcs' + '_buffer' + '.shp'
                                        
                                        # final check for '-' because arcpy does not allow shapefile names with '-' for buffer analysis apparently.. 
                                        if '-' in SHAPEFILE_PCS_BUFFER:
                                            SHAPEFILE_PCS_BUFFER = SHAPEFILE_PCS_BUFFER.replace('-','_')

                                        # technically, this should be called after the shapefile with gcs is created, \n\
                                        # but we need to re-name Portland to respective state \n\
                                        #grab Landsat SRID                                                                                
                                        # charlotte in charlottetown == True; disambiguation fix 
#                                        if  filtered_image_spatial_ref[i][0] == 'charlottetown': 
#                                            spatial_ref = arcpy.SpatialReference(32620)
#                                            
#                                        else:
                                        spatial_ref = arcpy.SpatialReference(filtered_image_spatial_ref[i][1])

                                       
                                        ## below clause is no longer considered; revised lookup table for disambiguation of portland as 'portland_or' & 'portland_me'
                                        # add clause to remain Portland Oregon vs. Portland Maine... 
#                                        if 'Portland_OR' in GCS_fc and filtered_image_spatial_ref[i][0] == 'portland_or':
#                                            # Portland Oregon; spatial_ref.factorycode; spatial_ref.PCSCode 
#                                            spatial_ref = arcpy.SpatialReference(filtered_image_spatial_ref[i][1])
#                                            
#                                            print 'test'
#                                            SHAPEFILE_GCS = GCS_fc + '_gcs.shp'
#                                            SHAPEFILE_PCS = PCS_fc + '_pcs' + '.shp'
#                                            SHAPEFILE_PCS_BUFFER = PCS_fc + '_pcs' + '_buffer' + '.shp'
#                                            
#                                        if 'Portland_ME' in GCS_fc and filtered_image_spatial_ref[i][0] == 'portland_me':    
#                                            
#                                            print 'test'
#                                            # Portland Maine;  spatial_ref.factorycode; spatial_ref.PCSCode 
#                                            spatial_ref = arcpy.SpatialReference(filtered_image_spatial_ref[i][1])
#                                            
#                                            print spatial_ref.PCSCode
#
#                                            SHAPEFILE_GCS = GCS_fc + '_gcs.shp'
#                                            SHAPEFILE_PCS = PCS_fc  + '_pcs' + '.shp'
#                                            SHAPEFILE_PCS_BUFFER = PCS_fc  + '_pcs' + '_buffer' + '.shp'

                                        arcpy.MakeFeatureLayer_management(EVENT_LAYER, EVENT_LAYER_TO_TEMP_SHAPEFILE) 
                                        arcpy.CopyFeatures_management(EVENT_LAYER_TO_TEMP_SHAPEFILE, os.path.join(dir_shp_gcs, SHAPEFILE_GCS)) # save shapefile with gcs 

                                        arcpy.env.workspace = dir_shp_gcs # direct workspace for gcs shapefile input 
                                          
                                        # project and buffer 
                                        buffer_zone = buffer_radius
                                        arcpy.Project_management(SHAPEFILE_GCS, os.path.join(dir_shp_pcs, SHAPEFILE_PCS), spatial_ref)
                                        
                                        arcpy.env.workspace = dir_shp_pcs # re-direct workspace for pcs shapefile input 
                                        
                                        arcpy.Buffer_analysis(SHAPEFILE_PCS, os.path.join(dir_shp_pcs_buffer,SHAPEFILE_PCS_BUFFER), buffer_zone) #default meter distance   
        
        
                                        arcpy.Delete_management(EVENT_LAYER)
                                        arcpy.Delete_management(EVENT_LAYER_TO_TEMP_SHAPEFILE)   
                                    
                                    # non-alpha and cities with 2 words         
#                                    elif filtered_image_spatial_ref[i][0] in filenames[file][:-14].replace('_', ' ').lower() or filtered_image_spatial_ref[i][0] in filenames[file][:-14].replace('-', ' ').lower():
                                    elif filtered_image_spatial_ref[i][0] == filenames[file][:-4].replace('_', ' ').lower() or filtered_image_spatial_ref[i][0] == filenames[file][:-4].replace('-', ' ').lower():

                                        # for csv with non-alphabetic char(s) in name AND WITH CITY NAMES WITH 2 WORDS
#                                        if not filenames[file][:-14].split('_')[0].isalpha() or not filenames[file][:-14].split('_')[1].isalpha():
                                        if not filenames[file][:-4].split('_')[0].isalpha() or not filenames[file][:-4].split('_')[1].isalpha():
                                            
                                            # if non-alpha char(s) in name, replace with '_' (arcpy does not accept non-alpha char)
                                            GCS_fc = re.sub('[^0-9a-zA-Z]+', '_', (filenames[file]))
                                            
                                            # python logic will not honor declared variables GCS_fc or PCS_fc above, so we must re-initialize
                                            # final check for '-' because arcpy does not allow shapefile names with '-' ...
#                                            GCS_fc = filenames[file][:-14].replace('-','_')
#                                            PCS_fc = filenames[file][:-14].replace('-','_')
                                            GCS_fc = filenames[file][:-4].replace('-','_')
                                            PCS_fc = filenames[file][:-4].replace('-','_')
                                            
                                            print GCS_fc
                                            
    #                                            WKID_EPSG = str(filtered_image_spatial_ref[i][1]) # WKID = well-known id; EPSG = European Petroleum Survey Group
                                                                       

                                            EVENT_LAYER_TO_TEMP_SHAPEFILE = GCS_fc + '_station' + '_gcs.shp'
                                            SHAPEFILE_GCS = GCS_fc + '_station' + '_gcs' + '.shp'
    #                                            SHAPEFILE_PCS = PCS_fc  + '_' + WKID_EPSG + '_pcs' + '.shp'
    #                                            SHAPEFILE_PCS_BUFFER = PCS_fc + '_' + WKID_EPSG + '_buffer' + '.shp'
                                            SHAPEFILE_PCS = PCS_fc + '_station' + '_pcs' + '.shp'
                                            SHAPEFILE_PCS_BUFFER = PCS_fc + '_station' + '_pcs' + '_buffer' + '.shp'
                 
                    
                                            arcpy.MakeFeatureLayer_management(EVENT_LAYER, EVENT_LAYER_TO_TEMP_SHAPEFILE)   
                                            
                                            
                                            arcpy.CopyFeatures_management(EVENT_LAYER_TO_TEMP_SHAPEFILE, os.path.join(dir_shp_gcs, SHAPEFILE_GCS)) # save shapefile with gcs 
                                            
                                            #grab Landsat SRID
                                            spatial_ref = arcpy.SpatialReference(filtered_image_spatial_ref[i][1])
                                            
                                            arcpy.env.workspace = dir_shp_gcs # direct workspace for gcs shapefile input 
                                            
    #                                            raster_code = str(filtered_image_spatial_ref[i][1])
                                            
                                            # project and buffer 
                                            buffer_zone = buffer_radius
                                            arcpy.Project_management(SHAPEFILE_GCS, os.path.join(dir_shp_pcs, SHAPEFILE_PCS), spatial_ref)
                                            
                                            arcpy.env.workspace = dir_shp_pcs # re-direct workspace for pcs shapefile input 
                                            arcpy.Buffer_analysis(SHAPEFILE_PCS, os.path.join(dir_shp_pcs_buffer,SHAPEFILE_PCS_BUFFER), buffer_zone) #default meter distance   
            
                                            arcpy.Delete_management(EVENT_LAYER)
                                            arcpy.Delete_management(EVENT_LAYER_TO_TEMP_SHAPEFILE)   
                                            

                            except Exception as e:
                                print( "Error: {0} on line {1}".format(e, sys.exc_info()[-1].tb_lineno) )
                                pass  
                                

def main():
    
    raster_file = raster(r'G:\Landsat_Download', 'sr_composite.tif')
#    dir_csv = r'D:\GLUE Datasets\popMeans_allCities'
#    dir_csv = r'D:\GLUE Datasets\woodstock_test'
    dir_csv = r'D:\GLUE Datasets\popMeans_allCities'
    
    dir_shp_gcs = r'D:\GLUE Results\Process_Transect\gcs'
    dir_shp_pcs = r'D:\GLUE Results\Process_Transect\pcs'
    dir_shp_pcs_buffer = r'D:\GLUE Results\Process_Transect\pcs_buffer'
    buffer_radius = 100
    city_landsat_list = lookup_table(r'D:\Python scripts - FINAL')  
    csv_to_shp(dir_csv, raster_file, dir_shp_gcs, dir_shp_pcs, dir_shp_pcs_buffer, buffer_radius, city_landsat_list)
    
    
if __name__ == "__main__":
    
    main()

