# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 13:23:23 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15

# step 1:
mosaic DEM if necessary
    - rename derivative raster as <city>_dem_gcs

# step 2:
    - using landsat code, project raster as <city>_dem_pcs 
    
logic only works if no other tifs in the directory to be processed BUT the original files 


# using gdal; call in sub_process 
https://gis.stackexchange.com/questions/236746/calling-gdal-merge-into-python-script

"""    
import os, sys
from shutil import copyfile

try: 
    directory = r'D:\Python scripts - FINAL'
    os.chdir(directory)
    
    from _00_glue_utils import raster 
    from _00_glue_utils import lookup_table
    
except ImportError as IE:
    print (IE)
    print ("These functions requires arcpy to run")  
    sys.exit(1)  

finally:
    try: 
        import arcpy 
        
    except ImportError as IE:
        print (IE)
        print ("These functions requires arcpy to run")  
        sys.exit(1)  
   
         
def delete():
    '''
    Description:
        delete processed dem files. 
    '''
    import os
    for root, dirnames, filenames in os.walk(directory):
        for file in range(len(filenames)):
            if 'dem' in filenames[file]:
                os.remove(os.path.join(root,filenames[file]))
     
        
def process_dem(directory):
    '''
    Description:
        process SRTM or ASTER GDEM for: 
            (1) rename to common standardized naming convention                           
            (2) mosaic and rename to standardized naming convention
                              
        Please use process_dem_pcs() function to project renamed \n\
        DEM to 2-D coord ref sys (UTM) for feature extraction

    Args:
        directory (str):
            
    Returns:
        No returns or exchanges.
    '''    
    images = []
    GCSE = arcpy.SpatialReference(4326)
    
    arcpy.env.overwriteOutput = True
    
    # modify for Landsat? SRTM should benefit from this 
    # set nodata limit to avoid issues
    #    arcpy.env.nodata = "MAXIMUM"
    
    for root, dirnames, filenames in os.walk(directory):        
        
        
        # if file exists already, skip 
        if os.path.isfile(root + '\\' + os.path.basename(root) + '_dem_gcs.tif'):        
            print os.path.basename(root) + '_dem_gcs.tif' + ' ...already exists and was not processed.'
               
        else:   
            # detect SRTM-1
            for i in filenames:
                
                # detect ASTER v.2
                if i.startswith('ASTGTM2') and i.endswith('dem.tif'):
                    images.append(i)

                # detect ASTER v.3
                if i.startswith('ASTGTMV003') and i.endswith('dem.tif'):
                    images.append(i)
                        
                # detect SRTM-1
                elif i.endswith('v3.tif'):
                    images.append(i)
    
            if len(images) > 0:
                
                # if mosaick not necessary, copy file and rename to <city>_dem_gcs.tif 
                if len(images) == 1:
                    # if mosaick not necessary, copy file and rename to <city>_dem_gcs.tif
                    copyfile(root + '\\' + images[0], root + '\\' + os.path.basename(root) + '_dem_gcs.tif')                
                    print os.path.basename(root)
                    print images
                    
                # if mosaic necessary, mosaic and rename to <city>_dem_gcs.tif
                elif len(images) > 1:
        
                    left_bracket_remove = str(images).replace('[','')
                    right_bracket_remove = str(left_bracket_remove).replace(']','')
                    quote_remove = str(right_bracket_remove).replace("'",'')
                    add_semi_colon = str(quote_remove).replace(', ',';')
        #            final_image_list = str_quotes(add_semi_colon)
                    
                    print add_semi_colon
                    print os.path.basename(root)
                    print images
                    
                    arcpy.env.workspace = root     
                    
                    # process; 16 bit unsigned (-32,768 to 32,767) for SRTM or ASTER
                    arcpy.MosaicToNewRaster_management(add_semi_colon, root, os.path.basename(root) + '_dem_gcs.tif', 
                                                    GCSE,'16_BIT_SIGNED','', '1', 'MAXIMUM','')
                           
                images = []            


def batch_process_dem():
    '''
    Require to subset data into subfolder in directory; arcpy cannot handle and does not execute; freezes
    '''
    try:        
        directory = ['D:\\GLUE Datasets\\DEM\\africa',
                     'D:\\GLUE Datasets\\DEM\\asia',
                     'D:\\GLUE Datasets\\DEM\\middle east',
                     'D:\\GLUE Datasets\\DEM\\oceania',
                     'D:\\GLUE Datasets\\DEM\\europe',
                     'D:\\GLUE Datasets\\DEM\\south america',
                     'D:\\GLUE Datasets\\DEM\\usa',
                     'D:\\GLUE Datasets\\DEM\\canada']        
        
        for folder in range(len(directory)):
            print directory[folder]
            process_dem(directory[folder])
            
    except Exception as e:
        print e
    
batch_process_dem()



def process_dem_pcs(directory, cities, raster_file):
    '''
    Description:
        Converts SRTM-1 or ASTER GDEM from default WGS 1984 datum (geographic coordinate system using lat/lon) to UTM  coordinate system (projected coordindate system). 
        The logic is set-up to use a hard-coded list with a city name and corresponding Landsat path/row (e.g., ['Canada_ON','Toronto','018030']). This function will
        pull the UTM coordinate system from a Landat image's path/row that corresponds to the location of a city and project accordingly. 
        
    Args:
        directory (str): specify directory of SRTM-1 and/or ASTER GDEM rasters for processing  (e.g., r'D:\GLUE Datasets\SRTM - Working')
        cities (list of lists): specify the hard-coded list of lists with cities and corresponding Landsat row/path  (e.g., cities(r'D:\Python scripts - WORKING'))
        raster_file (str): specify directory of composited Landsat images e.g., r'D:\Scenes_L2 - Processed\image_composites'
        
    Returns:
        No returns or exchanges
        
    '''    
    arcpy.env.overwriteOutput = True
    get_image_spatial_ref = []
    
    
    # step 1. get all image spatial ref and filter to a single entry for each city for processing   
    for root, dirnames, filenames in os.walk(directory):
        
        for file in range(len(filenames)):
            if filenames[file].endswith('dem_gcs.tif'):
                
                for i in range(len(cities)):
                    
                    # standardized for name match 
                    if ' ' in cities[i][1]:
                        cities[i][1] = cities[i][1].replace(' ' , '_')
                        
                    if cities[i][1].lower() in filenames[file].lower():
                        
                        # get Landsat image 
                        for j in range(len(raster_file)): 
                            
                            #splice for code 
                            a = os.path.basename(raster_file[j])[:10]
                            b = os.path.basename(raster_file[j])[16:]
                        
                            raster_code = os.path.basename(raster_file[j].split(a)[-1].split(b)[0])
                            
                            #compare city name to raster code; will only compare to composited image
                            if raster_code in cities[i]:
                                
                                inRas = arcpy.Raster(raster_file[j])
                                
                                # arcpy.ProjectRaster_management DOES NOT ACCEPT str, 
                                # therefore do not use spatialReference.name, but instead spatialReference.factorycode (corresponds to EPSG code)
                                spatial_ref = arcpy.Describe(inRas).spatialReference
                                
                                # specify EPSG equivalent
                                ### If an Esri well-known ID is below 32767, it corresponds to the EPSG ID. WKIDs that are 32767 or above are Esri-defined.
                                print ('%s, %s' % (cities[i][1], spatial_ref.factorycode))
                                
                                get_image_spatial_ref.append((cities[i][1].lower(), spatial_ref.factorycode))
                                
                                # Using set theory, only 1 unique occurrence allowed
                                filtered_image_spatial_ref = sorted(set(get_image_spatial_ref))
                                
    # step 2. for every city + spatial ref in list, compare to dem images, and use that spatial ref for projecting to projected coordinate system (e.g., 2d coord sys)
    for root, dirnames, filenames in os.walk(directory):
        
        for file in range(len(filenames)):
            if filenames[file].endswith('dem_gcs.tif'):
                              
                for i in range(len(filtered_image_spatial_ref)):
                    
                    if filtered_image_spatial_ref[i][0].lower() in filenames[file].lower():                   
                        print filenames[file]

                        # if file exists already, skip 
                        if os.path.isfile(root + '\\' + os.path.basename(root) + '_dem_pcs.tif'):        
                            print os.path.basename(root) + '_dem_pcs.tif' + ' ...already exists and was not processed.'
                            
                        # else process...     
                        else:
#                            print filenames[file], filtered_image_spatial_ref[i][1], root 
                        
                            out_name = filenames[file][:-7] + 'pcs.tif'
                            spatia_ref_final = arcpy.SpatialReference(filtered_image_spatial_ref[i][1])
                            
                            arcpy.env.workspace = root
                            print root, filenames[file], out_name, spatia_ref_final

                            # https://pdfs.semanticscholar.org/e5b9/c673d89ffffd8588a55a8643d33a1005f858.pdf
                            # BILINEAR interpolation OK for continuous data, but paper found CUBIC convolution to have better accuracy (RMSE error)
                            arcpy.ProjectRaster_management(filenames[file], out_name, spatia_ref_final, 'CUBIC') 


def batch_process_dem_to_pcs():
    '''
    Require to subset data into subfolder in directory; arcpy cannot handle and does not execute; freezes
    '''

    try:
        raster_file = raster(r'G:\Landsat_Download','sr_composite.tif')
        cities_lookup = lookup_table(r'D:\Python scripts - FINAL')  
                
        directory = ['D:\\GLUE Datasets\\DEM\\africa',
                     'D:\\GLUE Datasets\\DEM\\asia',
                     'D:\\GLUE Datasets\\DEM\\middle east',
                     'D:\\GLUE Datasets\\DEM\\oceania',
                     'D:\\GLUE Datasets\\DEM\\europe',
                     'D:\\GLUE Datasets\\DEM\\south america',
                     'D:\\GLUE Datasets\\DEM\\usa',
                     'D:\\GLUE Datasets\\DEM\\canada']

##        directory = ['D:\\GLUE Datasets\\DEM\\europe']  
        
        for folder in range(len(directory)):
            print directory[folder]
            
            process_dem_pcs(directory[folder], cities_lookup, raster_file)
            
    except Exception as e:
        print e


if __name__ == '__main__':
    
    batch_process_dem_to_pcs()
 

