# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 15:04:38 2019

@author: Alexander Tong

Developed and tested with Python 2.7.15
"""

# select out revised csv
import os, sys, shutil

try: 
    import arcpy 
    
except ImportError as IE:
    print (IE)
    print ("These functions requires arcpy to run")  
    sys.exit(1)   

def select_revised_csv_for_shp_conversion_and_extraction(parent_dir,copy_to_dir,lst_of_GLUE_csv):
    '''
    Description:
        Select out revised GLUE city csv from Dropbox or wherever current GLUE city \n\
        csv fil(s) are stored for conversion to shapefile and for subsequent feature extraction. 
        
    Args:
        parent_dir (str): specify GLUE city directory with all csv file
                                e.g., r'D:\GLUE Datasets\popMeans_allCities'
        copy_to_dir (str): specify new directory to copy out select GLUE city csv file(s)
                                e.g., r'D:\GLUE Datasets\popMeans_allCities_select'
        lst_of_GLUE_csv (list): specify the GLUE city csv to copy out
                                e.g., ['atlantic_city','charlotte']
                                
    Returns:
        pass
    '''    
    for i in os.listdir(parent_dir):
        
        for j in lst_of_GLUE_csv:
            
            if i[:-4].lower() in j:
                
                if not os.path.exists(os.path.join(copy_to_dir,i)):
                    
                    shutil.copyfile(os.path.join(parent_dir,i), os.path.join(copy_to_dir,i))
                
    #            os.remove(os.path.join(directory,i))


def delete_old_shp_to_be_replaced_by_new_shp_for_GLUE_city(shp_dir_gcs, shp_dir_pcs, lst_of_GLUE_csv):
    '''
    Description:
        Delete old shapefile to be replaced by revised shapefile for GLUE citie(s).
        You must re-run _04_process_csv_to_shp.py to generate the replacement shp(s)  
    
    Args:
        shp_dir_gcs (str): specify directory of shapefiles in GCS 
                           e.g., r'D:\GLUE Results\Process_Transect\gcs' 
        shp_dir_pcs (str): specify directory of shapefiles in PCS 
                           e.g., r'D:\GLUE Results\Process_Transect\pcs'         
        lst_of_GLUE_csv (list): specify the GLUE city csv to copy out
                                e.g., ['atlantic_city','charlotte'] 
                                
                            
    Returns:
        Deleted selected shapefiles
    '''
    arcpy.env.workspace = shp_dir_gcs 
    
    fclist = arcpy.ListFeatureClasses()
    
    for fc in fclist:
        for j in lst_of_GLUE_csv:
            if fc[:-8].lower() in j:
                arcpy.Delete_management(fc)
        
    arcpy.env.workspace = shp_dir_pcs 
    
    fclist = arcpy.ListFeatureClasses()
    
    for fc in fclist:
        for j in lst_of_GLUE_csv:
            if fc[:-8].lower() in j:
                arcpy.Delete_management(fc)
                
                
def delete_ancillary_files_for_csv_outputs(directory,lst_of_GLUE_csv):
    '''
    !!!WARNING: MUST REMOVE ALL ARCMAP GENERATED ANCILLARY FILES OR ELSE ARCPY LOGIC BREAKS!!!
    !!!WARNING: ALSO REMOVE schema.ini AND info folder
    
    Description:
        Remove all traces of csv outputs and ancillary files from the creation of csv outputs \n\
        using Esri Arcmap. This includes the 'info' folder and a schema.ini file. This is \n\
        extremely important as you cannot overwrite outputs without deleting all \n\
        previous traces of the output.
    
    Args:
        directory (str): specify directory to delete ALL files.
                         e.g., r'D:\GLUE Results'   
        lst_of_GLUE_csv (list): specify the GLUE city csv to copy out
                                e.g., ['atlantic_city','charlotte']
            
    Returns:
        pass
    '''       
    for i in os.listdir(directory):    
        
        if i.startswith('Extracted'):  
            
            # go into sub-folder; remove all ancillary files for GLUE city
            for root, dirnames, filenames in os.walk(os.path.join(directory,i)):          
                for files in range(len(filenames)):          
                    
                    if filenames[files]:    
                        
                        for j in lst_of_GLUE_csv:                        
                            if j in filenames[files].lower():           
                                pass
                            
                                print os.path.join( root, filenames[files])
                                os.remove(  os.path.join( root, filenames[files]) )
                    
                    # schema.ini 
                    if filenames[files] == 'schema.ini':
                        
                        print filenames[files]
                        os.remove(  os.path.join( root, filenames[files]) )
                        
                # info 
                if 'info' in root:
                    print root 
                    shutil.rmtree(root)
    
  
                           
## create revised output folders for easier deployment to dropbox ...
#
#def ig_f(dir, files):
#    '''
#    https://stackoverflow.com/questions/15663695/shutil-copytree-without-files
#    '''
#    return [f for f in files if os.path.isfile(os.path.join(dir, f))]
#    
#for i in os.listdir(directory):    
#     
#    if i.startswith('Extracted'):  
#        
#        shutil.copytree(os.path.join(directory,i), os.path.join(directory,i + '_select_revised'), ignore=ig_f)
#        
#        
#
## remove folder
#for i in os.listdir(directory):    
#    
#    if i.endswith('_select_revised'):  
#        print i 
#       # Remove (delete) the directory path. Only works when the directory is empty, 
#       # otherwise, OSError is raised.
##        os.rmdir( os.path.join (directory, i) )
#       
#       # recursively deletes a directory, even if it has contents.
#        shutil.rmtree(os.path.join(directory, i))
#
#        
        
if __name__ == '__main__':
    
    lst_of_GLUE_csv = ['atlantic_city','charlotte','halifax','little_rock','mexico_city',
                  'morelia','ottawa','stockholm','newcastle']
    
    parent_dir = r'D:\GLUE Datasets\popMeans_allCities'
    copy_to_dir = r'D:\GLUE Datasets\popMeans_allCities_select'

#    select_revised_csv_for_shp_conversion_and_extraction(parent_dir,copy_to_dir,lst_of_GLUE_csv)      
     
    
    shp_dir_gcs = r'D:\GLUE Results\Process_Transect\gcs'
    shp_dir_pcs = r'D:\GLUE Results\Process_Transect\pcs'
#    delete_old_shp_to_be_replaced_by_new_shp_for_GLUE_city(shp_dir_gcs, shp_dir_pcs, lst_of_GLUE_csv)
       
    
    directory = r'D:\GLUE Results'   
#    delete_ancillary_files_for_csv_outputs(directory,lst_of_GLUE_csv)
            
            
        
        
        
        
        
        
        
        