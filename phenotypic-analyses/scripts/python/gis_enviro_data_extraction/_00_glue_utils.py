# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 18:32:56 2019

@author: Alexander Tong 

Developed and tested with Python 2.7.15
"""
import os

def raster(directory, ext):
    '''
    Description: retrieve list of raster images for processing
    
        
    Args:
        directory (str): directory of raster files
        
        ext (str): specify raster with wildcard and/or file extension 
                   (e.g., 'pcs.tif',
                          'NDVI.tif',
                          'NDSI.tif',
                          'LST_Celsius.tif'
                          'wc2.0_30s_tavg_')
                                                                          
    Returns:
        list of raster(s) 
        
    '''
#    import os
    
    raster = []
    for root, dirnames, filenames in os.walk(directory):
        for file in range(len(filenames)):
            
            if filenames[file].endswith(ext):
                raster.append(os.path.join(root,filenames[file]))
                
            elif filenames[file].startswith(ext):
                raster.append(os.path.join(root,filenames[file]))
                
    return raster

def raster_for_glue_maps_summer(directory, ext):
    '''
    Description: retrieve list of raster images for processing
    
        
    Args:
        directory (str): directory of raster files
        
        ext (str): specify raster with wildcard and/or file extension 
                (e.g., 'pcs.tif',
                        'NDVI.tif',
                        'NDSI.tif',
                        'LST_Celsius.tif'
                        'wc2.0_30s_tavg_')
                                                                          
    Returns:
        list of raster(s) 
        
    '''
#    import os
    
    raster = []
    for root, dirnames, filenames in os.walk(directory):
        if root.endswith('image_composites_summer'): # testing; change to lst when surfaces are computed
            for file in range(len(filenames)):
                
                if filenames[file].endswith(ext):
                    raster.append(os.path.join(root,filenames[file]))
                    
                elif filenames[file].startswith(ext):
                    raster.append(os.path.join(root,filenames[file]))               
                
    return raster

def raster_for_glue_maps_winter(directory, ext):
    '''
    Description: retrieve list of raster images for processing
    
        
    Args:
        directory (str): directory of raster files
        
        ext (str): specify raster with wildcard and/or file extension 
                (e.g., 'pcs.tif',
                        'NDVI.tif',
                        'NDSI.tif',
                        'LST_Celsius.tif'
                        'wc2.0_30s_tavg_')
                                                                          
    Returns:
        list of raster(s) 
        
    '''
#    import os
    
    raster = []
    for root, dirnames, filenames in os.walk(directory):
        if root.endswith('image_composites_winter'): # testing; change to lst when surfaces are computed
            for file in range(len(filenames)):
                
                if filenames[file].endswith(ext):
                    raster.append(os.path.join(root,filenames[file]))
                    
                elif filenames[file].startswith(ext):
                    raster.append(os.path.join(root,filenames[file]))
                               
    return raster

def raster_for_glue_maps_lst_summer(directory, ext):
    '''
    Description: retrieve list of raster images for processing
    
        
    Args:
        directory (str): directory of raster files
        
        ext (str): specify raster with wildcard and/or file extension 
                (e.g., 'pcs.tif',
                        'NDVI.tif',
                        'NDSI.tif',
                        'LST_Celsius.tif'
                        'wc2.0_30s_tavg_')
                                                                          
    Returns:
        list of raster(s) 
        
    '''
#    import os
    
    raster = []
    for root, dirnames, filenames in os.walk(directory):
        if root.endswith('image_lst_summer'): # testing; change to lst when surfaces are computed
            for file in range(len(filenames)):
                
                if filenames[file].endswith(ext):
                    raster.append(os.path.join(root,filenames[file]))
                    
                elif filenames[file].startswith(ext):
                    raster.append(os.path.join(root,filenames[file]))
                              
    return raster

def raster_for_glue_maps_lst_winter(directory, ext):
    '''
    Description: retrieve list of raster images for processing
    
        
    Args:
        directory (str): directory of raster files
        
        ext (str): specify raster with wildcard and/or file extension 
                (e.g., 'pcs.tif',
                        'NDVI.tif',
                        'NDSI.tif',
                        'LST_Celsius.tif'
                        'wc2.0_30s_tavg_')
                                                                          
    Returns:
        list of raster(s) 
        
    '''
#    import os
    
    raster = []
    for root, dirnames, filenames in os.walk(directory):
        if root.endswith('image_lst_winter'): # testing; change to lst when surfaces are computed
            for file in range(len(filenames)):
                
                if filenames[file].endswith(ext):
                    raster.append(os.path.join(root,filenames[file]))
                    
                elif filenames[file].startswith(ext):
                    raster.append(os.path.join(root,filenames[file]))
                              
    return raster


def lookup_table(directory):
    '''
    Description: 
        retrieve lookup table (lst of lsts) of GLUE city to Landsat code for processing
        
    Args:
        directory (str): specify directory of _00_lookup_table.py  
        
    Returns:
        list of list 
    '''
#    import os
    
    try:
        os.chdir(directory)
        
        from _01_lookup_table import lookup_table

        return lookup_table() 
    
    except:
        print ('could not locate lookup table') 
        

    
def lookup_table_for_glue_maps(directory):
    '''
    Description: 
        retrieve lookup table (lst of lsts) of GLUE city to Landsat code for processing
        
    Args:
        directory (str): specify directory of _00_lookup_table.py  
        
    Returns:
        list of list 
    '''
#    import os
    
    try:
        os.chdir(directory)
        
        from _01_lookup_table_for_glue_maps import lookup_table

        return lookup_table() 
    
    except:
        print ('could not locate lookup table') 
