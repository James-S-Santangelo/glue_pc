# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 11:12:10 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15
"""
import os, sys, re

# include this and arcpy licensing ad nauseum 
if sys.version_info[0] != 2:
    print("This script requires Python version 2.xx")
    sys.exit(1)
    
try:
    import arcpy
    
except Exception as e:
    print (e)
    print ("This script requires arcpy to run")


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def global_aridity_rescale(folder, outpath):
    '''
    Description: 
        Rescale Aridity Index (AI) product values from CGIAR. Raster values represent \n\
        annual average aridity over the 1950-2000 period. AI is modelled using the data \n\
        available from the WorldClim Global Climate Data (Hijmans et al. 2005) as \n\
        input parameters. The WorldClim, based on a high number of climate observations \n\
        and SRTM topographical data, is a high-resolution global geo-database \n\
        (30 arc seconds or ~ 1km at equator) of monthly average data (1950-2000) \n\
        for the following climatic parameters: precipitation, mean, minimum and maximum \n\
        temperature. This set of parameters is insufficient to fully parameterize physical \n\
        radiation-based PET equations (i.e. the FAO-PM), though can parameterize \n\
        simpler temperature-based PET equations. \n\
        
        From Read Me Doc:
        
        The Aridity Index values reported within the Global-Aridity geodataset \n\
        have been multiplied by a factor of 10,000 to derive and distribute the \n\
        data as integers (with 4 decimal accuracy). This multiplier has been used \n\
        to increase the precision of the variable values without using decimals \n\
        (real or floating values are less efficient in terms of computing time \n\
        and space compared to integer values). Global-Aridity values need to be \n\
        multiplied for 0.0001 to retrieve the values in the correct units.
            
    Dependencies:
        arcpy 
        
    Args:
        folder (str): specify folder of aridity index raster for processing (ai_yr).
        outpath (str): specify outpath location for rescaled aridity index raster.
        
    Returns:
        No returns or exchanges.
    '''
    arcpy.env.workspace = folder 
   
    RasList = arcpy.ListRasters()
    Ras_to_be_processed = []
    
    for Ras in RasList:
        Ras_to_be_processed.append(Ras)
    
    
    for i in Ras_to_be_processed:
        if 'ai_yr' in Ras_to_be_processed:
            aridity_index_rescale = arcpy.Raster(Ras_to_be_processed[0])*0.0001
        
        #save
        arcpy.CopyRaster_management(aridity_index_rescale, os.path.join(outpath, str(i) + '_rescaled.tif'))


def convert_adf_to_tif_format(root_directory):
    '''
    Description:
        Convert CGIAR datasets .adf to .tif format for processing. Ersi native grid format is useless.
        
    **NOTE: outpath folder is hard-coded for 3 different output folders in this function. Change as necessary. 
        
    **WARNING: WHEN DEALING WITH .adf, you do not specify the raster; you only specify the folder it is contained in... 
                e.g., r'D:\GLUE Datasets\PET_he_annual\PET_he_annual\pet_he_yr' 
                
                     instead of:       
                         
                     arcpy.env.workspace = r'D:\GLUE Datasets\PET_he_annual\PET_he_annual\pet_he_yr'
                     
                     RasList = arcpy.ListRasters()
                     Ras_to_be_processed = []
            
                     for Ras in RasList:
                         Ras_to_be_processed.append(Ras) 
                    
                     Arcpy.Raster(Ras_to_be_processed[0])
                     
                     
    Args:
        root_directory (str): specify root (parent) directory of all CGIAR rasters to be processed 
        
    Returns:
        No returns or exchanges.
    '''
    
    import os 
    import arcpy
        
    PET_HE_MONTHLY_directory = []
    PET_HE_YR_directory = []
    AI_directory = []
    

    # PET_HE_MONTHLY
    for root, dirs, files in os.walk(root_directory):
        if 'PET_he_monthly' in root:
            PET_HE_MONTHLY_directory.append(root)  
            
    PET_dir = []
    
    PET_dir_count = 12
    
    # get only pet_he_int subdirs 
    while PET_dir_count > 0:
        
        # use ele at index 0 to get root folder for recursive parsing
        for root, dirs, files in os.walk(PET_HE_MONTHLY_directory[0]): 
            if 'pet_he_' + str(PET_dir_count) in root: 
                PET_dir.append(root)
                PET_dir.sort(key=natural_keys)
                PET_dir_count -= 1
    
    for folder in PET_dir:
        
        Ras_to_be_processed = []
        
        arcpy.env.workspace = folder
        RasList = arcpy.ListRasters()
        
        for Ras in RasList:
            Ras_to_be_processed.append(Ras) 
                
            if os.path.isfile(r'D:\GLUE Datasets\PET_he_monthly\PET_he_monthly_tif' + '\\' + Ras_to_be_processed[0] + '.tif'):
                print ('{0} has already been converted to tif!'.format(Ras_to_be_processed[0])) 
                    
            else:
                outpath = r'D:\GLUE Datasets\PET_he_monthly\PET_he_monthly_tif'
                arcpy.RasterToOtherFormat_conversion(folder, outpath, 'TIFF') 


    # PET_HE_YR
    for root, dirs, files in os.walk(root_directory):
        if 'pet_he_yr' in root:
            PET_HE_YR_directory.append(root)        
               
    arcpy.env.workspace = PET_HE_YR_directory[0]        
    RasList = arcpy.ListRasters()
    Ras_to_be_processed = [] 
    
    RasList = arcpy.ListRasters()
    
    for Ras in RasList:
        Ras_to_be_processed.append(Ras)     
            
    if os.path.isfile(r'D:\GLUE Datasets\PET_he_annual\PET_he_annual_tif' + '\\' + Ras_to_be_processed[0] + '.tif'):
        print ('{0} has already been converted to tif!'.format(Ras_to_be_processed[0]))    
        
    else:
        outpath = r'D:\GLUE Datasets\PET_he_annual\PET_he_annual_tif'
        arcpy.RasterToOtherFormat_conversion(PET_HE_YR_directory[0], outpath, 'TIFF') 


    # AI_YR
    for root, dirs, files in os.walk(root_directory):
        if root.endswith('ai_yr'):
            AI_directory.append(root)
            
    arcpy.env.workspace = AI_directory[0]   
    RasList = arcpy.ListRasters()
    Ras_to_be_processed = [] 
    
    RasList = arcpy.ListRasters()        
    
    for Ras in RasList:
        Ras_to_be_processed.append(Ras)    
    
    if os.path.isfile(r'D:\GLUE Datasets\AI_annual\AI_annual_tif' + '\\' + Ras_to_be_processed[0] + '.tif'):
        print ('{0} has already been converted to tif!'.format(Ras_to_be_processed[0]))      
           
    else: 
        outpath = r'D:\GLUE Datasets\AI_annual\AI_annual_tif'
        arcpy.RasterToOtherFormat_conversion(AI_directory[0], outpath, 'TIFF') 
        
        
if __name__ == "__main__":    
    
#    folder = r'D:\GLUE Datasets\AI_annual\AI_annual\ai_yr'
#    outpath = r'D:\GLUE Datasets\AI_annual\AI_annual\ai_yr_rescaled'
    
    folder = r'D:\GLUE Datasets\AI_annual\AI_annual_tif'
    outpath = r'D:\GLUE Datasets\AI_annual\AI_annual_tif'
    global_aridity_rescale(folder, outpath)

    root_directory = r'D:\GLUE Datasets'
    convert_adf_to_tif_format(root_directory)  

