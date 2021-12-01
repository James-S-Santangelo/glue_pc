# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 16:01:53 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15
"""
import os, sys

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
    

        
def index(directory, outpath, index):  
    '''
    ## removed argument outpath_unconstrained
    # not necessary at the moment
    
   **WARNING** logic is dependent on default Landsat naming convention and must be strictly adhered to.
   
    Description: 
        This function will recursively go into a directory and compute a <spectral index> product\n\
        from composited surface reflectance Landsat datasets. 
    
    arcpy is forcing creation of temporary layers in memory to process. 
    
            Normalized Difference Vegetation Index
            NDVI = NIR - Red / NIR + Red
 
            Landsat 5/7:(Band 4 - Band 3)/(Band 4 + Band 3)
            Landsat 8: (Band 5 - Band 4)/(Band 5 + Band 4) 
            
            
            Normalized Difference Snow Index
            NDSI = green - SWIR / green + SWIR

            Landsat 5/7:(Band 2 - Band 5/(Band 2 + Band 5) 
            Landsat 8: (Band 3 - Band 6/(Band 3 + Band 6) 
    
    
            Normalized Difference Built-up Index
            NDBI = NIR - SWIR / NIR + SWIR
            
            Landsat 5/7: (Band 5 - Band 4)/(Band 5 + Band 4) 
            Landsat 8: (Band 6 - Band 5)/(Band 6 + Band 5) 
    
    Args:
        directory (str): input directory to be parsed. e.g., C:\\Users\\EvoEco\\Desktop\\test_img\\'
        outpath (str): specify outpath. e.g., 'G:\\Landsat_Download\\canada\\image_composites_winter_ndvi' 
        index (str): specify spectral index to be calculated. e.g., (NDVI, NDSI)
            
    Returns:
        No returns 
    
    $ to be implemented: optimize code; reduce amount of lines    
    $ to be implemented: if outpath does not exist, create folder, else nothing
    
    '''
#    arcpy.CheckOutExtension("spatial")
#    directory = r'G:\Landsat_Download'
    arcpy.env.overwriteOutput = True
    
    for root, dirnames, filenames in os.walk(directory):
        for file in range(len(filenames)):
            if filenames[file].endswith('sr_composite.tif'):
                
                if index == 'NDVI':
                    
                    if os.path.isfile(os.path.join(outpath, filenames[file][:41] + '_NDVI' + '.tif')):
                        print filenames[file][:41] + '_NDVI' + '.tif' + ' already exists!'
                        
                    else: 
                        # Landsat 5 and 7
                        if 'LE07' in filenames[file] or 'LT05' in filenames[file]:
                            print filenames[file]
    
                            arcpy.env.workspace = root 
                            print filenames[file]
                            Band_3 = arcpy.Raster(filenames[file] + '\Band_3')*0.0001  # apply scale factor
                            Band_4 = arcpy.Raster(filenames[file] + '\Band_4')*0.0001  # apply scale facto
                            
                            arcpy.MakeRasterLayer_management(Band_3,'Red_out')
                            arcpy.MakeRasterLayer_management(Band_4,'NIR_out') 
                        
                        # Landsat 8
                        elif 'LC08' in filenames[file]:
                            
                            arcpy.env.workspace = root 
                            print filenames[file]
                            Band_4 = arcpy.Raster(filenames[file] + '\Band_4')*0.0001  # apply scale factor
                            Band_5 = arcpy.Raster(filenames[file] + '\Band_5')*0.0001  # apply scale factor
                            
                            arcpy.MakeRasterLayer_management(Band_4,'Red_out')
                            arcpy.MakeRasterLayer_management(Band_5,'NIR_out') 
                    
                    
                        Num = arcpy.sa.Float(Raster('NIR_out') - Raster('Red_out'))
                        Denom = arcpy.sa.Float(Raster('NIR_out') + Raster('Red_out'))
                        arcpy.MakeRasterLayer_management(Num,'Num_out')
                        arcpy.MakeRasterLayer_management(Denom,'Denom_out') 
                        
                        NDVI = arcpy.sa.Divide('Num_out', 'Denom_out')
                        arcpy.MakeRasterLayer_management(NDVI,'NDVI_out')                         
    
                        # constrained -1 to 1 
                        NDVI_rescaled = SetNull('NDVI_out', 'NDVI_out','VALUE < -1 OR VALUE > 1' )
        
                        # fill artifacts with NoData; only true for isolated areas; else large NoData areas are likely cloud, ice/snow/water, edge of water
                        neighborhood = NbrRectangle(3, 3, "CELL")
                        infill_ndvi = arcpy.sa.Con(IsNull(NDVI_rescaled), FocalStatistics(NDVI_rescaled, neighborhood, "MEAN"), NDVI_rescaled)
                        
                        OUTPUT_NDVI = os.path.join(outpath, filenames[file][:41] + '_NDVI' + '.tif')
                        
                        arcpy.CopyRaster_management(infill_ndvi, OUTPUT_NDVI)
                                           
                        arcpy.Delete_management('NIR_out')
                        arcpy.Delete_management('Red_out')
                        arcpy.Delete_management('Num_out')    
                        arcpy.Delete_management('Denom_out')
                        
                
                if index == 'NDSI':
                    if os.path.isfile(os.path.join(outpath, filenames[file][:41] + '_NDSI' + '.tif')):
                        print filenames[file][:41] + '_NDSI' + '.tif' + ' already exists!'
                        
                    else: 
                        # Landsat 5 and 7
                        if 'LE07' in filenames[file] or 'LT05' in filenames[file]:
                            
                            arcpy.env.workspace = root 
                            print filenames[file]
                            Band_2 = arcpy.Raster(filenames[file] + '\Band_2') 
                            Band_5 = arcpy.Raster(filenames[file] + '\Band_5') 
                            
                            arcpy.MakeRasterLayer_management(Band_2,'Green_out')
                            arcpy.MakeRasterLayer_management(Band_5,'SWIR_out') 
                            
                        # Landsat 8  
                        elif 'LC08' in filenames[file]:
                            
                            arcpy.env.workspace = root 
                            print filenames[file]
                            Band_3 = arcpy.Raster(filenames[file] + '\Band_3') 
                            Band_6 = arcpy.Raster(filenames[file] + '\Band_6') 
                            
                            arcpy.MakeRasterLayer_management(Band_3,'Green_out')
                            arcpy.MakeRasterLayer_management(Band_6,'SWIR_out') 
                                
                            
                        Num = arcpy.sa.Float(Raster('Green_out') - Raster('SWIR_out'))
                        Denom = arcpy.sa.Float(Raster('Green_out') + Raster('SWIR_out'))
                        arcpy.MakeRasterLayer_management(Num,'Num_out')
                        arcpy.MakeRasterLayer_management(Denom,'Denom_out') 
                        
                        NDSI = arcpy.sa.Divide('Num_out', 'Denom_out')
                        arcpy.MakeRasterLayer_management(NDSI,'NDSI_out')                         
    
                        # fix below: add moving window to fill voids 
                        # constrained -1 to 1  
                        NDSI_rescaled = SetNull('NDSI_out', 'NDSI_out','VALUE < -1 OR VALUE > 1' )                    
                        
                        neighborhood = NbrRectangle(3, 3, "CELL")
                        infill_ndsi = arcpy.sa.Con(IsNull(NDSI_rescaled), FocalStatistics(NDSI_rescaled, neighborhood, "MEAN"), NDSI_rescaled)
                        
                        OUTPUT_NDSI = os.path.join(outpath, filenames[file][:41] + '_NDSI' + '.tif')
                        arcpy.CopyRaster_management(infill_ndsi, OUTPUT_NDSI)
                                                       
                        arcpy.Delete_management('NIR_out')
                        arcpy.Delete_management('Red_out')
                        arcpy.Delete_management('Num_out')    
                        arcpy.Delete_management('Denom_out')
                      
                    
def batch_process_index():
    '''
    Description:
        Require to subset data into subfolder in directory; arcpy cannot handle too many files at once; \\
        does not execute and freezes
    '''
    try:  
        ## below for NDVI
        # asia
#        directory = [r'G:\Landsat_Download\asia\image_composites_summer', 
#                     r'G:\Landsat_Download\asia\image_composites_summer_ndvi',
#                     r'G:\Landsat_Download\asia\image_composites_winter',
#                     r'G:\Landsat_Download\asia\image_composites_winter_ndvi']
        
        # oceania
#        directory = [r'G:\Landsat_Download\oceania\image_composites_summer', 
#                     r'G:\Landsat_Download\oceania\image_composites_summer_ndvi',
#                     r'G:\Landsat_Download\oceania\image_composites_winter',
#                     r'G:\Landsat_Download\oceania\image_composites_winter_ndvi']
        
        # canada and usa 
        directory = [
###                     r'G:\Landsat_Download\canada\image_composites_summer', 
###                     r'G:\Landsat_Download\canada\image_composites_summer_ndvi',
###                     r'G:\Landsat_Download\canada\image_composites_winter',
###                     r'G:\Landsat_Download\canada\image_composites_winter_ndvi',
##                     r'G:\Landsat_Download\usa\image_composites_summer', 
##                     r'G:\Landsat_Download\usa\image_composites_summer_ndvi',
                     r'G:\Landsat_Download\usa\image_composites_winter',
                     r'G:\Landsat_Download\usa\image_composites_winter_ndvi'] 
                 
#        # europe_et_al
#        directory = [r'G:\Landsat_Download\europe_et_al\image_composites_summer',
#                     r'G:\Landsat_Download\europe_et_al\image_composites_summer_ndvi',
#                     r'G:\Landsat_Download\europe_et_al\image_composites_winter',
#                     r'G:\Landsat_Download\europe_et_al\image_composites_winter_ndvi']       
            
      # south_america
#        directory = [r'G:\Landsat_Download\south_america\image_composites_summer',
#                     r'G:\Landsat_Download\south_america\image_composites_summer_ndvi',
#                     r'G:\Landsat_Download\south_america\image_composites_winter',
#                     r'G:\Landsat_Download\south_america\image_composites_winter_ndvi']    
      
        for folder in range(len(directory)):
#            
            if folder%2 == 0:
                print directory[folder]
                images = directory[folder]
                
            else:                
                outfolder = directory[folder] 
                
                print images, outfolder
            
                index(images, outfolder, 'NDVI')

      
        ## below for NDSI
        # asia
#        directory = [r'G:\Landsat_Download\asia\image_composites_winter',
#                     r'G:\Landsat_Download\asia\image_composites_winter_ndsi']
        # oceania
##        directory = [ r'G:\Landsat_Download\oceania\image_composites_winter',
##                     r'G:\Landsat_Download\oceania\image_composites_winter_ndsi']
        
        # canada and usa 
##        directory = [
##                     r'G:\Landsat_Download\canada\image_composites_winter',
##                     r'G:\Landsat_Download\canada\image_composites_winter_ndsi',
##                     r'G:\Landsat_Download\usa\image_composites_winter',
##                     r'G:\Landsat_Download\usa\image_composites_winter_ndsi'] 
                 
        # europe_et_al
#        directory = [r'G:\Landsat_Download\europe_et_al\image_composites_winter',
#                     r'G:\Landsat_Download\europe_et_al\image_composites_winter_ndsi']       
            
      # south_america
#        directory = [r'G:\Landsat_Download\south_america\image_composites_winter',
#                     r'G:\Landsat_Download\south_america\image_composites_winter_ndsi']    
        
                
##        for folder in range(len(directory)):
##            
##            if folder%2 == 0:
###                print directory[folder]
##                images = directory[folder]
##                
##            else:                
##                outfolder = directory[folder] 
##                
##                print images, outfolder
##            
##                index(images, outfolder, 'NDSI')
###            
    except Exception as e:
        print e
    
batch_process_index()

                    

#if __name__ == "__main__":                
                        
    #directory = 'C:\\Users\\EvoEco\\Desktop\\test_img\\'    
    #directory = 'C:\\Users\\EvoEco\\Desktop\\Scenes_SR\\test_ndvi_from_sr\\'
    #directory = 'C:\\Users\\Alex\\Desktop\\test_ndvi\\'
    #
    #directory = r'C:\Users\EvoEco\Desktop\Scenes_SR\test_ndvi_from_sr'
    
    ### FOR NDVI
    #directory = r'E:\Scenes_L2\image_composites' 
    ##outpath = r'E:\Scenes_L2\ndvi' 
    #outpath_unconstrained = r'E:\Scenes_L2\ndvi_unconstrained' 
    #outpath_constrained =  r'E:\Scenes_L2\ndvi_constrained'
    #
    #
    #NDVI(directory,outpath_constrained,outpath_unconstrained,8,'NDVI')
    
    
    ### FOR NDBI
    #directory = r'E:\Scenes_L2\image_composites' 
    ##outpath = r'E:\Scenes_L2\ndvi' 
    #outpath_unconstrained = r'E:\Scenes_L2\ndbi_unconstrained' 
    #outpath_constrained =  r'E:\Scenes_L2\ndbi_constrained'
    #
    #
    #NDVI(directory,outpath_unconstrained,outpath_constrained,8,'NDBI')


    ### FOR IBI
    #directory = r'C:\Users\EvoEco\Desktop\test_cal\image_composites' 
    #outpath = r'E:\Scenes_L2\ndvi' 
    #outpath_unconstrained = r'C:\Users\EvoEco\Desktop\test_cal\ndbi_binary' 
    #outpath_constrained =  r'C:\Users\EvoEco\Desktop\test_cal\ndbi_binary' 
    
    
#    index(directory,outpath_unconstrained,outpath_constrained,8,'NDBI')


#-----------------------------------------------------------------------------#
### FOR LST VALIDATION; 20 TORONTO SCENES
#-----------------------------------------------------------------------------#
#    directory = r'D:\validate_landsat_lst_toronto\espa-alexander.tong@mail.utoronto.ca-09182018-093502-904\image_composites'
##    outpath = r'D:\image_ndvi' 
##    outpath_unconstrained = r'C:\Users\EvoEco\Desktop\test_cal\ndbi_binary' 
#    outpath_constrained = r'D:\validate_landsat_lst_toronto\espa-alexander.tong@mail.utoronto.ca-09182018-093502-904\image_ndvi'
#    
#    
#    index(directory,outpath_constrained,8,'NDVI')







