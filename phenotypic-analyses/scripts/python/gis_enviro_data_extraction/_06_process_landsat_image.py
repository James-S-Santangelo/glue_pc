# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 10:00:17 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15

## SET NO DATA FOR MULTISPECTRAL (MULTI-BAND) IMAGE
https://gis.stackexchange.com/questions/136014/setting-nodata-for-a-multiband-raster 

"""
import os, sys, re 

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
    
    
def atoi(text):
    '''
    https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside 
    '''
    return int(text) if text.isdigit() else text


def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]


def composite_image_recurve(directory, outfolder, *args):
    '''
    Description: this function will recursively go into a directory and composite \n\
    a multispectral image. If composite image already exists, skip. 
    
    #https://gist.github.com/whophil/2a999bcaf0ebfbd6e5c0d213fb38f489
    #https://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
    
    Args:
        directory (str): input directory to be parsed. e.g., 'C:\\Users\\EvoEco\\Desktop\\test_img\\'
        outfolder (str): out directory. Specify 'default' to be saved in same folder of input rasters, else specify outpath
        *args (str): currently supports 3 string arguments: 'band', 'sr', 'toa'. 
        
    Returns:
        No returns 
        
        $future implementation:
            replace arcpy with gdal -> arcpy is extremely temperamental. 
    ''' 
    
#    import os, arcpy  

    
    to_be_processed = []
    
    for arg in args: 
        if arg == 'band':
            to_be_processed.append(arg)
        elif arg == 'sr':
            to_be_processed.append(arg)
        elif arg == 'toa':
            to_be_processed.append(arg) 

    Original_count = 0
    SR_count = 0
    ToA_count = 0
    images = []
    
    for root, dirnames, filenames in os.walk(directory):
        
        # while still in directory, go through files, until exhaust... 
        # out of while loop, process... use a counter 
        Original_count = len(filenames) 
        SR_count = len(filenames) 
        ToA_count = len(filenames) 
        
        try:
            for i in to_be_processed:
                if 'band' in i:
                    print 'processing band'
                    while Original_count > 0:
                        for file in range(len(filenames)): 
                            if 'T1_b' in filenames[file]:
            #                    print filenames[file]
                                images.append(filenames[file])
                                images.sort(key=natural_keys)
                            
                            elif 'RT_b' in filenames[file]:
            #                    print filenames[file]
                                images.append(filenames[file])
                                images.sort(key=natural_keys)
                                
                            #excl. brightness temperature and band quality assessment bands                      
                            for i in range(len(images)):
                                if 'bt' in images[i]:
                                    del images[i]
                                elif 'bqa' in images[i]:
                                    del images[i]
                                
                            Original_count -= 1 
                            
                        left_bracket_remove = str(images).replace('[','')
                        right_bracket_remove = str(left_bracket_remove).replace(']','')
                        quote_remove = str(right_bracket_remove).replace("'",'')
                        add_semi_colon = str(quote_remove).replace(', ',';')
                        
                        images = []
                        
                        #if saving to same folder as input rasters 
                        if 'composite' not in add_semi_colon:
            #               directory = root 
                            arcpy.env.workspace = root 
                            
                            if outfolder == 'default':
                                print root
                                arcpy.CompositeBands_management(add_semi_colon, root + '//' + add_semi_colon[:41] + "composite.tif")
                                print add_semi_colon

                                # get composite image
                                composite_image = arcpy.Raster(root + '//' + add_semi_colon[:41] + "composite.tif")
                                
                                # pixels == 0 when border, else pixel == 1 is image 
                                mask = arcpy.sa.Con(1,0)

                                # set pixel border == 0 as NoData 
                                mask_SetNull = SetNull(mask, mask, "VALUE = 0")

                                # where border is NoData, excl. from final image 
                                Composited_cleaned = mask_SetNull * composite_image

                                #Save
                                arcpy.CopyRaster_management(Composited_cleaned, root + '//' + add_semi_colon[:41] + "composite_cleaned.tif")
                                
                            else:
                                print root
                                arcpy.CompositeBands_management(add_semi_colon, outfolder + '//' + add_semi_colon[:41] + "composite.tif")
                                print add_semi_colon

                                # get composite image
                                composite_image = arcpy.Raster(root + '//' + add_semi_colon[:41] + "composite.tif")
                                
                                # pixels == 0 when border, else pixel == 1 is image 
                                mask = arcpy.sa.Con(1,0)

                                # set pixel border == 0 as NoData 
                                mask_SetNull = SetNull(mask, mask, "VALUE = 0")

                                # where border is NoData, excl. from final image 
                                Composited_cleaned = mask_SetNull * composite_image

                                #Save
                                arcpy.CopyRaster_management(Composited_cleaned, root + '//' + add_semi_colon[:41] + "composite_cleaned.tif")
            #                print root
            #                print add_semi_colon
                        else:
                            print add_semi_colon[:40] + ' is already processed'
    
    #            print root 
    #            print add_semi_colon   
                elif 'sr' in i:
                    print 'processing sr'
                    while SR_count > 0:
                        for file in range(len(filenames)): 
                            if 'sr_band' in filenames[file]:
            #                    print filenames[file]
                                images.append(filenames[file])
                                images.sort(key=natural_keys)
                                
                            #excl. aerosol band and processed NDVI image # ADD TO THIS LIST FOR OTHER DERIVED IMAGE PRODUCTS
                            for i in range(len(images)):
                                if 'aerosol' in images[i]:
                                    del images[i]
                                elif 'ndvi' in images[i]:
                                    del images[i]
                                    
                            SR_count -= 1 
                            
                        left_bracket_remove = str(images).replace('[','')
                        right_bracket_remove = str(left_bracket_remove).replace(']','')
                        quote_remove = str(right_bracket_remove).replace("'",'')
                        add_semi_colon = str(quote_remove).replace(', ',';')
                        
                        images = []
                        
                        if 'composite' not in add_semi_colon:
            #               directory = root 
                            arcpy.env.workspace = root 
        
                            if outfolder == 'default':
                                arcpy.CompositeBands_management(add_semi_colon, add_semi_colon[:44] + "sr_composite.tif")
                                print add_semi_colon                        
                            else:
#                                arcpy.CompositeBands_management(add_semi_colon, outfolder + '//' + add_semi_colon[:44] + "composite.tif")
#                                print add_semi_colon 
                                
                                print root
                                arcpy.CompositeBands_management(add_semi_colon, outfolder + '//' + add_semi_colon[:41] + "sr_composite.tif")
                                print add_semi_colon
                                
#                                #think it is skipping because i did not make raster layers... 
#                                # get composite image
#                                composite_image = arcpy.Raster(root + '//' + add_semi_colon[:41] + "composite.tif")
#                                
#                                # pixels == 0 when border, else pixel == 1 is image 
#                                mask = arcpy.sa.Con(1,0)
#
#                                # set pixel border == 0 as NoData 
#                                mask_SetNull = SetNull(mask, mask, "VALUE = 0")
#
#                                # where border is NoData, excl. from final image 
#                                Composited_cleaned = mask_SetNull * composite_image
#
#                                #Save
#                                arcpy.CopyRaster_management(Composited_cleaned, root + '//' + add_semi_colon[:41] + "composite_cleaned.tif")
#                                
                                
                        else:
                            print add_semi_colon[:40] + ' is already processed'
                  
            #            print root 
            #            print add_semi_colon   
                elif 'toa' in i:
                    print 'processing toa'        
                    while ToA_count > 0:
                        for file in range(len(filenames)): 
                            if 'toa_band' in filenames[file]:
            #                    print filenames[file]
                                images.append(filenames[file])
                                images.sort(key=natural_keys)
                    
                            ToA_count -= 1 
                            
                        left_bracket_remove = str(images).replace('[','')
                        right_bracket_remove = str(left_bracket_remove).replace(']','')
                        quote_remove = str(right_bracket_remove).replace("'",'')
                        add_semi_colon = str(quote_remove).replace(', ',';')
                        
                        images = []
                        
                        if 'composite' not in add_semi_colon:
            #               directory = root 
                            arcpy.env.workspace = root 
        
                            if outfolder == 'default':
                                arcpy.CompositeBands_management(add_semi_colon, add_semi_colon[:45] + "toa_composite.tif")
                                print add_semi_colon
                            else:
                                arcpy.CompositeBands_management(add_semi_colon, outfolder + '//' + add_semi_colon[:45] + "toa_composite.tif")
                                print add_semi_colon
                        else:
                            print add_semi_colon[:40] + ' is already processed'
                #ignore folders with images that have been processed already 
        except:
            pass


#def raster(directory, ext):
#    '''
#    Description: retrieve list of raster images for processing
#            
#    Args:
#        directory (str): directory of raster files
#        ext (str): specify raster file extension (e.g., '.tif')
#        
#    Returns:
#        list of raster(s).       
#    '''
#    import os
#    
#    raster = []
#    for root, dirnames, filenames in os.walk(directory):
#        for file in range(len(filenames)):
#            if filenames[file].endswith(ext):
#                raster.append(filenames[file])
#            
#    return raster
#
#raster_file = raster(r'D:\Scenes_L2 - Processed\image_composites', 'sr_composite.tif')


def splitword(string):  
    '''
    Description:
        splits a string in half of total length for processing (used in process_landsat_composite function)
        
        e.g., '123456'  -->  ('123' , '456')
        
        
    Args:
        string (str): string value to be split in half
            
    Returns:
        tuple containing split string as half of total length 
        
    Source: 
        https://stackoverflow.com/questions/22108306/how-to-split-a-string-into-two-parts    
    '''
    split = -((-len(string))//2)                   
    return string[:split], string[split:]


def process_landsat_composite(directory):
    '''
    Description:
        Selectively mosaic Landsat images based on standard Landsat naming convention. \n\ 
        The logic parses the Landsat path, row and acquisition date for mosaicking such that: \n\
            1) like images share same path and acquisition date and will only  \n\
                process if like images have sequential rows (e.g., adjacent) 
        
        The logic of this function assumes that images with same path and acquisition  \n\
        and sequential rows are at most 3 rows, therefore 3 images (worst case scenario),  \n\
        but most likely 2 rows and therefore 2 images (best case scenario) to be mosaicked. 
        Note that images with different paths are not captured on the same date.
        
        The logic for mosaicking is based on the purpose of the GLUE project;  \n\
        retrieval of images that centre on the area of a city per global basis.  \n\
        Given that the Landsat satellite trajectory follows a pre-defined trajectory  \n\
        with image swaths taken on the basis of a standardized path/row system,   \n\
        some images do not entirely overlap a city area and therefore require an  \n\
        additional image to complete the picture. Additionally, because weather systems  \n\
        are often localized and may obstruct ground features (e.g., city area),  \n\
        having clear and usable images with the same path, acquisition and sequential rows is  \n\
        highly unlikely to be beyond several images (worst case scenario).
        
        Future development should focus on mosaicking Landsat images such that  \n\
        only 2 images can be mosaicked (i.e., at max, 2 images that covers a city of interest);  \n\
        this can be achieved by using a hard-coded list of lists...
        
        ....

        Images mosaicked are given the following naming convention:
            e.g., 
            
            # image 1 + image 2
                LC08_L1TP_018029_20180611_20180615_01_T1_sr_composite.tif + LC08_L1TP_018030_20180611_20180615_01_T1_sr_composite.tif
            
            # image 1 path/row + image2 path/row
                LC08_L1TP_018029_ + 018030 + _20180611_20180615_01_T1_sr_composite.tif
            
            # output (original prefix + path/row + path_row + original suffix)
                LC08_L1TP_018029_018030_20180611_20180615_01_T1_sr_composite.tif 
                     
        Please use <to be made> function to project renamed DEM to 2-D coord ref sys (UTM) for feature extraction 
                  
    Args:
        directory (str): specify directory of composited images for mosaicking. (e.g., r'D:\Scenes_L2 - Processed\image_composites')
        
    Returns:
        No returns or exchanges. 
        
    $ to be implemented: replace arcpy method with gdal_merge utils via sub_process cmd 
    '''
    import os, arcpy
    from itertools import groupby 
    from operator import itemgetter
        
    arcpy.env.overwriteOutput = True

    # modify for Landsat? SRTM should benefit from this 
    # set nodata limit to avoid issues
    #    arcpy.env.nodata = "MAXIMUM"
    
    images = []
    to_mosaic = []
    to_mosaic_2nd_pass = []
    image_count_list = []

    # Recurse through directory 
    # step 0. get images 
    # any image that requires mosaicking will use same root... non-issue here
    for root, dirnames, filenames in os.walk(directory): 
        for file in range(len(filenames)):
            if filenames[file].endswith('sr_composite.tif'):  
#                print filenames[file]
                image_count_list.append(filenames[file]) 
                    
        image_count = len(image_count_list)
        while image_count > 0:
            
            for file in range(len(filenames)):
                if filenames[file].endswith('sr_composite.tif'):  
                    images.append(filenames[file])  
                    
                image_count -= 1
                
    print directory
                               
    # step 1. split landsat file naming convention into individual elements for processing
    split_images = sorted([x.split('_') for x in images])

    # step 2. split row and path into separate elements and replace combined row and path element in list 
    for i in range(len(split_images)):
        row =  splitword(split_images[i][2])[0]  # target row at index [2] 
        path = splitword(split_images[i][2])[1]  # target path
        
        split_images[i].insert(2,row)
        split_images[i].insert(3,path)
        
        del split_images[i][4]
    
    # step 3. group individual images by landsat path and date for processing 
    # first sort by acquisition date; this will set up the grouping 
    # itemgetter does not work as intended, instead sort by the integer value of landsat path and date  
    split_images.sort(key=lambda x: int(x[4])) # VERY IMPORTANT FOR SORTING ELEMENTS AT SAME INDEX THAT REPEAT TOO MANY TIMES 
      
    group_images_by_path_and_date = [list(g) for key,g in groupby(split_images, key = lambda x: (int(x[2]),(int(x[4]))) )]
#    group_images_by_path_and_date = [list(g) for key,g in groupby(split_images, key = itemgetter(2,4))]  # Python itemgetter is finicky and does not work if too many elements are the same in a single position; will not sort with groupby 
    
#    for i in range(len(group_images_by_path_and_date)):
#        if len(group_images_by_path_and_date[i]) > 1:
#            print group_images_by_path_and_date[i]

    # step 4. after sorting and grouping image files, need to re-combine landsat path and row for processing 
    #         whilst also deleting the separated landsat path and row elements in list 
    for i in range(len(group_images_by_path_and_date)):
        for j in range(len(group_images_by_path_and_date[i])):
    #        print test[i][j]
            row =  group_images_by_path_and_date[i][j][2]
            path = group_images_by_path_and_date[i][j][3]
    #        print row, path
    #        print ''.join((row, path))
            
            group_images_by_path_and_date[i][j].insert(2,''.join((row, path)))
            
            # delete intermediate step (landsat path)
            del group_images_by_path_and_date[i][j][3]
    
    
    for i in range(len(group_images_by_path_and_date)):
        for j in range(len(group_images_by_path_and_date[i])):
            
             # delete intermediate step (landsat row)
            del group_images_by_path_and_date[i][j][3]
            
    # step 5. re-combine all individual string elements of split landsat file naming convention   
    for i in range(len(group_images_by_path_and_date)):
        for j in range(len(group_images_by_path_and_date[i])):
            group_images_by_path_and_date[i][j] =  '_'.join(group_images_by_path_and_date[i][j])   
            
    
    # step 6. check image files for sequential landsat row if mosaicking is required (e.g., localizes the area) 
    temp_seq_list = []
    temp_seq_not_processed_list = []
    not_processed = [] 
    for i in range(len(group_images_by_path_and_date)):
        
        len_sublist = len(group_images_by_path_and_date[i])
        
        if len(group_images_by_path_and_date[i]) > 1:
            
            while len_sublist > 0:
                
                for j in range(len(group_images_by_path_and_date[i])):
                    
                    a = group_images_by_path_and_date[i][j][:13] # slice front for Landsat row 
                    b = group_images_by_path_and_date[i][j][16:] # slice back  for Landsat row
                    c = group_images_by_path_and_date[i][j][:17] # slice front for Landsat date
                    d = group_images_by_path_and_date[i][j][25:] # slice back for Landsat date
                    
                    # get Landsat row to determine if sequential 
                    temp_seq_list.append(int(group_images_by_path_and_date[i][j].split(a)[-1].split(b)[0]))
                    temp_seq_not_processed_list.append(int(group_images_by_path_and_date[i][j].split(a)[-1].split(b)[0]))
#                    print group_images_by_path_and_date[i][j].split(a)[-1].split(b)[0]   
#                    print group_images_by_path_and_date[i][j]   

                    len_sublist -= 1
            
   
            # check for sequential landsat rows for mosaicking (e.g., adjacent images); if not, remove values 
            while (sorted(temp_seq_list) == range(min(temp_seq_list), max(temp_seq_list)+1)) != True:
                temp_seq_list.pop()
                
            # check for image files not processed (e.g., some may be sequential but further down in list)
            for row in temp_seq_not_processed_list:
                if row not in temp_seq_list:
                    not_processed.append(row)                             
            
            # step 7.a. using stored sequential landsat rows list (e.g., temp_seq_list), compare with master list of image files 
            #         for processing and append to final list for processing
            len_temp_seq_list = len(temp_seq_list)
            
            while len_temp_seq_list > 0:
                
                # this first conditional is a provisionary measure; may be simplified and removed... 
                # if length of sublist is not same as seq list... 
                if len(group_images_by_path_and_date[i]) != len_temp_seq_list:
                    difference = len(group_images_by_path_and_date[i]) - len_temp_seq_list
#                    print difference
                    for j in range(len(group_images_by_path_and_date[i])-difference):
#                        print group_images_by_path_and_date[i][j]
                        for k in range(len(temp_seq_list)):
                            if str(temp_seq_list[k]) in group_images_by_path_and_date[i][j].split(a)[-1].split(b)[0]:
#                                print str(temp_seq_list[k])
#                                print group_images_by_path_and_date[i][j]
                    
                                to_mosaic.append(group_images_by_path_and_date[i][j])
                                                                            
                                len_temp_seq_list -= 1
                                
                # if length of sublist is same as seq list...                
                elif len(group_images_by_path_and_date[i]) == len_temp_seq_list:                        
                    for j in range(len(group_images_by_path_and_date[i])):
                        for k in temp_seq_list:
                            if str(k) in group_images_by_path_and_date[i][j].split(a)[-1].split(b)[0]:
                                
#                                            print group_images_by_path_and_date[i][j]
                                
                                to_mosaic.append(group_images_by_path_and_date[i][j])
                                
                                len_temp_seq_list -= 1
                                                        
            # step 7.b. using stored not processed list (e.g., not_processed), compare with master list of image files 
            #           for processing and append to list for further processing 
            for l in range(len(group_images_by_path_and_date[i])):
#                print group_images_by_path_and_date[i][l]
                for m in range(len(not_processed)):
                    if str(not_processed[m]) in group_images_by_path_and_date[i][l].split(a)[-1].split(b)[0]:
                        to_mosaic_2nd_pass.append(group_images_by_path_and_date[i][l])
                        print to_mosaic_2nd_pass

            # step 8.a. first pass for process images for mosaicking (e.g., to_mosaic), only if > 1 images 
            if len(to_mosaic) > 1:   
                                                     
                to_mosaic_final = sorted(list(set(to_mosaic)))
                
                left_bracket_remove = str(to_mosaic_final).replace('[','')
                right_bracket_remove = str(left_bracket_remove).replace(']','')
                quote_remove = str(right_bracket_remove).replace("'",'')
                add_semi_colon = str(quote_remove).replace(', ',';')  
                 
                # compute output name 
                output_name_split = []
                
                for i in range(len(to_mosaic)):
       
                    a = to_mosaic[i][:10] # slice front for Landsat path/row 
                    b = to_mosaic[i][16:] # slice back  for Landsat path/row
                    
#                    print to_mosaic[i][:17] + ' ' + to_mosaic[i].split(a)[-1].split(b)[0] 
                    
                    
                    if i == 0: # get first part of name
                        output_name_split.append(to_mosaic[i][:17])
                    elif i > 0: # add path row of subsequent images to be mosaicked
                        output_name_split.append(to_mosaic[i].split(a)[-1].split(b)[0] + '_')
                        
                        if i == len(to_mosaic)-1:
                            output_name_split.append(to_mosaic[i][17:])
                    
                output_name = ''.join(output_name_split)
                
                print output_name
#                print to_mosaic
                print add_semi_colon                        
                print '...'
                
                # if exists, skip 
                if os.path.isfile(root + '\\'+  output_name): 
                    print '{0} already processed... '.format(output_name)
                    
                else:
                    arcpy.env.workspace = root 
                    if 'LC08' in add_semi_colon:
                        arcpy.MosaicToNewRaster_management(add_semi_colon, root, output_name, 
                                                           '','16_BIT_SIGNED','', '7', 'MAXIMUM','')
                        
                    elif 'LE07' in add_semi_colon or 'LT05' in add_semi_colon: 
                        arcpy.MosaicToNewRaster_management(add_semi_colon, root, output_name, 
                                                           '','16_BIT_SIGNED','', '6', 'MAXIMUM','')
                    
            # step 8.b. second pass for processing images for mosaicking (e.g., to_mosaic_2nd_pass)
            if len(to_mosaic_2nd_pass) > 0:      
                
                to_mosaic_final_2nd_pass = sorted(list(set(to_mosaic_2nd_pass)))
                
                to_be_processed_2nd_pass = []
                
                # split remaining landsat rows into lists (e.g., if sequential values found, will be put into own list, else separate list)
                # https://stackoverflow.com/questions/3149440/python-splitting-list-based-on-missing-numbers-in-a-sequence 
                for k, g in groupby(enumerate(not_processed), lambda (i, x): i-x):
                    val =  map(itemgetter(1), g)
                    to_be_processed_2nd_pass.append(val)

                # convert str landsat row to numeric for element matching
                for i in range(len(to_be_processed_2nd_pass)):
                    if len(to_be_processed_2nd_pass[i]) > 1:
                        for j in range(len(to_be_processed_2nd_pass[i])):
                            to_be_processed_2nd_pass[i][j] = str(to_be_processed_2nd_pass[i][j])
#                            print to_be_processed_2nd_pass[i]
                         
                # step 8.b.a process for image files with sequential landsat rows 
                index = 0  # used for testing purposes, ignore
                
                for i in range(len(to_be_processed_2nd_pass)):  
                                                                           
                    for j in range(len(to_be_processed_2nd_pass[i])):                                    
 
                        if isinstance(to_be_processed_2nd_pass[i][j], str) == True:

                            for k in range(len(to_mosaic_final_2nd_pass)):
                                
                                a2 = to_mosaic_final_2nd_pass[k][:13] # slice front for Landsat row 
                                b2 = to_mosaic_final_2nd_pass[k][16:] # slice back  for Landsat row

#                                print to_mosaic_final_2nd_pass[k].split(a2)[-1].split(b2)[0]
#                                print to_be_processed_2nd_pass[j][k]
                                if to_be_processed_2nd_pass[i][j] in to_mosaic_final_2nd_pass[k].split(a2)[-1].split(b2)[0]:
                                   
#                                    print to_mosaic_final_2nd_pass[k]
#                                    print to_be_processed_2nd_pass[i][j]
#                                    print '...'
                                    #why not just insert at index at sublist??, and delete all elemnts with len 1
                                    to_be_processed_2nd_pass[i].append(to_mosaic_final_2nd_pass[k])
#                                    print index 
                                    index += 1  # used for testing purposes, ignore
                                                                              
                    index = 0  # used for testing purposes, ignore
                    
            
            # step 8.b.b only if image files that are not processed, then proceed to mosaicking for 2nd pass    
            output_name_to_be_processed_2nd_pass = [] 
            if len(not_processed) > 0:
                
                for x in range(len(to_be_processed_2nd_pass)):
                    len_split = len(to_be_processed_2nd_pass[x])//2 
                    if len(to_be_processed_2nd_pass[x]) > 1:
                        print to_be_processed_2nd_pass[x][len_split:]
                        
                        left_bracket_remove_2nd_pass = str(to_be_processed_2nd_pass[x][len_split:]).replace('[','')
                        right_bracket_remove_2nd_pass = str(left_bracket_remove_2nd_pass).replace(']','')
                        quote_remove_2nd_pass = str(right_bracket_remove_2nd_pass).replace("'",'')
                        add_semi_colon_2nd_pass = str(quote_remove_2nd_pass).replace(', ',';')
                        
                        output_name_to_be_processed_2nd_pass = to_be_processed_2nd_pass[x][len_split:]
                        
                        # compute output name 
                        output_name_split_2nd_pass = []
                        
                        for i in range(len(output_name_to_be_processed_2nd_pass)):
               
                            a = output_name_to_be_processed_2nd_pass[i][:10] # slice front for Landsat path/row 
                            b = output_name_to_be_processed_2nd_pass[i][16:] # slice back  for Landsat path/row
                            
#                            print output_name_to_be_processed_2nd_pass[i][:17] + ' ' + output_name_to_be_processed_2nd_pass[i].split(a)[-1].split(b)[0] 
                            
                            if i == 0:
                                output_name_split_2nd_pass.append(output_name_to_be_processed_2nd_pass[i][:17])
                            elif i > 0:
                                output_name_split_2nd_pass.append(output_name_to_be_processed_2nd_pass[i].split(a)[-1].split(b)[0] + '_')
                                
                                if i == len(output_name_to_be_processed_2nd_pass)-1:
                                    output_name_split_2nd_pass.append(output_name_to_be_processed_2nd_pass[i][17:])
                            
                        output_name_2nd_pass = ''.join(output_name_split_2nd_pass)
                        
                        print output_name_2nd_pass
#                        
                        print 'not to be left behind...'
                        print add_semi_colon_2nd_pass
                        print '......'               
                        
                        # if exists, skip 
                        if os.path.isfile(root + '\\'+  output_name_2nd_pass): 
                            print '{0} already processed... '.format(output_name_2nd_pass)
                            
                        else:
                            arcpy.env.workspace = root
    #                        # landsat 8
                            if 'LC08' in add_semi_colon_2nd_pass:
                                arcpy.MosaicToNewRaster_management(add_semi_colon_2nd_pass, root, output_name_2nd_pass, 
                                                '','16_BIT_SIGNED','', '7', 'MAXIMUM','')
                            # landsat 5/7    
                            elif 'LE07' in add_semi_colon or 'LT05' in add_semi_colon:
                                arcpy.MosaicToNewRaster_management(add_semi_colon_2nd_pass, root, output_name_2nd_pass, 
                                                '','16_BIT_SIGNED','', '6', 'MAXIMUM','')   
                                                     
            image_count_list = []            
            temp_seq_list = []
            temp_seq_not_processed_list = []
            to_mosaic = []
            to_mosaic_2nd_pass = []
            not_processed = []                        


def batch_process_mosaic():
#    import os 
#    directory = [r'G:\Landsat_Download\asia\image_composites_summer',
#                 r'G:\Landsat_Download\asia\image_composites_winter',
#                 r'G:\Landsat_Download\oceania\image_composites_summer',
#                 r'G:\Landsat_Download\oceania\image_composites_winter',
#                 r'G:\Landsat_Download\canada\image_composites_summer',
#                 r'G:\Landsat_Download\canada\image_composites_winter',
#                 r'G:\Landsat_Download\usa\image_composites_summer',
#                 r'G:\Landsat_Download\usa\image_composites_winter',
#                 r'G:\Landsat_Download\europe_et_al\image_composites_summer',
#                 r'G:\Landsat_Download\europe_et_al\image_composites_winter',      
#                 r'G:\Landsat_Download\south_america\image_composites_summer',
#                 r'G:\Landsat_Download\south_america\image_composites_winter'] 
        
    directory = [r'G:\Landsat_Download\usa\image_composites_summer',
                 r'G:\Landsat_Download\usa\image_composites_winter',] 
    
    for folder in directory:
        process_landsat_composite(folder)
         
#        # need to rename 'composite' to 'sr_composite'; accidentally did not fix output name for 'sr' image composite exe
#        for folder in directory:
#            for root, dirnames, filenames in os.walk(folder):
##            process_landsat_composite(folder)
#                
#                for file in range(len(filenames)):
#                    if 'sr_composite' in filenames[file]:
#                        print filenames[file] + ' exists!'
#                    else:
#                        os.rename(os.path.join(root,filenames[file]), os.path.join(root,filenames[file].replace('composite', 'sr_composite')))
#                    
#batch_process_mosaic()


def batch_process_composite():
    '''
    Description:
        Require to subset data into subfolder in directory; arcpy cannot handle too many files at once; \\
        does not execute and freezes
    '''
    try:  
        # asia
#        directory = [r'G:\Landsat_Download\asia\espa-alexander.tong@mail.utoronto.ca-10292018-105823-178', 
#                     r'G:\Landsat_Download\asia\image_composites_summer',
#                     r'G:\Landsat_Download\asia\espa-alexander.tong@mail.utoronto.ca-10292018-105905-194',
#                     r'G:\Landsat_Download\asia\image_composites_winter']
        # oceania
#        directory = [r'G:\Landsat_Download\oceania\espa-alexander.tong@mail.utoronto.ca-10292018-113819-365', 
#                     r'G:\Landsat_Download\oceania\image_composites_summer',
#                     r'G:\Landsat_Download\oceania\espa-alexander.tong@mail.utoronto.ca-10292018-114004-367',
#                     r'G:\Landsat_Download\oceania\image_composites_winter']
        
        # canada and usa 
##        directory = [r'G:\Landsat_Download\canada\espa-alexander.tong@mail.utoronto.ca-10292018-105600-988', 
##                     r'G:\Landsat_Download\canada\image_composites_summer',
##                     r'G:\Landsat_Download\canada\espa-alexander.tong@mail.utoronto.ca-10292018-105724-706',
##                     r'G:\Landsat_Download\canada\image_composites_winter',
##                     r'G:\Landsat_Download\usa\espa-alexander.tong@mail.utoronto.ca-10292018-113443-256', 
##                     r'G:\Landsat_Download\usa\image_composites_summer',
##                     r'G:\Landsat_Download\usa\espa-alexander.tong@mail.utoronto.ca-10292018-113734-294',
##                     r'G:\Landsat_Download\usa\image_composites_winter'] 
         
         # addendum canada winter scene
#        directory = [r'G:\Landsat_Download\canada\espa-alexander.tong@mail.utoronto.ca-11092018-101009-387',
#                     r'G:\Landsat_Download\canada\image_composites_winter']
                 
#        # europe_et_al
#        directory = [r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-10292018-110638-866',
#                     r'G:\Landsat_Download\europe_et_al\image_composites_summer',
#                     r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-10292018-113302-004',
#                     r'G:\Landsat_Download\europe_et_al\image_composites_winter']       
#            
       # addendum europe_et_al winter scenes  
#        directory = [r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-11072018-104859-769',
#                     r'G:\Landsat_Download\europe_et_al\image_composites_winter',
#                     r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-11072018-122905-326',
#                     r'G:\Landsat_Download\europe_et_al\image_composites_winter']
 
      # south_america
#        directory = [r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-10292018-120915-536',
#                     r'G:\Landsat_Download\south_america\image_composites_summer',
#                     r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-10292018-120942-303',
#                     r'G:\Landsat_Download\south_america\image_composites_winter']     
      
       # south_america addendum 
#        directory = [r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-11192018-131412-079',
#                     r'G:\Landsat_Download\south_america\image_composites_summer',
#                     r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-11192018-131512-240',
#                     r'G:\Landsat_Download\south_america\image_composites_winter']          
       
       
        # south_america addendum 2
#        directory = [r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-03042019-161814-357',
#                     r'G:\Landsat_Download\south_america\image_composites_summer',
#                     r'G:\Landsat_Download\south_america\espa-alexander.tong@mail.utoronto.ca-03042019-161848-728',
#                     r'G:\Landsat_Download\south_america\image_composites_winter']    
      
        # oceania addendum
#        directory = [r'G:\Landsat_Download\oceania\espa-alexander.tong@mail.utoronto.ca-03042019-161727-973', 
#                     r'G:\Landsat_Download\oceania\image_composites_summer',
#                     r'G:\Landsat_Download\oceania\espa-alexander.tong@mail.utoronto.ca-03042019-161752-773',
#                     r'G:\Landsat_Download\oceania\image_composites_winter']
        
#        # europe_et_al addedum 2
##        directory = [r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-03192019-134034-964',
##                     r'G:\Landsat_Download\europe_et_al\image_composites_summer',
##                     r'G:\Landsat_Download\europe_et_al\espa-alexander.tong@mail.utoronto.ca-03192019-134103-601',
##                     r'G:\Landsat_Download\europe_et_al\image_composites_winter']

        # Fairbanks
##        directory = [r'G:\Landsat_Download\usa\espa-james.santangelo37@gmail.com-07232020-082513-144',
##                     r'G:\Landsat_Download\usa\image_composites_summer',
##                     r'G:\Landsat_Download\usa\espa-james.santangelo37@gmail.com-07232020-082414-856',
##                     r'G:\Landsat_Download\usa\image_composites_winter']

        # Cobourg addendum
##        directory = [r'G:\Landsat_Download\canada\espa-cobourg-addendum',
##                     r'G:\Landsat_Download\canada\image_composites_summer']

        # Anchorage addendum
        directory = [r'G:\Landsat_Download\usa\espa-Anchorage-summer-addendum',
                     r'G:\Landsat_Download\usa\image_composites_summer',
                     r'G:\Landsat_Download\usa\espa-Anchorage-winter-addendum',
                     r'G:\Landsat_Download\usa\image_composites_winter']
        
        for folder in range(len(directory)):
            if folder%2 == 0:
#                print directory[folder]
                images = directory[folder]
            else:                
                outfolder = directory[folder] 
                
                print images, outfolder
            
                composite_image_recurve(images,outfolder,'sr')
            
    except Exception as e:
        print e
    
batch_process_composite()

#
#if __name__ == "__main__":
#
##    directory = r'D:\validate_landsat_lst_toronto\test\LC080180302013071501T1 - SC20180919132357'
##    outfolder = r'D:\validate_landsat_lst_toronto\test\image_composites'
##    
##    directory = r'G:\Landsat_Download\asia\espa-alexander.tong@mail.utoronto.ca-10292018-105823-178'
##    outfolder = r'G:\Landsat_Download\asia\image_composites_summer' 
#
#    composite_image_recurve(directory,outfolder,'sr')                


