# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 15:24:24 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15

Many more method calls available in ArcPro 

This workflow was tested and developed using Python 2.7.15 and ArcMap 10.6.1
 
**WARNING**
To use Data Driven Pages, you must manually enable the option in the ArcMap GUI \n\
Menu ... Customize > Toolbars > Data Driven Pages 

# need to use a mask layer with page definition to select polygon
- copy the so-called index layer twice (i.e., your shapefile being used) and call \n\
  the duplicate 'mask'
- using Data Driven Pages, select your original index layer 
- right-click properties on 'mask' > Definition Query > Page Definition... > 
  > select same attribute field as original index layer > Do not match 
  
  Source: https://support.esri.com/en/technical-article/000011376 

SQL query to filter for select countries... 



Data Driven Pages Workflow... 
1. Use forloop to activate data driven pages
2. get field row from shapefile you want to associate for the data driven pages for parsing (e.g., will display that row in the DataFrame)
3. set custom extents for each field row 
4. using city_to_landsat_code list, match to the field row of shapefile 
5. maintain the current page mxd.dataDrivenPages.currentPageID = pageNum so it does not move to each field row value
6. using transect point shapefiles, match to the city_to_landsat_code list and get the extent of the transect point shapefiles
7. when city of transect point shapefile == city of city_to_landsat_code list, then match the landsat image row/path and get the extent of the landsat image
8. now determine if the city of transect point shapefile extent falls within the extent of the landsat 
9. if they match, create the city of transect point shapefile and landsat image as layers 
10. set symbology for city of transect point shapefile and add layer to the respective dataframe in the map template for parsing
11. set symbology for landat image and add layer to the respective dataframe in the map template for parsing
12. to set a zoom level of the transect + landsat image, create buffer around transect, add symbology to make it transparent, add as a layer and delete 
13. need to add shapefile as point of interest (POI) for world map inset here
14. if the city of transect point shapefile matches the shapefile as POI for world map inset, add shapefile as layer to world map, then begin to set symbology:
    - set dynamic text elements
    - set symbology for legend
    - export map
    - remove layers from dataframe
    - delete shapefiles, temporary layers
    - refresh table of contents and dataframes... 


"""
try: 
    import os, sys, arcpy
    
    directory = r'D:\Python scripts - FINAL'
    os.chdir(directory)
    
    from _00_glue_utils import raster_for_glue_maps_summer 
    from _00_glue_utils import raster_for_glue_maps_winter 
    from _00_glue_utils import lookup_table_for_glue_maps
    
except ImportError as IE:
    print (IE)
    print ("These functions requires arcpy to run")  
    sys.exit(1)    

# include this and arcpy licensing ad nauseum 
if sys.version_info[0] != 2:
    print("This script requires Python version 2.xx")
    sys.exit(1)
    
     
# step 1. select out countries for GLUE cities 
def select_countries():
    """
    Description: Create shapefile of only select countries and cities that are within shapefile with all countries and cities. 
    
    Args: hard-coded within function 

    """
    import arcpy
    
    in_dir = r'D:\GLUE Datasets\Natural_Earth\ne_10m_admin_0_countries'
    out_dir = r'D:\GLUE Results\GLUE_maps\select_countries'
    
    # COPY OUT ROW BY SELECT COUNTRIES 
    where_clause =  """"ADMIN"  IN ('Argentina', 'Australia', 'Belgium', 'Bolivia', 'Brazil', 'Canada', 'Chile', 'China', 'Colombia', 
                 'Ecuador', 'Finland', 'France', 'Germany', 'Greece', 'Iran', 'Italy', 'Japan', 'Mexico', 
                 'Netherlands', 'New Zealand', 'Norway', 'Papua New Guinea', 'Poland', 'Portugal', 'Russia', 'South Africa', 
                 'Spain', 'Sweden', 'Switzerland', 'United Kingdom', 'United States of America')""" 
    arcpy.Select_analysis(in_dir + '\\' + 'ne_10m_admin_0_countries.shp', out_dir + '\\' + 'select_countries.shp', where_clause)


#select_countries()

# step 2. specify custom extents for countries; behaviour of default data driven pages extent is inconsistent. 
def get_country_extent():
    '''
     Description:
         uses default extent for most countries, but always sets custom extents \n\
         for countries where tehe default extent is not representative of the geographical \n\
         area. 
    '''
    import arcpy 
    
    count = 0 
    country_extent = []
    
       
    arcpy.env.workspace = r'D:\GLUE Datasets\Natural_Earth\select_countries'
    datasets = arcpy.ListFeatureClasses()
    with arcpy.da.SearchCursor(datasets[0], ['SHAPE@','ADMIN']) as cursor:
        for row in cursor:
            extent = row[0].extent
          
            west = extent.XMin
            south = extent.YMin
            east = extent.XMax
            north = extent.YMax  
            
            # specify custom extent 
            
            if row[1] == 'Chile':
    
                west = -95.958628711 
                south = -63.583789545 
                east = -35.663086824 
                north = -12.436652538 
                
                country_extent.append([row[1],north,east,south,west])
                print row[1], north, east, south, west 
            
            elif row[1] == 'Bolivia':
    
                west = -89.735845329 
                south = -40.552596703 
                east = -29.440303441
                north = 10.594540305  
    
                country_extent.append([row[1],north,east,south,west])
                print row[1], north, east, south, west 
        
            elif row[1] == 'Argentina':
    
                west = -95.958628711 
                south = -63.583789545 
                east = -35.663086824 
                north = -12.436652538 
                
                country_extent.append([row[1],north,east,south,west])
                print row[1], north, east, south, west 
                
            elif row[1] == 'China':
    
                west = 59.35383622 
                south = -3.813839007 
                east = 149.797149051 
                north = 72.906866504 
    
                country_extent.append([row[1],north,east,south,west])
                print row[1], north, east, south, west             
            
            elif row[1] == 'France':
    
                west = -10.102177231
                south = 36.437326647
                east = 13.537705666
                north =  56.490423303
    
                country_extent.append([row[1],north,east,south,west])
                print row[1], north, east, south, west 
    
            elif row[1] == 'South Africa':
                
                west = -15.177652879 
                south = -50.355953875 
                east = 60.191774481 
                north =  13.577967385  
                
                country_extent.append([row[1],north,east,south,west])
                count += 1 
                print row[1], north, east, south, west 
           
            elif row[1] == 'Brazil':
                
                west = -89.735845329 
                south = -40.552596703 
                east = -29.440303441
                north = 10.594540305   
                
                country_extent.append([row[1],north,east,south,west])
                count += 1 
                print row[1], north, east, south, west             
            
            elif row[1] == 'Russia':
                
#                west = 13.822938491
#                south = -12.946966758 
#                east = 179.635678682
#                north =  127.707660013  
                
                west = -5.563033337
                south = -66.546588497 
                east = 178.144833154
                north = 89.286246973 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1 
                print row[1], north, east, south, west 
                
            elif row[1] == 'Germany':
                
                west = -2.851696721 
                south = 40.988034197 
                east = 21.266520034 
                north =  61.446889 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1 
                print row[1], north, east, south, west                     
            
            elif row[1] == 'Norway':
                
                west = 2.389784802 
                south = 53.520252439 
                east = 40.25544563 
                north =  85.640705819 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west  
            
            elif row[1] == 'Sweden':
                
                west = 1.883745844 
                south = 49.454898651 
                east = 32.031516788 
                north =  75.028467155  
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west  
    
            elif row[1] == 'Finland':
                
                west = 10.401094417 
                south = 53.430189011 
                east = 40.548865361 
                north =  79.003757515 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west             
            
            elif row[1] == 'Belgium':
                
                west = -4.195096904 
                south = 43.412718885 
                east = 13.893565662 
                north =  58.756859988 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west   
             
            elif row[1] == 'Spain':
                
                west = -22.336728039 
                south = 22.135224202 
                east = 15.347985641 
                north =  54.102184832 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west              
                
            elif row[1] == 'Poland':
                
                west = 3.326223045 
                south = 40.520918661 
                east = 33.473993989 
                north =  66.094487165 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west        
            
            elif row[1] == 'United Kingdom':
                
                west = -16.223266159
                south = 43.951432338
                east = 7.894950596 
                north = 64.410287141 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west 
            
            elif row[1] == 'Greece':
                
                west = 15.122595096 
                south = 31.589589288 
                east = 33.211257663 
                north = 46.933730391 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west 
            
            elif row[1] == 'Italy':
                
                west = -0.246624867
                south = 32.343401006 
                east = 23.871591888
                north = 52.802255809 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west             
            
            elif row[1] == 'Switzerland':
                
                west = -0.837967318 
                south = 39.21421247 
                east = 17.250695249 
                north = 54.558353572 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west              
            
            elif row[1] == 'Iran':
                
                west = 35.133310067  
                south = 17.06507072
                east = 72.818023747 
                north = 49.03203135 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1   
                print row[1], north, east, south, west 
                
            elif row[1] == 'Netherlands':
                
                west = -4.133273924  
                south = 43.27282228 
                east = 13.955388643 
                north = 58.616963382 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1    
                print row[1], north, east, south, west 
             
            elif row[1] == 'Ecuador':
                
                west = -94.274730185 
                south = -11.465979534 
                east = -56.590016505   
                north = 20.500981096 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1    
                print row[1], north, east, south, west        
                
            elif row[1] == 'Colombia':
                
                west = -94.274730185 
                south = -11.465979534 
                east = -56.590016505   
                north = 20.500981096 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1    
                print row[1], north, east, south, west           
                
            elif row[1] == 'Portugal':
                
                west = -32.109404764
                south = 22.572467989 
                east = 10.097474558 
                north = 58.375463894 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1    
                print row[1], north, east, south, west                 
    
            elif row[1] == 'Canada':
                
#                west = -151.956149309 
#                south = 12.693329888 
#                east = -46.438951006 
#                north = 102.200819652 
                
                west = -152.054184636 
                south = -4.826926291 
                east = -41.567986521 
                north = 88.894577645 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1 
                
            elif row[1] == 'United States of America':
                
#                west = -167.631577235 
#                south = -0.02829203 
#                east = -60.629406842 
#                north =  90.738860861 
                
                west = -173.816144904 
                south = -11.426297544 
                east = -61.685647936 
                north =  83.690006504 
                
                country_extent.append([row[1],north,east,south,west])
                count += 1 
                print row[1], north, east, south, west 
                          
            elif row[1] == 'Papua New Guinea':
                
                west = 132.211843757 
                south = -16.573188732 
                east = 159.344837607 
                north = 6.443022921  
                
                country_extent.append([row[1],north,east,south,west])
                count += 1 
                print row[1], north, east, south, west 
                            
            elif row[1] == 'Australia':
                
                west = 96.66558671
                south = -57.075587446
                east = 171.017725967
                north =  5.995394825
                
                country_extent.append([row[1],north,east,south,west])
                count += 1 
                print row[1], north, east, south, west 
                
            elif row[1] == 'New Zealand':
                
                west = 130.54680753
                south = -62.945036883 
                east = 180.177692786 
                north =  -20.846296916
                
                country_extent.append([row[1],north,east,south,west])
                count += 1    
                print row[1], north, east, south, west  
                
            elif row[1] == 'Japan':
                
                west = 112.952540954
                south = 19.950018671
                east = 158.17419737
                north =  58.310371427
                
                country_extent.append([row[1],north,east,south,west])
                count += 1    
                print row[1], north, east, south, west    
                
            else:
                country_extent.append([row[1],north,east,south,west])
                count += 1             
                print row[1], north, east, south, west 
            
    return country_extent


def get_shapefiles(directory):
    '''
    Description:
        pass
    Args:
        directory (str): specify directory of shapefiles 
    Returns:
        list of shapefiles 
    '''
    import os
    
    fclist = []
    
    for i in os.listdir(directory):
        if i.endswith('.shp'):
            fclist.append(os.path.join(directory,i))        
        
    return fclist 


def refine_cities(cities_list):
    '''
    Description:
        refine country name to match attribute table country name from world map 
        
    refine country name
    refine city name... 
    '''
    for i in range(len(cities_list)):
        
        # refine country names (if applicable)
        if '_' in cities_list[i][0]:
            target_country = cities_list[i][0].split('_')[0]
            
            if target_country == 'USA':
                target_country = 'United States of America'
                cities_list[i][0] = cities_list[i][0].replace(cities_list[i][0],target_country)
                
            else:
                cities_list[i][0] = cities_list[i][0].replace(cities_list[i][0],target_country)
                
        elif cities_list[i][0] == 'UK':
            cities_list[i][0] = cities_list[i][0].replace(cities_list[i][0],'United Kingdom')
                
        # replace whitespace with underscores 
        if ' ' in cities_list[i][2]:  
            cities_list[i][2] = cities_list[i][2].replace(' ','_')

        # replace hyphens with underscores 
        elif '-' in cities_list[i][2]:  
            cities_list[i][2] = cities_list[i][2].replace('-','_')
            
        # replace periods with underscores     
        elif '.' in cities_list[i][2]:  
            cities_list[i][2] = cities_list[i][2].replace('.','')
            cities_list[i][2] = cities_list[i][2].replace(' ','_')        
            
    return cities_list


'''
    !!! TO USE CORRECTLY, MUST ENSURE INPUT, OUTPUT DIRECTORIES ARE CONSISTENT 
    
    Description:
        Using Data Driven Pages (DDP) we can automate map production. 
        This workflow processes ONLY natural color composite Landsat images (summer/winter) \n\
        This workflow uses the point shapefiles for each GLUE city 
        (e.g., r'D:\GLUE Results\Process_Transect\pcs')
    
        For summer images, set the following arguments:
            directory = r'D:\GLUE Results\Process_Transect\pcs'
            raster_file = raster_for_glue_maps_summer(r'G:\Landsat_Download','tif')    
            outpath = r'D:\GLUE Results\GLUE_maps\output_summer'
            processed_maps_dir = r'D:\GLUE Results\GLUE_maps\output_summer'
            
        For winter images, set the following arguments:
            directory = r'D:\GLUE Results\Process_Transect\pcs'
            raster_file = raster_for_glue_maps_winter(r'G:\Landsat_Download','tif')
            outpath = r'D:\GLUE Results\GLUE_maps\output_winter'
            processed_maps_dir = r'D:\GLUE Results\GLUE_maps\output_winter'        
    
        Use _17_process_filter_outputs.py to filter out exact Landsat image used for each GLUE city.
    
'''
    
# get shapefiles
directory = r'D:\GLUE Results\Process_Transect\pcs'
#directory = r'D:\GLUE Results\Process_Transect\pcs_select'
shapefiles = get_shapefiles(directory)


cities_list = lookup_table_for_glue_maps(r'D:\Python scripts - FINAL')
cities_list_refined = refine_cities(cities_list)


raster_file = raster_for_glue_maps_summer(r'G:\Landsat_Download','tif')
#raster_file = raster_for_glue_maps_winter(r'G:\Landsat_Download','tif')


### !!!WARNING: this shapefile is based off of the Natural Earth Dataset originally \n\
### found in 'D:\GLUE Datasets\Natural_Earth', but copied over to 'D:\GLUE Results\Process_City_Location'
### this shapefile is processed from _05_process_select_landsat_image.py and the string name
### of glue cities are modified for string matching (e.g., replace spaces with underscores)
world_cities = r'D:\GLUE Results\GLUE_maps\glue_cities\select_cities_mod_for_str_matching_ddp.shp'

outpath = r'D:\GLUE Results\GLUE_maps\output_summer'
#outpath = r'D:\GLUE Results\GLUE_maps\output_winter'
#outpath = r'D:\GLUE Results\GLUE_maps\output_select_summer'
#outpath = r'D:\GLUE Results\GLUE_maps\output_select_winter'


country_extent = get_country_extent()

            

arcpy.env.overwriteOutput = True

# get mxd file 
mxd =  arcpy.mapping.MapDocument(r'D:\GLUE Results\GLUE_maps\map_template_final.mxd')


#-------------------------------DataFrames Start------------------------------#
# get dataframes 
#dataframes = arcpy.mapping.ListDataFrames(mxd)
#df = []
#for frame in range(len(dataframes)):
#    df.append(dataframes[frame])

# get dataframes   
df1 = arcpy.mapping.ListDataFrames(mxd, "Country_DataFrame")[0]
df2 = arcpy.mapping.ListDataFrames(mxd, "Coerce_Project")[0]
df3 = arcpy.mapping.ListDataFrames(mxd, "Get_Raster_PCS")[0]
legend = arcpy.mapping.ListLayoutElements(mxd,"LEGEND_ELEMENT")[0]
Legend_Style = arcpy.mapping.ListStyleItems(r'C:\Users\EvoEcoLab\AppData\Roaming\ESRI\Desktop10.6\ArcMap\EvoEcoLab.style', "Legend Items", "Legend_GLUE_Project")[0]
#--------------------------------DataFrames End-------------------------------#


#-----------------------------RGB Display Start-------------------------------#
# Modify display gun channels for Landsat Images (RGB) to Natural Colour
LC08_RGB = arcpy.mapping.Layer(r'D:\GLUE Results\GLUE_maps\LC08_RGB.lyr')
LE07_RGB = arcpy.mapping.Layer(r'D:\GLUE Results\GLUE_maps\LE07_RGB.lyr')
LT05_RGB = arcpy.mapping.Layer(r'D:\GLUE Results\GLUE_maps\LT05_RGB.lyr')
LC08_RGB_Gamma_2 = arcpy.mapping.Layer(r'D:\GLUE Results\GLUE_maps\LC08_RGB_Gamma_2.lyr')
LC08_RGB_Gamma_3 = arcpy.mapping.Layer(r'D:\GLUE Results\GLUE_maps\LC08_RGB_Gamma_3.lyr')
LE07_RGB_Gamma_2 = arcpy.mapping.Layer(r'D:\GLUE Results\GLUE_maps\LE07_RGB_Gamma_2.lyr')
LT05_RGB_Gamma_2 = arcpy.mapping.Layer(r'D:\GLUE Results\GLUE_maps\LT05_RGB_Gamma_2.lyr')

# select Landsat images to default gamma stretch of 1 
LC08_RGB_Gamma_1_cities = ['Angus','Barrie','Bradford','Charlottetown','Cobourg','Elmira','Everett',
                           'Halifax','Helsinki','Moncton','New_Tecumseth','North_Bay','Port_Hope',
                           'Rochester','Whitchurch_Stouffville']

# select Landsat images to gamma stretch of 3
LC08_RGB_Gamma_3_cities = ['Athens','Atlanta','Charlotte','Christchurch','Curitiba','Bern','Bogota','Glasgow',
                           'Kyoto','Loja','Louisville','Medellin','Mexico_City','Minneapolis',
                           'Palmerston_North','Paris','Quito','Raleigh','Reading','Sapporo','Seattle',
                           'Sioux_Falls','Tacoma','Uruapan','Zurich',]

#------------------------------RGB Display End---------------------------------#
processed_glue_maps = []
for i in os.listdir(outpath):
    if i.endswith('.png'):
        processed_glue_maps.append(os.path.basename(i)[:-30])
            
# start parsing pages
for pageNum in range(1,mxd.dataDrivenPages.pageCount + 1): 
    print pageNum
    mxd.dataDrivenPages.currentPageID = pageNum
    
    # **WARNING; set the layer you want using the data driven pages tool and then this will be USED TO DRIVE THE PAGES
    # "ADMIN" is field containing country names in the attribute table 
    fieldValue = mxd.dataDrivenPages.pageRow.getValue("ADMIN") 

#    try: 
    # set custom extents for each country using dataframe extent via country_extent()
    for country in range(len(country_extent)): 
        if country_extent[country][0] == fieldValue: 
            
            # get country extent
            df1.extent = arcpy.Extent(country_extent[country][4],country_extent[country][3], country_extent[country][2], country_extent[country][1])    
    
            ## df2
            ### deploy logic to match getValue("ADMIN"); country to country match at index cities_to_landsat_code[0]
            ## match country name from master list to country name from world shapefile
            # get country to country name 
            for i in range(len(cities_list_refined)):
                if cities_list_refined[i][0] == str(fieldValue):
                    
                    # maintain pageNum until forloop end for each country
                    mxd.dataDrivenPages.currentPageID = pageNum
                    
                    # get shapefile (transect) to city name match 
                    for j in range(len(shapefiles)):
                        
                        # !!! checks if image already processed based on GLUE city name; 
                        # if multiple outputs for each GLUE city, need remove all images for that GLUE city
                        if os.path.basename(shapefiles[j][:-8]) in processed_glue_maps:
                            pass
                        
                        else:
                            # [:-11] if _buffer.shp; [:-8] if _pcs.shp or _gcs.shp
                            if os.path.basename(shapefiles[j][:-8]) == cities_list_refined[i][2]:
                                                            
                                # grab spatial ref and extent info of transect points 
                                fc_describe = arcpy.Describe(shapefiles[j])
                                fc_spatial_ref = fc_describe.spatialReference.name
                                fc_extent = fc_describe.extent
                                
                                # modify to only match with 'lst.tif' for instance; 
                                # get raster code match from city name match that matches shapefile
                                for k in range(len(raster_file)): 
                                        
                                    # splice for code 
                                    a = os.path.basename(raster_file[k])[:10] 
                                    b = os.path.basename(raster_file[k])[16:] 
                                    
                                    raster_code = os.path.basename(raster_file[k].split(a)[-1].split(b)[0])
                    
                                    # compare city name to raster code 
                                    if raster_code in cities_list_refined[i][3]:
    #                                    print fieldValue, os.path.basename(shapefiles[j][:-11]), cities_list_refined[i][2]
    
                                        # grab spatial ref and extent info of raster 
                                        ras_describe = arcpy.Describe(raster_file[k])
                                        ras_spatial_ref = ras_describe.spatialReference.name
                                        ras_extent = ras_describe.extent
                                    
                                        print os.path.basename(raster_file[k])
                                        print os.path.basename(shapefiles[j][:-8]), raster_code                                        
                                        print fc_spatial_ref, ras_spatial_ref
    
                                        
                                        ### once shapefile matches Landsat image, we can then start adding stylized elements to the map template ###
                                        # do disjoint to find proper match... (e.g., if they fall within same extent)
                                        if not fc_extent.disjoint(ras_extent):
                                            
                                            transect_layer = arcpy.mapping.Layer(shapefiles[j]) # set transect vector as layer                             
                                            image_layer = arcpy.mapping.Layer(raster_file[k])  # set Landsat raster as a layer
    
                                            # apply common symbol to all transect points 
                                            arcpy.ApplySymbologyFromLayer_management (transect_layer, os.path.join(r'D:\GLUE Results\GLUE_maps', 'transect_points_symbol.lyr'))
                                            # add transect
                                            arcpy.mapping.AddLayer(df2, transect_layer, "AUTO_ARRANGE")
                                            
    
                                            # 1. update RGB symbology of Landsat Raster and add Landsat image
                                            # 2. apply image enhancement if image looks too dark (e.g., gamma set 2.0 FOR ALL AS DEFAULT)
                                            if 'LC08' in os.path.basename(raster_file[k]):
                                                # gamma stretch of 3 
                                                if os.path.basename(shapefiles[j][:-8]) in LC08_RGB_Gamma_3_cities:
                                                    arcpy.mapping.UpdateLayer(df2, image_layer, LC08_RGB_Gamma_3, True) 
                                                # gamma stretch of 1 (default)
                                                elif os.path.basename(shapefiles[j][:-8]) in LC08_RGB_Gamma_1_cities:
                                                    arcpy.mapping.UpdateLayer(df2, image_layer, LC08_RGB, True) 
                                                else:
                                                    arcpy.mapping.UpdateLayer(df2, image_layer, LC08_RGB_Gamma_2, True) 
                                                       
                                            elif 'LE07' in os.path.basename(raster_file[k]):
                                                arcpy.mapping.UpdateLayer(df2, image_layer, LE07_RGB_Gamma_2, True) 
                                                
                                            elif 'LT05' in os.path.basename(raster_file[k]):
                                                arcpy.mapping.UpdateLayer(df2, image_layer, LT05_RGB_Gamma_2, True) 
    
                                            # add Landsat image
                                            arcpy.mapping.AddLayer(df2, image_layer, "AUTO_ARRANGE") 
                                            
                                            # use a temporary buffer layer as the extent instead of transect extent for zoom  
                                            TRANSECT_SHAPEFILE = shapefiles[j]
                                            TEMP_TRANSECT_SHAPEFILE = os.path.basename(shapefiles[j][:-4]) + '_layer'
                                            TEMP_TRANSECT_SHAPEFILE_BUFFER = os.path.basename(shapefiles[j][:-4]) + '_buffer_layer'
                                            buffer_zone = 10000
                                            buffer_directory = r'D:\GLUE Results\GLUE_maps\temp_buffer'
                                            
                                            arcpy.MakeFeatureLayer_management(TRANSECT_SHAPEFILE, TEMP_TRANSECT_SHAPEFILE)
                                            
                                            arcpy.Buffer_analysis(TEMP_TRANSECT_SHAPEFILE, os.path.join(buffer_directory, TEMP_TRANSECT_SHAPEFILE_BUFFER), buffer_zone) 
                                            
                                            buffer_layer = arcpy.mapping.Layer(os.path.join(buffer_directory, TEMP_TRANSECT_SHAPEFILE_BUFFER + '.shp'))
                                            
                                            
                                            arcpy.ApplySymbologyFromLayer_management (buffer_layer, os.path.join(r'D:\GLUE Results\GLUE_maps', 'buffer_hide.lyr'))
                                            
                                            arcpy.mapping.AddLayer(df2, buffer_layer, "BOTTOM") 
                                            
                                            # change extent to buffer (10 km) zoom level for transect 
                                            df2.extent = buffer_layer.getExtent()
                                                                                             
                                            # add Landsat image for coordinate system text display in map
                                            arcpy.mapping.AddLayer(df3, image_layer, "BOTTOM") 
            
                                            # add shapefile as point of interest to world map for each country, for each city
                                            with arcpy.da.SearchCursor(world_cities, ["OID@",'nameascii','adm0name','adm1name']) as cursor:
                                                for row in cursor:
                                                    city_oid = []
                                                    
                                                    # get OBJECTID 
                                                    # if space in name, replace
                                                    if ' ' in row[1]:
                                                        city_name = row[1].replace(' ', '_')
                                                        city_oid.append(row[0])
                #                                        if os.path.basename(shapefiles[j][:-11]) in city_name:
                #                                        print city_name
                
                                                    # get OBJECTID 
                                                    else:
                                                        city_name = row[1]
                                                        city_oid.append(row[0])
                #                                        if os.path.basename(shapefiles[j][:-11]) in city_name:
                #                                            print city_name
                
                                                    # once you get OBJECTID, match city to city name
                                                    # **WARNING: may want to replace '==' with 'in'... city names should match however... 
                                                    # **Amendment: DO NOT USE 'in' as charlotte vs. charlottetownl instead have to rely on perfect string match 
                                                    if os.path.basename(shapefiles[j][:-8]).lower() == city_name.lower():
    
                                                        #  set up SQL-eque query to get OBJECTID (OID) match
                                                        # https://gis.stackexchange.com/questions/213930/creating-expression-by-using-object-id-in-arcpy
                                                        oid_fieldname = arcpy.ListFields(world_cities,"","OID")[0].name 
                                                        where_clause = "{:s} = {:d}".format(oid_fieldname,row[0])
                                                        
                                                        # https://gis.stackexchange.com/questions/145188/getting-index-of-active-data-frame-in-arcpy
            #                                            for index, item in enumerate(arcpy.mapping.ListDataFrames(mxd)):
            #                                                print (index, item)
            #                                                if item.name == 'Country_DataFrame':
            
                                                        # make temporary layer to pass as layer into dataframe 
                                                        arcpy.MakeFeatureLayer_management(world_cities,'select_city', where_clause)
                                                        
                                                        
                                                        #-------------df1 symbology start-----------#
                                                        layer = arcpy.mapping.Layer("select_city")
                                                        
            #                                            layer.labelClasses[0].expression = LabelExpression
            #                                            layer.showLabels = True
            
                                                        # change to paper pin symbol 
                                                        arcpy.ApplySymbologyFromLayer_management (layer, os.path.join(r'D:\GLUE Results\GLUE_maps', 'select_cities_symbol_5.lyr'))
                                                        
                                                        arcpy.mapping.AddLayer(df1, layer)
                                                        
                                                        # re-initialize country extent
                                                        df1.extent = arcpy.Extent(country_extent[country][4],country_extent[country][3], country_extent[country][2], country_extent[country][1])    
                                                        #-------------df1 symbology end-----------#
                                                        
    
                                                        #-------dynamic text elements start-------#
                                                        # dynamic city, country names                                 
                                                        # https://gis.stackexchange.com/questions/129682/having-title-which-includes-dynamic-text-of-layer-name 
                                                        Title_Element = arcpy.mapping.ListLayoutElements(mxd, "TEXT_ELEMENT", 'Title')[0]
                                                        
                                                        if os.path.basename(shapefiles[j][:-8]).replace('_', ' ') == 'Whitchurch Stouffville':
                                                            Title_Element.text ='Whitchurch-Stouffville' + ', ' + cities_list_refined[i][1] + ', ' + cities_list_refined[i][0]
                                                            
                                                        elif os.path.basename(shapefiles[j][:-8]) == 'Halle':
                                                            Title_Element.text = 'Halle (Saale)' + ', ' + cities_list_refined[i][1] + ', ' + cities_list_refined[i][0]
                                                            
                                                        elif os.path.basename(shapefiles[j][:-8]).replace('_', ' ') == 'Portland OR':
                                                            Title_Element.text = 'Portland' + ', ' + cities_list_refined[i][1] + ', ' + cities_list_refined[i][0]
                                                            
                                                        elif os.path.basename(shapefiles[j][:-8]).replace('_', ' ') == 'Portland ME':
                                                            Title_Element.text = 'Portland' + ', ' + cities_list_refined[i][1] + ', ' + cities_list_refined[i][0]    
                                                            
                                                        elif os.path.basename(shapefiles[j][:-8]) == 'St_Louis':
                                                            Title_Element.text = 'St. Louis' + ', ' + cities_list_refined[i][1] + ', ' + cities_list_refined[i][0]  
                                                            
                                                        elif os.path.basename(shapefiles[j][:-8]) == 'St_Albert':
                                                            Title_Element.text = 'St. Albert' + ', ' + cities_list_refined[i][1] + ', ' + cities_list_refined[i][0]    
                                                            
                                                        elif os.path.basename(shapefiles[j][:-8]) == 'St_Thomas':
                                                            Title_Element.text = 'St. Thomas' + ', ' + cities_list_refined[i][1] + ', ' + cities_list_refined[i][0]  
                                                            
                                                        elif os.path.basename(shapefiles[j][:-8]) == 'St_John':
                                                            Title_Element.text = 'St. John' + ', ' + cities_list_refined[i][1] + ', ' + cities_list_refined[i][0]  
    
                                                        elif os.path.basename(shapefiles[j][:-8]) == 'St_Johns':
                                                            Title_Element.text = "St. John's" + ', ' + cities_list_refined[i][1] + ', ' + cities_list_refined[i][0]  
                                                            
                                                        else: 
                                                            Title_Element.text = os.path.basename(shapefiles[j][:-8]).replace('_', ' ') + ', ' + cities_list_refined[i][1] + ', ' + cities_list_refined[i][0]   
                                                        
                                                        # dynamic UTM coord for Landsat images
                                                        # when the dataframe refreshes, because it is specified with utm coord to begin with \n\
                                                        # it will not refresh to another UTM zone; may be fixed in ArcPro??
                                                        CRS_Element = arcpy.mapping.ListLayoutElements(mxd, "TEXT_ELEMENT", 'CRS')[0]
                                                        CRS_Element.text = ('Coordinate System: {0} \nProjection: Transverse Mercator \nDatum: WGS 1984 \nUnits: Metre'.format(str(ras_spatial_ref).replace('_', ' '))) 
                                                       
                                                        # dynamic image and image source
                                                        # text box code: <dyn type="document" property="name"/>"
                                                        # Element Name: Image Source
                                                        # http://desktop.arcgis.com/en/arcmap/10.3/map/page-layouts/working-with-dynamic-text.htm
                                                        Data_Source_Element = arcpy.mapping.ListLayoutElements(mxd, "TEXT_ELEMENT", 'Image Source')[0] 
                                                        
                                                        # fix for toronto
                                                        if os.path.basename(raster_file[k]) == 'LC08_L1TP_018029_018030_20180611_20180615_01_T1_sr_composite.tif':
                                                            Data_Source_Element.text = ('Landsat Image: LC08_L1TP_018029_20180611_20180615_01_T1; LC08_L1TP_018030_20180611_20180615_01_T1 \nSource: U.S. Geological Survey')
                                                            
                                                        # fix for new york 
                                                        elif os.path.basename(raster_file[k]) == 'LC08_L1TP_014031_014032_014033_20170730_20170811_01_T1_sr_composite.tif':
                                                            
                                                            Data_Source_Element.text = ('Landsat Image: LC08_L1TP_014031_20170730_20170811_01_T1; LC08_L1TP_014032_20170730_20170811_01_T1 \nSource: U.S. Geological Survey')
                                                            
                                                            
                                                        else:
                                                            Data_Source_Element.text = ('Landsat Image: {0} \nSource: U.S. Geological Survey'.format(os.path.basename(raster_file[k])[:-17]))
                                              
                                                        #-------dynamic text elements end-------#
                                                        
                                                        
                                                        # legend
                                                        for lyr in legend.listLegendItemLayers():
                                                            #only keep transect points 
                                                            if lyr.name != os.path.basename(shapefiles[j][:-4]):
                                                                legend.removeItem(lyr)
                                                            else:
                                                                # re-name layer 
                                                                lyr.name = 'Transect Point'
                                                                # apply style to layer
                                                                legend.updateItem(lyr, Legend_Style) 
    
                                                        
                                                        # save map as <x> format 
                                                        print "Exporting page {0} of {1} for {2}".format(str(mxd.dataDrivenPages.currentPageID), str(mxd.dataDrivenPages.pageCount), str(city_name))
                                                        outname = os.path.join(outpath,os.path.basename(shapefiles[j][:-8]) + '_' + os.path.basename(raster_file[k])[:25])
    
    #                                                    arcpy.mapping.ExportToPNG(mxd, outname)                                                  
                                                        arcpy.mapping.ExportToPNG(mxd, outname, resolution=300)
                                                            
                                                        # remove df1 layers
                                                        for lyr in arcpy.mapping.ListLayers(mxd, "", df1):
                                                            if lyr.name == 'select_city':
                                                                arcpy.mapping.RemoveLayer(df1, lyr)        
                                                            
                                                        # remove df2 layers
                                                        for lyr in arcpy.mapping.ListLayers(mxd, "", df2):
                                                            arcpy.mapping.RemoveLayer(df2, lyr)                              
                                                        
                                                        # remove df3 layers
                                                        for lyr in arcpy.mapping.ListLayers(mxd, "", df3):
                                                            arcpy.mapping.RemoveLayer(df3, lyr)
                                                        
                                                        # cursor deletion is causing an error; should not, but is... 
    #                                                    del cursor
                                                        
                                                        arcpy.Delete_management(TEMP_TRANSECT_SHAPEFILE)
                                                        arcpy.Delete_management(os.path.join(buffer_directory, TEMP_TRANSECT_SHAPEFILE_BUFFER + '.shp'))
                                            
                                                        arcpy.Delete_management('select_city')
                                                        arcpy.Delete_management(transect_layer)
            
                                                        # Refresh things
                                                        arcpy.RefreshActiveView()
                                                        arcpy.RefreshTOC()
                                                        
    # print out each country template... need's fixing...
#    finally:
#        for country in range(len(country_extent)): 
#            if country_extent[country][0] == fieldValue: 
#                
#                df1.extent = arcpy.Extent(country_extent[country][4],country_extent[country][3], country_extent[country][2], country_extent[country][1])                                                       
#                print "Exporting page {0} of {1} for {2}".format(str(mxd.dataDrivenPages.currentPageID), str(mxd.dataDrivenPages.pageCount), fieldValue)
#                
#                # re-intialize the extent 
#                #                                mxd.activeView = df1
#                # mod to match based on country to country match... 
##                df1.extent = arcpy.Extent(country_extent[country][4],country_extent[country][3], country_extent[country][2], country_extent[country][1]) 
#                
#                Title_Element = arcpy.mapping.ListLayoutElements(mxd, "TEXT_ELEMENT", 'Title')[0]
#                Title_Element.text = fieldValue
#                
#                arcpy.mapping.ExportToPNG(mxd,os.path.join(outpath,fieldValue + ".jpg"))   
#            
#            # Refresh things
#            arcpy.RefreshActiveView()
#            arcpy.RefreshTOC()
                    
#            country_extent.pop(0)
    
