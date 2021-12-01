# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:38:15 2018

@author: Alexander Tong

Developed and tested with Python 2.7.15
"""

def lookup_table():
    
    GLUE_city_to_Landsat = [['Argentina','Buenos Aires','225084'],              
                            ['Australia','Armidale','089081'],
                            ['Australia','Armidale','090081'],
                            ['Australia','Canberra','090085'],
                            ['Australia','Melbourne','093086'],
                            ['Australia','Newcastle','089083'],
                            ['Australia','Parramatta','089083'],
                            ['Australia','Perth','112082'],
                            ['Australia','Sydney','089083'],
                            ['Australia','Hobart','090090'],  
                            ['Australia','Wollongong','089084'],
                            ['Belgium','Antwerp','198024'],   
                            ['Belgium','Brussels','198024'],#
                            ['Belgium','Brussels','198025'],#
                            ['Belgium','Brussels','199024'],#
                            ['Belgium','Brussels','199025'],#              
                            ['Belgium','Ghent','199024'],
                            ['Belgium','Leuven','198025'],               
                            ['Bolivia','La Paz','001071'],#
                            ['Bolivia','La Paz','001072'],#            
                            ['Brazil','Curitiba','220078'],#              
                            ['Brazil','Santa Maria','223080'],
                            ['Brazil','Santa Maria','223081'],                      
                            ['Canada_AB','Calgary','042025'],
                            ['Canada_AB','Edmonton','043023'],
                            ['Canada_AB','St. Albert','043023'], # remove maybe
                            ['Canada_AB','St Albert','043023'],              
                            ['Canada_BC','Vancouver','047026'],           
                            ['Canada_BC','Victoria','047026'], 
                            ['Canada_MB','Dauphin','033024'],
                            ['Canada_MB','Winnipeg','030025'],
                            ['Canada_MB','Winnipeg','031025'],              
                            ['Canada_NB','Fredericton','010028'],
                            ['Canada_NB','Moncton','009028'],             
                            ['Canada_NB','Saint John','009029'],#
                            ['Canada_NB','Saint John','010028'],#
                            ['Canada_NB','Saint John','010029'],#
                            ['Canada_NF',"Saint John's",'001027'],#
                            ['Canada_NF',"Saint John's",'002027'],#    
                            ['Canada_NF','Saint Johns','001027'],#
                            ['Canada_NF','Saint Johns','002027'],#                
                            ['Canada_NF',"St. John's",'001027'],#
                            ['Canada_NF',"St. John's",'002027'],# 
                            ['Canada_NF','St Johns','001027'],#
                            ['Canada_NF','St Johns','002027'],#  
                            ['Canada_NS','Halifax','008029'],              
                            ['Canada_ON','Acton','018030'],
                            ['Canada_ON','Angus','018029'],             
                            ['Canada_ON','Barrie','018029'],
                            ['Canada_ON','Bradford','018029'],
                            ['Canada_ON','Brantford','018030'],
                            ['Canada_ON','Cobourg','017029'],
                            ['Canada_ON','Elmira','019030'],              
                            ['Canada_ON','Everett','018029'],                
                            ['Canada_ON','Fergus','019030'],              
                            ['Canada_ON','Georgetown','018030'],
                            ['Canada_ON','Guelph','018030'],
                            ['Canada_ON','Kingston','016029'],              
                            ['Canada_ON','London','019030'],              
                            ['Canada_ON','New Tecumseth','018029'],              
                            ['Canada_ON','North Bay','018028'],
                            ['Canada_ON','Orangeville','018030'],          
                            ['Canada_ON','Ottawa','015029'],
                            ['Canada_ON','Port Hope','017029'],
                            ['Canada_ON','St Thomas','019030'],                 
                            ['Canada_ON','Stratford','018030'],
                            ['Canada_ON','Stratford','019030'],
                            ['Canada_ON','Sudbury','019028'], 
                            ['Canada_ON','Timmins','020026'],              
                            ['Canada_ON','Toronto','018029'],
                            ['Canada_ON','Toronto','018030'],              
                            ['Canada_ON','Waterloo','018030'],
                            ['Canada_ON','Waterloo','019030'],               
                            ['Canada_ON','Whitchurch-Stouffville','018029'],               
                            ['Canada_ON','Whitchurch Stouffville','018029'], 
                            ['Canada_ON','Woodstock','018030'],
                            ['Canada_ON','Woodstock','019030'],             
                            ['Canada_PEI','Charlottetown','008028'],    
                            ['Canada_QC','Montreal','014028'],
                            ['Canada_QC','Montreal','015028'],
                            ['Canada_QC','Quebec City','013028'],
                            ['Chile','Punta Arenas','228097'],
                            ['Chile','Punta Arenas','229097'],            
                            ['Chile','Rancagua','233084'],              
                            ['Chile','Santiago','233083'],              
                            ['Chile','Talca','233085'],       
                            ['Chile','Temuco','233087'],#
                            ['China','Beijing','123032'],             
                            ['China','Chengdu','129039'],#
                            ['China','Chengdu','130039'],#            
                            ['China','Kunming','129043'],             
#                            ['China','Lanzhou','130035'],#
                            ['China','Lanzhou','131035'],#              
                            ['China','Shanghai','118038'],              
                            ['China','Wuhan','122038'],#
                            ['China','Wuhan','122039'],#
                            ['China','Wuhan','123038'],#
                            ['China','Wuhan','123039'],#             
                            ['Colombia','Bogota','007056'],#
                            ['Colombia','Bogota','007057'],#
                            ['Colombia','Bogota','008056'],#
                            ['Colombia','Bogota','008057'],#
                            ['Colombia','Bucaramanga','007055'],#
                            ['Colombia','Bucaramanga','008055'],#
                            ['Colombia','Medellin','009055'],#
                            ['Colombia','Medellin','009056'],#
                            ['Colombia','Santa Marta','008052'],#
                            ['Colombia','Santa Marta','008053'],#
                            ['Colombia','Santa Marta','009052'],#
                            ['Colombia','Santa Marta','009053'],#              
                            ['Czech Republic','Prague','191025'],#
                            ['Czech Republic','Prague','192025'],#              
                            ['Ecuador','Loja','010063'],       
                            ['Ecuador','Latacunga','010061'], #
                            ['Ecuador','Quito','010060'],
                            ['Finland','Helsinki','188018'],                    
                            ['France','Montpellier','196030'],#
                            ['France','Montpellier','197029'],#
                            ['France','Montpellier','197030'],#              
                            ['France','Paris','199026'],
                            ['France','Rennes','201027'],  
                            ['France','Rennes','202026'],
                            ['Germany','Berlin','193023'],#
                            ['Germany','Berlin','193024'],#
                            ['Germany','Cologne','197024'],             
                            ['Germany','Frankfurt','196025'],  
#                            ['Germany','Halle (Saale)','193024'],
#                            ['Germany','Halle (Saale)','194024'],                  
                            ['Germany','Halle','193024'],
                            ['Germany','Halle','194024'],     
                            ['Germany','Landshut','192026'],
#                            ['Germany','Landshut','193026'],
                            ['Germany','Munich','193026'],              
                            ['Germany','Münster','196024'],
                            ['Germany','Münster','197024'],             
                            ['Germany','Munster','196024'],
                            ['Germany','Munster','197024'],               
                            ['Germany','Reutlingen','194026'],
                            ['Germany','Reutlingen','195026'],              
                            ['Germany','Stuttgart','194026'],
                            ['Germany','Stuttgart','195026'],
                            ['Greece','Thessaloniki','184032'],
                            ['Iran', 'Tehran', '164035'],              
                            ['Italy','Milan','193028'],#
                            ['Italy','Milan','193029'],#
                            ['Italy','Milan','194028'],#              
                            ['Japan','Hiroshima','112036'],
                            ['Japan','Kyoto','110036'],
                            ['Japan','Sapporo','108030'],
                            ['Japan','Tokyo','107035'],
#                            ['Mexico','Guadalajara (Jalisco)','029046'],
                            ['Mexico','Mexico City','026047'],
                            ['Mexico','Morelia','027046'],            
#                            ['Mexico','Puebla','025047'],#
#                            ['Mexico','Puebla','026046'],#
#                            ['Mexico','Puebla','026047'],#   
                            ['Mexico','Toluca', '026047'],
                            ['Mexico','Uruapan','028047'],#  
                            ['Mexico','Xalapa','025046'],             
                            ['Netherlands','Amsterdam','198023'],#
                            ['Netherlands','Amsterdam','198024'],#              
                            ['Netherlands','Leiden','198024'],
                            ['Netherlands','Leiden','199024'],               
                            ['New Zealand','Christchurch','073090'],#
                            ['New Zealand','Christchurch','074090'],#              
                            ['New Zealand','Palmerston North','072088'],              
                            ['Norway','Trondheim','199016'],
                            ['Papua New Guinea','Wewak','098063'],#
                            ['Papua New Guinea','Wewak','099062'],#
                            ['Papua New Guinea','Wewak','099063'],#             
                            ['Poland','Warsaw','187024'],
                            ['Poland','Warsaw','188024'],            
                            ['Portugal', 'Almada','204033'],
                            ['Portugal','Lisbon','204033'],              
                            ['Russia','Moscow','178021'],#
                            ['Russia','Moscow','179021'],#                
#                            ['South Africa','Cape Town','175083'],
                            ['South Africa','Cape Town','175084'],
                            ['South Africa','Pietermaritzburg','168080'],#
                            ['South Africa','Pietermaritzburg','168081'],#          
                            ['Sweden','Linköping','193019'],
                            ['Sweden','Linköping','194019'],
                            ['Sweden','Linkoping','193019'],
                            ['Sweden','Linkoping','194019'],             
                            ['Sweden','Malmö','194021'],
                            ['Sweden','Malmo','194021'],              
#                            ['Sweden','Norrtälje','192018'],
#                            ['Sweden','Norrtälje','193018'],
                            ['Sweden','Norrtalje','192018'],
#                            ['Sweden','Norrtalje','193018'],              
                            ['Sweden','Stockholm','192019'],              
                            ['Sweden','Uppsala','193018'],#              
                            ['Switzerland','Bern','195027'],
                            ['Switzerland','Zurich','194027'],
                            ['Switzerland','Zurich','195027'],
                            ['UK','Birmingham','202023'],#
                            ['UK','Birmingham','202024'],#
                            ['UK','Birmingham','203023'],#              
                            ['UK','Brighton','201025'],
                            ['UK','Cambridge','201024'],              
                            ['UK','Glasgow','205021'],#
                            ['UK','Glasgow','206021'],#              
                            ['UK','Manchester','203023'],
                            ['UK','Manchester','204023'],
                            ['UK','Reading','202024'],
                            #['USA_AK','Anchorage','068017'],
                            ['USA_AK','Anchorage','069017'],
                            #['USA_AK','Anchorage','070017'],
                            #['USA_AK','Fairbanks','068015'],
                            #['USA_AK','Fairbanks','069014'],
                            #['USA_AK','Fairbanks','069015'],
                            ['USA_AK','Fairbanks','070014'],
                            #['USA_AK','Fairbanks','070015'],            
                            ['USA_AR','Little Rock','024036'],            
                            ['USA_CA','La Verne','040036'],# this is claremont; might need to rename 
                            ['USA_CA','Pomona','040036'],# this is claremont; might need to rename               
                            ['USA_CA','Los Angeles','041036'],
                            ['USA_CA','Sacramento','044033'],#                
                            ['USA_CA','Santa Cruz','044034'],#
                            ['USA_CO','Denver','033033'],
                            ['USA_CO','Fort Collins','034032'],#
                            ['USA_CT','New Haven','013031'],
                            ['USA_DC','Washington','015033'],#
                            ['USA_FL','Jacksonville','016039'],#
                            ['USA_FL','Jacksonville','017039'],#
                            ['USA_FL','Tampa','017040'],              
                            ['USA_GA','Athens','018037'],
                            ['USA_GA','Atlanta','019036'],#
                            ['USA_GA','Atlanta','019037'],#
                            ['USA_IL','Chicago','023031'],
                            ['USA_IL','Rockford','024031'],
                            ['USA_IN','Indianapolis','021032'],
                            ['USA_KY','Louisville','021033'],
                            ['USA_LA','Baton Rouge','023039'],
                            ['USA_LA','Lake Charles','024039'],
                            ['USA_MA','Boston','012031'],#
                            ['USA_MA','Northampton','013031'],
                            ['USA_MD','Baltimore','015033'],#
                            ['USA_ME','Augusta','012029'],
#                            ['USA_ME','Portland','012030'], # disambiguation; error when process_csv_to_shp
                            ['USA_ME','Portland_ME','012030'],
                            ['USA_ME','Waterville','011029'],
                            ['USA_MI','Detroit','020030'],#
                            ['USA_MI','Detroit','020031'],#
                            ['USA_MI','Kalamazoo','021031'],#
                            ['USA_MI','Kalamazoo','022030'],#
                            ['USA_MI','Kalamazoo','022031'],#
                            ['USA_MI','Lansing','021030'],#
                            ['USA_MN','Minneapolis','027028'],#
                            ['USA_MN','Minneapolis','027029'],#            
                            ['USA_MO','Ellisville','024033'],#
                            ['USA_MO','St. Louis','024033'], # remove maybe
                            ['USA_MO','St Louis','024033'],              
                            ['USA_MS','Starkville','022037'],
                            ['USA_NC','Charlotte','017035'],#
#                            ['USA_NC','Charlotte','017036'],#
                            ['USA_NC','Greenville','015035'],#
                            ['USA_NC','Raleigh','016035'],
                            ['USA_NJ','Atlantic City','014033'],
                            ['USA_NJ','Freehold','014032'],
                            ['USA_NM','Albuquerque','034036'],
                            ['USA_NY','Hempstead','013032'],
                            ['USA_NY','New York','014031'],
                            ['USA_NY','Rochester','016030'],
                            ['USA_OH','Cincinnati','020033'],#
                            ['USA_OH','Cleveland','018031'],#
                            ['USA_OH','Cleveland','019031'],#
#                            ['USA_OR','Portland','046028'], # disambiguation; error when process_csv_to_shp
                            ['USA_OR','Portland_OR','046028'], 
                            ['USA_OR','Salem','046029'],#
                            ['USA_PA','Philadelphia','014032'],#
                            ['USA_PA','Pittsburgh','017032'],#
                            ['USA_RI','Providence','012031'],
                            ['USA_SC','Charleston','016037'],#
                            ['USA_SD','Sioux Falls','029030'],        
                            ['USA_TN','Memphis','023036'],
                            ['USA_TX','Houston','025040'],
                            ['USA_VA','Norfolk','014034'],#
                            ['USA_VA','Norfolk','014035'],#
                            ['USA_VT','Burlington','014029'],
                            ['USA_WA','Seattle','046027'],
                            ['USA_WA','Tacoma','046027'],
                            ['USA_WI','Madison','024030']]
        
    return GLUE_city_to_Landsat

if __name__ == "__main__":
    
    lookup_table()
    
