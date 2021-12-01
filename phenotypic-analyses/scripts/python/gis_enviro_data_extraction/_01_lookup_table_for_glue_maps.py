# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:00:00 2019

@author: Alexander Tong

Developed and tested with Python 2.7.15
"""

def lookup_table():
    
    GLUE_city_to_Landsat = [['Argentina','Buenos Aires Province','Buenos Aires','225084'],              
                            ['Australia','New South Wales','Armidale','089081'],
                            ['Australia','New South Wales','Armidale','090081'],
                            ['Australia','Australian Capital Territory','Canberra','090085'],
                            ['Australia','Victoria','Melbourne','093086'],
                            ['Australia','New South Wales','Newcastle','089083'],
                            ['Australia','New South Wales','Parramatta','089083'],
                            ['Australia','Western Australia','Perth','112082'],
                            ['Australia','New South Wales','Sydney','089083'],
                            ['Australia','Tasmania','Hobart','090090'],  
                            ['Australia','New South Wales', 'Wollongong','089084'],
                            ['Belgium','Antwerp Province','Antwerp','198024'],   
                            ['Belgium','Brussels-Capital Region','Brussels','198024'],#
                            ['Belgium','Brussels-Capital Region','Brussels','198025'],#
                            ['Belgium','Brussels-Capital Region','Brussels','199024'],#
                            ['Belgium','Brussels-Capital Region','Brussels','199025'],#              
                            ['Belgium','East Flanders','Ghent','199024'],
                            ['Belgium','Flemish Brabant','Leuven','198025'],               
                            ['Bolivia','La Paz Department','La Paz','001071'],#
                            ['Bolivia','La Paz Department','La Paz','001072'],#            
                            ['Brazil','Parana','Curitiba','220078'],#              
                            ['Brazil','Rio Grande do Sul','Santa Maria','223080'],
                            ['Brazil','Rio Grande do Sul','Santa Maria','223081'],                      
                            ['Canada_AB','Alberta','Calgary','042025'],
                            ['Canada_AB','Alberta','Edmonton','043023'],
                            ['Canada_AB','Alberta','St. Albert','043023'], # remove maybe
                            ['Canada_AB','Alberta','St Albert','043023'],              
                            ['Canada_BC','British Columbia','Vancouver','047026'],           
                            ['Canada_BC','British Columbia','Victoria','047026'], 
                            ['Canada_MB','Manitoba','Dauphin','033024'],
                            ['Canada_MB','Manitoba','Winnipeg','030025'],
                            ['Canada_MB','Manitoba','Winnipeg','031025'],              
                            ['Canada_NB','New Brunswick','Fredericton','010028'],
                            ['Canada_NB','New Brunswick','Moncton','009028'],             
                            ['Canada_NB','New Brunswick','Saint John','009029'],#
                            ['Canada_NB','New Brunswick','Saint John','010028'],#
                            ['Canada_NB','New Brunswick','Saint John','010029'],#
                            ['Canada_NF','Newfoundland & Labrador',"Saint John's",'001027'],#
                            ['Canada_NF','Newfoundland & Labrador',"Saint John's",'002027'],#    
                            ['Canada_NF','Newfoundland & Labrador','Saint Johns','001027'],#
                            ['Canada_NF','Newfoundland & Labrador','Saint Johns','002027'],#                
                            ['Canada_NF','Newfoundland & Labrador',"St. John's",'001027'],#
                            ['Canada_NF','Newfoundland & Labrador',"St. John's",'002027'],# 
                            ['Canada_NF','Newfoundland & Labrador','St Johns','001027'],#
                            ['Canada_NF','Newfoundland & Labrador','St Johns','002027'],#  
                            ['Canada_NS','Nova Scotia','Halifax','008029'],              
                            ['Canada_ON','Ontario','Acton','018030'],
                            ['Canada_ON','Ontario','Angus','018029'],             
                            ['Canada_ON','Ontario','Barrie','018029'],
                            ['Canada_ON','Ontario','Bradford','018029'],
                            ['Canada_ON','Ontario','Brantford','018030'],
                            ['Canada_ON','Ontario','Cobourg','017029'],
                            ['Canada_ON','Ontario','Elmira','019030'],              
                            ['Canada_ON','Ontario','Everett','018029'],                
                            ['Canada_ON','Ontario','Fergus','019030'],              
                            ['Canada_ON','Ontario','Georgetown','018030'],
                            ['Canada_ON','Ontario','Guelph','018030'],
                            ['Canada_ON','Ontario','Kingston','016029'],              
                            ['Canada_ON','Ontario','London','019030'],              
                            ['Canada_ON','Ontario','New Tecumseth','018029'],              
                            ['Canada_ON','Ontario','North Bay','018028'],
                            ['Canada_ON','Ontario','Orangeville','018030'],          
                            ['Canada_ON','Ontario','Ottawa','015029'],
                            ['Canada_ON','Ontario','Port Hope','017029'],
                            ['Canada_ON','Ontario','St Thomas','019030'],                 
                            ['Canada_ON','Ontario','Stratford','018030'],
                            ['Canada_ON','Ontario','Stratford','019030'],
                            ['Canada_ON','Ontario','Sudbury','019028'], 
                            ['Canada_ON','Ontario','Timmins','020026'],              
                            ['Canada_ON','Ontario','Toronto','018029'],
                            ['Canada_ON','Ontario','Toronto','018030'],              
                            ['Canada_ON','Ontario','Waterloo','018030'],
                            ['Canada_ON','Ontario','Waterloo','019030'],               
                            ['Canada_ON','Ontario','Whitchurch-Stouffville','018029'],               
                            ['Canada_ON','Ontario','Whitchurch Stouffville','018029'], 
                            ['Canada_ON','Ontario','Woodstock','018030'],
                            ['Canada_ON','Ontario','Woodstock','019030'],             
                            ['Canada_PEI','Prince Edward Island','Charlottetown','008028'],    
                            ['Canada_QC','Quebec','Montreal','014028'],
                            ['Canada_QC','Quebec','Montreal','015028'],
                            ['Canada_QC','Quebec','Quebec City','013028'],
                            ['Chile','Magallanes Region','Punta Arenas','228097'],
                            ['Chile','Magallanes Region','Punta Arenas','229097'],            
                            ['Chile',"O'Higgins Region",'Rancagua','233084'],              
                            ['Chile','Santiago Metropolitan Region','Santiago','233083'],              
                            ['Chile','Maule','Talca','233085'],       
                            ['Chile','Araucania Region','Temuco','233087'],#
                            ['China','Beijing Municipality','Beijing','123032'],             
                            ['China','Sichuan','Chengdu','129039'],#
                            ['China','Sichuan','Chengdu','130039'],#            
                            ['China','Yunnan','Kunming','129043'],             
#                            ['China','Lanzhou','130035'],#
                            ['China','Gansu','Lanzhou','131035'],#              
                            ['China','Shanghai Municipality','Shanghai','118038'],              
                            ['China','Hubei','Wuhan','122038'],#
                            ['China','Hubei','Wuhan','122039'],#
                            ['China','Hubei','Wuhan','123038'],#
                            ['China','Hubei','Wuhan','123039'],#             
                            ['Colombia','District Capital','Bogota','007056'],#
                            ['Colombia','District Capital','Bogota','007057'],#
                            ['Colombia','District Capital','Bogota','008056'],#
                            ['Colombia','District Capital','Bogota','008057'],#
                            ['Colombia','Santander','Bucaramanga','007055'],#
                            ['Colombia','Santander','Bucaramanga','008055'],#
                            ['Colombia','Antioquia','Medellin','009055'],#
                            ['Colombia','Antioquia','Medellin','009056'],#
                            ['Colombia','Magdalena','Santa Marta','008052'],#
                            ['Colombia','Magdalena','Santa Marta','008053'],#
                            ['Colombia','Magdalena','Santa Marta','009052'],#
                            ['Colombia','Magdalena','Santa Marta','009053'],#              
                            ['Czech Republic','Prague','191025'],#
                            ['Czech Republic','Prague','192025'],#              
                            ['Ecuador','Loja Province','Loja','010063'],       
                            ['Ecuador','Cotopaxi','Latacunga','010061'], #
                            ['Ecuador','Pichincha','Quito','010060'],
                            ['Finland','Helsinki-Uusimaa Region','Helsinki','188018'],                    
                            ['France','Occitanie','Montpellier','196030'],#
                            ['France','Occitanie','Montpellier','197029'],#
                            ['France','Occitanie','Montpellier','197030'],#              
                            ['France','Paris Region','Paris','199026'],
                            ['France','Bretagne','Rennes','201027'],  
                            ['France','Bretagne','Rennes','202026'],
                            ['Germany','Berlin','Berlin','193023'],#
                            ['Germany','Berlin','Berlin','193024'],#
                            ['Germany','Nordrhein-Westfalen','Cologne','197024'],             
                            ['Germany','Hessen','Frankfurt','196025'],  
#                            ['Germany','Halle (Saale)','193024'],
#                            ['Germany','Halle (Saale)','194024'],                  
                            ['Germany','Sachsen-Anhalt','Halle','193024'],
                            ['Germany','Sachsen-Anhalt','Halle','194024'],
                            ['Germany','Bayern','Landshut','192026'],
#                            ['Germany','Bayern','Landshut','193026'],
                            ['Germany','Bayern','Munich','193026'],              
                            ['Germany','North Rhine-Westphalia','Münster','196024'],
                            ['Germany','North Rhine-Westphalia','Münster','197024'],             
                            ['Germany','North Rhine-Westphalia','Munster','196024'],
                            ['Germany','North Rhine-Westphalia','Munster','197024'],               
                            ['Germany','Baden-Wurttemberg','Reutlingen','194026'],
                            ['Germany','Baden-Wurttemberg','Reutlingen','195026'],              
                            ['Germany','Baden-Wurttemberg','Stuttgart','194026'],
                            ['Germany','Baden-Wurttemberg','Stuttgart','195026'],
                            ['Greece','Central Macedonia','Thessaloniki','184032'],
                            ['Iran','Tehran Province','Tehran', '164035'],              
                            ['Italy','Lombardia','Milan','193028'],#
                            ['Italy','Lombardia','Milan','193029'],#
                            ['Italy','Lombardia','Milan','194028'],#              
                            ['Japan','Hiroshima Prefecture','Hiroshima','112036'],
                            ['Japan','Kyoto Prefecture','Kyoto','110036'],
                            ['Japan','Hokkaido Prefecture','Sapporo','108030'],
                            ['Japan','Tokyo Metropolis','Tokyo','107035'],
#                            ['Mexico','Guadalajara (Jalisco)','029046'],
                            ['Mexico','Federal District','Mexico City','026047'],
                            ['Mexico','Michoacan','Morelia','027046'],            
#                            ['Mexico','Puebla','025047'],#
#                            ['Mexico','Puebla','026046'],#
#                            ['Mexico','Puebla','026047'],#   
                            ['Mexico','Mexico State','Toluca', '026047'],
                            ['Mexico','Michoacan','Uruapan','028047'],#  
                            ['Mexico','Veracruz','Xalapa','025046'],             
                            ['Netherlands','North Holland','Amsterdam','198023'],#
                            ['Netherlands','North Holland','Amsterdam','198024'],#              
                            ['Netherlands','South Holland','Leiden','198024'],
                            ['Netherlands','South Holland','Leiden','199024'],               
                            ['New Zealand','Canterbury','Christchurch','073090'],#
                            ['New Zealand','Canterbury','Christchurch','074090'],#              
                            ['New Zealand','Manawatu-Wanganui','Palmerston North','072088'],              
                            ['Norway','Trondelag','Trondheim','199016'],
                            ['Papua New Guinea','Wewak','098063'],#
                            ['Papua New Guinea','Wewak','099062'],#
                            ['Papua New Guinea','Wewak','099063'],#             
                            ['Poland','Masovian','Warsaw','187024'],
                            ['Poland','Masovian','Warsaw','188024'],            
                            ['Portugal','Lisboa Region','Almada','204033'],
                            ['Portugal','Lisboa Region','Lisbon','204033'],              
                            ['Russia','Central Federal District','Moscow','178021'],#
                            ['Russia','Central Federal District','Moscow','179021'],#         
#                            ['South Africa','Western Cape','Cape Town','175083'],
                            ['South Africa','Western Cape','Cape Town','175084'],
                            ['South Africa','KwaZulu-Natal','Pietermaritzburg','168080'],#
                            ['South Africa','KwaZulu-Natal','Pietermaritzburg','168081'],#          
                            ['Sweden','East Gothland','Linköping','193019'],
                            ['Sweden','East Gothland','Linköping','194019'],
                            ['Sweden','East Gothland','Linkoping','193019'],
                            ['Sweden','East Gothland','Linkoping','194019'],             
                            ['Sweden','Scania','Malmö','194021'],
                            ['Sweden','Scania','Malmo','194021'],              
#                            ['Sweden','Uppland','Norrtälje','192018'],
#                            ['Sweden','Uppland','Norrtälje','193018'],
                            ['Sweden','Uppland','Norrtalje','192018'],
#                            ['Sweden','Uppland','Norrtalje','193018'],              
                            ['Sweden','Uppland','Stockholm','192019'],              
                            ['Sweden','Uppland','Uppsala','193018'],#              
                            ['Switzerland','Canton of Bern','Bern','195027'],
                            ['Switzerland','Canton of Zurich','Zurich','194027'],
                            ['Switzerland','Canton of Zurich','Zurich','195027'],
                            ['UK','West Midlands','Birmingham','202023'],#
                            ['UK','West Midlands','Birmingham','202024'],#
                            ['UK','West Midlands','Birmingham','203023'],#              
                            ['UK','Brighton and Hove','Brighton','201025'],
                            ['UK','Cambridgeshire','Cambridge','201024'],              
                            ['UK','Scotland','Glasgow','205021'],#
                            ['UK','Scotland','Glasgow','206021'],#              
                            ['UK','Greater Manchester','Manchester','203023'],
                            ['UK','Greater Manchester','Manchester','204023'],
                            ['UK','Oxfordshire','Reading','202024'],
                            ['USA_AK','Alaska','Anchorage','068017'],#
                            ['USA_AK','Alaska','Anchorage','069017'],#
                            ['USA_AK','Alaska','Anchorage','070017'],#
                            ['USA_AK','Alaska','Fairbanks','068015'],#
                            ['USA_AK','Alaska','Fairbanks','069014'],#
                            ['USA_AK','Alaska','Fairbanks','069015'],#
                            ['USA_AK','Alaska','Fairbanks','070014'],#
                            ['USA_AK','Alaska','Fairbanks','070015'],#             
                            ['USA_AR','Arkansas','Little Rock','024036'],            
                            ['USA_CA','California','La Verne','040036'],# this is claremont; might need to rename 
                            ['USA_CA','California','Pomona','040036'],# this is claremont; might need to rename               
                            ['USA_CA','California','Los Angeles','041036'],
                            ['USA_CA','California','Sacramento','044033'],#                
                            ['USA_CA','California','Santa Cruz','044034'],#
                            ['USA_CO','Colorado','Denver','033033'],
                            ['USA_CO','Colorado','Fort Collins','034032'],#
                            ['USA_CT','Connecticut','New Haven','013031'],
                            ['USA_DC','District of Columbia','Washington','015033'],#
                            ['USA_FL','Florida','Jacksonville','016039'],#
                            ['USA_FL','Florida','Jacksonville','017039'],#
                            ['USA_FL','Florida','Tampa','017040'],              
                            ['USA_GA','Georgia','Athens','018037'],
                            ['USA_GA','Georgia','Atlanta','019036'],#
                            ['USA_GA','Georgia','Atlanta','019037'],#
                            ['USA_IL','Illinois','Chicago','023031'],
                            ['USA_IL','Illinois','Rockford','024031'],
                            ['USA_IN','Indiana','Indianapolis','021032'],
                            ['USA_KY','Kentucky','Louisville','021033'],
                            ['USA_LA','Louisiana','Baton Rouge','023039'],
                            ['USA_LA','Louisiana','Lake Charles','024039'],
                            ['USA_MA','Massachusetts','Boston','012031'],#
                            ['USA_MA','Massachusetts','Northampton','013031'],
                            ['USA_MD','Maryland','Baltimore','015033'],#
                            ['USA_ME','Maine','Augusta','012029'],
#                            ['USA_ME','Portland','012030'], # disambiguation; error when process_csv_to_shp
                            ['USA_ME','Maine','Portland_ME','012030'],
                            ['USA_ME','Maine','Waterville','011029'],
                            ['USA_MI','Michigan','Detroit','020030'],#
                            ['USA_MI','Michigan','Detroit','020031'],#
                            ['USA_MI','Michigan','Kalamazoo','021031'],#
                            ['USA_MI','Michigan','Kalamazoo','022030'],#
                            ['USA_MI','Michigan','Kalamazoo','022031'],#
                            ['USA_MI','Michigan','Lansing','021030'],#
                            ['USA_MN','Minnesota','Minneapolis','027028'],#
                            ['USA_MN','Minnesota','Minneapolis','027029'],#            
                            ['USA_MO','Missouri','Ellisville','024033'],#
                            ['USA_MO','Missouri','St. Louis','024033'], # remove maybe
                            ['USA_MO','Missouri','St Louis','024033'],              
                            ['USA_MS','Mississippi','Starkville','022037'],
                            ['USA_NC','North Carolina','Charlotte','017035'],#
#                            ['USA_NC','Charlotte','017036'],#
                            ['USA_NC','North Carolina','Greenville','015035'],#
                            ['USA_NC','North Carolina','Raleigh','016035'],
                            ['USA_NJ','New Jersey','Atlantic City','014033'],
                            ['USA_NJ','New Jersey','Freehold','014032'],
                            ['USA_NM','New Mexico','Albuquerque','034036'],
                            ['USA_NY','New York','Hempstead','013032'],
                            ['USA_NY','New York','New York','014031'],
                            ['USA_NY','New York','Rochester','016030'],
                            ['USA_OH','Ohio','Cincinnati','020033'],#
                            ['USA_OH','Ohio','Cleveland','018031'],#
                            ['USA_OH','Ohio','Cleveland','019031'],#
#                            ['USA_OR','Portland','046028'], # disambiguation; error when process_csv_to_shp
                            ['USA_OR','Oregon','Portland_OR','046028'], 
                            ['USA_OR','Oregon','Salem','046029'],#
                            ['USA_PA','Pennsylvania','Philadelphia','014032'],#
                            ['USA_PA','Pennsylvania','Pittsburgh','017032'],#
                            ['USA_RI','Rhode Island','Providence','012031'],
                            ['USA_SC','South Carolina','Charleston','016037'],#
                            ['USA_SD','South Dakota','Sioux Falls','029030'],        
                            ['USA_TN','Tennessee','Memphis','023036'],
                            ['USA_TX','Texas','Houston','025040'],
                            ['USA_VA','Virginia','Norfolk','014034'],#
                            ['USA_VA','Virginia','Norfolk','014035'],#
                            ['USA_VT','Vermont','Burlington','014029'],
                            ['USA_WA','Washington','Seattle','046027'],
                            ['USA_WA','Washington','Tacoma','046027'],
                            ['USA_WI','Wisconsin','Madison','024030']]
        
    return GLUE_city_to_Landsat

if __name__ == "__main__":
    
    lookup_table()
    
