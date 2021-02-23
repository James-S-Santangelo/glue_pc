# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 16:32:37 2018

@author: Alex

Developed and tested with Python 2.7.15
"""


def join_johnson_csv(directory, outpath):
    '''
    Description:
        Received csv with HCN data but without the lat/lon coord for each transect point. \n\

        Therefore need to join:
            HCN data == '20_cities_cyano_poplevel.csv'
            Lat/Lon data == 'Johnson et al_ON clover_lat_long_for Tong.csv'

        Primary/Foreign Key Join on 'Site'

        **NOTE: Both csv for processing are hard-coded into logic...

    Args:
        directory (str): specify inpath directory/folder
        outpath (str): specify outpath directory/folder

    Output:
         '20_cities_cyano_poplevel_lat_lon.csv'
    '''
    import os
    import pandas as pd

    for i in os.listdir(directory):
        if '20_cities_cyano_poplevel' in i:
            # print (i)
            df_hcn = pd.read_csv(os.path.join(directory, i))
            # print (df_hcn.dtypes)
        if 'Johnson et al_ON clover_lat_long_for Tong' in i:
            # print (i)
            df_coord = pd.read_csv(os.path.join(directory, i))
            # print (df_coord.dtypes)

        # all same merge, different ways of verbalizing it
#        merged_df = df_hcn.merge(df_coord, how = 'inner', on = ['Site','Site'])
#        merged_df = df_hcn.merge(df_coord, left_on='Site', right_on='Site', how='outer', validate = 'one_to_one')
#        merged_df = df_hcn.merge(df_coord, left_index ='Site', right_index ='Site', how='outer', validate = 'one_to_one')

        # keep only specified columns from join (SQL left join)
        # https://stackoverflow.com/questions/17978133/python-pandas-merge-only-certain-columns

    # James: Moved indentation one to the left (i.e. outside of for loop)
    # to prevent "local variable 'df_hcn' referenced before assignment" error.
    merged_df = pd.merge(
        df_hcn, df_coord.iloc[:, 3:5], how='left', left_on=df_hcn['Site'], right_on=df_coord['Site'])

    # save
    merged_df.to_csv(os.path.join(
        outpath, '20_cities_cyano_poplevel_lat_lon' + '.csv'), sep=',', index=False)


def split_cities(directory, csv, outpath):
    '''
    Description:
        purpose built to process a single csv at a time. This function will handle \n\
        cases associated with:
            '20_cities_cyano_poplevel_lat_lon.csv'
            'Santangelo_AllCities_AllPlants.csv'

    Args:
        directory (str): specify directory filepath
        csv (str): specify csv filepath
        outpath (str): specify filepath for output csv files

    Returns:
        pass
    '''
    import os
    import re
    import pandas as pd

    for i in os.listdir(directory):
        if csv in i:
            # print (i)

            df = pd.read_csv(os.path.join(directory, i))

            # step 1. group each dataframe by city name
            # satisfy for Marc Johnson's dataset
            if '20_cities_cyano_poplevel_lat_lon' in i:
                # clean up and process Marc Johnson's csv
                # remove first col 'key_0' and second col 'Order'
                # print(df)
                # James: Change indexing to drop 0:1 rather than 0:2, which
                # was resulting in a shift of the 'site' column
                df = df.drop(df.columns[0:1], axis=1)
                # print(df)

                # rename col 'cy' to 'HCN_Result'
                df = df.rename(index=str, columns={"cy": "HCN_Result"})

                # Make sure there are no spaces in city/file names
                df['city'] = df['city'].str.replace(" ", "_")

                # process sites names; separate city acronym from population num (e.g., Ac 1)
                sites = df[df.columns[2]].values.tolist()
                # split city code; we just want to the population number...
                # https://stackoverflow.com/questions/3340081/product-code-looks-like-abcd2343-what-to-split-by-letters-and-numbers
                for i in range(len(sites)):
                    # print(sites[i])
                    # remove spaces
                    #            sites[i] = sites[i].replace(' ', '')
                    # separate alphabet from numeric
                    #            sites[i] = re.split('(\d+)',sites[i])

                    #            print (re.split('\d*\D+', sites[i])) # return number with empty string in preceding index
                    sites[i] = re.findall('\d+', sites[i])  # return only number
                    # print(sites[i])
                    sites[i] = int(sites[i][0])

        #        # remove empty string (deprecated usage from:  sites[i] = re.split('(\d+)',sites[i]))
        #        for i in range(len(sites)):
        #            del (sites[i][-1])

                # insert population number
                df.insert(loc=2, column='population', value=sites)
                # remove city acronym
                df = df.drop(df.columns[3], axis=1)
                df = df.drop(df.columns[0], axis=1)
                # group subsequent dataframes by city name
                grouped = df.groupby(df.columns[0])
            # satisfy for James Santangelo's dataset
            else:
                # group subsequent dataframes by city name
                grouped = df.groupby(df.columns[0])

            # for group by city name in column at index 0, process...
            for i in grouped.groups:

                individual_group = grouped.get_group(i)

                try:
                    # fix for Marc Johnson's naming convention of latitude and longitude headers
                    lat_index_search = individual_group.columns.str.contains(
                        'Lat')  # returns boolean array
                    lon_index_search = individual_group.columns.str.contains(
                        'Long')  # returns boolean array

                    for j in range(len(lat_index_search)):
                        if lat_index_search[j] == True:
                            lat_index = j

                            # rename col header by index
                            individual_group = individual_group.rename(
                                columns={individual_group.columns[lat_index]: 'population_latitude'})

    #                        print (individual_group.columns[lat_index])

                    for j in range(len(lon_index_search)):
                        if lon_index_search[j] == True:
                            lon_index = j

                            # rename col header by index
                            individual_group = individual_group.rename(
                                columns={individual_group.columns[lon_index]: 'population_longitude'})

    #                        print (individual_group.columns[lon_index])

                except:
                    # fix for Jame Santangelo's naming convention of latitude and longitude headers
                    lat_index_search = individual_group.columns.str.contains(
                        'Lat.pop')  # returns boolean array
                    lon_index_search = individual_group.columns.str.contains(
                        'Long.pop')  # returns boolean array

                    for j in range(len(lat_index_search)):
                        if lat_index_search[j] == True:
                            lat_index = j

                            # rename col header by index
                            individual_group = individual_group.rename(
                                columns={individual_group.columns[lat_index]: 'population_latitude'})

    #                        print (individual_group.columns[lat_index])

                    for j in range(len(lon_index_search)):
                        if lon_index_search[j] == True:
                            lon_index = j

                            # rename col header by index
                            individual_group = individual_group.rename(
                                columns={individual_group.columns[lon_index]: 'population_longitude'})

    #                        print (individual_group.columns[lon_index])

                # fix naming convention from Jame Santangelo's aggregated val
                if individual_group.iloc[0, 0] == 'NewYork':
                    individual_group.iloc[:, 0] = 'New_York'

                    # save
                    individual_group.to_csv(os.path.join(
                        outpath, individual_group.iloc[0, 0] + '.csv'), sep=',', index=False)

                # fix naming convention from Jame Santangelo's aggregated val
                elif individual_group.iloc[0, 0] == 'Washington D.C.':
                    individual_group.iloc[:, 0] = 'Washington'

                    # save
                    individual_group.to_csv(os.path.join(
                        outpath, individual_group.iloc[0, 0] + '.csv'), sep=',', index=False)

                else:
                    # save
                    individual_group.to_csv(os.path.join(
                        outpath, individual_group.iloc[0, 0] + '.csv'), sep=',', index=False)


#            df = pd.read_csv(os.path.join(directory,i))
#            grouped = df.groupby(df.columns[0])
#
#            # same as below
#            split_list = []
#            for i in grouped.groups:
#                individual_group = grouped.get_group(i)
#
#                individual_group.to_csv(os.path.join(outpath, individual_group.iloc[0,0] + '.csv'), sep=',', index=False)

    #            print (individual_group)
    #            split_list.append(individual_group)

            # same as above
#            split_list = [grouped.get_group(x) for x in grouped.groups]


if __name__ == "__main__":
    import os
    import time

    fn = os.path.join(os.path.dirname(__file__), '..', '..',)
    os.chdir(fn)
    print 'The current working directory is ' + os.getcwd()

    print "Merging Marc's cities with their lat/long coordinates."
    time.sleep(3)
    directory = "data/raw/mtjj_jss"
    outpath = "data/raw/mtjj_jss"
    join_johnson_csv(directory, outpath)

    print "Splitting Marc's cities"
    time.sleep(3)
    csv = '20_cities_cyano_poplevel_lat_lon'
    outpath = "data/raw/mtjj_jss/mtjj_popMeans/"
    split_cities(directory, csv, outpath)

    print "Splitting James' cities"
    csv = 'Santangelo_AllCities_AllPlants'
    outpath = "Data/raw/mtjj_jss/jss_splitCities/"
    split_cities(directory, csv, outpath)
