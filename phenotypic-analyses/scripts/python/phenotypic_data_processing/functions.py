# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 12:18:26 2018
Modified by Jameson Thursday Jan. 10, 2018

@author: Alex
@contributor: James S Santangelo

Developed and tested with Python 2.7.15

There are many Python 2 specific syntax that are different Python 3 usage; \n\
(e.g., handling of unicode)

***ONLY IN PYTHON 2xx IS IT REQUIRED FOR UNICODE READ SUCH THAT u'string'


Functions used for cleaning/processing CSV submitted by GLUE clollaborators.

Baton_Rouge.csv ... lon val 91.13364 should be changed to -91.13364
Providence.csv ... lon val -7142947 should be changed to -71.42947
Madison.csv ... lat val 113.08156 should be changed to 43.08156
"""

import pandas as pd
import numpy as np
import os
import traceback
import re


def create_directory(directory):
    """Creates directory if it does not exist

    Args:
        directory ('str'): Directory to be created

    Returns:
        None
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def negative(x):
    """Return the negative of input

    Args:
        x (numeber): Integer of float

    Returns:
        -x (number): Negative of input
    """
    return -x


def swap_val(a, i, j):
    '''
    Description:
        Swap a value at given index

    Args:
        a (list): specify list
        i (list[index]): specify index of element to be swapped in list
        j (list[index]): specify index of element to be swapped in list

    Returns:
        list with swapped elements at specified index.
    '''
    temp = a[i]
    a[i] = a[j]
    a[j] = temp

    return a


def DMS_to_DD(string):
    '''
    Description:
        convert degrees, minutes, seconds (DMS) notation to decimal degrees (DD) notation.

        degrees = degrees
        minutes = minutes/60
        seconds = seconds/3600

    Args:
       string (str): DMS notation from row in csv file

    Returns:
        decimal degrees notation (float)
    '''
    # break up DMS for transformation
    # dms = re.split('[°\'\"]+', string)
    dms = re.split('[^a-zA-Z\d\s:\.]+', string)

    degrees = float(dms[0])
    minutes = float(dms[1]) / 60
    seconds = float(dms[2]) / 3600

    try:
        # factor for latitude southern hemisphere
        if 'S' in dms[-1]:
            dd = (degrees + minutes + seconds) * -1
        # factor for longitude west of GMT
        elif 'W' in dms[-1]:
            dd = (degrees + minutes + seconds) * -1
#
        else:
            dd = (degrees + minutes + seconds)

    except Exception as E:
        print (E)
        pass

    return dd


def rename_col(df, abbrv, col_name):
    """Renames columns, if present, adds column, if absent

    Args:
        df (pandas dataframe): specify dataframe for processing
        abbrv (str): Substring found in column to be renamed
        col_name (str): New column name

    Returns:
        Dataframe columns (list): List with new column names to replace existing names.
    """

    # Find columns in dataframe with names containing substrings.
    if abbrv == 'pop':
        # Only required for abbrv == 'pop' since 'pop' is in additional column names
        columns = [col for col in df.columns if abbrv in col and 'latit' not in col and 'long' not in col]
    else:
        columns = [col for col in df.columns if abbrv in col]

#     print columns
    # Get column name and index. Rename column
    if columns:
        col = columns[0]
        index = df.columns.get_loc(col)
        df[df.columns[index]]
        df = df.rename(columns={df.columns[index]: col_name})
    # If no column found containing substring, add column
    else:
        df[col_name] = np.nan

    return df.columns


def fix_columns(df):
    """Converts headers to lover case, renames and orders columns

    Args:
        df (pandas dataframe): specify dataframe for processing

    Retruns:
        dataframe (pandas dataframe): Modified dataframe with reordered and renamed columns
    """

    # Colvert column names to lowercase
    df.columns = df.columns.str.lower()

    # Rename columns
    df.columns = rename_col(df, 'city', 'city')
    df.columns = rename_col(df, 'pop', 'population')
    df.columns = rename_col(df, 'plant', 'plant')
    df.columns = rename_col(df, 'latit', 'population_latitude')
    df.columns = rename_col(df, 'long', 'population_longitude')
    df.columns = rename_col(df, 'habitat', 'habitat_type')
    df.columns = rename_col(df, 'hcn', 'hcn_result')
    df.columns = rename_col(df, 'plate', 'plate_no')

    # Columns to keep
    standard_headers = ["city", "population", "plant", "plate_no",
                        "population_latitude", "population_longitude",
                        "habitat_type", "hcn_result"]

    # Subset df to keep only specific columns. Also reorders.
    df = df[standard_headers]

    return df


def lat_lon_conditionals(dataframe, index):
    '''
    Description:
        processes longitude coordinate value (decimal degrees) from string to float value

        By Python definition '\W == [^a-zA-Z0-9_], which excludes all numbers, letters and _

    Args:
        dataframe (pandas dataframe): specify dataframe for processing
        dataframe_index (int): specify dataframe index of longitude values to be processed

    Returns:
        processed longitude coordinate value (decimal degrees) from string to float value
    '''
    def check_str(x):
        if isinstance(x, str):
            return DMS_to_DD(x)
        else:
            return x

    try:
        df = dataframe

        # Create subsetted dataframe with no NANs
        no_nan = df[df[df.columns[index]].notnull()]

        # If all rows in lat/long column of subsetted data match regex patter
        # Works for all datasets in Degrees, Minutes, Seconds, inccluding Lanzhou
        if all(no_nan[no_nan.columns[index]].str.contains('^(?=.*°)(?=.*\')(?=.*\")|(?=.*_)', regex=True)):
            # print "Coordinates are in degrees, minutes, seconds. Converting to decimal degrees"
            # print ""
            # print no_nan[no_nan.columns[index]].head()
            # Only convert if row value is string. If NaN, leave as is.
            df[df.columns[index]] = df[df.columns[index]].apply(check_str)

        # Convert to string, strip, replace commas and alpha numeric characters.
        df[df.columns[index]] = df[df.columns[index]].astype(str)
        df[df.columns[index]].str.strip()
        df[df.columns[index]] = df[df.columns[index]].str.replace(',', '.')
        df[df.columns[index]] = df[df.columns[index]].str.replace(r'[^a-zA-Z\d\s:\.\-]', '')

        # Remove north and east identifiers in columns, if present
        if any(df[df.columns[index]].str.contains('^(?=.*N)|(?=.*E)', regex=True)):
            # print df[df.columns[index]].head()
            df[df.columns[index]] = df[df.columns[index]].str.replace('[A-Za-z]', '')
            df[df.columns[index]] = pd.to_numeric(df[df.columns[index]])

        # Remove south and west identifiers in columns, if present. Convert to negative.
        elif any(df[df.columns[index]].str.contains('^(?=.*S)|(?=.*W)', regex=True)):
            # print df[df.columns[index]].head()
            df[df.columns[index]] = df[df.columns[index]].str.replace('[A-Za-z]', '')
            df[df.columns[index]] = pd.to_numeric(df[df.columns[index]]).apply(negative)

    except Exception as exc:
        print df.loc[0][0]
        print traceback.format_exc()
        print exc
        pass

    # Round columns to 5 digits
    df[df.columns[index]] = pd.to_numeric(df[df.columns[index]], errors='coerce')
    # print(df[df.columns[index]].dtype)
    df[df.columns[index]] = df[df.columns[index]].round(7)

    return df[df.columns[index]]


def remove_special_quirks(df, file, lat_index, lon_index):
    """Handles specific edge cases found in certain datasets

    Args:
        df (pandas dataframe): specify dataframe for processing
        file (str): Name of file

    Returns:
        dataframe (pandas dataframe): Modified dataframe with quirks removed
    """

    # Remove quotes from Paris dataset
    if 'Paris' in file:
        df = df.replace('""', np.nan)

    # Remove pound sign from plate_no column in Landshut
    if 'Landshut' in file:
        df['plate_no'] = df['plate_no'].str.replace('#', '')
        df['plate_no'] = pd.to_numeric(df['plate_no'], errors='coerce')

    # Split Population column for Atlantic City into Population and Plant
    if 'Atlantic' in file:
        new = df["population"].astype(str).str.split('.', n=1, expand=True)
        df["population"] = new[0]
        df["plant"] = new[1]
        df['population'] = pd.to_numeric(df['population'], errors='coerce')
        df['plant'] = pd.to_numeric(df['plant'], errors='coerce')

        # Swap lat longs for Los Angeles and La Verne
    if u'Los_Angeles.csv' in file or u'La_Verne.csv' in file:
        df_list = list(df.columns.values)
        df = df.reindex(columns=swap_val(df_list, lat_index, lon_index))
        df.rename(columns={df.columns[lat_index]: df.columns[lon_index],
                           df.columns[lon_index]: df.columns[lat_index]}, inplace=True)
        df[df.columns[lon_index]] = df[df.columns[lon_index]].apply(negative)

    # Make lat longs negative for Burlington
    if u'Burlington.csv' in file:
        df[df.columns[lon_index]] = df[df.columns[lon_index]].apply(negative)

    # Change HCN columns for Quito and Lake Charles to numeric
    if u'Lake_Charles.csv' in file or u'Quito.csv' in file:
        df = df.replace({'hcn_result': {'Positive': 1, 'Negative': 0}})
        df['hcn_result'] = pd.to_numeric(df['hcn_result'])

    if u'Quito.csv' in file:
        df['plant'] = df['plant'].str[-2:]
        df['plant'] = pd.to_numeric(df['plant'], errors='coerce')

    # Change HCN columns for St Louis
    if u'Louis.csv' in file:
        df = df.replace({'hcn_result': {'positive': 1, 'negative': 0}})
        df['hcn_result'] = pd.to_numeric(df['hcn_result'])

    # Change NewHaven longitude to negative
    if u'NewHaven.csv' in file:
        df[df.columns[lon_index]] = df[df.columns[lon_index]].apply(negative)

    # Change Buenos Aires lat longs to negative
    if u'Buenos_Aires.csv' in file:
        df[df.columns[lon_index]] = df[df.columns[lon_index]].apply(negative)
        df[df.columns[lat_index]] = df[df.columns[lat_index]].apply(negative)

    if u'Fairbanks.csv' in file:
        df[df.columns[lon_index]] = df[df.columns[lon_index]].apply(negative)

    return df


def standardize_habitat_type(df):
    """Standardize habitat type factor levels

    Args:
        df (pandas dataframe): specify dataframe for processing

    Return:
        df (pandas dataframe): New dataframe with renamed factor levels
    """

    # Rename factor levels based on first letter of row value in habitat type column
    try:
        df.loc[df['habitat_type'].str[0].str.lower() == 'u', 'habitat_type'] = 'Urban'
        df.loc[df['habitat_type'].str[0].str.lower() == 'r', 'habitat_type'] = 'Rural'
        df.loc[df['habitat_type'].str[0].str.lower() == 'p', 'habitat_type'] = 'Periurban'
        df.loc[df['habitat_type'].str[0].str.lower() == 's', 'habitat_type'] = 'Suburban'

    # Continue even for empty columns. Print exception
    except Exception as E:
        print E
        pass

    return df


def dataframe_city_name(infile, dataframe):
    '''
    Description:
        Standardize string name for column 'city' for processed csv

    Args:
        infile (csv): specify csv to be processed
        dataframe (pandas DataFrame): specify dataframe

    Returns:
        processed cleaned up city name (str) for column in pandas dataframe
    '''
    df = dataframe
    file = infile

    # csv has no name in index 0 (col 1); fix to generalize to other cases?
    if file == 'Atlantic_City.csv':
        # print (file)
        df['city'] = file[: -4].title()

    elif file == 'NewHaven.csv':
        # print (file)
        df['city'] = "New_Haven"

    elif file == 'FortCollins.csv':
        # print (file)
        df['city'] = "Fort_Collins"

    # re-name city name
    elif '_' in file[: -4]:
        df.iloc[:, 0] = file[: -4].title()

    elif '-' in file[: -4]:
        df.iloc[:, 0] = file[: -4].replace('-', '_').title()

    elif ' ' in file[: -4]:
        df.iloc[:, 0] = file[: -4].replace(' ', '_').title()

    elif '.' in file[: -4]:
        df.iloc[:, 0] = file[: -4].replace('.', '_').title()

    elif 'Quebec' in file[: -4]:
        df.iloc[:, 0] = file[: -4].title() + '_City'

    else:
        df.iloc[:, 0] = file[0].upper() + file[1: -4].lower()

    return df.iloc[:, 0]


def clean_csv_outname(infile):
    '''
    Description:
        remove hyphens, space, periods and any addendums to city names

    Args:
        infile (csv): specify csv to be processed

    Returns:
        processed outname (city name) for csv files
    '''
    file = infile

    if '-' in file[: -4]:
        outname = file[: -4].replace('-', '_')
    elif ' ' in file[: -4]:
        outname = file[: -4].replace(' ', '_')
    elif '.' in file[: -4]:
        outname = file[: -4].replace('.', '_')
    elif 'Quebec' in file[: -4]:
        outname = file[: -4].title() + '_City'
    else:
        outname = file[: -4]

    return str(outname.title() + ".csv")


def process_csv(inpath, outpath):
    '''
    '''
    for file in os.listdir(inpath):
        # process csv
        if file.endswith('.csv'):
            print "------------------------------------------"
            print "NOW PROCESSING " + file
            print ""

            # Open file as read only
            with open(os.path.join(inpath, file), 'r') as f:
                first_line = f.readline()

                semiColon_count = first_line.count(';')
                comma_count = first_line.count(',')

                # print semiColon_count, comma_count
#                 check for ',' delimiter
                if comma_count > semiColon_count:
                    # print "Comma"
                    df = pd.read_csv(inpath + '/' + file, sep=",", skip_blank_lines=True)

                # check for ';' delimiter
                elif semiColon_count > comma_count:
                    # print "semi-colon"
                    df = pd.read_csv(inpath + '/' + file, sep=";", skip_blank_lines=True)

                # print "Dataframe read. Standardizing column headers"
                df = fix_columns(df)

                # print "Standardizing city name in 'city' column"
                dataframe_city_name(file, df)

                # print "Processing latitude and longitude columns"
                # get index of latitude and longitude columns
                lat_col_name = [col for col in df.columns if 'latit' in col][0]
                lat_index = df.columns.get_loc(lat_col_name)
                df = df.rename(columns={df.columns[lat_index]: 'population_latitude'})
                # print lat_index

                lon_col_name = [col for col in df.columns if 'long' in col][0]
                lon_index = df.columns.get_loc(lon_col_name)
                df = df.rename(columns={df.columns[lon_index]: 'population_longitude'})
                # print lon_index

                # print df[['population_latitude', 'population_longitude']].head()
                # print df[df.columns[lat_index]].dtype
                if df[df.columns[lat_index]].dtype != np.number:
                    # Process latitude if it's not already a number
                    df[df.columns[lat_index]] = lat_lon_conditionals(df, lat_index)
                else:
                    # If it is a number, round.
                    df[df.columns[lat_index]] = df[df.columns[lat_index]].round(7)

                if df[df.columns[lon_index]].dtype != np.number:
                    df[df.columns[lon_index]] = lat_lon_conditionals(df, lon_index)
                else:
                    df[df.columns[lon_index]] = df[df.columns[lon_index]].round(7)
                # print df[["population_latitude", "population_longitude"]].head()

                # print "Handling specific edge cases"
                df = remove_special_quirks(df, file, lat_index, lon_index)

                # print "Forward filling rows, if necessary"
                fill_cols = ['city', 'population', 'plant', 'plate_no', 'population_latitude',
                             'population_longitude', 'habitat_type']
                df[fill_cols] = df[fill_cols].ffill()

                # print "Standardizing factor levels for habitat type"
                df = standardize_habitat_type(df)

                # Remove non-binary numeric values from HCN column
                df["hcn_result"] = df["hcn_result"].replace('[^0-1]', np.nan, regex=True)

                try:

                    # Creat name of file to write to disk
                    outname = clean_csv_outname(file)
                    print "Writing: " + outname

#                     # index = False; remove pandas auto index from output
                    df.to_csv(os.path.join(outpath, outname), sep=',', index=False)

                except Exception as E:
                    print E
                    pass
