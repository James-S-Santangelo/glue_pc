# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 12:18:26 2018
Modified by Jameson Thursday Jan. 10, 2018

@author: Alexander Tong
@contributor: James S Santangelo

Developed and tested with Python 2.7.15

There are many Python 2 specific syntax that are different Python 3 usage; \n\
(e.g., handling of unicode)

***ONLY IN PYTHON 2xx IS IT REQUIRED FOR UNICODE READ SUCH THAT u'string'


Original CSVs that should be modified with changes explained:

Baton_Rouge.csv ... lon val 91.13364 should be changed to -91.13364
Providence.csv ... lon val -7142947 should be changed to -71.42947
Madison.csv ... lat val 113.08156 should be changed to 43.08156
"""

import pandas as pd
import numpy as np
import os
import traceback
import re
import subprocess
import time


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
        if all(no_nan[no_nan.columns[index]].str.contains('^((?=.*°)(?=.*\')(?=.*\")|(?=.*_))', regex=True)):
            print "Coordinates are in degrees, minutes, seconds. Converting to decimal degrees"
            print ""

            # Only convert if row value is string. If NaN, leave as is.
            df[df.columns[index]] = df[df.columns[index]].apply(check_str)

        # Convert to string, strip, replace commas and alpha numeric characters.
        df[df.columns[index]] = df[df.columns[index]].astype(str)
        df[df.columns[index]].str.strip()
        df[df.columns[index]] = df[df.columns[index]].str.replace(',', '.')
        df[df.columns[index]] = df[df.columns[index]].str.replace(r'[^a-zA-Z\d\s:\.\-]', '')

        # Remove north and east identifiers in columns, if present
        if any(no_nan[no_nan.columns[index]].str.contains('^(?=.*(N|E))', regex=True)):
            df[df.columns[index]] = df[df.columns[index]].str.replace('[A-Za-z]', '')
            df[df.columns[index]] = pd.to_numeric(df[df.columns[index]])

        # Remove south and west identifiers in columns, if present. Convert to negative.
        elif any(no_nan[no_nan.columns[index]].str.contains('^(?=.*(S|W))', regex=True)):
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
    df[df.columns[index]] = df[df.columns[index]].round(5)

    return df[df.columns[index]]


def remove_special_quirks(df, file):
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
        indices = df['habitat_type'].str[0].str.lower() == 'p'
        df['habitat_type'][indices] = "Periurban"
        indices = df['habitat_type'].str[0].str.lower() == 's'
        df['habitat_type'][indices] = "Suburban"
        indices = df['habitat_type'].str[0].str.lower() == 'r'
        df['habitat_type'][indices] = "Rural"
        indices = df['habitat_type'].str[0].str.lower() == 'u'
        df['habitat_type'][indices] = "Urban"

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

                # check for ',' delimiter
                if ',' in first_line:

                    df = pd.read_csv(inpath + '/' + file, sep=",", skip_blank_lines=True)

                # check for ';' delimiter
                elif ';' in first_line:

                    df = pd.read_csv(inpath + '/' + file, sep=";", skip_blank_lines=True)

                print "Dataframe read. Standardizing column headers"
                df = fix_columns(df)

                print "Standardizing city name in 'city' column"
                dataframe_city_name(file, df)

                print "Processing latitude and longitude columns"
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
                    df[df.columns[lat_index]] = df[df.columns[lat_index]].round(5)

                if df[df.columns[lon_index]].dtype != np.number:
                    df[df.columns[lon_index]] = lat_lon_conditionals(df, lon_index)
                else:
                    df[df.columns[lon_index]] = df[df.columns[lon_index]].round(5)
                # print df[["population_latitude", "population_longitude"]].head()

#                  # wrap around another function if this conditional becomes reoccurrent...
                if u'Los_Angeles.csv' in file or u'La_Verne.csv' in file:
                    print "Swapping lat longs for " + file
                    # swap by col at index
                    df_list = list(df.columns.values)

                    # swap lat and lon
                    df = df.reindex(columns=swap_val(df_list, lat_index, lon_index))

                    # rename lon, lat back to correct lat, lon
                    df.rename(columns={df.columns[lat_index]: df.columns[lon_index],
                                       df.columns[lon_index]: df.columns[lat_index]}, inplace=True)

                    # make longitude values negative (west of GMT)
                    df[df.columns[lon_index]] = df[df.columns[lon_index]].apply(negative)

                if u'Burlington.csv' in file:

                    print "Changing longitude to negative for " + file
                    # make longitude values negative (west of GMT)
                    df[df.columns[lon_index]] = df[df.columns[lon_index]].apply(negative)

                if u'Lake_Charles.csv' in file:

                    print "Fixing HCN column for " + file
                    df = df.replace({'hcn_result': {'Positive': 1, 'Negative': 0}})
                    df['hcn_result'] = pd.to_numeric(df['hcn_result'])

                if u'Louis.csv' in file:

                    print "Fixing HCN column for " + file

                    df = df.replace({'hcn_result': {'positive': 1, 'negative': 0}})
                    df['hcn_result'] = pd.to_numeric(df['hcn_result'])

                if u'NewHaven.csv' in file:

                    print "Fixing Longitude for " + file

                    df[df.columns[lon_index]] = df[df.columns[lon_index]].apply(negative)

                print "Handling specific edge cases"
                df = remove_special_quirks(df, file)

                print "Forward filling rows, if necessary"
                fill_cols = ['city', 'population', 'plant', 'plate_no', 'population_latitude',
                             'population_longitude', 'habitat_type']
                df[fill_cols] = df[fill_cols].ffill()

                print "Standardizing factor levels for habitat type"
                df = standardize_habitat_type(df)

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


if __name__ == "__main__":

    inpaths = ['GLUE_Datasets/', 'data/raw/mtjj_jss/jss_splitCities/']
    outpath = 'data/clean/individualPlant_allCities/'

    fn = os.path.join(os.path.dirname(__file__), '..', '..',)
    os.chdir(fn)
    print 'The current working directory is ' + os.getcwd()

    print "Creating output directories"
    utf8_directory = "data/raw/utf8_encoded/"
    create_directory(utf8_directory)
    create_directory(outpath)

    for path in inpaths:

        # Coversion to UTF8 done using the 'iconv' command-line utility
        print "Converting all CSVs from " + path + " to UTF8 format"
        time.sleep(1)
        utf8_command = "sh scripts/shell/convert_to_utf8.sh " + path + " " + utf8_directory
        subprocess.call(utf8_command, shell=True)
        print "Done converting CSVs to UTF8. Converted datasets at in " + utf8_directory

    print "Changing permission to allow overwrite, if necessary"
    permissions_command = "chmod ugo+w " + outpath + "*.csv"
    subprocess.call(permissions_command, shell=True)

    print "Will now begin cleaning and standardizing all UTF8-encoded CSVs"
    time.sleep(5)
    process_csv(utf8_directory, outpath)
    print "Done processing CSVs from " + path
    time.sleep(1)

    print "Changing permissions of datasets to read only."
    print time.sleep(1)

    permissions_command = "chmod ugo-w " + outpath + "*.csv"
    subprocess.call(permissions_command, shell=True)
