import functions
import os
import subprocess
import time
import pandas as pd


def replace_hcn_calls(clean_calls_path, dirtyCalls_path, outpath):

    for file in os.listdir(dirtyCalls_path):

        if file.endswith('.csv'):

            city_name = file.split('.')[0]
            # print(city_name)

            df_original = pd.read_csv(dirtyCalls_path  + '/' + file, sep=",", skip_blank_lines=True)
            # df_original['plant'] = df_original['plant'].astype(str)
            # df_original['population'] = df_original['population'].astype(str)

            # print(df_original.plant.dtype, df_original.population.dtype)

            clean_call_file = [clean_file for clean_file in os.listdir(clean_calls_path) if clean_file.endswith('.csv') and city_name in clean_file]

            if clean_call_file:

                clean_file = clean_call_file[0]

                df_clean = pd.read_csv(clean_calls_path + '/' + clean_file, sep=",", skip_blank_lines=True)

                df_clean = df_clean[['population', 'plant', 'hcn_result']]
                # df_clean['plant'] = df_clean['plant'].astype(str)
                # df_clean['population'] = df_clean['population'].astype(str)

                # print(df_clean.plant.dtype, df_clean.population.dtype)

                new_df = pd.merge(df_original, df_clean, on = ['population', 'plant'], how = 'left')
                new_df['hcn_result_x'] = new_df['hcn_result_y']
                del new_df['hcn_result_y']
                new_df.rename(columns={'hcn_result_x': 'hcn_result'}, inplace=True)

                # Drop duplicates required to handle Bogota dplication of population 40 plant 19
                new_df.drop_duplicates(keep='first', inplace=True)
                # Creat name of file to write to disk
                new_df.to_csv(os.path.join(outpath, file), sep=',', index=False)

            else:

                df_original.to_csv(os.path.join(outpath, file), sep=',', index=False)


if __name__ == "__main__":

    inpaths = ['/Users/jamessantangelo/Sync/Academia/Doctorate_PhD/Projects/in-progress/GLUE_Global.Urban.Clines.in.Cyanogenesis/data/collaborator_submissions/GLUE_Datasets/', 'data/raw/mtjj_jss/jss_splitCities/']
    # inpaths = ['test_df/']
    dirtyCalls_path = 'data/raw/collabDatasets_dirtyCalls/'
    # outpath = 'test_clean_df/'

    fn = os.path.join(os.path.dirname(__file__), '..', '..',)
    os.chdir(fn)
    print 'The current working directory is ' + os.getcwd()
    time.sleep(1)

    print "Creating output directories"
    time.sleep(1)
    utf8_directory = "data/raw/utf8_encoded/"
    functions.create_directory(utf8_directory)
    functions.create_directory(dirtyCalls_path)

    for path in inpaths:

        # Coversion to UTF8 done using the 'iconv' command-line utility
        print "Converting all CSVs from " + path + " to UTF8 format"
        time.sleep(1)
        utf8_command = "sh scripts/shell/convert_to_utf8.sh " + path + " " + utf8_directory
        subprocess.call(utf8_command, shell=True)
        print "Done converting CSVs to UTF8. Converted datasets in " + utf8_directory

    print "Changing permission to allow overwrite, if necessary"
    time.sleep(1)
    permissions_command = "chmod ugo+w " + dirtyCalls_path + "*.csv"
    subprocess.call(permissions_command, shell=True)

    print "Will now begin cleaning and standardizing all UTF8-encoded CSVs"
    time.sleep(1)
    functions.process_csv(utf8_directory, dirtyCalls_path)
    print "Done processing CSVs from " + path
    time.sleep(1)

    print "Replacing HCN calls, if necessary"
    clean_calls_path = "/Users/jamessantangelo/Sync/Academia/Doctorate_PhD/Projects/in-progress/GLUE_Global.Urban.Clines.in.Cyanogenesis/data/collaborator_submissions/GLUE_Datasets-clean/standardized/"
    outpath = "data/clean/individualPlant_allCities/"
    functions.create_directory(outpath)

    replace_hcn_calls(clean_calls_path, dirtyCalls_path, outpath)
