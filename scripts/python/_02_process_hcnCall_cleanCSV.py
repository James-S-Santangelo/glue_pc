
import functions
import os
import subprocess
import time


if __name__ == "__main__":

    inpaths = ['/Users/jamessantangelo/Sync/Academia/Doctorate_PhD/Projects/in-progress/GLUE_Global.Urban.Clines.in.Cyanogenesis/data/collaborator_submissions/GLUE_Datasets-clean/']
    # inpaths = ['test_df/']
    # outpath = 'GLUE_Datasets-clean/standardized/'
    outpath = '/Users/jamessantangelo/Sync/Academia/Doctorate_PhD/Projects/in-progress/GLUE_Global.Urban.Clines.in.Cyanogenesis/data/collaborator_submissions/GLUE_Datasets-clean/standardized/'

    fn = os.path.join(os.path.dirname(__file__), '..', '..',)
    os.chdir(fn)
    print 'The current working directory is ' + os.getcwd()
    time.sleep(1)

    print "Creating output directories"
    time.sleep(1)
    # utf8_directory = "data/raw/utf8_encoded/"
    # utf8_directory = "GLUE_Datasets-clean/utf8_encoded/"
    utf8_directory = "/Users/jamessantangelo/Sync/Academia/Doctorate_PhD/Projects/in-progress/GLUE_Global.Urban.Clines.in.Cyanogenesis/data/collaborator_submissions/GLUE_Datasets-clean/utf8_encoded/"
    functions.create_directory(utf8_directory)
    functions.create_directory(outpath)

    for path in inpaths:

        # Conversion to UTF8 done using the 'iconv' command-line utility
        print "Converting all CSVs from " + path + " to UTF8 format"
        time.sleep(1)
        utf8_command = "sh scripts/shell/convert_to_utf8.sh " + path + " " + utf8_directory
        subprocess.call(utf8_command, shell=True)
        print "Done converting CSVs to UTF8. Converted datasets in " + utf8_directory

    print "Changing permission to allow overwrite, if necessary"
    time.sleep(1)
    permissions_command = "chmod ugo+w " + outpath + "*.csv"
    subprocess.call(permissions_command, shell=True)

    print "Will now begin cleaning and standardizing all UTF8-encoded CSVs"
    time.sleep(1)
    functions.process_csv(utf8_directory, outpath)
    print "Done processing CSVs from " + path
    time.sleep(1)

    print "Changing permissions of datasets to read only."
    time.sleep(1)

    permissions_command = "chmod ugo-w " + outpath + "*.csv"
    subprocess.call(permissions_command, shell=True)
