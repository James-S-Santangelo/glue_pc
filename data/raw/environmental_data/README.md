The folders in this directory contain the environmental data submitted by Alex Tong on February 10th, 2019.

The monthlyPET data is not currently being used in the analyses (decided at the meeting February 14 between James and Marc)

At the same meeting, we decided to include human population density, city size, and human population size as predictors in some of the analyses. These data have not yet been collected and were not part of Alex's contract.

The data need to be processed (e.g. forward filling) before merging with the population-mean datasets.

Check population indices of the following cities. They must have duplicates.

1. Atlantic City: 25 (fixed in raw individualPlant in GLUE_Datasets/)
2. Charlotte: 8 (fixed in raw individualPlant data in raw/data/mtjj_jss_data)
3. Halifax: 28 (fixed in raw individualPlant in GLUE_Datasets/)
4. Little Rock: 20 (fixed in raw individualPlant in GLUE_Datasets/)
5. Mexico City: 18 (fixed in raw individualPlant in GLUE_Datasets/)
6. Morelia: 17 (fixed in raw individualPlant in GLUE_Datasets/)
7. Ottawa: 22 (fixed in raw individualPlant in GLUE_Datasets/)
8. Stockholm: 43 (fixed in raw individualPlant in GLUE_Datasets/)
9. Newcastle: 11 and 37 (fixed in raw individualPlant in GLUE_Datasets/)

The duplicates were present in the original population mean datasets sent to Alex.

I regenerated the population mean datasets on March 3, 2019.

The environmental data for these cities will have to be regenerated so that duplucate-free population-mena datasets can be created.

Note two LST datasets for city of Toronto due to size of transect requiring multiple landstat images.

Data were re-downloaded on March 27 due to updates by Alex.

March 27 update: Regenerated population-mean datasets with environmental data. Some of the cities listed above still have duplicates, but only a single duplicate rather than many. Issue seems to be on Alex's end since origin pop-mean datasets are not duplicated.
