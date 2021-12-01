This directory contains raw data in the following folders:

- [city_data](./city_data): City characteristics (e.g., Age, Size, etc.)
- [environmental_data](./environmental_data): One subdirectory for each of the 9 environmental variables we analysed. With the exception of GMIS and HII, extraction was done using [custom Python scripts](../scripts/python). GMIS was extracted using [this script](../../scripts/r/data-extraction/gmis_extraction.R), and HII using [this script](../scripts/r/data-extraction/hii_extraction.R)
- [mtjj_jss](./mtjj_jss): Population-mean datasets from [Johnson _et al._ 2018](https://royalsocietypublishing.org/doi/10.1098/rspb.2018.1019) and individual-plant-level data from [Santangelo _et al._ 2020](https://onlinelibrary.wiley.com/doi/10.1002/evl3.163)

In addition, this directory contains the following 2 files:

- [latLongs_cityCenters](./latLongs_cityCenters.csv): Latitude and Longitudes for city centres for each city
- [johnson_2018_plantsperPop.csv](./johnson_2018_plantsperPop.csv): Number of plants per population from [Johnson _et al._ 2018](https://royalsocietypublishing.org/doi/10.1098/rspb.2018.1019)
- [Johnson et al_20_city_clover 07.12.16_RAW.csv](./Johnson&#32et&#32al_20_city_clover&#3207.12.16_RAW.csv): Raw indnividual-plant data from [Johnson _et al._ 2018](https://royalsocietypublishing.org/doi/10.1098/rspb.2018.1019)