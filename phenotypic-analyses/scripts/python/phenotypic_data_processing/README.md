Python scripts to clean and standardize individual-palnt phenotype data.

Scripts need to be run in the order in which they are numbered. 

1. Split and clean individual-plant data collected by [Santangelo _et al._ 2020](https://onlinelibrary.wiley.com/doi/10.1002/evl3.163) and [Johnson _et al._ 2018](https://royalsocietypublishing.org/doi/10.1098/rspb.2018.1019)
2. Clean and standardize individual-plant datasets submitted by collaborators
3. Clean up HCN calls in collaborator datasets (if necessary). 

__NOTE__: These scripts are a mess and not general enough that they can be runn by anyone else without some decent effort. They are provided for transparency are are reasonably well documented to provide an overview of functionality. Briefly, the scripts perform some basic data cleaning and processing such as ensuring all datasets are UTF8-encoded, standardizing column names and types, ensuring consistent missing data encoding, etc. These scripts shouldn't need to be run; we will distribute the cleaned and standardized datasets with the manuscript to faciliitate downstream use.
