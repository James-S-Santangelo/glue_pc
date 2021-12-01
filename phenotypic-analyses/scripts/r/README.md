R scripts used to analyse environmental and phenotypic data

[main.R](./main.R) calls all other scripts in the correct order and can more or less be run from top to bottom. Here is a brief description of the scripts in each subfolder. Detailed documentation can be found in each of the individual scripts.

- [analyses](./analyses): Perform the analyses presented in the paper. There include (1) Analyzing the environmental data, (2) Analyzing the phenotypic data (i.e., clines), (3) Predicting the clines from the environmental data
- [data-extraction](./data-extraction): Used to extract population-level % impervious surface and Human Influence Index (HII).
- [data-processing](./data-processing): Used to summarise and clean raw data for use in downstream analyses.
- [figures-tables](./figures-tables): Generate figures and tables for main text and supplement
- [misc](./misc): Miscellaneous scripts -- functions and descriptive statistics.