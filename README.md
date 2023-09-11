# ML-ethanol-optimization

## Description 
This repository represents the machine learning (ML) scripts used in the publication 
"Optimizing Ethanol Production in Saccharomyces cerevisiae at Ambient and Elevated Temperatures through Machine Learning-Guided Combinatorial Promoter Modifications" 
(https://pubs.acs.org/doi/full/10.1021/acssynbio.3c00199)

Overall, we used XGBoost to train the ML model for predicting ethanol production using ethanol and promoter strengths as input.
We first searched for suitable ML algorithms using the caret package, and XGBoost was used for comprehensive model training.

All analyses were performed under R version 4.1.1

## Data availability
Raw data required for analysis are deposited in this repository and can be found in supplementary materials in the publication.

## R Packages used in the analysis
* [caret](https://topepo.github.io/caret/)
* [xgboost](https://xgboost.readthedocs.io/en/stable/)
* [ggplot2](https://ggplot2.tidyverse.org/)
* [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)
* [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
* [dplyr](https://github.com/tidyverse/dplyr)
* [readxl](https://readxl.tidyverse.org/)
