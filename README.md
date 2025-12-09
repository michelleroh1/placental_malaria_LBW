# placental_malaria_LBW

## Overview

This repository provides Stata and R code for the manuscript, "Quantifying the contributions of active and past placental malaria infection to low birthweight in a high transmission area of Uganda". 

Datasets needed to reproduce results are available at https://github.com/Grantdorsey/DPSP-study-data-files:   
      
For any additional questions, please reach out to Michelle Roh (rohmi@ohsu.edu). 

## To run analysis:

  1. Download the following raw data files from https://github.com/Grantdorsey/DPSP-study-data-files:   
       a. `DPSP enrollment analysis database_FINAL.dta`   
       b. `DPSP individual level analysis database_FINAL.dta`   
       c. `DPSP delivery analysis database_FINAL.dta`  
  
  2. Install the following packages:   
       a. In Stata:
           <pre> ```ssc install hotdeckvar``` </pre>
       b. In R:
           <pre> ```install.packages(c('openxlsx', 'readstata13', 'dplyr', 'tidyr', 'ggplot2’, ’tibble’, ‘purrr’, ‘boot’, ‘Evalue’, ‘logistf’))``` </pre>   
  3. In each Stata do-file or R script and update the project path directory under “User inputs”.
  4. Code is numbered numerically. You must run Stata do-files (1-3) to generate the full analytic database prior to running analyses.
