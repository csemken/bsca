# Overview

This repository contains the code needed to reproduce the results of “Analysis of multiple specifications: statistical validity and a Bayesian proposal” (Semken and Rossell 2021).

Preprint (pdf): https://osf.io/cahyq

Online supplement (pdf): https://osf.io/7fsqu

# Run R code

To reproduce the analysis, first follow the instructions in `code/bsca_preanalysis.Rmd` to construct the datasets. Second, run the code by completing the following steps: 

1. Set the environment variable `BSCA_PATH` to the absolute path of the project folder.
2. Run `renv/activate.R` and `renv::restore()` using R 4.0. This will install the right versions of all dependencies.
3. Follow the instructions in and compile `code/bsca_preanalysis.Rmd`, `code/bsca_supplement.Rmd` and `code/bma_simulation.Rmd` to reproduce the analysis.

On a bash shell run from the project directory the analysis can be reproduced using the following command:
```bash
export BSCA_PATH=$(pwd) 
Rscript "renv/activate.R"
Rscript -e "renv::restore()"
Rscript -e "setwd('code'); library(rmarkdown); render('bsca_preanalysis.Rmd'); render('bsca_supplement.Rmd'); render('bma_simulation.Rmd')"
```
