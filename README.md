# Overview

This repository contains the code needed to reproduce the results of “Analysis of multiple specifications: statistical validity and a Bayesian proposal” (Semken and Rossell 2021).

Preprint (pdf): https://osf.io/cahyq

Online supplement (pdf): https://osf.io/7fsqu

To use BSCA, download `code/functions.R` and follow the instructions in the online supplement.

To reproduce the analysis, first follow the instructions in `code/bsca_preanalysis.Rmd` to construct the dataset and perform a pre-analysis.  Then obtain the main results and create the online supplement by running `code/bsca_supplement.Rmd`.  Make sure to have the below requirements installed.

# Requirements

Main software:
- R 4.0.3
- Pandoc 2.10.1
- texlive 2019 (with packages float, placeins, cleveref)
- mombf R package compiled from latest [source code](https://github.com/davidrusi/mombf) (e.g. using `devtools::install_github("davidrusi/mombf")`)

R packages from CRAN:
```
BAS 1.5.5
tidyverse 1.3.0
cowplot 1.1.0
gridGraphics 0.5.0
kableExtra 1.2.1
specr 0.2.1
knitr 1.30
rmarkdown 2.5
```

