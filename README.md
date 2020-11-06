# Overview

This repository contains the code needed to reproduce the results of “Bayesian Specification Curve Analysis” (Semken and Rossell 2020).

To reproduce the analysis, first follow the instructions in `code/bsca_preanalysis.Rmd` to construct the dataset and perform a pre-analysis.  Then obtain the main results and create the supplementary information by compiling `code/bsca_supplement.Rmd`.  Make sure to have have the below requirements installed.

# Requirements

Main software:
- R 4.0.3
- Pandoc 2.10.1
- texlive 2019 (with packages float, placeins, cleveref)

R packages:
```
mombf 2.4.0
BAS 1.5.5
tidyverse 1.3.0
cowplot 1.1.0
gridGraphics 0.5.0
kableExtra 1.2.1
specr 0.2.1
knitr 1.30
rmarkdown 2.5
```

