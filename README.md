regionReport
===============

Generate HTML reports for a set of regions such as those from 
[derfinder](https://github.com/lcolladotor/derfinder) results.

Check the documentation for _derfinderReport()_ for an example on how to create 
the necessary input files and generating the HTML report for __derfinder__ 
results. Or use:

```R
example('derfinderReport', 'regionReport', ask=FALSE)
```

You may also want to consult the [vignette](http://lcolladotor.github.io/regionReport/) for `regionReport`.

# Installation instructions

Get R 3.1.1 or newer from [CRAN](http://cran.r-project.org/) as well as pandoc 
1.12.3 or newer from [here](http://johnmacfarlane.net/pandoc/installing.html).

```R
## If needed
install.packages('devtools')

## Pre-requisites from CRAN
install.packages(c('knitr', 'ggplot2', 'gridExtra', 'data.table',
    'knitcitations', 'RColorBrewer', 'mgcv', 'rmarkdown', 'xtable',
    'knitrBootstrap'))
    
## You might need to install the Cairo package, but that will depend on the device 
## you use for the plots
install.packages('Cairo')

## Pre-requisites from Bioconductor
source('http://bioconductor.org/biocLite.R')
biocLite(c('IRanges', 'GenomicRanges', 'biovizBase', 'ggbio', 
    'TxDb.Hsapiens.UCSC.hg19.knownGene', 'GenomeInfoDb'))

## GitHub dependencies
library('devtools')
install_github('ramnathv/rCharts')
## If you use the GitHub version of knitrBoostrap you get better looking reports
# install_github('jimhester/knitrBootstrap')

## For BioC-devel use:
install_github('lcolladotor/derfinder@master')
install_github('lcolladotor/derfinderPlot')

## regionReport itself
install_github('lcolladotor/regionReport')
```

# Citation

Below is the citation output from using `citation('regionReport')` in R. 
Please run this yourself to check for any updates on how to cite 
__regionReport__.

---

To cite package __regionReport__ in publications use:

Leonardo Collado-Torres, Andrew Jaffe and Jeffrey Leek (2014). regionReport: Generate HTML reports for exploring a set of regions. R package version 0.99.0. https://github.com/lcolladotor/regionReport

A BibTeX entry for LaTeX users is

@Manual{,
    title = {regionReport: Generate HTML reports for exploring a set of regions},
    author = {Leonardo Collado-Torres and Andrew Jaffe and Jeffrey Leek},
    year = {2014},
    note = {R package version 0.99.0},
    url = {https://github.com/lcolladotor/regionReport},
}

# Travis CI

This package is automatically tested thanks to [Travis CI](travis-ci.org) and [r-travis](https://github.com/craigcitro/r-travis). If you want to add this to your own package use:

```R
## Use devtools to create the .travis.yml file
library('devtools')
use_travis('yourPackage')

## Read https://github.com/craigcitro/r-travis/wiki to configure .travis.yml appropriately

## Add a status image by following the info at http://docs.travis-ci.com/user/status-images/
```
