regionReport
===============

Generate HTML reports for a set of regions such as those from 
[derfinder](https://github.com/lcolladotor/derfinder) results.

Check the documentation for _generateReport()_ for an example on how to create 
the necessary input files and generating the HTML report for __derfinder__ 
results. Or use:

```R
example('generateReport', 'regionReport', ask=FALSE)
```

For a full example on how to use __derfinder__ and __regionReport__ check 
https://github.com/lcolladotor/derfinderExample.

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

## For BioC-release use:
#install_github('lcolladotor/derfinder@release')

## regionReport itself
install_github('lcolladotor/regionReport')
```

# Citation

Below is the citation output from using `citation('regionReport')` in R. 
Please run this yourself to check for any updates on how to cite 
__regionReport__.

---

To cite package __regionReport__ in publications use:

Leonardo Collado-Torres, Andrew Jaffe and Jeffrey Leek (2014). regionReport: 
Generate HTML reports for derfinder (https://github.com/lcolladotor/derfinder) 
results. R package version 0.0.18. 
https://github.com/lcolladotor/regionReport

A BibTeX entry for LaTeX users is

@Manual{,
    title = {regionReport: Generate HTML reports for derfinder
        (https://github.com/lcolladotor/derfinder) results},
    author = {Leonardo Collado-Torres and Andrew Jaffe and Jeffrey Leek},
    year = {2014},
    note = {R package version 0.0.18},
    url = {https://github.com/lcolladotor/regionReport},
}
