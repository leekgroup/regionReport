regionReport [![Build Status](https://travis-ci.org/leekgroup/regionReport.svg?branch=master)](https://travis-ci.org/leekgroup/regionReport)
===============

Generate HTML reports for a set of regions such as those from 
[derfinder](https://github.com/lcolladotor/derfinder) results or any other pipeline that defines a set of genomic regions.

Check the documentation for _derfinderReport()_ for an example on how to create 
the necessary input files and generating the HTML report for __derfinder__ 
results. Or use:

```R
example('derfinderReport', 'regionReport', ask=FALSE)
```

Similarly, check _renderReport()_ for an example of a general report, or use:

```R
example('renderReport', 'regionReport', ask=FALSE)
```


You may also want to consult the [vignette](http://leekgroup.github.io/regionReport/) for `regionReport` as well as the supplementary website [regionReportSupp](http://leekgroup.github.io/regionReportSupp/).

# Installation instructions

Get R 3.2.0 from [CRAN](http://cran.r-project.org/) as well as pandoc 
1.12.3 or newer from [here](http://johnmacfarlane.net/pandoc/installing.html).

```R
## From Bioconductor
source('http://bioconductor.org/biocLite.R')
biocLite('regionReport')
```

# Vignette

The vignette for this package can be viewed [here](http://leekgroup.github.io/regionReport/) or via [Bioconductor's website](http://www.bioconductor.org/packages/devel/bioc/html/regionReport.html).


# Citation

Below is the citation output from using `citation('regionReport')` in R. 
Please run this yourself to check for any updates on how to cite 
__regionReport__.


To cite package __regionReport__ in publications use:

Collado-Torres L, Jaffe AE and Leek JT (2015). _regionReport: Generate HTML reports for exploring a set of regions_.
https://github.com/leekgroup/regionReport - R package version 1.3.0, <URL:
http://www.bioconductor.org/packages/release/bioc/html/regionReport.html>.


A BibTeX entry for LaTeX users is

@Manual{,
    title = {regionReport: Generate HTML reports for exploring a set of regions},
    author = {Leonardo Collado-Torres and Andrew E. Jaffe and Jeffrey T. Leek},
    year = {2015},
    url = {http://www.bioconductor.org/packages/release/bioc/html/regionReport.html},
    note = {https://github.com/leekgroup/regionReport - R package version 1.3.0},
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

Testing on R-devel for Bioc-devel is feasible thanks to [r-builder](https://github.com/metacran/r-builder).
