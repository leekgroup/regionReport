derfinderReport
===============

Generate HTML reports for [derfinder](https://github.com/lcolladotor/derfinder) 
results.

Check the documentation for _generateReport()_ for an example on how to create 
the necessary input files and generating the HTML report for __derfinder__ 
results. Or use:

```S
example("generateReport", "derfinderReport", ask=FALSE)
```

For a full example on how to use __derfinder__ and __derfinderReport__ check 
https://github.com/lcolladotor/derfinderExample.

# Installation instructions

Get R 3.1.0 or newer from [CRAN](http://cran.r-project.org/) as well as pandoc 
1.12.3 or newer from [here](http://johnmacfarlane.net/pandoc/installing.html).

```S
## If needed
install.packages("devtools")

## Note that knitr 1.5.27 is currently not available via CRAN
## Following https://github.com/yihui/knitr
install.packages('knitr', repos = c('http://rforge.net', 'http://cran.rstudio.org'),
    type = 'source')

## Pre-requisites from CRAN
install.packages(c("xtable", "ggplot2", "gridExtra", "data.table", "knitcitations",
    "RColorBrewer", "mgcv"))
    
## You might need to install the Cairo package, but that will depend on the device 
## you use for the plots
install.packages("Cairo")

## Pre-requisites from Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(c("IRanges", "GenomicRanges", "biovizBase", "ggbio", 
    "TxDb.Hsapiens.UCSC.hg19.knownGene"))

## GitHub dependencies
library("devtools")
install_github("rstudio/rmarkdown")
install_github("jimhester/knitrBootstrap")
install_github("ramnathv/rCharts", ref="dev")
install_github("lcolladotor/derfinder")

## derfinderReport itself
library("devtools")
install_github("lcolladotor/derfinderReport")
```

# Citation

Below is the citation output from using `citation("derfinderReport")` in R. 
Please run this yourself to check for any updates on how to cite 
__derfinderReport__.

---

To cite package __derfinderReport__ in publications use:

Leonardo Collado-Torres, Andrew Jaffe and Jeffrey Leek (2014). derfinderReport: 
Generate HTML reports for derfinder (https://github.com/lcolladotor/derfinder) 
results. R package version 0.0.16. 
https://github.com/lcolladotor/derfinderReport

A BibTeX entry for LaTeX users is

@Manual{,
    title = {derfinderReport: Generate HTML reports for derfinder
        (https://github.com/lcolladotor/derfinder) results},
    author = {Leonardo Collado-Torres and Andrew Jaffe and Jeffrey Leek},
    year = {2014},
    note = {R package version 0.0.16},
    url = {https://github.com/lcolladotor/derfinderReport},
}
