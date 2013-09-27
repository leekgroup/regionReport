derfinderReport
===============

Generate HTML reports for [derfinder](https://github.com/lcolladotor/derfinder) results.

For a full example on how to use __derfinderReport__ check https://github.com/lcolladotor/derfinderExample.

# Installation instructions

Get R 3.0.1 or newer from [CRAN](http://cran.r-project.org/).

```S
## If needed
install.packages("devtools")

## Pre-requisites from CRAN
install.packages(c("knitr", "ggplot2", "gridExtra", "data.table", "knitr", "knitcitations",
	"xtable", "RColorBrewer"))

## Pre-requisites from Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(c("IRanges", "GenomicRanges", "biovizBase", "ggbio", "TxDb.Hsapiens.UCSC.hg19.knownGene"))

## GitHub dependencies
library(devtools)
install_github("rCharts", "ramnathv", ref="dev")
install_github(username="jimhester", repo="knitrBootstrap")
install_github("derfinder", "lcolladotor")

## derfinderReport itself
library(devtools)
install_github("derfinderReport", "lcolladotor")
```

# Citation

Below is the citation output from using `citation("derfinderReport")` in R. Please run this yourself to check for any updates on how to cite __derfinderReport__.

---

To cite package __derfinderReport__ in publications use:

Leonardo Collado-Torres, Andrew Jaffe and Jeffrey Leek (2013). derfinderReport: Generate HTML reports for derfinder (https://github.com/lcolladotor/derfinder) results. R package version 0.0.1. https://github.com/lcolladotor/derfinderReport

A BibTeX entry for LaTeX users is

@Manual{, title = {derfinderReport: Generate HTML reports for derfinder (https://github.com/lcolladotor/derfinder) results}, author = {Leonardo Collado-Torres and Andrew Jaffe and Jeffrey Leek}, year = {2013}, note = {R package version 0.0.1}, url = {https://github.com/lcolladotor/derfinderReport}, }
