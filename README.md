<a href="http://www.bioconductor.org/packages/release/bioc/html/regionReport.html#since"><img border="0" src="http://www.bioconductor.org/shields/years-in-bioc/regionReport.svg" title="How long since the package was first in a released Bioconductor version (or is it in devel only)."></a> <a href="https://bioconductor.org/packages/stats/bioc/regionReport/"><img border="0" src="http://www.bioconductor.org/shields/downloads/regionReport.svg" title="Percentile (top 5/20/50% or 'available') of downloads over last 6 full months. Comparison is done across all package categories (software, annotation, experiment)."></a> <a href="https://support.bioconductor.org/t/regionReport/"><img border="0" src="http://www.bioconductor.org/shields/posts/regionReport.svg" title="Support site activity, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts."></a> <a href="http://www.bioconductor.org/packages/release/bioc/html/regionReport.html#svn_source"><img border="0" src="http://www.bioconductor.org/shields/commits/bioc/regionReport.svg" title="average Subversion commits (to the devel branch) per month for the last 6 months"></a>

Status: Travis CI [![Build Status](https://travis-ci.org/leekgroup/regionReport.svg?branch=master)](https://travis-ci.org/leekgroup/regionReport),
Bioc-release <a href="http://www.bioconductor.org/packages/release/bioc/html/regionReport.html#archives"><img border="0" src="http://www.bioconductor.org/shields/availability/release/regionReport.svg" title="Whether the package is available on all platforms; click for details."></a> <a href="http://bioconductor.org/checkResults/release/bioc-LATEST/regionReport/"><img border="0" src="http://www.bioconductor.org/shields/build/release/bioc/regionReport.svg" title="build results; click for full report"></a>,
Bioc-devel <a href="http://www.bioconductor.org/packages/devel/bioc/html/regionReport.html#archives"><img border="0" src="http://www.bioconductor.org/shields/availability/devel/regionReport.svg" title="Whether the package is available on all platforms; click for details."></a> <a href="http://bioconductor.org/checkResults/devel/bioc-LATEST/regionReport/"><img border="0" src="http://www.bioconductor.org/shields/build/devel/bioc/regionReport.svg" title="build results; click for full report"></a>.

Bioc-release <a href="https://bioconductor.org/developers/how-to/unitTesting-guidelines/#coverage"><img border="0" src="http://www.bioconductor.org/shields/coverage/release/regionReport.svg" title="Test coverage percentage, or 'unknown'"></a>, Bioc-devel <a href="https://bioconductor.org/developers/how-to/unitTesting-guidelines/#coverage"><img border="0" src="http://www.bioconductor.org/shields/coverage/devel/regionReport.svg" title="Test coverage percentage, or 'unknown'"></a>

regionReport
===============

Generate HTML reports for a set of regions such as those from 
[derfinder](http://www.bioconductor.org/packages/derfinder) results or any other pipeline that defines a set of genomic regions.

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

For DESeq2 or edgeR results check _DESeq2Report()_ and _edgeReport()_.


You may also want to consult the [vignette](http://www.bioconductor.org/packages/regionReport) for `regionReport` as well as the supplementary website [regionReportSupp](http://leekgroup.github.io/regionReportSupp/).

# Installation instructions

Get R 3.3.x from [CRAN](http://cran.r-project.org/) as well as pandoc 
1.17.0.3 or newer from [here](http://johnmacfarlane.net/pandoc/installing.html).

```R
## From Bioconductor
source('http://bioconductor.org/biocLite.R')
biocLite('regionReport')
```

# Vignette

The vignette for this package can be viewed via [Bioconductor's website](http://www.bioconductor.org/packages/regionReport).


# Citation

Below is the citation output from using `citation('regionReport')` in R. 
Please run this yourself to check for any updates on how to cite 
__regionReport__.


To cite package __regionReport__ in publications use:

Collado-Torres L, Jaffe AE and Leek JT (2016). “regionReport: Interactive reports for region-level and feature-level genomic analyses [version2; referees: 2 approved, 1 approved with reservations].” _F1000Research_, *4*, pp. 105. doi: 10.12688/f1000research.6379.2 (URL:
http://doi.org/10.12688/f1000research.6379.2), <URL: http://f1000research.com/articles/4-105/v2>.


A BibTeX entry for LaTeX users is

@Article{,
    author = {Leonardo Collado-Torres and Andrew E. Jaffe and Jeffrey T. Leek},
    title = {regionReport: Interactive reports for region-level and feature-level genomic analyses [version2; referees: 2 approved, 1 approved with reservations]},
    journal = {F1000Research},
    year = {2016},
    doi = {10.12688/f1000research.6379.2},
    url = {http://f1000research.com/articles/4-105/v2},
    volume = {4},
    pages = {105},
}

# Testing

Testing on Bioc-devel is feasible thanks to [R Travis](http://docs.travis-ci.com/user/languages/r/) as well as Bioconductor's nightly build.
