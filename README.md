<a href="http://www.bioconductor.org/packages/release/bioc/html/regionReport.html#since"><img border="0" src="http://www.bioconductor.org/shields/years-in-bioc/regionReport.svg" title="How long since the package was first in a released Bioconductor version (or is it in devel only)."></a> <a href="http://bioconductor.org/packages/stats/bioc/regionReport.html"><img border="0" src="http://www.bioconductor.org/shields/downloads/regionReport.svg" title="Percentile (top 5/20/50% or 'available') of downloads over last 6 full months. Comparison is done across all package categories (software, annotation, experiment)."></a> <a href="https://support.bioconductor.org/t/regionReport/"><img border="0" src="http://www.bioconductor.org/shields/posts/regionReport.svg" title="Support site activity, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts."></a> <a href="http://www.bioconductor.org/packages/release/bioc/html/regionReport.html#svn_source"><img border="0" src="http://www.bioconductor.org/shields/commits/bioc/regionReport.svg" title="average Subversion commits (to the devel branch) per month for the last 6 months"></a>

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


You may also want to consult the [vignette](http://www.bioconductor.org/packages/regionReport) for `regionReport` as well as the supplementary website [regionReportSupp](http://leekgroup.github.io/regionReportSupp/).

# Installation instructions

Get R 3.2.2 from [CRAN](http://cran.r-project.org/) as well as pandoc 
1.12.3 or newer from [here](http://johnmacfarlane.net/pandoc/installing.html).

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

Collado-Torres L, Jaffe AE and Leek JT (2015). _regionReport: Generate HTML reports for exploring a set of regions_.
https://github.com/leekgroup/regionReport - R package version 1.5.1, <URL:
http://www.bioconductor.org/packages/release/bioc/html/regionReport.html>.


A BibTeX entry for LaTeX users is

@Manual{,
    title = {regionReport: Generate HTML reports for exploring a set of regions},
    author = {Leonardo Collado-Torres and Andrew E. Jaffe and Jeffrey T. Leek},
    year = {2015},
    url = {http://www.bioconductor.org/packages/release/bioc/html/regionReport.html},
    note = {https://github.com/leekgroup/regionReport - R package version 1.5.1},
}

# Testing

Testing on Bioc-devel is feasible thanks to [r-builder](https://github.com/metacran/r-builder) as well as Bioconductor's nightly build.
