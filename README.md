
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regionReport

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/regionReport.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/regionReport)
[![R build
status](https://github.com/leekgroup/regionReport/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/leekgroup/regionReport/actions)
<!-- badges: end -->

Generate HTML reports for a set of regions such as those from
*[derfinder](https://bioconductor.org/packages/3.11/derfinder)* results
or any other pipeline that defines a set of genomic regions.

Check the documentation for `derfinderReport()` for an example on how to
create the necessary input files and generating the HTML report for
*[derfinder](https://bioconductor.org/packages/3.11/derfinder)* results.
Or use:

``` r
example('derfinderReport', 'regionReport', ask=FALSE)
```

Similarly, check `renderReport()` for an example of a general report, or
use:

``` r
example('renderReport', 'regionReport', ask=FALSE)
```

For *[DESeq2](https://bioconductor.org/packages/3.11/DESeq2)* or
*[edgeR](https://bioconductor.org/packages/3.11/edgeR)* results check
`DESeq2Report()` and `edgeReport()`.

## Documentation

For more information about `derfinderPlot` check the vignettes [through
Bioconductor](http://bioconductor.org/packages/regionReport) or at the
[documentation website](http://leekgroup.github.io/regionReport).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `regionReport` using
from [Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("regionReport")
```

## Citation

Below is the citation output from using `citation('regionReport')` in R.
Please run this yourself to check for any updates on how to cite
**regionReport**.

``` r
print(citation('regionReport'), bibtex = TRUE)
#> 
#> Collado-Torres L, Jaffe AE, Leek JT (2016). "regionReport: Interactive
#> reports for region-level and feature-level genomic analyses [version2;
#> referees: 2 approved, 1 approved with reservations]." _F1000Research_,
#> *4*, 105. doi: 10.12688/f1000research.6379.2 (URL:
#> https://doi.org/10.12688/f1000research.6379.2), <URL:
#> http://f1000research.com/articles/4-105/v2>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     author = {Leonardo Collado-Torres and Andrew E. Jaffe and Jeffrey T. Leek},
#>     title = {regionReport: Interactive reports for region-level and feature-level genomic analyses [version2; referees: 2 approved, 1 approved with reservations]},
#>     journal = {F1000Research},
#>     year = {2016},
#>     doi = {10.12688/f1000research.6379.2},
#>     url = {http://f1000research.com/articles/4-105/v2},
#>     volume = {4},
#>     pages = {105},
#>   }
#> 
#> Collado-Torres L, Nellore A, Frazee AC, Wilks C, Love MI, Langmead B,
#> Irizarry RA, Leek JT, Jaffe AE (2017). "Flexible expressed region
#> analysis for RNA-seq with derfinder." _Nucl. Acids Res._. doi:
#> 10.1093/nar/gkw852 (URL: https://doi.org/10.1093/nar/gkw852), <URL:
#> http://nar.oxfordjournals.org/content/early/2016/09/29/nar.gkw852>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Flexible expressed region analysis for RNA-seq with derfinder},
#>     author = {Leonardo Collado-Torres and Abhinav Nellore and Alyssa C. Frazee and Christopher Wilks and Michael I. Love and Ben Langmead and Rafael A. Irizarry and Jeffrey T. Leek and Andrew E. Jaffe},
#>     year = {2017},
#>     journal = {Nucl. Acids Res.},
#>     doi = {10.1093/nar/gkw852},
#>     url = {http://nar.oxfordjournals.org/content/early/2016/09/29/nar.gkw852},
#>   }
#> 
#> Collado-Torres L, Jaffe AE, Leek JT (2017). _regionReport: Generate
#> HTML or PDF reports for a set of genomic regions or DESeq2/edgeR
#> results_. doi: 10.18129/B9.bioc.regionReport (URL:
#> https://doi.org/10.18129/B9.bioc.regionReport),
#> https://github.com/leekgroup/regionReport - R package version 1.21.7,
#> <URL: http://www.bioconductor.org/packages/regionReport>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {regionReport: Generate HTML or PDF reports for a set of genomic regions or DESeq2/edgeR results},
#>     author = {Leonardo Collado-Torres and Andrew E. Jaffe and Jeffrey T. Leek},
#>     year = {2017},
#>     url = {http://www.bioconductor.org/packages/regionReport},
#>     note = {https://github.com/leekgroup/regionReport - R package version 1.21.7},
#>     doi = {10.18129/B9.bioc.regionReport},
#>   }
```

Please note that the `regionReport` was only made possible thanks to
many other R and bioinformatics software authors, which are cited either
in the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the regionReport project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Development tools

  - Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*,
    *[sysreqs](https://github.com/r-hub/sysreqs)* and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.11/BiocCheck)*.
  - Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
  - The [documentation
    website](http://leekgroup.github.io/derfinderPlot) is automatically
    updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
  - The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
  - The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.
