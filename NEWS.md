# regionReport 1.21.8

BUG FIXES

* In an attempt to fix the bug I reported at
https://stat.ethz.ch/pipermail/bioc-devel/2020-April/016645.html and in
preparation to reply to this thread, I'm experimenting with suggesting that
`magick` is installed, at which point `BiocStyle` will try to crop. Another
option might be to use `crop = NULL` on all Rmd files.

# regionReport 1.21.5

BUG FIXES

* The example in `renderReport()` was failing because the first time the
function was called, it ran `derfinder::makeGenomicState()` which in turn
uses `GenomicFeatures::isActiveSeq()`. `derfinder` version 1.21.5 fixed this
bug and is thus required by `regionReport` now.

# regionReport 1.21.4

SIGNIFICANT USER-VISIBLE CHANGES

* Add links to the example reports created when deploying the documentation
website with `pkgdown::deploy_to_branch()`. This eliminates the need for
the `vignetttes/Makefile` and the fake vignettes. It should also avoid
confusing `pkgdown`.

# regionReport 1.21.3

SIGNIFICANT USER-VISIBLE CHANGES

* Documentation website is now available at
http://leekgroup.github.io/derfinderPlot/. It gets updated with every commit on
the master branch (bioc-devel) using GitHub Actions and pkgdown.

# regionReport 1.21.2

SIGNIFICANT USER-VISIBLE CHANGES

* Now use `GenomeIndoDb::getChromInfoFromUCSC()` instead of
`data(hg19Ideogram, package = 'biovizBase')` to get the hg19 chromosome lengths.

# regionReport 1.19.2

SIGNIFICANT USER-VISIBLE CHANGES

* Added a `NEWS.md` file to track changes to the package.

# regionReport 1.19.2


SIGNIFICANT USER-VISIBLE CHANGES

* Renamed `load_install()` as `load_check()` as this function now only checks
that the package(s) was installed and returns an error if missing. The
error shows the user how to install the package(s) they are missing
instead of installing them automatically. This complies with Marcel
Ramos' request at https://github.com/leekgroup/recount/issues/14.


# regionReport 1.17.2


NEW FEATURES

* Add ORCID's following changes at
http://bioconductor.org/developers/package-guidelines/#description

# regionReport 1.15.4


BUG FIXES

* Fix a bug in the order that was reported and fixed by @bounlu at
https://github.com/leekgroup/regionReport/pull/9/


# regionReport 1.15.3


BUG FIXES

* Fixed an issue with `DESeq2Exploration.Rmd` that affected both
`DESeq2Report()` and `edgeReport()`. This should also fix the recount
bioc-release (3.7) and bioc-devel (3.8) branches.

* Fixed a `NAMESPACE` issue with `rmarkdown::html_document` and
`BiocStyle::html_document`


# regionReport 1.15.2


SIGNIFICANT USER-VISIBLE CHANGES

* Use BiocManager

# regionReport 1.15.1


BUG FIXES

* Fix namespace issue in relation to `BiocStyle::html_document2`

# regionReport 1.11.6


SIGNIFICANT USER-VISIBLE CHANGES

* Changed the default style to `BiocStyle::html_document` to mirror recent
changes in `BiocStyle`.


# regionReport 1.11.4


SIGNIFICANT USER-VISIBLE CHANGES

* Vignette now uses the new `BiocStyle::html_document` that was recently
released.


# regionReport 1.11.2


BUG FIXES

* Fixed the citations.
* Fixed `DESeq2Report()` for `limma`-results so that it will properly cite
`limma`.

# regionReport 1.11.1


SIGNIFICANT USER-VISIBLE CHANGES

* `DESeq2Report()` can now be used with other software if their results are
 made to look like `DESeq2` results. For example, with `limma`-voom results.

BUG FIXES

* Made the `DESeq2Report()` more robust in case` rlog()` fails initially.


# regionReport 1.9.1


SIGNIFICANT USER-VISIBLE CHANGES

* Changed the default style to `BiocStyle::html_document2`.

# regionReport 1.7.10


SIGNIFICANT USER-VISIBLE CHANGES

* Help pages now document advanced arguments.

# regionReport 1.7.2


SIGNIFICANT USER-VISIBLE CHANGES

* Dropped defunct functions.


# regionReport 1.5.48


BUG FIXES

* Fixed a bug in `derfinderReport()` for a case when there are significant
regions but not all regions have finite areas.

# regionReport 1.5.43


SIGNIFICANT USER-VISIBLE CHANGES

* `edgeReport()` now includes two edgeR specific plots: one showing the BCV
and another showing a 2-dim MDS. Also added more `edgeR` citations that I
missed earlier: thank you Gordon Smyth!

# regionReport 1.5.33


NEW FEATURES

* Added the function `edgeReport()` for creating HTML or PDF reports based on
edgeR results. Together with `DESeq2Report()` now `regionReport` supports
the two most used packages for RNA-seq feature-level analysis.


# regionReport 1.5.19


NEW FEATURES

* Added the templates `templatePvalueHistogram` and `templateHistogram`
to be used with `renderReport()` if you prefer histogram plots instead
of density plots.

# regionReport 1.5.12


NEW FEATURES

* Added the function `DESeq2Report()` for creating HTML or PDF reports based
on `DESeq2` results. This should also be useful to explore derfinder
results created from the expressed regions-level approach.

SIGNIFICANT USER-VISIBLE CHANGES

* Added a `digits` argument to control how to round some numerical variables
in all type of reports.
* Added a `theme` argument to allow setting the `ggplot2` theme for the plots.

BUG FIXES

* Improved the PDF versions of all reports by hiding code and shortening
tables. Also added a warning to switch the device to `pdf` for PDF
output since it looks better than the default `png`.

# regionReport 1.5.6


SIGNIFICANT USER-VISIBLE CHANGES

* Switched to using `rmarkdown` instead of `knitrBootstrap` as the default
engine for creating the reports.

# regionReport 1.3.8


SIGNIFICANT USER-VISIBLE CHANGES

* `renderReport()` and `derfinderReport()` now show Manhattan plots for p-value
variables (p-value, q-value, FWER adjusted p-value).


# regionReport 1.3.7


NEW FEATURES

* `renderReport()` now has the `densityTemplates` argument via which users can
customize the density plots for the p-value variables and the continuous
variables. This addresses one of David Robinson's requests at
http://f1000research.com/articles/4-105/v1



# regionReport 1.3.6


NEW FEATURES

* Added a vignette with an example report from bumphunter results.

# regionReport 1.3.5


SIGNIFICANT USER-VISIBLE CHANGES

* Merged pull request https://github.com/leekgroup/regionReport/pull/7
* Added `template` argument to `renderReport()` and `derfinderReport()` to
customize the `knitr` template used
* Wrapped code that works in a temporary directory in `with_wd()` function,
which evaluates in the directory but returns to the original directory
in the case of a user interrupt or error (with `on.exit()`)



# regionReport 1.3.4


NEW FEATURES

* Reports now have a link to the BibTeX file used for the references. This
addresses http://f1000research.com/articles/4-105/v1#reflist Karthik
Ram's bullet point number 4.

# regionReport 1.3.3


NEW FEATURES

* Now uses `derfinderPlot::vennRegions()` to show venn diagram of genomic
states. Requires `derfinderPlot` version 1.3.2 or greater.
* `derfinderReport()` now has a `significantVar` argument that allows users to
choose between determining significant regions by P-values, FDR adjusted
P-values, or FWER adjusted P-values (if FWER adjusted P-values are
absent, then FDR adjusted P-values are used instead, with a warning).


# regionReport 1.3.2


SIGNIFICANT USER-VISIBLE CHANGES

* Deprecated functions with underscores in their names in favor of
camelCase functions. This was done to simplify the package.

# regionReport 1.3.1


BUG FIXES

* Fixed `renderReport()` and `derfinderReport()` so they'll open the correct URL
when `interactive() == TRUE` and the user has knitrBootstrap version 0.9.0
installed instead of the latest GitHub version.


# regionReport 1.1.9


NEW FEATURES

* Introduced `renderReport()` which creates a simple exploratory report for
any set of genomic regions. It allows the user to further customize the
report by using a child file.
* You can now use the `output_format` advanced parameter on both
`renderReport()` and `derfinderReport()` to output a PDF file instead
of an HTML file. The interactive tables are lost and only the top 20
rows are shown.

# regionReport 1.1.8


SIGNIFICANT USER-VISIBLE CHANGES

* Adapted to work with `bumphunter` >= 1.7.6


# regionReport 1.1.7


NEW FEATURES

* Users can now control `output_format` and `clean` options from
`rmarkdown::render()` when running `derfinderReport()`


# regionReport 1.1.3


BUG FIXES

* Adapted `derfinderReport()` to `derfinder` 1.1.5


# regionReport 0.99.0


NEW FEATURES

* Preparing to submit to Bioconductor.

SIGNIFICANT USER-VISIBLE CHANGES

* Updated the vignette and the package to work with recent versions of
the packages this package depends on.
* Renamed the package from `derfinderReport` to `regionReport` and
`generateReport()` to `derfinderReport()`. In the future we will add
another report for a general GRanges object.
* Simplified `derfinderReport()`'s call by using advanced arguments.
* Added Travis integration.


# regionReport 0.0.18


* Now `derfinderReport()` has a `chrsStyle` argument to match changes in
`derfinder` version 0.0.60. `chrsStyle` is set to `UCSC` by default.

# regionReport 0.0.17


SIGNIFICANT USER-VISIBLE CHANGES

* Made more robust for cases where there is a small number of significant
DERs: need at least 3 observations by chr for the chr to be included in
the density plots.

# regionReport 0.0.16


SIGNIFICANT USER-VISIBLE CHANGES

* MA-style plots now use the scaling factor.
* Using a GAM smoother instead of loess for MA-style plots. Helps for cases
with many regions.

# regionReport 0.0.13


SIGNIFICANT USER-VISIBLE CHANGES

* Added a vignette

# regionReport 0.0.12


BUG FIXES

* complying with `BiocCheck` version 1.0.0

# regionReport 0.0.11


SIGNIFICANT USER-VISIBLE CHANGES

* `genomicState` data moved to `derfinder` 0.0.53

# regionReport 0.0.8


SIGNIFICANT USER-VISIBLE CHANGES

* Now requires `knitrBootstrap` 1.0.0
* Matches `derfinder` version 0.0.49

# regionReport 0.0.3


SIGNIFICANT USER-VISIBLE CHANGES

* Matches `derfinder` version 0.0.34

# regionReport 0.0.1


NEW FEATURES

* Migrated from `derfinder`
