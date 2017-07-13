#' Generate a HTML/PDF report exploring a set of genomic regions
#'
#' This function generates a HTML report with quality checks, genome location
#' exploration, and an interactive table with the results. Other output formats
#' are possible such as PDF but lose the interactivity. Users can easily append
#' to the report by providing a R Markdown file to \code{customCode}, or can
#' customize the entire template by providing an R Markdown file to
#' \code{template}.
#'
#' @param regions The set of genomic regions of interest as a \code{GRanges}
#' object. All sequence lengths must be provided.
#' @param project The title of the project.
#' @param pvalueVars The names of the variables with values between 0 and 1 to
#' plot density values by chromosome and a table for commonly used cutoffs.
#' Most commonly used to explore p-value distributions. If a named
#' character vector is provided, the names are used in the plot titles.
#' @param densityVars The names of variables to use for making density plots
#' by chromosome. Commonly used to explore scores and other variables given
#' by region.  If a named character vector is provided, the names are used in 
#' the plot titles.
#' @param significantVar A \code{logical} variable differentiating statistically
#' significant regions from the rest. When provided, both types of regions
#' are compared against each other to see differences in width, location, etc.
#' @param annotation The output from \link[bumphunter]{matchGenes} used on
#' \code{regions}. Note that this can take time for a large set of regions
#' so it's better to pre-compute this information and save it.
#' @param nBestRegions The number of regions to include in the interactive 
#' table.
#' @param customCode An absolute path to a child R Markdown file with code to be
#' evaluated before the reproducibility section. Its useful for users who want
#' to customize the report by adding conclusions derived from the data and/or
#' further quality checks and plots.
#' @param outdir The name of output directory.
#' @param output The name of output HTML file (without the html extension).
#' @param browse If \code{TRUE} the HTML report is opened in your browser once 
#' it's completed.
#' @param txdb Specify the transcription database to use for identifying the
#' closest genes via \link[bumphunter]{matchGenes}. If \code{NULL} it will
#' use TxDb.Hsapiens.UCSC.hg19.knownGene by default.
#' @param device The graphical device used when knitting. See more at 
#' http://yihui.name/knitr/options (\code{dev} argument).
#' @param densityTemplates A list of length 3 with templates for the p-value 
#' density plots (variables from \code{pvalueVars}), the continuous 
#' variables density plots (variables from \code{densityVars}), and Manhattan 
#' plots for the p-value variables (\code{pvalueVars}). These templates 
#' are processed by \link[whisker]{whisker.render}. Check the default templates 
#' for more information. The \code{densityTemplates} argument is available for 
#' those users interested in customizing these plots. For example, to show 
#' histograms instead of density plots use \code{templatePvalueHistogram} and
#' \code{templateHistogram} instead of \code{templatePvalueDensity} and
#' \code{templateDensity} respectively.
#' @param template Template file to use for the report. If not provided, will
#' use the default file found in regionExploration/regionExploration.Rmd
#' within the package source.
#' @param theme A ggplot2 \link[ggplot2]{theme} to use for the plots made with
#' ggplot2.
#' @param digits The number of digits to round to in the interactive table of
#' the top \code{nBestRegions}. Note that p-values and adjusted p-values won't 
#' be rounded.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{overviewParams }{ A two element list with \code{base_size} and 
#' \code{areaRel} that control the text size for the genomic overview plots.}
#' \item{output_format }{ Either \code{html_document}, \code{pdf_document} or
#' \code{knitrBootstrap::bootstrap_document} unless you modify the YAML
#' template.}
#' \item{clean }{ Logical, whether to clean the results or not. Passed to
#' \link[rmarkdown]{render}.}
#' }
#'
#' @return An HTML report with a basic exploration for the given set of
#' genomic regions.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @importFrom knitcitations cleanbib cite_options write.bibtex read.bibtex
#' citep bibliography
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom utils browseURL citation packageVersion
#' @importFrom rmarkdown render html_document
#' @importFrom GenomicRanges mcols 'mcols<-'
#' @importFrom knitrBootstrap knit_bootstrap
#' @importFrom BiocStyle html_document2
#' @importFrom methods is
#'
#' @details
#' Set \code{output_format} to \code{'knitrBootstrap::bootstrap_document'} or 
#' \code{'pdf_document'} if you want a HTML report styled by knitrBootstrap or
#' a PDF report respectively. If using knitrBootstrap, we recommend the version
#' available only via GitHub at https://github.com/jimhester/knitrBootstrap
#' which has nicer features than the current version available via CRAN. You can
#' also set the \code{output_format} to \code{'html_document'} for a HTML
#' report styled by rmarkdown. The default is set to 
#' \code{'BiocStyle::html_document2'}.
#'
#' If you modify the YAML front matter of \code{template}, you can use other 
#' values for \code{output_format}.
#'
#' The HTML report styled with knitrBootstrap can be smaller in size than the
#' \code{'html_document'} report.
#'
#' @examples
#'
#' ## Load derfinder for an example set of regions
#' library('derfinder')
#' regions <- genomeRegions$regions
#'
#' ## Assign chr length
#' library('GenomicRanges')
#' seqlengths(regions) <- c('chr21' = 48129895)
#'
#' ## The output will be saved in the 'renderReport-example' directory
#' dir.create('renderReport-example', showWarnings = FALSE, recursive = TRUE)
#'
#' ## Generate the HTML report
#' report <- renderReport(regions, 'Example run', pvalueVars = c(
#'     'Q-values' = 'qvalues', 'P-values' = 'pvalues'), densityVars = c(
#'     'Area' = 'area', 'Mean coverage' = 'meanCoverage'), 
#'     significantVar = regions$qvalues <= 0.05, nBestRegions = 20,
#'     outdir = 'renderReport-example')
#'
#' if(interactive()) {
#'     ## Browse the report
#'     browseURL(report)
#' }
#'
#' \dontrun{
#' ## Note that you can run the example using:
#' example('renderReport', 'regionReport', ask=FALSE)
#' }
#'
#' ## Check the default templates. For users interested in customizing these 
#' ## plots.
#' ## For p-value variables:
#' cat(templatePvalueDensity)
#'
#' ## For continous variables:
#' cat(templateDensity)
#'
#' ## For Manhattan plots
#' cat(templateManhattan)
#'


renderReport <- function(regions, project = "", 
    pvalueVars = c('P-values' = 'pval'),
    densityVars = NULL, significantVar = mcols(regions)$pval <= 0.05,
    annotation = NULL, nBestRegions = 500, customCode = NULL,
    outdir = 'regionExploration', output = 'regionExploration',
    browse = interactive(), txdb = NULL, device = 'png',
    densityTemplates = list(Pvalue = templatePvalueDensity,
        Common = templateDensity, Manhattan = templateManhattan),
    template = NULL, theme = NULL, digits = 2, ...) {
    ## Save start time for getting the total processing time
    startTime <- Sys.time()
    
    
    ## Check inputs
    stopifnot(is(regions, 'GRanges'))
    stopifnot(!any(is.na(seqlengths(regions))))
    stopifnot(is.list(densityTemplates) & length(densityTemplates) == 3 & all(c('Pvalue', 'Common', 'Manhattan') %in% names(densityTemplates)))
    hasCustomCode <- !is.null(customCode)
    if(hasCustomCode) stopifnot(length(customCode) == 1)
    if(!is.null(annotation)) stopifnot(nrow(annotation) == length(regions))
    if(!is.null(theme)) stopifnot(is(theme, c('theme', 'gg')))
    
    # @param overviewParams A two element list with \code{base_size} and 
    # \code{areaRel} that control the text size for the genomic overview plots.
    overviewParams <- .advanced_argument('overviewParam', list(base_size = 10,
            areaRel = 5), ...)
    
    
    ## Check that overviewParams is correctly specified
    stopifnot(sum(names(overviewParams) %in% c('base_size', 'areaRel')) == 2)
    
    ## Are there p-value vars?
    hasPvalueVars <- length(pvalueVars) > 0
    if(hasPvalueVars) stopifnot(is.character(pvalueVars))
    ## Fix p-value variable names
    colnames(mcols(regions))[which(colnames(mcols(regions)) %in% pvalueVars)] <- make.names(pvalueVars)
    pvalueVars[seq_len(length(pvalueVars))] <- make.names(pvalueVars)
        
    ## Are there density vars?
    hasDensityVars <- length(densityVars) > 0
    if(hasDensityVars) stopifnot(is.character(densityVars))
    ## Fix density variable names
    colnames(mcols(regions))[which(colnames(mcols(regions)) %in% densityVars)] <- make.names(densityVars)
    densityVars[seq_len(length(densityVars))] <- make.names(densityVars)
        
    ## Are there significant regions?
    hasSignificant <- length(significantVar) > 0
    if(hasSignificant) {
        stopifnot(length(significantVar) == length(regions))
        stopifnot(is.logical(significantVar))
        hasSignificant <- any(significantVar)
        if(!hasSignificant) warning("There are no statistically significant regions in this set.")
    }   
    
    ## Template for pvalueVars
    if(hasPvalueVars) {
        templatePvalueDensityInUse <- densityTemplates$Pvalue
        templateManhattanInUse <- densityTemplates$Manhattan
    }
    
    ## Template for densityVars
    if(hasDensityVars) templateDensityInUse <- densityTemplates$Common
    
    
    ## Create outdir
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    workingDir <- getwd()
    
    ## Locate Rmd if one is not provided
    if (is.null(template)) {
        templateNull <- TRUE
        template <- system.file(
            file.path('regionExploration', 'regionExploration.Rmd'),
            package = 'regionReport', mustWork = TRUE
        )
    } else {
        templateNull <- FALSE
    }
    
    ## Load knitcitations with a clean bibliography
    cleanbib()
    cite_options(hyperlink = 'to.doc', citation_format = 'text', style = 'html')
    # Note links won't show for now due to the following issue
    # https://github.com/cboettig/knitcitations/issues/63
    
    ## Install suggested packages that are needed for citation to work
    for(pkg in c('bumphunter', 'derfinderPlot', 'DT', 'ggbio', 'ggplot2', 'whisker')) load_install(pkg)
    
    ## Write bibliography information
    bib <- c(
        knitcitations = citation('knitcitations'), 
        regionReport = citation('regionReport')[1],
        derfinderPlot = citation('derfinderPlot')[1],
        DT = citation('DT'), 
        ggbio = citation('ggbio'),
        ggplot2 = citation('ggplot2'),
        knitr = citation('knitr')[3],
        rmarkdown = citation('rmarkdown'),
        whisker = citation('whisker'),
        bumphunter = citation('bumphunter')[1],
        derfinder = citation('derfinder')[1])
    write.bibtex(bib, file = file.path(outdir, paste0(output, '.bib')))
    
    ## Save the call
    theCall <- match.call()
    
    ## Generate report
    ## Perform code within the output directory.
    tmpdir <- getwd()
    with_wd(outdir, {
        file.copy(template, to = paste0(output, '.Rmd'))
    
        ## Output format
        output_format <- .advanced_argument('output_format',
            'BiocStyle::html_document2', ...)
        outputIsHTML <- output_format %in% c('html_document',
            'knitrBootstrap::bootstrap_document', 'BiocStyle::html_document2')
        if(!outputIsHTML) {
            if(device == 'png') warning("You might want to switch the 'device' argument from 'png' to 'pdf' for better quality plots.")
        }
    
        ## Check knitrBoostrap version
        knitrBootstrapFlag <- packageVersion('knitrBootstrap') < '1.0.0'
            if(knitrBootstrapFlag & output_format == 'knitrBootstrap::bootstrap_document') {
            ## CRAN version
            tmp <- knit_bootstrap(paste0(output, '.Rmd'), chooser = c('boot',
                'code'), show_code = TRUE)
            res <- file.path(tmpdir, outdir, paste0(output, '.html'))
            unlink(paste0(output, '.md'))
        } else {
            res <- render(paste0(output, '.Rmd'), output_format,
                clean = .advanced_argument('clean', TRUE, ...))
        }
        if(templateNull) file.remove(paste0(output, '.Rmd'))
    
        ## Open
        if (browse) browseURL(res)
    })
    
    ## Finish
    return(invisible(res))
}


#' @rdname renderReport
#' @export
templatePvalueDensity <- "
## {{{densityVarName}}}

```{r pval-density-{{{varName}}}, fig.width=10, fig.height=10}
p1{{{varName}}} <- ggplot(regions.df.plot, aes(x={{{varName}}}, colour=seqnames)) +
    geom_line(stat='density') + xlim(0, 1) +
    labs(title='Density of {{{densityVarName}}}') + xlab('{{{densityVarName}}}') +
    scale_colour_discrete(limits=chrs) + theme(legend.title=element_blank())
p1{{{varName}}}
```


```{r 'pval-summary-{{{varName}}}'}
summary(mcols(regions)[['{{{varName}}}']])
```


This is the numerical summary of the distribution of the {{{densityVarName}}}.

```{r pval-tableSummary-{{{varName}}}, results='asis'}
{{{varName}}}table <- lapply(c(1e-04, 0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 1), function(x) {
    data.frame('Cut' = x, 'Count' = sum(mcols(regions)[['{{{varName}}}']] <= x))
})
{{{varName}}}table <- do.call(rbind, {{{varName}}}table)
if(outputIsHTML) {
    kable({{{varName}}}table, format = 'markdown', align = c('c', 'c'))
} else {
    kable({{{varName}}}table)
}
```

This table shows the number of regions with {{{densityVarName}}} less or equal than some commonly used cutoff values. 

"

#' @rdname renderReport
#' @export
templateDensity <- "
## {{{densityVarName}}}

```{r density-{{{varName}}}, fig.width=14, fig.height=14, eval=hasSignificant, echo=hasSignificant}
xrange <- range(regions.df.plot[, '{{{varName}}}']) * c(0.95, 1.05)
p3a{{{varName}}} <- ggplot(regions.df.plot[is.finite(regions.df.plot[, '{{{varName}}}']), ], aes(x={{{varName}}}, colour=seqnames)) +
    geom_line(stat='density') + labs(title='Density of {{{densityVarName}}}') +
    xlab('{{{densityVarName}}}') + scale_colour_discrete(limits=chrs) +
    xlim(xrange) + theme(legend.title=element_blank())
p3b{{{varName}}} <- ggplot(regions.df.sig[is.finite(regions.df.sig[, '{{{varName}}}']), ], aes(x={{{varName}}}, colour=seqnames)) +
    geom_line(stat='density') +
    labs(title='Density of {{{densityVarName}}} (significant only)') +
    xlab('{{{densityVarName}}}') + scale_colour_discrete(limits=chrs) +
    xlim(xrange) + theme(legend.title=element_blank())
grid.arrange(p3a{{{varName}}}, p3b{{{varName}}})
```

```{r density-solo-{{{varName}}}, fig.width=10, fig.height=10, eval=!hasSignificant, echo=!hasSignificant}
p3a{{{varName}}} <- ggplot(regions.df.plot[is.finite(regions.df.plot[, '{{{varName}}}']), ], aes(x={{{varName}}}, colour=seqnames)) +
    geom_line(stat='density') + labs(title='Density of {{{densityVarName}}}') +
    xlab('{{{densityVarName}}}') + scale_colour_discrete(limits=chrs) +
    theme(legend.title=element_blank())
p3a{{{varName}}}
```

This plot shows the density of the {{{densityVarName}}} for all regions. `r ifelse(hasSignificant, 'The bottom panel is restricted to significant regions.', '')`

"

#' @rdname renderReport
#' @export
templateManhattan <- "
## Manhattan {{{densityVarName}}}

```{r manhattan-{{{varName}}}, fig.width=10, fig.height=10, message = FALSE}

regions.manhattan <- regions
mcols(regions.manhattan)[['{{{varName}}}']] <- - log(mcols(regions.manhattan)[['{{{varName}}}']], base = 10)
pMan{{{varName}}} <- plotGrandLinear(regions.manhattan, aes(y = {{{varName}}}, colour = seqnames)) + theme(axis.text.x=element_text(angle=-90, hjust=0)) + ylab('-log10 {{{densityVarName}}}')
pMan{{{varName}}}
rm(regions.manhattan)
```

This is a Manhattan plot for the {{{densityVarName}}} for all regions. A single dot is shown for each region, where higher values in the y-axis mean that the {{{densityVarName}}} are closer to zero.

"

#' @rdname renderReport
#' @export
templatePvalueHistogram <- "
## {{{densityVarName}}}

```{r histPval-{{{varName}}}, fig.width=10, fig.height=10}
p1{{{varName}}} <- ggplot(regions.df.plot, aes(x={{{varName}}}, colour=seqnames)) +
    geom_histogram(bins = 50, alpha=.5, position='identity') +
    xlim(c(0, 1.0005)) +
    labs(title='Histogram of {{{densityVarName}}}') + 
    xlab('{{{densityVarName}}}') +
    scale_colour_discrete(limits=chrs) + theme(legend.title=element_blank())
p1{{{varName}}}
```

This plot shows the distribution of {{{densityVarName}}} in a histogram. It might be skewed right or left, or flat as shown in the [Wikipedia examples](https://en.wikipedia.org/wiki/Histogram#Examples). The shape depends on the percent of regions that are differentially expressed. For further information on how to interpret a histogram of p-values check [David Robinson's post on this topic](http://varianceexplained.org/statistics/interpreting-pvalue-histogram/).


```{r 'summary-{{{varName}}}'}
summary(mcols(regions)[['{{{varName}}}']])
```


This is the numerical summary of the distribution of the {{{densityVarName}}}.

```{r tableSummary-{{{varName}}}, results='asis'}
{{{varName}}}table <- lapply(c(1e-04, 0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 1), function(x) {
    data.frame('Cut' = x, 'Count' = sum(mcols(regions)[['{{{varName}}}']] <= x))
})
{{{varName}}}table <- do.call(rbind, {{{varName}}}table)
if(outputIsHTML) {
    kable({{{varName}}}table, format = 'markdown', align = c('c', 'c'))
} else {
    kable({{{varName}}}table)
}
```

This table shows the number of regions with {{{densityVarName}}} less or equal than some commonly used cutoff values. 

"

#' @rdname renderReport
#' @export
templateHistogram <- "

## {{{densityVarName}}}

```{r histogram-{{{varName}}}, fig.width=14, fig.height=14, eval=hasSignificant, echo=hasSignificant}
xrange <- range(regions.df.plot[, '{{{varName}}}'])
p3a{{{varName}}} <- ggplot(regions.df.plot[is.finite(regions.df.plot[, '{{{varName}}}']), ], aes(x={{{varName}}}, fill=seqnames)) +
    geom_histogram(alpha=.5, position='identity', bins = 50) +
    labs(title='Histogram of {{{densityVarName}}}') +
    xlab('{{{densityVarName}}}') +
    xlim(xrange) + theme(legend.title=element_blank())
p3b{{{varName}}} <- ggplot(regions.df.sig[is.finite(regions.df.sig[, '{{{varName}}}']), ], aes(x={{{varName}}}, fill=seqnames)) +
    geom_histogram(alpha=.5, position='identity', bins = 50) +
    labs(title='Histogram of {{{densityVarName}}} (significant only)') +
    xlab('{{{densityVarName}}}') +
    xlim(xrange) + theme(legend.title=element_blank())
grid.arrange(p3a{{{varName}}}, p3b{{{varName}}})
```

```{r histogram-solo-{{{varName}}}, fig.width=10, fig.height=10, eval=!hasSignificant, echo=!hasSignificant}
p3a{{{varName}}} <- ggplot(regions.df.plot[is.finite(regions.df.plot[, '{{{varName}}}']), ], aes(x={{{varName}}}, fill=seqnames)) +
    geom_histogram(alpha=.5, position='identity', bins = 50) +
    labs(title='Histogram of {{{densityVarName}}}') +
    xlab('{{{densityVarName}}}') +
    theme(legend.title=element_blank())
p3a{{{varName}}}
```

This plot shows the distribution of {{{densityVarName}}} in a histogram for all regions. `r ifelse(hasSignificant, 'The bottom panel is restricted to significant regions.', '')`

"
