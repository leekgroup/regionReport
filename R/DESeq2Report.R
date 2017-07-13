#' Generate a HTML/PDF report exploring DESeq2 results
#'
#' This function generates a HTML report with exploratory data analysis plots
#' for DESeq2 results created with \link[DESeq2]{DESeq}. Other output formats
#' are possible such as PDF but lose the interactivity. Users can easily append
#' to the report by providing a R Markdown file to \code{customCode}, or can
#' customize the entire template by providing an R Markdown file to
#' \code{template}.
#'
#' @param dds A \link[DESeq2]{DESeqDataSet} object with the results from
#' running \link[DESeq2]{DESeq}.
#' @param project The title of the project.
#' @param intgroup interesting groups: a character vector of names in 
#' \code{colData(x)} to use for grouping. This parameter is passed to functions
#' such as \link[DESeq2]{plotPCA}.
#' @param colors vector of colors used in heatmap. If \code{NULL}, then a
#' a default set of colors will be used. This argument is passed to
#' \link[pheatmap]{pheatmap}.
#' @param res A \link[DESeq2]{DESeqResults} object. If \code{NULL}, then
#' \link[DESeq2]{results} will be used on \code{dds} with default parameters.
#' @param nBest The number of features to include in the interactive 
#' table. Features are ordered by their adjusted p-values.
#' @param nBestFeatures The number of best features to make plots of their 
#' counts. We recommend a small number, say 20.
#' @param customCode An absolute path to a child R Markdown file with code to be
#' evaluated before the reproducibility section. Its useful for users who want
#' to customize the report by adding conclusions derived from the data and/or
#' further quality checks and plots.
#' @param outdir The name of output directory.
#' @param output The name of output HTML file (without the html extension).
#' @param browse If \code{TRUE} the HTML report is opened in your browser once 
#' it's completed.
#' @param device The graphical device used when knitting. See more at 
#' http://yihui.name/knitr/options (\code{dev} argument).
#' @param template Template file to use for the report. If not provided, will
#' use the default file found in DESeq2Exploration/DESeq2Exploration.Rmd
#' within the package source.
#' @param searchURL A url used for searching the name of the features in
#' the web. By default \code{http://www.ncbi.nlm.nih.gov/gene/?term=} is used
#' which is the recommended option when features are genes. It's only used
#' when the output is a HTML file.
#' @param theme A ggplot2 \link[ggplot2]{theme} to use for the plots made with
#' ggplot2.
#' @param digits The number of digits to round to in the interactive table of
#' the top \code{nBestFeatures}. Note that p-values and adjusted p-values won't 
#' be rounded.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{software }{ The name of the package used for performing the 
#' differential expression analysis. Either \code{DESeq2} or \code{edgeR}.}
#' \item{dge }{ A \link[edgeR]{DGEList} object. \code{NULL} by default and only 
#' used by \link{edgeReport}.}
#' \item{theCall }{ The function call. \code{NULL} by default and only used by
#' \link{edgeReport}.}
#' \item{output_format }{ Either \code{html_document}, \code{pdf_document} or
#' \code{knitrBootstrap::bootstrap_document} unless you modify the YAML
#' template.}
#' \item{clean }{ Logical, whether to clean the results or not. Passed to
#' \link[rmarkdown]{render}.}
#' }
#'
#' @return An HTML report with a basic exploration for the given set of DESeq2
#' results.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom DESeq2 results
#' @importFrom knitcitations cleanbib cite_options write.bibtex read.bibtex
#' citep bibliography
#' @importFrom utils browseURL citation packageVersion
#' @importFrom rmarkdown render html_document
#' @importFrom GenomicRanges mcols
#' @importFrom knitrBootstrap knit_bootstrap
#' @importFrom RefManageR BibEntry
#' @importFrom BiocStyle html_document2
#' @import knitr
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
#' ## Load example data from the pasilla package
#' library('pasilla')
#' library('DESeq')
#' library('DESeq2')
#'
#' ## Create DESeqDataSet object from the pasilla package
#' data('pasillaGenes')
#' countData <- counts(pasillaGenes)
#' colData <- pData(pasillaGenes)[, c('condition', 'type')]
#' dds <- DESeqDataSetFromMatrix(countData = countData,
#'    colData = colData,
#'    design = ~ condition)
#' dds <- DESeq(dds)
#'
#' ## The output will be saved in the 'DESeq2Report-example' directory
#' dir.create('DESeq2Report-example', showWarnings = FALSE, recursive = TRUE)
#'
#' ## Generate the HTML report
#' report <- DESeq2Report(dds, 'DESeq2-example', c('condition', 'type'),
#'     outdir = 'DESeq2Report-example')
#'
#' if(interactive()) {
#'     ## Browse the report
#'     browseURL(report)
#' }
#'
#' \dontrun{
#' ## Note that you can run the example using:
#' example('DESeq2Report', 'regionReport', ask=FALSE)
#' }
#'
#'


DESeq2Report <- function(dds, project = "", intgroup, colors = NULL, res = NULL,
    nBest = 500, nBestFeatures = 20, customCode = NULL,
    outdir = 'DESeq2Exploration', output = 'DESeq2Exploration',
    browse = interactive(), device = 'png', template = NULL, 
    searchURL = 'http://www.ncbi.nlm.nih.gov/gene/?term=', theme = NULL,
    digits = 2, ...) {
    ## Save start time for getting the total processing time
    startTime <- Sys.time()
    
    
    ## Check inputs
    stopifnot(is(dds, 'DESeqDataSet'))
    if (!"results" %in% mcols(mcols(dds))$type)
        stop("couldn't find results. you should first run DESeq()")
    if (!all(intgroup %in% names(colData(dds)))) 
        stop("all variables in 'intgroup' must be columns of colData")
    if(is.null(res)) {
        ## Run results with default parameters
        res <- results(dds)
    } else {
        stopifnot(is(res, 'DESeqResults'))
        stopifnot(identical(nrow(res), nrow(dds)))
    }
    stopifnot(is.null(searchURL) | length(searchURL) == 1)
    if(!is.null(theme)) stopifnot(is(theme, c('theme', 'gg')))
        
# @param software The name of the package used for performing the differential
# expression analysis. Either \code{DESeq2}, \code{edgeR} or \code{other}.
    software <- .advanced_argument('software', 'DESeq2', ...)
    stopifnot(software %in% c('DESeq2', 'edgeR', 'other'))
    isEdgeR <- software == 'edgeR'
    
# @param dge A \link[edgeR]{DGEList} object.
    dge <- .advanced_argument('dge', NULL, ...)
    if(isEdgeR) stopifnot(is(dge, 'DGEList'))
    
    ## Is there custom code?
    hasCustomCode <- !is.null(customCode)
    if(hasCustomCode) stopifnot(length(customCode) == 1)

    
    ## Create outdir
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    workingDir <- getwd()
    
    ## Locate Rmd if one is not provided
    if (is.null(template)) {
        templateNull <- TRUE
        template <- system.file(
            file.path('DESeq2Exploration', 'DESeq2Exploration.Rmd'),
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
    pkgs <- c('DT', 'ggplot2', 'pheatmap', 'RColorBrewer')
    if(isEdgeR) pkgs <- c(pkgs, 'edgeR')
    for(pkg in pkgs) load_install(pkg)
    
    ## Write bibliography information
    write.bibtex(c(
        knitcitations = citation('knitcitations'), 
        regionReport = citation('regionReport')[1],
        DT = citation('DT'), 
        ggplot2 = citation('ggplot2'),
        knitr = citation('knitr')[3],
        rmarkdown = citation('rmarkdown'),
        pheatmap = citation('pheatmap'),
        RColorBrewer = citation('RColorBrewer'),
        DESeq2 = citation('DESeq2'),
        if(isEdgeR) edgeR1 = citation('edgeR')[1] else NULL,
        if(isEdgeR) edgeR2 = citation('edgeR')[2] else NULL,
        if(isEdgeR) edgeR5 = citation('edgeR')[5] else NULL,
        if(isEdgeR) edgeR6 = RefManageR::BibEntry('inbook', key = 'edgeR6', author = 'Chen, Yunshun and Lun, Aaron T. L. and Smyth, Gordon K.', title = 'Differential expression analysis of complex RNA-seq experiments using edgeR', booktitle = 'Statistical Analysis of Next Generation Sequencing Data', year = 2014, editor = 'Datta, Somnath and Nettleton, Dan', publisher = 'Springer', location = 'New York', pages = '51-74') else NULL
        ), file = file.path(outdir, paste0(output, '.bib'))
    )
    bib <- read.bibtex(file.path(outdir, paste0(output, '.bib')))
    
    ## Assign short names
    bib.names <- c('knitcitations', 'regionReport', 'DT', 'ggplot2', 'knitr',
        'rmarkdown', 'pheatmap', 'RColorBrewer', 'DESeq2')
    if(isEdgeR) bib.names <- c(bib.names, 'edgeR1', 'edgeR2', 'edgeR5', 'edgeR6')
    names(bib) <- bib.names 
    
    ## Save the call
    theCall <- .advanced_argument('theCall', NULL, ...)
    if(!is(theCall, 'call')) theCall <- match.call()
    
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
