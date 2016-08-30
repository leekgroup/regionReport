#' Generate a HTML/PDF report exploring edgeR results
#'
#' This function generates a HTML report with exploratory data analysis plots
#' for edgeR results created. Other output formats are possible such as PDF 
#' reports but they lose the interactivity. Users can easily append
#' to the report by providing a R Markdown file to \code{customCode}, or can
#' customize the entire template by providing an R Markdown file to
#' \code{template}.
#'
#' @param dge A \link[edgeR]{DGEList} object.
#' @param object A \link[edgeR:DGEList-class]{DGEExact} or 
#' \link[edgeR:DGELRT-class]{DGELRT} object that contains p-values stored in
#' \code{object$table$PValue}.
#' @param pAdjustMethod the method to use for adjusting p-values, see 
#' \link[stats]{p.adjust}. This argument will be passed to 
#' \link[DESeq2]{results}.
#' @param alpha the significance cutoff used for optimizing the independent 
#' filtering (by default 0.1). If the adjusted p-value cutoff (FDR) will be a 
#' value other than 0.1, alpha should be set to that value. This argument will
#' be passed to \link[DESeq2]{results}.
#' @param independentFiltering logical, whether independent filtering should be 
#' applied automatically. By default it's set to \code{FALSE} in contrast with
#' the default used in \link[DESeq2]{results} to match \code{edgeR}'s behavior.
#' @param filter the vector of filter statistics over which the independent filtering 
#' will be optimized. By default the logCPM will be used if
#' \code{independentFiltering} is set to \code{TRUE}. It can also be a length
#' 1 character vector specifying one of the column names of \code{object$table}.
#' @param theta the quantiles at which to assess the number of rejections from independent filtering. This argument is passed \link[DESeq2]{results}.
#' @param filterFun an optional custom function as described in 
#' \link[DESeq2]{results}.
#' @inheritParams DESeq2Report
#'
#' @return An HTML report with a basic exploration for the given set of edgeR
#' results.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @importFrom S4Vectors DataFrame metadata 'metadata<-'
#' @importFrom DEFormats as.DESeqDataSet
#' @importFrom DESeq2 DESeqResults
#' @importFrom GenomicRanges mcols 'mcols<-'
#' @importFrom methods is
#'
#' @details
#' Set \code{output_format} to \code{'knitrBootstrap::bootstrap_document'} or 
#' \code{'pdf_document'} if you want a HTML report styled by knitrBootstrap or
#' a PDF report respectively. If using knitrBootstrap, we recommend the version
#' available only via GitHub at https://github.com/jimhester/knitrBootstrap
#' which has nicer features than the current version available via CRAN.
#'
#' If you modify the YAML front matter of \code{template}, you can use other 
#' values for \code{output_format}.
#'
#' This report is similar to the one created by \link{DESeq2Report} with two
#' additional plots exclusive for edgeR results. We designed the reports to be
#' very similar intentionally and use the Bioconductor package DEFormats to
#' achieve this goal.
#'
#' @examples
#'
#' ## Create example data using DEFormats
#' library('DEFormats')
#' set.seed(20160407)
#' counts <- simulateRnaSeqData()
#' group <- rep(c("A", "B"), each = 3)
#' 
#' ## Create DGEList object
#' library('edgeR')
#' dge <- DGEList(counts, group = group)
#'
#' ## Perform DE analysis with edgeR
#' design <- model.matrix( ~ group)
#' dge <- estimateDisp(dge, design)
#' fit <- glmFit(dge, design)
#' lrt <- glmLRT(fit, coef = 2)
#'
#' ## The output will be saved in the 'edgeReport-example' directory
#' dir.create('edgeReport-example', showWarnings = FALSE, recursive = TRUE)
#'
#' ## Generate the HTML report
#' report <- edgeReport(dge, lrt, project = 'edgeR-example', intgroup = 'group',
#'     outdir = 'edgeReport-example')
#'
#' if(interactive()) {
#'     ## Browse the report
#'     browseURL(report)
#' }
#'
#' \dontrun{
#' ## Note that you can run the example using:
#' example('edgeReport', 'regionReport', ask=FALSE)
#' }
#'
#'


edgeReport <- function(dge, object, project = "", intgroup, colors = NULL,
    pAdjustMethod = 'BH', alpha = 0.1, independentFiltering = FALSE, filter,
    theta, filterFun, nBest = 500, nBestFeatures = 20, customCode = NULL,
    outdir = 'edgeRexploration', output = 'edgeRexploration',
    browse = interactive(), device = 'png', template = NULL, 
    searchURL = 'http://www.ncbi.nlm.nih.gov/gene/?term=', theme = NULL,
    digits = 2, ...) {
    
    ## Check inputs
    stopifnot(is(dge, 'DGEList'))
    stopifnot(is(object, 'DGEExact') | is(object, 'DGELRT'))
    if (is.null(object$table)) stop("Need to run exactTest or glmLRT first")
    stopifnot(nrow(dge$counts) == nrow(object$table))
    stopifnot(rownames(dge$counts) == rownames(object$table))
    stopifnot(all(c('PValue', 'logCPM', 'logFC') %in% colnames(object$table)))
    
    ## Transform to DESeq2 format
    dds <- as.DESeqDataSet(dge)
    
    ## Fake having created results
    mcols(dds) <- DataFrame(pvalue = object$table$PValue)
    mcols(mcols(dds)) <- DataFrame(type = 'results',
        description = 'manual incomplete conversion from edgeR to DESeq2')
    
    ## Fix column names
    colnames(object$table)[colnames(object$table) == 'PValue'] <- 'pvalue'
    colnames(object$table)[colnames(object$table) == 'logFC'] <- 'log2FoldChange'
    colnames(object$table)[colnames(object$table) == 'logCPM'] <- 'baseMean'
    
    ## Define res input
    deseqRes <- DESeqResults(DataFrame(object$table))
    
    ## Use as default the logCPM (since it's renamed to baseMean)
    if(!missing(filter)) {
        if(is.character(filter) & length(filter) == 1) {
            stopifnot(filter %in% colnames(object$table))
            filter <- object$table[, filter]
        }
    }
    
    ## Obtain results
    if(missing(filterFun)) {
        deseqRes <- DESeq2:::pvalueAdjustment(deseqRes, 
            independentFiltering = independentFiltering, filter = filter,
            theta = theta, alpha = alpha, pAdjustMethod = pAdjustMethod)
    } else {
        deseqRes <- filterFun(deseqRes, filter = filter, alpha = alpha,
            pAdjustMethod = pAdjustMethod)
    }
    ## Save alpha value
    metadata(deseqRes)[["alpha"]] <- alpha
    
    ## Call used
    theCall <- match.call()
    
    ## Create report
    DESeq2Report(dds, project = project, intgroup = intgroup, colors = colors,
        res = deseqRes, nBest = nBest, nBestFeatures = nBestFeatures,
        customCode = customCode, outdir = outdir, output = output,
        browse = browse, device = device, template = template,
        searchURL = searchURL, theme = theme, digits = digits,
        software = 'edgeR', theCall = theCall, dge = dge, ...)
    
}
