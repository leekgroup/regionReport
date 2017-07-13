#' Generate a HTML/PDF report exploring the basic results from derfinder
#'
#' This function generates a HTML report exploring the basic results from 
#' single base-level approach derfinder analysis results 
#' (www.bioconductor.org/packages/derfinder). The HTML report itself 
#' is generated using rmarkdown (http://rmarkdown.rstudio.com/). It works best 
#' after using \link[derfinder]{mergeResults}.
#' 
#' @param prefix The main data directory path where 
#' \link[derfinder]{mergeResults} was run. It should be the same as 
#' \code{mergeResults(prefix)}.
#' @param outdir The name of output directory relative to \code{prefix}.
#' @param output The name of output HTML file (without the html extension).
#' @param project The title of the project.
#' @param browse If \code{TRUE} the HTML report is opened in your browser once 
#' it's completed.
#' @param nBestRegions The number of region plots to make, ordered by area.
#' @param makeBestClusters If \code{TRUE}, \link[derfinderPlot]{plotCluster} is 
#' used on the \code{nBestClusters} regions by area. Note that these plots take 
#' some time to make.
#' @param nBestClusters The number of region cluster plots to make by taking 
#' the \code{nBestClusters} regions ranked by area of the cluster.
#' @param fullCov A list where each element is the result from 
#' \link[derfinder]{loadCoverage} used with \code{cutoff=NULL}. Can be 
#' generated using \link[derfinder]{fullCoverage}.
#' @param hg19 If \code{TRUE} then the reference is assumed to be hg19 and 
#' chromosome lengths as well as the default transcription database 
#' (TxDb.Hsapiens.UCSC.hg19.knownGene) will be used.
#' @param p.ideos A list where each element is the result of 
#' \link[ggbio]{plotIdeogram}. If it's \code{NULL} and \code{hg19=TRUE} then 
#' they are created for the hg19 human reference.
#' @param txdb Specify the transcription database to use for making the plots 
#' for the top regions by area. If \code{NULL} and \code{hg19=TRUE} then 
#' TxDb.Hsapiens.UCSC.hg19.knownGene is used.
#' @param device The graphical device used when knitting. See more at 
#' http://yihui.name/knitr/options (\code{dev} argument).
#' @param significantVar A character variable specifying whether to use the
#' p-values, the FDR adjusted p-values or the FWER adjusted p-values to
#' determine significance. Has to be either \code{'pvalue'}, \code{'qvalue'}
#' or \code{'fwer'}.
#' @param customCode An absolute path to a child R Markdown file with code to be
#' evaluated before the reproducibility section. Its useful for users who want
#' to customize the report by adding conclusions derived from the data and/or
#' further quality checks and plots.
#' @param template Template file to use for the report. If not provided, will
#' use the default file found in basicExploration/basicExploration.Rmd
#' within the package source.
#' @param theme A ggplot2 \link[ggplot2]{theme} to use for the plots made with
#' ggplot2.
#' @param digits The number of digits to round to in the interactive table of
#' the top \code{nBestRegions}. Note that p-values and adjusted p-values won't 
#' be rounded.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{chrsStyle }{ The naming style of the chromosomes. By default, UCSC. 
#' See \link[GenomeInfoDb]{seqlevelsStyle}.}
#' \item{species }{ Species name. See \link[derfinder]{extendedMapSeqlevels}  
#' for more information.}
#' \item{currentStyle }{ Current naming style used. See 
#' \link[derfinder]{extendedMapSeqlevels} for more information.}
#' \item{fullRegions }{ Part of the output of \link[derfinder]{mergeResults}. 
#' Specify it only if you have already loaded it in memory.}
#' \item{fullNullSummary }{ Part of the output of 
#' \link[derfinder]{mergeResults}. Specify it only if you have already loaded 
#' it in memory.}
#' \item{fullAnnotatedRegions }{ Part of the output of 
#' \link[derfinder]{mergeResults}. Specify it only if you have already loaded 
#' it in memory.}
#' \item{optionsStats }{ Part of the output of \link[derfinder]{analyzeChr}. 
#' Specify it only if you have already loaded it in memory.}
#' \item{optionsMerge }{ Part of the output of \link[derfinder]{mergeResults}. 
#' Specify it only if you have already loaded it in memory.}
#' \item{overviewParams }{ A two element list with \code{base_size} and 
#' \code{areaRel} that control the text size for the genomic overview plots.}
#' \item{output_format }{ Either \code{html_document}, \code{pdf_document} or
#' \code{knitrBootstrap::bootstrap_document} unless you modify the YAML
#' template.}
#' \item{clean }{ Logical, whether to clean the results or not. Passed to
#' \link[rmarkdown]{render}.}
#' }
#' Passed to \link[derfinder]{extendedMapSeqlevels}.
#'
#' @return An HTML report with a basic exploration of the derfinder results.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link[derfinder]{mergeResults}, \link[derfinder]{analyzeChr}, 
#' \link[derfinder]{fullCoverage}
#' @export
#'
#' @importFrom derfinder extendedMapSeqlevels
#' @importFrom knitcitations cleanbib cite_options write.bibtex read.bibtex
#' citep bibliography
#' @importFrom GenomeInfoDb seqlevels renameSeqlevels
#' @importFrom utils browseURL citation packageVersion
#' @importFrom rmarkdown render html_document
#' @importFrom GenomicRanges mcols
#' @importFrom knitrBootstrap knit_bootstrap
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
#' ## Load derfinder
#' library('derfinder')
#'
#' ## The output will be saved in the 'derfinderReport-example' directory
#' dir.create('derfinderReport-example', showWarnings = FALSE, recursive = TRUE)
#' 
#' ## For convenience, the derfinder output has been pre-computed
#' file.copy(system.file(file.path('extdata', 'chr21'), package='derfinder', 
#'     mustWork=TRUE), 'derfinderReport-example', recursive = TRUE)
#' 
#' \dontrun{
#' ## If you prefer, you can generate the output from derfinder
#' initialPath <- getwd()
#' setwd(file.path(initialPath, 'derfinderReport-example'))
#'
#' ## Collapse the coverage information
#' collapsedFull <- collapseFullCoverage(list(genomeData$coverage), 
#'     verbose=TRUE)
#' 
#' ## Calculate library size adjustments
#' sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5), nonzero=TRUE, 
#'     verbose=TRUE)
#' 
#' ## Build the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)
#'
#' ## Analyze chromosome 21
#' analyzeChr(chr='21', coverageInfo=genomeData, models=models, 
#'     cutoffFstat=1, cutoffType='manual', seeds=20140330, groupInfo=group, 
#'     mc.cores=1, writeOutput=TRUE, returnOutput=FALSE)
#'
#' ## Change the directory back to the original one
#' setwd(initialPath)
#' }
#'
#' ## Merge the results from the different chromosomes. In this case, there's 
#' ## only one: chr21
#' mergeResults(chrs = '21', prefix = 'derfinderReport-example', 
#'     genomicState = genomicState$fullGenome)
#'
#' ## Load the options used for calculating the statistics
#' load(file.path('derfinderReport-example', 'chr21', 'optionsStats.Rdata'))
#'
#' ## Generate the HTML report
#' report <- derfinderReport(prefix='derfinderReport-example', browse=FALSE, 
#'     nBestRegions=1, makeBestClusters=FALSE, 
#'     fullCov=list('21'=genomeDataRaw$coverage), optionsStats=optionsStats)
#'
#'
#' if(interactive()) {
#'     ## Browse the report
#'     browseURL(report)
#' }
#'
#' \dontrun{
#' ## Note that you can run the example using:
#' example('derfinderReport', 'regionReport', ask=FALSE)
#' }



derfinderReport <- function(prefix, outdir = 'basicExploration', 
    output = 'basicExploration', project = prefix, browse = interactive(),
    nBestRegions = 100, makeBestClusters = TRUE, nBestClusters = 2, 
    fullCov = NULL, hg19 = TRUE, p.ideos = NULL, txdb = NULL, 
    device = 'png', significantVar = 'qvalue', customCode = NULL,
    template = NULL, theme = NULL, digits = 2, ...) {
    
    stopifnot(length(significantVar) == 1)
    stopifnot(significantVar %in% c('pvalue', 'qvalue', 'fwer'))
    if(!is.null(theme)) stopifnot(is(theme, c('theme', 'gg')))
    
    hasCustomCode <- !is.null(customCode)
    if(hasCustomCode) stopifnot(length(customCode) == 1)
    
    ## Save start time for getting the total processing time
    startTime <- Sys.time()
    
    ## Advanced parameters
# @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
# \link[GenomeInfoDb]{seqlevelsStyle}.    
    chrsStyle <- .advanced_argument('chrsStyle', 'UCSC', ...)
    
# @param species Species name. See \link[derfinder]{extendedMapSeqlevels}  for 
# more information.
    species <- .advanced_argument('species', getOption('species', 'homo_sapiens'), ...)

# @param currentStyle Current naming style used. See 
# \link[derfinder]{extendedMapSeqlevels} for more information.
    currentStyle <- .advanced_argument('currentStyle', 'UCSC', ...)
    
    
# @param fullRegions Part of the output of \link[derfinder]{mergeResults}. 
# Specify it only if you have already loaded it in memory.
    fullRegions <- .advanced_argument('fullRegions', NULL, ...)
    
    
# @param fullNullSummary Part of the output of \link[derfinder]{mergeResults}. 
# Specify it only if you have already loaded it in memory.
    fullNullSummary <- .advanced_argument('fullNullSummary', NULL, ...)
    
    
# @param fullAnnotatedRegions Part of the output of 
# \link[derfinder]{mergeResults}. Specify it only if you have already loaded 
# it in memory.
    fullAnnotatedRegions <- .advanced_argument('fullAnnotatedRegions', NULL,
        ...)
    
    
# @param optionsStats Part of the output of \link[derfinder]{analyzeChr}. 
# Specify it only if you have already loaded it in memory.
    optionsStats <- .advanced_argument('optionsStats', NULL, ...)
    
    
# @param optionsMerge Part of the output of \link[derfinder]{mergeResults}. 
# Specify it only if you have already loaded it in memory.
    optionsMerge <- .advanced_argument('optionsMerge', NULL, ...)
    
    
# @param overviewParams A two element list with \code{base_size} and 
# \code{areaRel} that control the text size for the genomic overview plots.
    overviewParams <- .advanced_argument('overviewParam', list(base_size = 10,
        areaRel = 5), ...)
    
    
    ## Check that overviewParams is correctly specified
    stopifnot(sum(names(overviewParams) %in% c('base_size', 'areaRel')) == 2)
    
    ## Create outdir
    dir.create(file.path(prefix, outdir), showWarnings = FALSE,
        recursive = TRUE)
    workingDir <- file.path(getwd(), prefix)
    
    ## Locate Rmd if one is not provided
    if (is.null(template)) {
        templateNull <- TRUE
        template <- system.file(
            file.path('basicExploration', 'basicExploration.Rmd'),
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
    for(pkg in c('derfinderPlot', 'DT', 'ggbio', 'ggplot2')) load_install(pkg)
    
    ## Write bibliography information
    bib <- c(
        knitcitations = citation('knitcitations'), 
        derfinder = citation('derfinder')[1],
        derfinderPlot = citation('derfinderPlot')[1],
        regionReport = citation('regionReport')[1],
        DT = citation('DT'), 
        ggbio = citation('ggbio'),
        ggplot2 = citation('ggplot2'),
        knitr = citation('knitr')[3],
        rmarkdown = citation('rmarkdown'))
        
    write.bibtex(bib, file = file.path(prefix, outdir, paste0(output, '.bib')))
    
    ## Load files
    if (is.null(fullRegions)) 
        load(file.path(prefix, 'fullRegions.Rdata'))
    if (is.null(fullNullSummary)) 
        load(file.path(prefix, 'fullNullSummary.Rdata'))
    if (is.null(fullAnnotatedRegions)) 
        load(file.path(prefix, 'fullAnnotatedRegions.Rdata'))
    if (is.null(optionsStats)) 
        load(file.path(prefix, dir(prefix, pattern = 'chr')[1], 'optionsStats.Rdata'))
    if (is.null(optionsMerge)) 
        load(file.path(prefix, 'optionsMerge.Rdata'))
    
    ## Require fullCov
    if (makeBestClusters) {
        stopifnot(!is.null(fullCov))
        makeBestClusters <- !nBestClusters == 0
    } else {
        nBestClusters <- 0
    }
    
    ## Use UCSC names for homo_sapiens by default
    fullRegions <- renameSeqlevels(fullRegions, extendedMapSeqlevels(seqlevels(fullRegions), ...))
    names(fullCov) <- extendedMapSeqlevels(names(fullCov), ...)
    
    ##### Setup chunk options Are there any null regions? If not, then there 
    ##### won't be any p-values either.
    nullExist <- length(fullNullSummary) > 0
    fwerExist <- all(c('fwer', 'significantFWER') %in% colnames(mcols(fullRegions)))
    
    ## Were permutations used?
    seeds <- optionsStats$seeds
    usedPermutations <- length(optionsStats$nPermute) > 0 & !is.null(seeds)
    ## Are there significant regions?
    sigVar <- switch(significantVar, pvalue = 'significant', qvalue = 'significantQval', fwer = 'significantFWER')
    if(significantVar == 'fwer' & !fwerExist) {
        warning('There are no FWER adjusted P-values, will use FDR adjusted p-values instead.')
        significantVar <- 'significantQval'
    }
    pvalText <- switch(sigVar, significant = 'P-value', significantQval = 'FDR adjusted P-value', significantFWER = 'FWER adjusted P-value')
    pvalVar <- switch(sigVar, significant = 'pval', significantQval = 'qval', significantFWER = 'fwer')
    idx.sig <- as.logical(mcols(fullRegions)[[sigVar]])
    sigCut <- optionsMerge$significantCut[ifelse(sigVar == 'significantQval', 2, 1)]
    hasSig <- any(idx.sig)
    ## Are there regions with infite area?
    finite.area <- is.finite(fullRegions$area)
    hasArea <- any(idx.sig & finite.area)
    inf.area <- sum(!is.finite(fullRegions$area))
    
    ## Save the call
    theCall <- match.call()

    ## Generate report
    ## Perform code within the output directory.
    tmpdir <- getwd()
    with_wd(file.path(prefix, outdir), {
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
