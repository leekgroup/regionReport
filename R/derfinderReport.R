#' Generate a HTML report exploring the basic results from derfinder
#'
#' This function generates a HTML report exploring the basic results from
#' derfinder (https://github.com/lcolladotor/derfinder). The HTML report itself 
#' is generated using knitrBoostrap which uses knitr (http://yihui.name/knitr/) 
#' behind the scenes. It works best after using \link[derfinder]{mergeResults}.
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
#' @param makeBestClusters If \code{TRUE}, \link[derfinder]{plotCluster} is 
#' used on the \code{nBestClusters} regions by area. Note that these plots take 
#' some time to make.
#' @param nBestClusters The number of region cluster plots to make by taking 
#' the \code{nBestClusters} regions ranked by area of the cluster.
#' @param fullCov A list where each element is the result from 
#' \link[derfinder]{loadCoverage} used with \code{cutoff=NULL}. Can be 
#' generated using \link{fullCoverage}.
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
#' http://yihui.name/knitr/options (dev argument).
#' @param fullRegions Part of the output of \link[derfinder]{mergeResults}. 
#' Specify it only if you have already loaded it in memory.
#' @param fullNullSummary Part of the output of \link[derfinder]{mergeResults}. 
#' Specify it only if you have already loaded it in memory.
#' @param fullAnnotatedRegions Part of the output of 
#' \link[derfinder]{mergeResults}. Specify it only if you have already loaded 
#' it in memory.
#' @param optionsStats Part of the output of \link[derfinder]{analyzeChr}. 
#' Specify it only if you have already loaded it in memory.
#' @param optionsMerge Part of the output of \link[derfinder]{mergeResults}. 
#' Specify it only if you have already loaded it in memory.
#' @param overviewParams A two element list with \code{base_size} and 
#' \code{areaRel} that control the text size for the genomic overview plots.
#' @param chrsStyle The naming style of the chromosomes. By default, UCSC. See 
#' \link[GenomeInfoDb]{seqlevelsStyle}.
#'
#' @return An HTML report with a basic exploration of the results.
#'
#' @author Leonardo Collado-Torres
#' @seealso \link[derfinder]{mergeResults}, \link[derfinder]{analyzeChr}, 
#' \link[derfinder]{fullCoverage}
#' @export
#' @aliases derfinder_report
#'
#' @import derfinder
#' @import GenomicRanges
#' @import IRanges
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @importFrom ggbio plotIdeogram
#' @import rCharts
#' @import data.table
#' @import knitr
#' @import rmarkdown
#' @import knitrBootstrap
#' @import knitcitations
#' @import xtable
#' @import RColorBrewer
#' @import mgcv
#' @import GenomeInfoDb
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
#' file.copy(system.file(file.path("extdata", "chr21"), package="derfinder", 
#' mustWork=TRUE), 'derfinderReport-example', recursive=TRUE)
#' 
#' \dontrun{
#' ## If you prefer, you can generate the output from derfinder
#' initialPath <- getwd()
#' setwd(file.path(initialPath, 'derfinderReport-example'))
#'
#' ## Collapse the coverage information
#' collapsedFull <- collapseFullCoverage(list(genomeData$coverage), 
#' verbose=TRUE)
#' 
#' ## Calculate library size adjustments
#' sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5), nonzero=TRUE, 
#' verbose=TRUE)
#' 
#' ## Build the models
#' group <- genomeInfo$pop
#' adjustvars <- data.frame(genomeInfo$gender)
#' models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)
#'
#' ## Analyze chromosome 21
#' analyzeChr(chr='21', coverageInfo=genomeData, models=models, 
#' cutoffFstat=1, cutoffType='manual', seeds=20140330, groupInfo=group, 
#' mc.cores=1, writeOutput=TRUE, returnOutput=FALSE)
#'
#' ## Change the directory back to the original one
#' setwd(initialPath)
#' }
#'
#' ## Merge the results from the different chromosomes. In this case, there's 
#' ## only one: chr21
#' mergeResults(chrs='21', prefix='derfinderReport-example', 
#' genomicState=genomicState)
#'
#' ## Load the options used for calculating the statistics
#' load(file.path('derfinderReport-example', 'chr21', 'optionsStats.Rdata'))
#'
#' ## Generate the HTML report
#' report <- derfinderReport(prefix='derfinderReport-example', browse=FALSE, 
#' nBestRegions=1, makeBestClusters=FALSE, 
#' fullCov=list('21'=genomeDataRaw$coverage), optionsStats=optionsStats)
#'
#'
#' if(interactive()) {
#' ## Browse the report
#' browseURL(report)
#' }
#'
#' \dontrun{
#' ## Note that you can run the example using:
#' example('derfinderReport', 'regionReport', ask=FALSE)
#' }



derfinderReport <- function(prefix, outdir = "basicExploration", 
    output = "basicExploration", project = prefix, browse = interactive(),
    nBestRegions = 100, makeBestClusters = TRUE, nBestClusters = 2, 
    fullCov = NULL, hg19 = TRUE, p.ideos = NULL, txdb = NULL, 
    device = "CairoPNG", fullRegions = NULL, fullNullSummary = NULL,
    fullAnnotatedRegions = NULL, optionsStats = NULL, optionsMerge = NULL,
    overviewParams = list(base_size = 10, areaRel = 5), chrsStyle = "UCSC") {
    
    ## Save start time for getting the total processing time
    startTime <- Sys.time()
    
    ## Check that overviewParams is correctly specified
    stopifnot(sum(names(overviewParams) %in% c("base_size", "areaRel")) == 2)
    
    if (hg19) {
        library("biovizBase")
        library("TxDb.Hsapiens.UCSC.hg19.knownGene")
    }
    
    ## Create outdir
    dir.create(file.path(prefix, outdir), showWarnings = FALSE, recursive = TRUE)
    workingDir <- file.path(getwd(), prefix)
    
    ## Locate Rmd
    template <- system.file(file.path("basicExploration", "basicExploration.Rmd"), 
        package = "regionReport", mustWork = TRUE)
    
    ## Load knitcitations with a clean bibliography
    cleanbib()
    cite_options(tooltip = TRUE)
    
    ## Fix citep
    mycitep <- function(x, short = NULL) {
        res <- gsub("  </p>", "", citep(x))
        if (!is.null(short)) {
            res <- gsub("></a>", paste0(">", short, "</a>"), res)
        }
        return(res)
    }
    
    ## Write bibliography information
    write.bibtex(c(knitcitations = citation("knitcitations"), derfinder = citation("derfinder"), 
        regionReport = citation("regionReport"), knitrBootstrap = citation("knitrBootstrap"), 
        ggbio = citation("ggbio"), ggplot2 = citation("ggplot2"), rCharts = citation("rCharts"), 
        knitr = citation("knitr")[3], rmarkdown = citation("rmarkdown")), file = file.path(prefix, 
            outdir, "references.bib"))
    bib <- read.bibtex(file.path(prefix, outdir, "references.bib"))
    
    ## Load files
    if (is.null(fullRegions)) 
        load(file.path(prefix, "fullRegions.Rdata"))
    if (is.null(fullNullSummary)) 
        load(file.path(prefix, "fullNullSummary.Rdata"))
    if (is.null(fullAnnotatedRegions)) 
        load(file.path(prefix, "fullAnnotatedRegions.Rdata"))
    if (is.null(optionsStats)) 
        load(file.path(prefix, dir(prefix, pattern = "chr")[1], "optionsStats.Rdata"))
    if (is.null(optionsMerge)) 
        load(file.path(prefix, "optionsMerge.Rdata"))
    
    ## Require fullCov
    if (makeBestClusters) {
        stopifnot(!is.null(fullCov))
        makeBestClusters <- !nBestClusters == 0
    } else {
        nBestClusters <- 0
    }
    
    ##### Setup chunk options Are there any null regions? If not, then there 
    ##### won't be any p-values either.
    nullExist <- length(fullNullSummary) > 0
    ## Were permutations used?
    seeds <- optionsStats$seeds
    usedPermutations <- length(optionsStats$nPermute) > 0 & !is.null(seeds)
    ## Are there significant (by q-value) regions?
    idx.sig <- which(as.logical(fullRegions$significantQval))
    hasSig <- length(idx.sig) > 0
    ## Are there regions with infite area?
    finite.area <- which(is.finite(fullRegions$area))
    if (length(intersect(idx.sig, finite.area)) > 0) {
        hasArea <- TRUE
    } else {
        hasArea <- FALSE
    }
    inf.area <- sum(!is.finite(fullRegions$area))
    
    ## Save the call
    theCall <- match.call()
    
    ## knitrBoostrap chunk options
    opts_chunk$set(bootstrap.show.code = FALSE)
    
    ## Generate report
    tmpdir <- getwd()
    setwd(file.path(prefix, outdir))
    file.copy(template, to = paste0(output, ".Rmd"))
    res <- render(paste0(output, ".Rmd"))
    file.remove(paste0(output, ".Rmd"))
    
    ## Open
    if (browse) 
        browseURL(res)
    setwd(tmpdir)
    
    ## Finish
    return(invisible(res))
} 

#' @export
derfinder_report <- derfinderReport
