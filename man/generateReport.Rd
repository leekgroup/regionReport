\name{generateReport}
\alias{generateReport}
\title{Generate a HTML report exploring the basic results from derfinder}
\usage{
generateReport(prefix, outdir = "basicExploration",
  output = "basicExploration.html", project = prefix,
  browse = interactive(), nBestRegions = 100, makeBestClusters = TRUE,
  nBestClusters = 2, fullCov = NULL, hg19 = TRUE, p.ideos = NULL,
  txdb = NULL, installMissing = TRUE, device = "CairoPNG",
  fullRegions = NULL, fullNullSummary = NULL, fullAnnotatedRegions = NULL,
  optionsStats = NULL, optionsMerge = NULL,
  overviewParams = list(base_size = 10, areaRel = 5))
}
\arguments{
  \item{prefix}{The main data directory path where
  \link[derfinder]{mergeResults} was run. It should be the
  same as \code{mergeResults(prefix)}.}

  \item{outdir}{The name of output directory relative to
  \code{prefix}.}

  \item{output}{The name output HTML file.}

  \item{project}{The title of the project.}

  \item{browse}{If \code{TRUE} the HTML report is opened in
  your browser once it's completed.}

  \item{nBestRegions}{The number of region plots to make,
  ordered by area.}

  \item{makeBestClusters}{If \code{TRUE},
  \link[derfinder]{plotCluster} is used on the
  \code{nBestClusters} regions by area. Note that these
  plots take some time to make.}

  \item{nBestClusters}{The number of region cluster plots
  to make by taking the \code{nBestClusters} regions ranked
  by area of the cluster.}

  \item{fullCov}{A list where each element is the result
  from \link[derfinder]{loadCoverage} used with
  \code{cutoff=NULL}. The elements of the list should be
  named according to the chromosome number. Can be
  generated using \link{fullCoverage}.}

  \item{hg19}{If \code{TRUE} then the reference is assumed
  to be hg19 and chromosome lengths as well as the default
  transcription database
  (TxDb.Hsapiens.UCSC.hg19.knownGene) will be used.}

  \item{p.ideos}{A list where each element is the result of
  \link[ggbio]{plotIdeogram}. If it's \code{NULL} and
  \code{hg19=TRUE} then they are created for the hg19 human
  reference.}

  \item{txdb}{Specify the transcription database to use for
  making the plots for the top regions by area. If
  \code{NULL} and \code{hg19=TRUE} then
  TxDb.Hsapiens.UCSC.hg19.knownGene is used.}

  \item{installMissing}{If \code{TRUE} all missing required
  packages are installed. Note that some are development
  versions hosted in GitHub.}

  \item{device}{The graphical device used when knitting.
  See more at http://yihui.name/knitr/options (dev
  argument).}

  \item{fullRegions}{Part of the output of
  \link[derfinder]{mergeResults}. Specify it only if you
  have already loaded it in memory.}

  \item{fullNullSummary}{Part of the output of
  \link[derfinder]{mergeResults}. Specify it only if you
  have already loaded it in memory.}

  \item{fullAnnotatedRegions}{Part of the output of
  \link[derfinder]{mergeResults}. Specify it only if you
  have already loaded it in memory.}

  \item{optionsStats}{Part of the output of
  \link[derfinder]{analyzeChr}. Specify it only if you have
  already loaded it in memory.}

  \item{optionsMerge}{Part of the output of
  \link[derfinder]{mergeResults}. Specify it only if you
  have already loaded it in memory.}

  \item{overviewParams}{A two element list with
  \code{base_size} and \code{areaRel} that control the text
  size for the genomic overview plots.}
}
\value{
An HTML report with a basic exploration of the results.
}
\description{
This function generates a HTML report exploring the basic
results from derfinder
(https://github.com/lcolladotor/derfinder). The HTML report
itself is generated using knitrBoostrap
(https://github.com/jimhester/knitrBootstrap) which uses
knitr (http://yihui.name/knitr/) behind the scenes. It
works best after using \link[derfinder]{mergeResults}.
}
\examples{
## The output will be saved in the 'generateReport-example' directory
dir.create("generateReport-example", showWarnings = FALSE, recursive = TRUE)
initialPath <- getwd()
setwd(file.path(initialPath, "generateReport-example"))

## Generate output from derfinder

## Collapse the coverage information
collapsedFull <- collapseFullCoverage(list(genomeData$coverage), verbose=TRUE)

## Calculate library size adjustments
sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5), nonzero=TRUE, verbose=TRUE)

## Build the models
group <- genomeInfo$pop
adjustvars <- data.frame(genomeInfo$gender)
models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)

## Analyze chromosome 21
analyzeChr(chrnum="21", coverageInfo=genomeData, models=models, cutoffFstat=1, cutoffType="manual", seeds=20140330, groupInfo=group, mc.cores=1, writeOutput=TRUE, returnOutput=FALSE)

\dontrun{
## Create the GenomicState object for Hsapiens.UCSC.hg19.knownGene
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## Creating this GenomicState object takes around 8 min
GenomicState.Hsapiens.UCSC.hg19.knownGene <- makeGenomicState(txdb=txdb)

## For convinience, this object is already included in derfinderReport
identical(GenomicState.Hsapiens.UCSC.hg19.knownGene, genomicState)

## Rename to match code outside the don'r run section
genomicState <- GenomicState.Hsapiens.UCSC.hg19.knownGene
}

## Change the directory back to the original one
setwd(initialPath)

## Merge the results from the different chromosomes. In this case, there's only one: chr21
mergeResults(chrnums="21", prefix="generateReport-example", genomicState=genomicState)

## Load the options used for calculating the statistics
load(file.path("generateReport-example", "chr21", "optionsStats.Rdata"))

## Generate the HTML report
generateReport(prefix="generateReport-example", nBestRegions=15, makeBestClusters=TRUE, fullCov=list("21"=genomeDataRaw$coverage), optionsStats=optionsStats)
}
\author{
Leonardo Collado-Torres
}
\seealso{
\link[derfinder]{mergeResults},
\link[derfinder]{analyzeChr},
\link[derfinder]{fullCoverage}
}

