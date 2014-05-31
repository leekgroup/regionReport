## Based on https://github.com/hadley/ggplot2/blob/master/R/zzz.r
.onAttach <- function(...) {
    if (!interactive() || stats::runif(1) > 0.3) return()
        
    tips <- c(
        "Found a bug? Report it at https://github.com/lcolladotor/derfinderReport/issues",
        "Want to contribute a new feature? Fork derfinderReport at\n https://github.com/lcolladotor/derfinderReport/fork\nThen submit a pull request =)",
        paste("Find out what's changed in derfinderReport with\n",
            "news(Version == \"", utils::packageVersion("derfinderReport"),
            "\", package = \"derfinderReport\")", sep = ""),
        "Use suppressPackageStartupMessages to eliminate package startup messages."
    )
    
    tip <- sample(tips, 1, prob = c(0.3, 0.1, 0.5, 0.1))
    packageStartupMessage(tip)
}
