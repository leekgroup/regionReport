## Based on https://github.com/hadley/ggplot2/blob/master/R/zzz.r
.onAttach <- function(...) {
    if (!interactive() || stats::runif(1) > 0.99) return()
        
    tips <- c(
        "Want to contribute a new feature? Fork regionReport at\n https://github.com/lcolladotor/regionReport/fork\nThen submit a pull request =)",
        paste("Find out what's changed in regionReport with\n",
            "news(Version == \"", utils::packageVersion("regionReport"),
            "\", package = \"regionReport\")", sep = ""),
        "Use suppressPackageStartupMessages to eliminate package startup messages."
    )
    
    tip <- sample(tips, 1, prob = c(0.1, 0.8, 0.1))
    packageStartupMessage(tip)
}
