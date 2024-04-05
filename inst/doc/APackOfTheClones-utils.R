## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

knitr::opts_chunk$set(echo = FALSE)
options(repos = c(CRAN = "http://cran.rstudio.com"))

quiet_load_all_CRAN <- function(...) {
  for (pkg in list(...)) {
    if (require(pkg, quietly = TRUE, character.only = TRUE)) next
    invisible(install.packages(
      pkg, quiet = TRUE, verbose = FALSE, character.only = TRUE
    ))
    suppressPackageStartupMessages(invisible(
      require(pkg, quietly = TRUE, character.only = TRUE)
    ))
  }
}

# load packages
quiet_load_all_CRAN("ggplot2", "Seurat", "magrittr", "APackOfTheClones")

# load data
pbmc <- get(data("combined_pbmc"))

## ----shared_clones, echo = TRUE-----------------------------------------------
getSharedClones(pbmc, clonecall = "aa")

## ----reduction_centroid, echo = TRUE------------------------------------------
head(getReductionCentroids(pbmc, "umap"))

## ----highlight, echo = TRUE---------------------------------------------------
# create the APackOfTheClones plot
apotc_plot <- pbmc %>%
    vizAPOTC(clonecall = "aa", show_labels = TRUE, verbose = FALSE)

# get the shared clonotypes
shared_clonotypes <- pbmc %>%
    getSharedClones(clonecall = "aa") %>%
    names()

# highlight the first 3 shared clones
apotc_plot %>%
    showCloneHighlight(shared_clonotypes[1:3])

