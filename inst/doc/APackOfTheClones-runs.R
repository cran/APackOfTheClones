## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

knitr::opts_chunk$set(echo = FALSE)
options(repos = c(CRAN = "http://cran.rstudio.com"))

# utility functions
head <- function(df) {
  knitr::kable(utils::head(df))
}

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
quiet_load_all_CRAN("ggplot2", "cowplot", "Seurat")

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages(library(APackOfTheClones))

## ----setup_seurat, eval = FALSE-----------------------------------------------
#  library(scRepertoire)
#  
#  pbmc <- scRepertoire::combineExpression(
#    scRepertoire::combineTCR(
#      get(data("contig_list", package = "scRepertoire")),
#      samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"),
#      removeNA = FALSE,
#      removeMulti = FALSE,
#      filterMulti = FALSE
#    ),
#    pbmc,
#    cloneCall = "gene",
#    proportion = TRUE
#  )
#  
#  print(pbmc)

## ----actual_print_pbmc, eval = TRUE, echo = TRUE, include = FALSE-------------
# TODO use a nicer looking dataset
pbmc <- get(data("combined_pbmc"))
print(pbmc)

## ----eval = FALSE-------------------------------------------------------------
#  # Here is the function ran with its default parameters
#  pbmc <- RunAPOTC(pbmc)
#  
#  #> Initializing APOTC run...
#  #> * Setting `clone_scale_factor` to 0.3
#  #> * id for this run: umap;CTstrict;_;_
#  #>
#  #> Packing clones into clusters
#  #> [==================================================] 100%
#  #>
#  #> repulsing all clusters | max iterations = 20
#  #> [==================================================] 100%
#  #>
#  #> Completed successfully, time elapsed: 0.155 seconds
#  #>

## ----runapotc_default, include = FALSE----------------------------------------
pbmc <- RunAPOTC(pbmc, verbose = FALSE)

## ----runapotc2----------------------------------------------------------------
pbmc <- RunAPOTC(
    pbmc, run_id = "sample17", orig.ident = c("P17B", "P17L"), verbose = FALSE
)

## ----apotcplot_subset_params, eval = FALSE------------------------------------
#  reduction_base = NULL,
#  clonecall = NULL,
#  ...,
#  extra_filter = NULL,

## ----apotcplot----------------------------------------------------------------
# Here, plots for samples 17, 18, and 19 as seen in the previous vignette are made, where
# `orig.ident` is a custom column in the example data with levels corresponding to sample ids:
# ("P17B" "P17L" "P18B" "P18L" "P19B" "P19L" "P20B" "P20L"). 

pbmc <- RunAPOTC(
  pbmc, run_id = "P17", orig.ident = c("P17B", "P17L"), verbose = FALSE
)

pbmc <- RunAPOTC(
  pbmc, run_id = "P18", orig.ident = c("P18B", "P18L"), verbose = FALSE
)

pbmc <- RunAPOTC(
  pbmc, run_id = "P19", orig.ident = c("P19B", "P19L"), verbose = FALSE
)

cowplot::plot_grid(
  vizAPOTC(pbmc, verbose = FALSE),
  APOTCPlot(pbmc, run_id = "P17"),
  APOTCPlot(pbmc, run_id = "P18"),
  APOTCPlot(pbmc), # run_id omitted as sample 19 was the latest run
  labels = c("all", "17", "18", "19")
)

