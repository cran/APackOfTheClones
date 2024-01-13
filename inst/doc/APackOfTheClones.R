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

## ----combineTCR, eval = FALSE-------------------------------------------------
#  library(scRepertoire)
#  
#  # load in the corresponding 6-sample TCR contigs from scRepertoire
#  contig_list <- get(data("contig_list", package = "scRepertoire"))
#  
#  # combine the TCR contigs into clones with custom samples
#  combined_contig_list <- scRepertoire::combineTCR(
#    contig_list,
#    samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"),
#    removeNA = FALSE,
#    removeMulti = FALSE,
#    filterMulti = FALSE
#  )

## ----combining, eval = FALSE--------------------------------------------------
#  # a seurat object corresponding to combined_contig_list named `pbmc` is loaded
#  pbmc <- scRepertoire::combineExpression(
#    combined_contig_list,
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

## ----umap, echo = TRUE--------------------------------------------------------
pbmc_umap_plot <- UMAPPlot(pbmc)
pbmc_umap_plot

## ----initial_vizapotc, echo = TRUE--------------------------------------------
default_apotc_plot <- vizAPOTC(pbmc, verbose = FALSE)
default_apotc_plot

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  reduction_base = "umap",
#  clonecall = "strict",

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  ...,
#  extra_filter = NULL,

## ----subsetting, echo = TRUE--------------------------------------------------
# `orig.ident` is a custom column in the example data with levels corresponding to sample ids:
# ("P17B" "P17L" "P18B" "P18L" "P19B" "P19L" "P20B" "P20L"). Here, it is subsetted
# by the keyword argument approach
subset_sample_17_plot <- vizAPOTC(
  pbmc, orig.ident = c("P17B", "P17L"), verbose = FALSE
)

# here, it is subsetted with `extra_filter` for sample 18 with dplyr syntax:
subset_sample_18_plot <- vizAPOTC(
  pbmc, extra_filter = "substr(orig.ident, 1, 3) == 'P18'", verbose = FALSE
)

# here, sample 19 is subsetted with both arguments to show that they work in conjunction
subset_sample_19_plot <- vizAPOTC(
  pbmc,
  orig.ident = "P19B",
  extra_filter = "orig.ident == 'P19L' | orig.ident == 'P19B'",
  verbose = FALSE
)

cowplot::plot_grid(
  default_apotc_plot,
  subset_sample_17_plot,
  subset_sample_18_plot,
  subset_sample_19_plot,
  labels = c("all", "17", "18", "19")
)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  repulse = TRUE,
#  repulsion_threshold = 1,
#  repulsion_strength = 1,
#  max_repulsion_iter = 10

## ----legend_params, echo = TRUE, eval = FALSE---------------------------------
#  add_size_legend = TRUE,
#  legend_sizes = "auto",
#  legend_position = "auto",
#  legend_buffer = 0.2,
#  legend_color = "#808080",
#  legend_spacing = "auto",
#  legend_label = "Clone sizes",
#  legend_text_size = 5,
#  add_legend_background = TRUE,

## ----other_params, echo = TRUE, eval = FALSE----------------------------------
#  order_clones = TRUE,
#  try_place = FALSE,
#  res = 360L,
#  linetype = "blank",
#  use_default_theme = TRUE,
#  retain_axis_scales = FALSE,
#  show_labels = FALSE,
#  label_size = 5,

## ----void_labelled_plot, echo = TRUE------------------------------------------
vizAPOTC(pbmc, show_labels = TRUE, use_default_theme = FALSE, verbose = FALSE)

