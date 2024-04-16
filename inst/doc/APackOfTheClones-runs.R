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
quiet_load_all_CRAN("ggplot2", "cowplot", "Seurat", "magrittr")

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages(library(APackOfTheClones))
# load data
pbmc <- get(data("combined_pbmc"))

## ----load_data, eval = TRUE, echo = FALSE, include = FALSE--------------------
pbmc <- get(data("combined_pbmc"))

## ----setup_seurat, echo = TRUE, eval = FALSE----------------------------------
#  library(scRepertoire)
#  
#  # A seurat object named `pbmc` is loaded with a corresponding `contig_list`
#  pbmc <- scRepertoire::combineExpression(
#    scRepertoire::combineTCR(
#      contig_list,
#      samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"),
#      removeNA = FALSE,
#      removeMulti = FALSE,
#      filterMulti = FALSE
#    ),
#    pbmc,
#    cloneCall = "gene",
#    proportion = TRUE
#  )

## ----actual_print_pbmc, eval = TRUE, echo = TRUE------------------------------
print(pbmc)

## ----echo = TRUE, eval = FALSE------------------------------------------------
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

## ----runapotc2, echo = TRUE---------------------------------------------------
pbmc <- RunAPOTC(
    pbmc, run_id = "sample17", orig.ident = c("P17B", "P17L"), verbose = FALSE
)

## ----apotcplot_subset_params, eval = FALSE------------------------------------
#  reduction_base = NULL,
#  clonecall = NULL,
#  ...,
#  extra_filter = NULL,
#  alt_ident = NULL

## ----apotcplot, echo = TRUE---------------------------------------------------
# Here, plots for samples 17 - 20 as seen in the previous vignette are made, where
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

pbmc <- RunAPOTC(
  pbmc, run_id = "P20", orig.ident = c("P20B", "P20L"), verbose = FALSE
)

cowplot::plot_grid(
  APOTCPlot(pbmc, run_id = "P17", retain_axis_scales = TRUE, add_size_legend = FALSE),
  APOTCPlot(pbmc, run_id = "P18", retain_axis_scales = TRUE, add_size_legend = FALSE),
  APOTCPlot(pbmc, run_id = "P19", retain_axis_scales = TRUE, add_size_legend = FALSE),
  APOTCPlot(pbmc, retain_axis_scales = TRUE, add_size_legend = FALSE), # defaults to latest
  labels = c("17", "18", "19", "20")
)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  new_rad_scale_factor = NULL,
#  new_clone_scale_factor = NULL,
#  relocate_cluster = NULL,
#  relocation_coord = NULL,
#  nudge_cluster = NULL,
#  nudge_vector = NULL,
#  recolor_cluster = NULL,
#  new_color = NULL,
#  rename_label = NULL,
#  new_label = NULL,
#  relocate_label = NULL,
#  label_relocation_coord = NULL,
#  nudge_label = NULL,
#  label_nudge_vector = NULL,
#  verbose = TRUE

## ----first_four_labeled, echo = TRUE------------------------------------------
# Do a run with just the first 4 seurat clusters, and rename labels
pbmc <- RunAPOTC(
    pbmc,
    run_id = "first_four",
    seurat_clusters = 1:4,
    verbose = FALSE
)

pbmc <- AdjustAPOTC(
    pbmc,
    run_id = "first_four",
    rename_label = 1:4,
    new_label = letters[1:4],
    verbose = FALSE
)

APOTCPlot(
    pbmc,
    run_id = "first_four",
    show_labels = TRUE,
    retain_axis_scales = TRUE
)

## ----repulse_again, echo = TRUE-----------------------------------------------
pbmc <- pbmc %>%
    RunAPOTC(run_id = "foo", verbose = FALSE) %>%
    AdjustAPOTC(
        run_id = "foo",
        repulse = TRUE,
        repulsion_threshold = 0.5,
        verbose = FALSE
    )

APOTCPlot(
    pbmc,
    show_labels = TRUE,
    retain_axis_scales = TRUE,
    add_size_legend = FALSE
)

