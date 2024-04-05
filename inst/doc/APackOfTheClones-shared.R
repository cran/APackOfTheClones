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
head(getSharedClones(pbmc, clonecall = "aa"))

## ----shared_clones_links, echo = TRUE-----------------------------------------
# get shared amino acid shared clones:
shared_clones_aa <- getSharedClones(pbmc, clonecall = "aa")

# generate the plot
vizAPOTC(
    pbmc,
    clonecall = "aa",
    show_shared = shared_clones_aa,
    verbose = FALSE
)

## ----shared_clones_links_from_1_blend, echo = TRUE----------------------------
vizAPOTC(
    pbmc,
    clonecall = "aa",
    show_shared = shared_clones_aa,
    only_link = 3, # only link clonotypes from cluster 3
    clone_link_color = "blend",
    clone_link_width = 2,
    clone_link_alpha = 0.9,
    show_labels = TRUE,
    verbose = FALSE
)

## ----shared_clones_links_highlight, eval = TRUE, echo = FALSE-----------------
# For convenience, do an APackOfTheClones Run first
pbmc <- RunAPOTC(pbmc, clonecall = "aa", verbose = FALSE)

# get shared amino acid shared clones for the last run -
# note that the run_id can be replaced with `clonecall = "aa"`
shared_clones_aa_top4 <- getSharedClones(
  pbmc,
  run_id = getLastApotcDataId(pbmc),
  top = 4
)

# generate the unhighlighted plot
linked_apotc_plot <- APOTCPlot(
    pbmc,
    show_shared = shared_clones_aa_top4,
    verbose = FALSE
)

# highlight the top 4 clones with the viridis palette
# also slightly dimming other clones
showCloneHighlight(
    linked_apotc_plot,
    clonotype = names(shared_clones_aa_top4),
    color_each = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF"),
    default_color = NULL,
    scale_bg = 0.95
)

