% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AdjustAPOTC.R
\name{AdjustAPOTC}
\alias{AdjustAPOTC}
\title{Adjust the parameters of the APackOfTheClones reduction in a seurat
object}
\usage{
AdjustAPOTC(
  seurat_obj,
  reduction_base = NULL,
  clonecall = NULL,
  ...,
  extra_filter = NULL,
  run_id = NULL,
  new_rad_scale_factor = NULL,
  new_clone_scale_factor = NULL,
  repulse = FALSE,
  repulsion_threshold = 1,
  repulsion_strength = 1,
  max_repulsion_iter = 10L,
  relocate_cluster = NULL,
  relocation_coord = NULL,
  nudge_cluster = NULL,
  nudge_vector = NULL,
  recolor_cluster = NULL,
  new_color = NULL,
  rename_label = NULL,
  new_label = NULL,
  relocate_label = NULL,
  label_relocation_coord = NULL,
  nudge_label = NULL,
  label_nudge_vector = NULL,
  verbose = TRUE
)
}
\arguments{
\item{seurat_obj}{The seurat object to be adjusted.}

\item{reduction_base}{character. The seurat reduction to base the clonal
expansion plotting on. Defaults to \code{'umap'} but can be any reduction present
within the reductions slot of the input seurat object, including custom ones.
If `'pca'``, the cluster coordinates will be based on PC1 and PC2.
However, generally APackOfTheClones is used for displaying UMAP and
occasionally t-SNE versions to intuitively highlight clonal expansion.}

\item{clonecall}{character. The column name in the seurat object metadata to
use. See \code{scRepertoire} documentation for more information about this
parameter that is central to both packages.}

\item{...}{additional "subsetting" keyword arguments indicating the rows
corresponding to elements in the seurat object metadata that should be
filtered by. E.g., \code{seurat_clusters = c(1, 9, 10)} will filter the cells to
those in the \code{seurat_clusters} column with any of the values 1, 9, and 10.
Unfortunately, column names in the seurat object metadata cannot
conflict with the keyword arguments. \emph{\strong{MAJOR NOTE}} if any subsetting
keyword arguments are a \emph{prefix} of any preceding argument names (e.g. a
column named \code{reduction} is a prefix of the \code{reduction_base} argument)
R will interpret it as the same argument unless \emph{both} arguments
are named. Additionally, this means any subsequent arguments \emph{must} be named.}

\item{extra_filter}{character. An additional string that should be formatted
\emph{exactly} like a statement one would pass into \link[dplyr:filter]{dplyr::filter} that does
\emph{additional} filtering to cells in the seurat object - on top of the other
keyword arguments - based on the metadata. This means that it will be
logically AND'ed with any keyword argument filters. This is a more flexible
alternative / addition to the filtering keyword arguments. For example, if
one wanted to filter by the length of the amino acid sequence of TCRs, one
could pass in something like \code{extra_filter = "nchar(CTaa) - 1 > 10"}. When
involving characters, ensure to enclose with single quotes.}

\item{run_id}{character. This will be the ID associated with the data of a
run, and will be used by other important functions like \link{APOTCPlot} and
\link{AdjustAPOTC}. Defaults to \code{NULL}, in which case the ID will be generated
in the following format:

\verb{reduction_base;clonecall;keyword_arguments;extra_filter}

where if keyword arguments and extra_filter are underscore characters if
there was no input for the \code{...} and \code{extra_filter} parameters.}

\item{new_rad_scale_factor}{a single numeric in (0, 1]. changes the
radius scaling factor of all circles.}

\item{new_clone_scale_factor}{a single numeric in (0, 1]. changes the
clone_scale_factor}

\item{repulse}{If \code{TRUE}, will attempt to push overlapping clusters away from
each other.}

\item{repulsion_threshold}{numeric. The radius that clonal circle clusters
overlap is acceptable when repulsing.}

\item{repulsion_strength}{numeric. The smaller the value the less the
clusters repulse each other per iteration, and vice versa.}

\item{max_repulsion_iter}{integer. The number of repulsion iterations.}

\item{relocate_cluster}{Numeric or Character. Indicates which
cluster(s) based on the index or label to relocate to new coordinates.}

\item{relocation_coord}{numeric of length two or a list of numerics of length
two of length of \code{relocate_cluster}. If its a list, indicates each coordinate
that the clusters in \code{relocate_cluster} should move to. If its just a
numeric, then will relocate all clusters in \code{relocate_cluster} to the input,
which is likely not desired behavior, so this should only be convenience
syntax if \code{relocate_cluster} has length 1.}

\item{nudge_cluster}{Numeric or Character. Indicates which
cluster(s) based on the index or label to "nudge"/translate their
coordinate(s) by.}

\item{nudge_vector}{numeric of length two or a list of numerics of length
two of length of \code{nudge_cluster}. If its a list, indicates each translation
vector (in other words, x-y coordinates) that the clusters in
\code{nudge_cluster} should be translate by. If its just a numeric, then will
translate all clusters in \code{nudge_cluster} by the input - which mostly is
syntactic sugar for translating a single cluster if the input of
\code{nudge_cluster} is of length 1.}

\item{recolor_cluster}{Numeric or Character. Indicates which
cluster(s) based on the index or label to change their color by.}

\item{new_color}{character of arbitrary length. Indicates the corresponding
new colors that selected clusters in \code{recolor_cluster} should be changed to.}

\item{rename_label}{Numeric or character. Indicates the index or name of
label(s) to be renamed.}

\item{new_label}{Character. Indicates the corresponding new label(s) that
selected label(s) in \code{rename_label} should be changed to.}

\item{relocate_label}{Numeric or character. Indicates the index or name of
label(s) to be relocated.}

\item{label_relocation_coord}{Numeric of length two or a list of numerics of
length two of length of \code{relocate_label}. If it's a list, indicates each
coordinate that the labels in \code{relocate_label} should move to. If it's just
a numeric, then will relocate all labels in \code{relocate_label} to the input,
which is likely not desired behavior, so this should only be convenience
syntax if \code{relocate_label} has length 1.}

\item{nudge_label}{Numeric or character. Indicates the index or name of
label(s) to be "nudged"/translated.}

\item{label_nudge_vector}{Numeric of length two or a list of numerics of
length two of length of \code{nudge_label}. If it's a list, indicates each
translation vector (in other words, x-y coordinates) that the labels in
\code{nudge_label} should be translated by. If it's just a numeric, then will
translate all labels in \code{nudge_label} by the input - which mostly is
syntactic sugar for translating a single label if the input of
\code{nudge_label} is of length 1.}

\item{verbose}{logical. Decides if visual cues are displayed to the R console
of the progress.}
}
\value{
The adjusted \code{seurat_obj}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

If the user is unsatisfied with the clonal expansion plot that
was generated from \code{RunAPOTC} and \code{APOTCPlot}, this function has a range of
arguments to modify the data and/or parameters of the visualization. Note
that some of the arguments may conflict with eachother.
}
\examples{
# do an APackOfTheClones run
pbmc <- RunAPOTC(get(data("combined_pbmc")), verbose = FALSE)

# adjust the rad_scale_factor, and nudge cluster 1 by x = 1, y = 1
pbmc <- AdjustAPOTC(
    pbmc,
    new_rad_scale_factor = 0.9,
    nudge_cluster = 1,
    nudge_vector = c(1, 1),
    verbose = FALSE
)

# plot the result
APOTCPlot(pbmc)

# rename some labels
pbmc <- AdjustAPOTC(
    pbmc, rename_label = c(2, 5), new_label = c("Foo", "Bar")
)

# perhaps multiple clusters need to be relocated and repulsed
pbmc <- AdjustAPOTC(
    pbmc,
    relocate_cluster = c("Foo", "C10"), # using labels works too
    relocation_coord = list(c(2, 3.5), c(0, 5)),
    repulse = TRUE,
    verbose = FALSE
)

# plot again to check results
APOTCPlot(pbmc, show_labels = TRUE, verbose = FALSE)

}
