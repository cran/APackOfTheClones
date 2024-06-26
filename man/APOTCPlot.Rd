% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/APOTCPlot.R
\name{APOTCPlot}
\alias{APOTCPlot}
\title{Various variations of visualizations of clonal expansion post-RunAPOTC}
\usage{
APOTCPlot(
  seurat_obj,
  reduction_base = NULL,
  clonecall = NULL,
  ...,
  extra_filter = NULL,
  run_id = NULL,
  show_shared = NULL,
  only_link = NULL,
  clone_link_width = "auto",
  clone_link_color = "black",
  clone_link_alpha = 0.5,
  res = 360L,
  linetype = "blank",
  use_default_theme = TRUE,
  retain_axis_scales = FALSE,
  alpha = 1,
  show_labels = FALSE,
  label_size = 5,
  add_size_legend = TRUE,
  legend_sizes = "auto",
  legend_position = "auto",
  legend_buffer = 0.2,
  legend_color = "#808080",
  legend_spacing = "auto",
  legend_label = "Clone sizes",
  legend_text_size = 5,
  add_legend_background = TRUE,
  add_legend_centerspace = 0,
  detail = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{seurat_obj}{A seurat object that has been integrated with clonotype
data and has had a valid run of \link{RunAPOTC}.}

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

\item{show_shared}{The output of \link{getSharedClones} can be inputted here,
and the resulting plot will overlay lines between clone circles if that
clonotype is common between clusters. Note that the input \emph{\strong{must}} be
generated from data in the correct \code{APackOfTheClones} run, and the behavior
is undefined otherwise and will likely error. The next 4 arguments allow for
aesthetic customization of these line links.}

\item{only_link}{Optional integer indicating to only display clone
links originating from this cluster if showing shared clones.}

\item{clone_link_width}{numeric. The width of the lines that connect shared
clones. Defaults to \code{"auto"} which will estimate a reasonable value depending
on circle sizes.}

\item{clone_link_color}{character. The color of the lines that connect shared
clones. Defaults to \code{"blend"} which will use the average colors of the two
connected clones. Else, any hex color or valid color string input will work,
and the corresponding color will be applied on all links.}

\item{clone_link_alpha}{numeric. The alpha of the lines that connect shared
clones.}

\item{res}{The number of points on the generated path per full circle. From
plot viewers, if circles seem slightly too pixelated, it is recommended to
first try to export the plot as an \code{.svg} before increasing \code{res} due to
increased plotting times from \link[ggforce:geom_circle]{ggforce::geom_circle}.}

\item{linetype}{The type of outline each circle should have. defaults to
\verb{"blank} meaning no outline. More information is in the function
documentation of \code{ggforce::geom_circle}.}

\item{use_default_theme}{logical that defaults to \code{TRUE}. If \code{TRUE},
the resulting plot will have the same theme as the seurat reference reduction
plot. Else, the plot will simply have a blank background.}

\item{retain_axis_scales}{If \code{TRUE}, approximately maintains the axis scales
of the original reduction plot. However, it will only attempt to extend the
axes and never shorten. Users are recommended to set this to \code{TRUE}
especially if working with subsetted versions of the clonal data to better
preserve the geometric relation to the original dimensional reduction.}

\item{alpha}{numeric. The alpha of the circles in (0, 1]. Defaults to 1.}

\item{show_labels}{If \code{TRUE}, will label each circle cluster at the centroid,
defaulting to "C0, C1, ...".}

\item{label_size}{The text size of labels if shown. Defaults to 5.}

\item{add_size_legend}{If \code{TRUE}, adds a legend to the plot visualizing the
relative sizes of clones. Note that it is simply an overlay and not a real
ggplot2 legend.}

\item{legend_sizes}{numeric vector. Indicates the circle sizes to be
displayed on the legend, and will always be sorted from smallest to greatest.
Defaults to \code{"auto"} which estimate a reasonable range of sizes to display.}

\item{legend_position}{character or numeric. Can be set to either
\code{"top_left"}, \code{"top_right"}, \code{"bottom_left"}, \code{"bottom_right"} and places the
legend roughly in the corresponding position. Otherwise, can be a numeric
vector of length 2 indicating the x and y position of the \emph{topmost (smallest)
circle} of the legend.}

\item{legend_buffer}{numeric. Indicates how much to "push" the legend towards
the center of the plot from the selected corner. If negative, will push away}

\item{legend_color}{character. Indicates the hex color of the circles
displayed on the legend. Defaults to the hex code for a gray tone}

\item{legend_spacing}{numeric. Indicates the horizontal distance between each
stacked circle on the size legend. Defaults to \code{"auto"} which will use an
estimated value depending on plot size}

\item{legend_label}{character. The title of the legend, which defaults to
\verb{"clone sizes}.}

\item{legend_text_size}{numeric. The text size of the letters and numbers on
the legend}

\item{add_legend_background}{logical. If \code{TRUE}, will add a border around the
legend and fill the background to be white, overlaying anything else.}

\item{add_legend_centerspace}{numeric. An additional amount of distance
changed between the circle sizes on the left side of the legend and the
numbers on the right. Useful to set to around 0.5 (or more / less) when there
are particularly large clone sizes that may cover the numbers.}

\item{detail}{logical. If \code{FALSE}, will only plot entire clusters as one
large circle, which may be useful in cases where there are a high number
of clones resulting in a large number of circles on the resulting ggplot,
which has increased plotting times, and certain aspects of the plot needs
to be finely adjusted with \link{AdjustAPOTC} or simply inspected. This should
not be set to \code{FALSE} for the actual clonal expansion plot.}

\item{verbose}{logical. Decides if visual cues are displayed to the R console
of the progress.}
}
\value{
A ggplot object of the APackOfTheClones clonal expansion plot of the
seurat object. There is an additional 10th element in the object named
\code{"APackOfTheClones"} used by other functions in this package and shouldn't
interfere with any other ggplot functionality. (As far as currently known)
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

Given a seurat object with an 'apotc' (APackOfTheClones) object
from running \link{RunAPOTC}, this function will read the information and return
a customizable ggplot2 object of the clonal expansion with a circle size
legend. If the user is unhappy about certain aspects of the plot, many
parameters can be adjusted with the \link{AdjustAPOTC} function.

The specific APackOfTheClones run to be plotted can be identified in two
ways: either by inputting the \code{run_id} associated with the run that was
either defined / auto-generated during \link{RunAPOTC}, or by inputting the
\code{reduction_base}, \code{clonecall}, \code{extra_filter} and any other keyword arguments
that corresponded to the run. Its heavily recommended to use the \code{run_id}.
If none of these parameters are inputted, the function defaults to returning
the plot of the latest run.
}
\examples{
data("combined_pbmc")

combined_pbmc <- RunAPOTC(
    combined_pbmc, run_id = "run1", verbose = FALSE
)

# plotting with default arguments will plot the latest "run1"
clonal_packing_plot <- APOTCPlot(combined_pbmc)

}
\seealso{
\link{AdjustAPOTC}
}
