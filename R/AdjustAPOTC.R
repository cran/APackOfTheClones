#' @title
#' Adjust the parameters of the APackOfTheClones reduction in a seurat
#' object
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' If the user is unsatisfied with the clonal expansion plot that
#' was generated from `RunAPOTC` and `APOTCPlot`, this function has a range of
#' arguments to modify the data and/or parameters of the visualization. Note
#' that some of the arguments may conflict with eachother.
#'
#' @inheritParams RunAPOTC
#'
#' @param seurat_obj The seurat object to be adjusted.
#' @param new_clone_scale_factor a single numeric in (0, 1]. changes the
#' clone_scale_factor
#' @param new_rad_scale_factor a single numeric in (0, 1]. changes the
#' radius scaling factor of all circles.
#' @param relocate_cluster numeric of arbitrary length. Indicates which
#' cluster(s) to relocate to new coordinates
#' @param relocation_coord numeric of length two or a list of numerics of length
#' two of length of `relocate_cluster`. If its a list, indicates each coordinate
#' that the clusters in `relocate_cluster` should move to. If its just a
#' numeric, then will relocate all clusters in `relocate_cluster` to the input,
#' which is likely not desired behavior, so this should only be convenience
#' syntax if `relocate_cluster` has length 1.
#' @param nudge_cluster numeric of arbitrary length. Indicates which
#' cluster(s) to "nudge"/translate their coordinate(s) by.
#' @param nudge_vector numeric of length two or a list of numerics of length
#' two of length of `nudge_cluster`. If its a list, indicates each translation
#' vector (in other words, x-y coordinates) that the clusters in
#' `nudge_cluster` should be translate by. If its just a numeric, then will
#' translate all clusters in `nudge_cluster` by the input - which mostly is
#' syntactic sugar for translating a single cluster if the input of
#' `nudge_cluster` is of length 1.
#' @param recolor_cluster numeric of arbitrary length. Indicates which
#' cluster(s) to change their color by.
#' @param new_color character of arbitrary length. Indicates the corresponding
#' new colors that selected clusters in `recolor_cluster` should be changed to.
#'
#' @return The adjusted `seurat_obj`
#' @export
#'
#' @examples
#' # do an APackOfTheClones run
#' pbmc <- RunAPOTC(get(data("combined_pbmc")), verbose = FALSE)
#'
#' # adjust the rad_scale_factor, and nudge cluster 1 by (1, 1)
#' pbmc <- AdjustAPOTC(
#'     pbmc,
#'     new_rad_scale_factor = 0.9,
#'     nudge_cluster = 1,
#'     nudge_vector = c(1, 1),
#'     verbose = FALSE
#' )
#'
#' # plot the result
#' APOTCPlot(pbmc)
#'
#' # perhaps multiple clusters need to be relocated and repulsed
#' pbmc <- AdjustAPOTC(
#'     pbmc,
#'     relocate_cluster = c(5, 10),
#'     relocation_coord = list(c(2, 3.5), c(0, 5)),
#'     repulse = TRUE,
#'     verbose = FALSE
#' )
#'
#' # plot again to check results
#' APOTCPlot(pbmc)
#'
AdjustAPOTC <- function(
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

	relocate_cluster = NULL, # can also be a vector
	relocation_coord = NULL, # vector or list of vectors

	nudge_cluster = NULL, # same as above
	nudge_vector = NULL,

	recolor_cluster = NULL, #same as above
	new_color = NULL,

	#interactive = FALSE,
	verbose = TRUE
) {
	varargs_list <- list(...)
	args <- hash::hash(as.list(environment()))
	AdjustAPOTC_error_handler(args)

	object_id <- infer_object_id_if_needed(args, varargs_list = varargs_list)
	apotc_obj <- getApotcData(seurat_obj, object_id)
	args <- hash::hash(as.list(environment()))

	# # TODO
	# if (interactive) {
	# 	return(runShinyAdjustAPOTC(args))
	# }

	if (should_change(new_clone_scale_factor)) {
		apotc_obj <- change_clone_scale(seurat_obj, args)
	}

	if (should_change(new_rad_scale_factor)) {
		apotc_obj <- change_rad_scale(apotc_obj, new_rad_scale_factor)
	}

	if (should_change(recolor_cluster)) {
		apotc_obj <- recolor_clusters(apotc_obj, recolor_cluster, new_color)
	}

	if (should_change(relocate_cluster)) {
		apotc_obj <- relocate_clusters(
			apotc_obj, relocate_cluster, relocation_coord
		)
	}

	if (should_change(nudge_cluster)) {
		apotc_obj <- nudge_clusters(apotc_obj, nudge_cluster, nudge_vector)
	}

	if (repulse) {
		apotc_obj <- repulseClusters(
			apotc_obj, repulsion_threshold, repulsion_strength,
			max_repulsion_iter, verbose
		)
	}

	seurat_obj <- setApotcData(seurat_obj, object_id, apotc_obj)
	seurat_obj
}

AdjustAPOTC_error_handler <- function(args) {
	# TODO - type check
	# TODO rest of errors

	check_apotc_identifiers(args)

	if (should_change(args$new_rad_scale_factor)) {
		if (args$new_rad_scale_factor < 0) {
			stop(call. = FALSE, "new_rad_scale_factor must be a positive number")
		}
	}

	if (should_change(args$new_clone_scale_factor)) {
		if (args$new_clone_scale_factor < 0) {
			stop(call. = FALSE, "new_rad_scale_factor must be a positive number")
		}
	}

	check_repulsion_params(args)

	# if (!is.null(args$relocation_coord) && (length(args$relocate_cluster) != length(args$relocation_coord))) {
	# return("length of relocate_cluster must be the same as the length of relocation_coord")
	# }
	# if (!is.null(args$nudge_vector) && (length(args$nudge_cluster) != length(args$nudge_vector))) {
	# return("length of nudge_cluster must be the same as the length of nudge_vector")
	# }
	# if ((args$relocate_cluster != -1) && (args$nudge_cluster != -1)) {
	# 	if (!is.null(intersect(args$relocation_cluster, args$nudge_cluster))) {
	# 		return("There are repeated elements in relocate_cluster and/or nudge_cluster")
	# 	}
	# }
}

# TODO it should be more efficient + safer to mathematically transform all vals
change_clone_scale <- function(seurat_obj, args) {

	if (args$verbose) {
		message("Repacking all clusters with new clone scale factor")
	}

	past_params <- find_seurat_command(
		seurat_obj = seurat_obj,
		func_name = "RunAPOTC",
		id = args$object_id
	)@params

	args$apotc_obj@clone_scale_factor <- args$new_clone_scale_factor

	circlepackClones(
		apotc_obj = args$apotc_obj,
		ORDER = past_params$order_clones,
		try_place = past_params$try_place,
		verbose = args$verbose
	)
}

change_rad_scale <- function(apotc_obj, new_factor) {

	old_factor <- get_rad_scale_factor(apotc_obj)
	conversion_num <- get_clone_scale_factor(apotc_obj) *
		(new_factor - old_factor)

	for (i in seq_len(get_num_clusters(apotc_obj))) {
		curr <- apotc_obj@clusters[[i]]
    	if (isnt_empty(curr)) {
    		apotc_obj@clusters[[i]]$rad <- curr$rad + conversion_num
    	}
	}

	apotc_obj@rad_scale_factor <- new_factor
	apotc_obj
}

recolor_clusters <- function(apotc_obj, recolor_cluster, new_color) {
	for (i in seq_along(recolor_cluster)) {
		apotc_obj@cluster_colors[recolor_cluster[i]] <- new_color[[i]]
	}
	apotc_obj
}

#FIXME (again :/)

relocate_clusters <- function(apotc_obj, relocate_cluster, relocation_coord) {

	if (is.numeric(relocation_coord)) {
		relocation_coord <- init_list(length(relocate_cluster), relocation_coord)
	}

	new_clusterlists <- get_clusterlists(apotc_obj)

	for (i in seq_along(relocate_cluster)) {
		cl_ind <- relocate_cluster[i]
		new_clusterlists[[cl_ind]] <- move_cluster(
	    	cluster = new_clusterlists[[cl_ind]],
	    	new_coord = relocation_coord[[i]]
	    )
	}

	setModifiedClusterlists(apotc_obj, new_clusterlists)
}

nudge_clusters <- function(apotc_obj, nudge_cluster, nudge_vector) {

	if (is.numeric(nudge_vector)) {
		nudge_vector <- init_list(length(nudge_cluster), nudge_vector)
	}

	relocate_clusters(
		apotc_obj,
		relocate_cluster = nudge_cluster,
		relocation_coord = operate_on_same_length_lists(
			func = add,
			l1 = nudge_vector,
			l2 = get_centroids(apotc_obj)[nudge_cluster]
		)
	)
}

# TODO label movement

# need functions for readjusting the apotc reduction for better visuals
# also possible to boot up a shiny window in the future?
