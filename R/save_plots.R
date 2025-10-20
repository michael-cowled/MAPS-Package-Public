#' Save Generated Plots
#'
#' Saves the barchart, cumulative histogram, and bubble chart to a 'Plots'
#' subdirectory within the specified output structure, respecting the 'unimelb'
#' user domain environment variable for file path destination. **Existing plot
#' files with the same name are explicitly deleted before saving.**
#'
#' @param plot_list A named list containing the plot objects: 'barchart_plot',
#'   'histogram_plot', and 'bubblechart_plot'. (Assumes ggplot objects).
#' @param folder Path to the local main output directory.
#' @param dataset.id Character string identifying the dataset (not used for plot names, but kept for function structure).
#'
#' @return Invisibly returns NULL. Prints messages on success or failure.
#' @export
#' @importFrom ggplot2 ggsave
save_plots <- function(plot_list, folder, dataset.id) {

  message("\nStarting Plot Saving...")

  # Determine the base directory for saving the plots
  if (Sys.getenv("USER_DOMAIN") == "unimelb") {
    # Using a distinct subfolder for unimelb context
    plot_dir <- paste0(folder, "/Plots_Unimelb_Data")
  } else {
    plot_dir <- paste0(folder, "/Plots")
  }

  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  # Define plots and file names
  plot_configs <- list(
    list(plot = plot_list$barchart_plot, name = "feature_count_barchart.png", width = 8, height = 6),
    list(plot = plot_list$histogram_plot, name = "cumulative_annotation_histogram.png", width = 8, height = 6),
    list(plot = plot_list$bubblechart_plot, name = "ms2query_class_bubblechart.png", width = 10, height = 7)
  )

  for (config in plot_configs) {
    file_path <- paste0(plot_dir, "/", config$name)

    # CHECK AND DELETE EXISTING FILE (Explicit pre-flight deletion)
    if (file.exists(file_path)) {
      message(paste("Checking for and deleting existing plot file:", config$name))
      tryCatch({
        # NOTE: Using file.remove() here to explicitly delete the old file
        file.remove(file_path)
      }, error = function(e) {
        warning(paste("Could not delete existing file:", file_path, ". Overwriting anyway."))
      })
    }

    p <- config$plot
    if (!is.null(p)) {
      tryCatch({
        ggplot2::ggsave(
          filename = file_path,
          plot = p,
          width = config$width,
          height = config$height,
          units = "in"
        )
        message(paste("Saved plot:", file_path))
      }, error = function(e) {
        warning(paste("Failed to save", config$name, "plot:", e$message))
      })
    }
  }

  invisible(NULL)
}
