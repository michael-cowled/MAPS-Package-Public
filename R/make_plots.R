#' Make Plots
#'
#' Executes a sequence of three plotting functions: barchart of total annotations,
#' cumulative histogram of unique annotations, and a bubble chart of metabolite classes
#' derived from MS2Query data. Plots are generated but not saved internally; they are
#' returned for subsequent saving by 'save_plots'.
#'
#' @param final_annotations A data frame containing processed feature annotations.
#'    Must include 'feature.ID', 'Samples', and 'compound.name'.
#' @param folder Path to the main output directory.
#' @param remove_sample Optional sample name to exclude.
#' @param ms2query_threshold MS2Query model prediction threshold (default 0.7).
#' @param noise_area_threshold Peak area threshold for quantification data (default 10000).
#'
#' @return A list containing all generated data frames and plot objects:
#'    'barchart_data', 'histogram_data', 'bubblechart_data',
#'    'barchart_plot', 'histogram_plot', and 'bubblechart_plot'.
#'
#' @export
#' @importFrom dplyr select rename mutate
make_plots <- function(
    final_annotations,
    folder,
    remove_sample = NULL,
    ms2query_threshold = 0.7,
    noise_area_threshold = 10000
) {

  # --- File Path Setup ---
  starburst.data <- final_annotations
  # 1. Check for Metadata file existence (case-sensitive)
  metadata_base_path <- "HGM/A - Analysis.csv"
  if (file.exists(metadata_base_path)) {
    metadata_file <- metadata_base_path
  } else {
    metadata_file <- NULL
    message("Metadata file 'HGM/A - Analysis.csv' not found. Barchart will use raw sample names.")
  }

  # 2. MS2Query file path
  ms2query_file <- paste0(folder, "/ms2query/ms2query.csv")

  # 3. Check for Quantification file existence (handling two common spellings)
  quant_file_base_1 <- paste0(folder, "/mzmine/DATA_iimn_gnps_quant.csv")
  quant_file_base_2 <- paste0(folder, "/mzmine/data_iimn_gnps_quant.csv") # lowercase d

  if (file.exists(quant_file_base_1)) {
    quant_file <- quant_file_base_1
  } else if (file.exists(quant_file_base_2)) {
    quant_file <- quant_file_base_2
  } else {
    quant_file <- NULL
    warning("Quantification file (DATA/data_iimn_gnps_quant.csv) not found in the mzmine folder. The bubble chart function will fail.")
  }

  # Initialize list to store results
  results <- list()

  # --- 1. BARCHART FUNCTION CALL ---

  message("\nStarting Barchart Generation (Feature Counts per Sample)...")
  # ASSUMPTION: generate_feature_barchart returns a list with $data and $plot
  barchart_results <- generate_feature_barchart(
    final_annotations = final_annotations,
    folder = folder,
    metadata.path = metadata_file
  )
  results$barchart_data <- barchart_results$data
  results$barchart_plot <- barchart_results$plot

  # --- 2. CUMULATIVE HISTOGRAM FUNCTION CALL ---

  message("\nStarting Cumulative Histogram Generation (Unique Annotations per Fraction)...")

  # Re-create the 'Annotations.with.samples' input
  Annotations.with.samples_input <- final_annotations %>%
    dplyr::select(feature.ID, Samples, compound.name) %>%
    dplyr::rename(annotation = compound.name) %>%
    dplyr::mutate(annotation = ifelse(is.na(annotation), 0, 1))

  # ASSUMPTION: generate_cumulative_histogram returns a list with $data and $plot
  histogram_results <- generate_cumulative_histogram(
    Annotations.with.samples = Annotations.with.samples_input,
    folder = folder,
    remove_sample = remove_sample # Pass NULL if not provided
  )
  results$histogram_data <- histogram_results$data
  results$histogram_plot <- histogram_results$plot

  # --- . STARBURST CHART FUNCTION CALL ---
  message("\nStarting Starburst Chart Generation (Class distribution for dataset)...")
  starburst_results <- generate_starburst(
    data = starburst.data
  )
  results$starburst_data <- starburst_results$data
  results$starburst_plot <- starburst_results$plot

  # --- 3. BUBBLE CHART FUNCTION CALL ---

  message("\nStarting Metabolite Class Bubble Chart Generation (MS2Query Data)...")

  if (is.null(quant_file) || !file.exists(ms2query_file)) {
    results$bubblechart_data <- NULL
    results$bubblechart_plot <- NULL
    warning("Skipping bubble chart generation because required files were not found.")
  } else {
    # ASSUMPTION: generate_ms2query_bubblechart returns a list with $data and $plot
    bubblechart_results <- generate_ms2query_bubblechart(
      ms2query.path = ms2query_file,
      quant.data.path = quant_file,
      folder = folder,
      threshold = ms2query_threshold,
      noise_area_threshold = noise_area_threshold
    )
    results$bubblechart_data <- bubblechart_results$data
    results$bubblechart_plot <- bubblechart_results$plot
  }



  # Return all data frames and plot objects
  return(results)
}
