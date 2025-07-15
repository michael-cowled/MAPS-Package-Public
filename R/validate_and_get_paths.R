#' Validate and Get Paths
#'
#' Check if the folder exists and verify the presence of expected data files.
#' Stops execution with informative errors if any file is missing.
#'
#' @param folder A character string specifying the folder path containing the data files.
#' @return A named list of valid file paths.
#' @export
validate_and_get_paths <- function(folder) {
  # Normalize folder path by replacing backslashes with forward slashes
  folder <- gsub("\\\\", "/", folder)

  # Check if folder exists
  if (!dir.exists(folder)) {
    stop(paste0("The folder '", folder, "' does not exist."))
  }

  # Define expected files with their relative paths inside the folder
  file_paths <- list(
    mzmine_data = file.path(folder, "ms1-and-ms2.csv"),
    mzmine_annotations = file.path(folder, "data_annotations.csv"),
    canopus_data = file.path(folder, "canopus_structure_summary.tsv"),
    csi_data = file.path(folder, "structure_identifications_top-100.tsv"),
    zodiac_data = file.path(folder, "formula_identifications.tsv"),
    ms2query_data = file.path(folder, "ms2query.csv"),
    cytoscape = file.path(folder, "cytoscape.csv")
  )

  # Corresponding error messages for missing files
  error_messages <- list(
    mzmine_data = "The file 'ms1-and-ms2.csv' is missing from the folder.",
    mzmine_annotations = "The file 'data_annotations.csv' is missing from the folder.",
    canopus_data = "The file 'canopus_structure_summary.tsv' is missing from the folder.",
    csi_data = "The file 'structure_identifications_top-100.tsv' is missing from the folder. Recompute for top K=100 hits.",
    zodiac_data = "The file 'formula_identifications.tsv' is missing from the folder.",
    ms2query_data = "The file 'ms2query.csv' is missing from the folder.",
    cytoscape = "The file 'cytoscape.csv' is missing from the folder."
  )

  # Check if each expected file exists; stop with message if missing
  for (name in names(file_paths)) {
    if (!file.exists(file_paths[[name]])) {
      stop(error_messages[[name]])
    }
  }

  # Return named list of file paths if all files exist
  return(file_paths)
}
