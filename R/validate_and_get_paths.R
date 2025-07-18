#' Validate and Get Paths
#'
#' Check if the folder exists and verify the presence of expected data files.
#' Stops execution with informative errors if any file is missing.
#'
#' @param folder A character string specifying the folder path containing the data files.
#' @return A named list of valid file paths.
#' @export
validate_and_get_paths <- function(folder) {
  # Normalize path slashes for consistency
  folder <- gsub("\\\\", "/", folder)

  # Correct subfolder paths for each file
  file_paths <- list(
    mzmine_data = paste0(folder, "/mzmine/ms1-and-ms2.csv"),
    mzmine_annotations = paste0(folder, "/mzmine/data_annotations.csv"),
    canopus_data = paste0(folder, "/sirius/canopus_structure_summary.tsv"),
    csi_data = paste0(folder, "/sirius/structure_identifications_top-100.tsv"),
    zodiac_data = paste0(folder, "/sirius/formula_identifications.tsv"),
    ms2query_data = paste0(folder, "/ms2query/ms2query.csv"),
    cytoscape = paste0(folder, "/gnps/cytoscape.csv")
  )

  error_messages <- list(
    mzmine_data = "The file 'ms1-and-ms2.csv' is missing from the mzmine folder.",
    mzmine_annotations = "The file 'data_annotations.csv' is missing from the mzmine folder.",
    canopus_data = "The file 'canopus_structure_summary.tsv' is missing from the sirius folder.",
    csi_data = "The file 'structure_identifications_top-100.tsv' is missing from the sirius folder. Recompute for top K=100 hits.",
    zodiac_data = "The file 'formula_identifications.tsv' is missing from the sirius folder.",
    ms2query_data = "The file 'ms2query.csv' is missing from the ms2query folder.",
    cytoscape = "The file 'cytoscape.csv' is missing from the gnps folder."
  )

  for (name in names(file_paths)) {
    if (!file.exists(file_paths[[name]])) {
      stop(error_messages[[name]])
    }
  }

  return(file_paths)
}
