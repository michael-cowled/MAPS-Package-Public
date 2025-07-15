#' Validate and Get Paths
#'
#' @param folder Folder path to check
#' @return A character vector of valid file paths
#' @export

## 12. validate_and_get_paths
#' Validate and Get Paths
#'
#' Check if the folder exists and return valid file paths
#'
#' @param folder A character string specifying the folder path
#' @return A character vector of file paths
#' @export
validate_and_get_paths <- function(folder) {
  # Clean up the folder path (replace backslashes with forward slashes)
  # Note the double backslashes in the pattern for escaping the backslash character itself
  folder <- gsub("\\\\", "/", folder)
  
  if (dir.exists(folder)) {
  } else {
    stop(paste0("The folder", folder, "does not exist."))  # Stop execution
  }
  
  # Create a named list of paths for easier iteration and return
  file_paths <- list(
    mzmine_data = mzmine.data,
    mzmine_annotations = mzmine.annotations,
    canopus_data = canopus.data,
    csi_data = csi.data,
    zodiac_data = zodiac.data,
    ms2query_data = ms2query.data,
    cytoscape = cytoscape
  )
  
  # Define specific error messages for each file
  error_messages <- list(
    mzmine_data = "The file or folder is missing for mzmine data (ms1-and-ms2.csv).",
    mzmine_annotations = "The file or folder is missing for mzmine data (data_annotations.csv).",
    canopus_data = "The file or folder is missing for canopus data (canopus_structure_summary.tsv).",
    csi_data = "The file or folder is missing for csi:fingerID data (structure_identifications_top-100.tsv). Recompute for top K=100 hits.",
    zodiac_data = "The file or folder is missing for zodiac data (formula_identifications.tsv).",
    ms2query_data = "The file or folder is missing for ms2query data (ms2query.csv).",
    cytoscape = "The file or folder is missing for cytoscape.csv."
  )
  
  # Check if each file exists
  for (name in names(file_paths)) {
    if (!file.exists(file_paths[[name]])) {
      stop(error_messages[[name]]) # Stop execution with the specific error message
    }
  }
  
  # If all checks pass, return the list of file paths
  return(file_paths)
}
