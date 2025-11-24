#' Validate and Get Paths
#'
#' Check if the folder exists and verify the presence of expected data files.
#' Prioritizes 'data_annotations_new.csv' over 'data_annotations.csv' for mzmine annotations.
#' Stops execution with informative errors if any required file is missing.
#'
#' @param folder A character string specifying the folder path containing the data files.
#' @param lv1 A logical value. If TRUE, includes paths and validates files for LV1 annotations (frag4, frag2, frag0, lv2_annotations).
#' @return A named list of valid file paths.
#' @export
validate_and_get_paths <- function(folder, lv1.subclasses = FALSE) {
  # Normalize path slashes for consistency
  folder <- gsub("\\\\", "/", folder)

  # --- 1. Handle Conditional File (mzmine_annotations) ---
  path_mzmine_annotations_final <- ""
  error_mzmine_annotations <- ""

  path_annotations_new <- paste0(folder, "/mzmine/data_annotations_frag6.csv")
  path_annotations_old <- paste0(folder, "/mzmine/data_annotations.csv")

  if (file.exists(path_annotations_new)) {
    # Prioritize new file if it exists
    path_mzmine_annotations_final <- path_annotations_new
    error_mzmine_annotations <- "The preferred file 'data_annotations_frag6.csv' was selected but is missing from the mzmine folder."
    message("Note: Using the preferred annotation file: 'data_annotations_frag6.csv'")
  } else {
    # Fallback to old file
    path_mzmine_annotations_final <- path_annotations_old
    error_mzmine_annotations <- paste(
      "The fallback file 'data_annotations.csv' is missing from the mzmine folder.",
      "The preferred file 'data_annotations_frag6.csv' was not found."
    )
    message("Note: Using the fallback annotation file: 'data_annotations.csv' (if it exists). 'data_annotations_frag6.csv' was not found.")
  }

  # --- 2. Assemble the REQUIRED file paths list ---
  file_paths <- list(
    mzmine_data = paste0(folder, "/mzmine/ms1-and-ms2.csv"),
    mzmine_annotations = path_mzmine_annotations_final, # Dynamically set
    canopus_data = paste0(folder, "/sirius/canopus_structure_summary.tsv"),
    csi_data = paste0(folder, "/sirius/structure_identifications_top-100.tsv"),
    zodiac_data = paste0(folder, "/sirius/formula_identifications.tsv"),
    ms2query_data = paste0(folder, "/ms2query/ms2query.csv"),
    cytoscape = paste0(folder, "/gnps/cytoscape.csv")
  )

  # --- 3. Assemble the REQUIRED error messages list ---
  error_messages <- list(
    mzmine_data = "The file 'ms1-and-ms2.csv' is missing from the mzmine folder.",
    mzmine_annotations = error_mzmine_annotations, # Dynamically set
    canopus_data = "The file 'canopus_structure_summary.tsv' is missing from the sirius folder.",
    csi_data = "The file 'structure_identifications_top-100.tsv' is missing from the sirius folder. Recompute for top K=100 hits.",
    zodiac_data = "The file 'formula_identifications.tsv' is missing from the sirius folder.",
    ms2query_data = "The file 'ms2query.csv' is missing from the ms2query folder.",
    cytoscape = "The file 'cytoscape.csv' is missing from the gnps folder."
  )

  # --- 4. Conditionally add LV1 files if lv1 == TRUE ---
  if (isTRUE(lv1.subclasses)) {

    # Conditionally add paths
    lv1_paths <- list(
      mzmine_annotations_4 = paste0(folder, "/mzmine/data_annotations_frag4.csv"),
      mzmine_annotations_2 = paste0(folder, "/mzmine/data_annotations_frag2.csv"),
      mzmine_annotations_0 = paste0(folder, "/mzmine/data_annotations_frag0.csv"),
      lv2_mzmine_annotations = paste0(folder, "/mzmine/lv2_annotations.csv")
    )
    file_paths <- c(file_paths, lv1_paths) # Append to main list

    # Conditionally add error messages
    lv1_errors <- list(
      mzmine_annotations_4 = "The file 'data_annotations_frag4.csv' is missing from the mzmine folder.",
      mzmine_annotations_2 = "The file 'data_annotations_frag2.csv' is missing from the mzmine folder.",
      mzmine_annotations_0 = "The file 'data_annotations_frag0.csv' is missing from the mzmine folder.",
      lv2_mzmine_annotations = "The file 'lv2_annotations.csv' is missing from the mzmine folder."
    )
    error_messages <- c(error_messages, lv1_errors) # Append to main list
  }


  # --- 5. Run the validation loop ---
  # The loop now runs over the combined list of REQUIRED + CONDITIONAL files
  for (name in names(file_paths)) {
    if (!file.exists(file_paths[[name]])) {
      stop(error_messages[[name]])
    }
  }

  return(file_paths)
}
