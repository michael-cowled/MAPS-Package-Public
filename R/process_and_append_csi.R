#' @title Process and Append CSI:FingerID Annotations
#' @description Cleans, standardizes, and appends CSI:FingerID annotations to an existing data frame.
#'
#' @param csi.data A data frame containing CSI:FingerID annotations.
#' @param existing_annotations A data frame of existing annotations (lv1, lv2, lv3).
#' @param csi.prob A numeric threshold for CSI:FingerID confidence scores.
#' @param cid_cache_df A data frame used for caching CID lookups.
#' @param lipids.file Path to a lipids file for standardization.
#' @param cid_database_path Path to the PubChem CID database.
#' @param compute_id_prob_func A function to compute identification probability (e.g., MAPS.Package::compute_id_prob).
#' @param deduplicate_data_func A function to remove duplicate annotations (e.g., MAPS.Package::deduplicate_data).
#' @param standardize_annotation_func A function to standardize annotations (e.g., MAPS.Package::standardize_annotation).
#'
#' @return A list containing the updated annotations data frame and the updated CID cache.
#' @export
process_and_append_csi <- function(
    csi.data,
    existing_annotations,
    csi.prob,
    cid_cache_df,
    lipids.file,
    cid_database_path,
    compute_id_prob_func,
    deduplicate_data_func,
    standardize_annotation_func
) {
  # Add missing columns before processing to prevent errors
  missing_cols <- colnames(existing_annotations)[!colnames(existing_annotations) %in% colnames(csi.data)]
  for (col in missing_cols) {
    csi.data[[col]] <- NA
  }

  # --- Data Cleaning and Standardization ---
  csi.data$confidence.score[csi.data$confidence.score == -Inf] <- 0
  csi.data$confidence.score <- as.numeric(csi.data$confidence.score)

  # Filter based on confidence probability
  csi.data <- csi.data %>%
    dplyr::filter(confidence.score >= csi.prob) %>%
    dplyr::filter(!grepl("PUBCHEM", compound.name, ignore.case = TRUE))

  # Compute ID probability
  csi.data <- compute_id_prob_func(csi.data, "confidence.score", csi.prob)

  # Deduplicate and standardize
  csi.data <- deduplicate_data_func(csi.data, compound.name, confidence.score)
  result <- standardize_annotation_func(
    csi.data,
    name_col = "compound.name",
    smiles_col = "smiles",
    cid_cache_df = cid_cache_df,
    lipids.file = lipids.file,
    cid_database_path = cid_database_path
  )

  csi.data <- result$data
  updated_cid_cache_df <- result$cache

  # Add missing columns after standardization
  if (!("Formula" %in% names(csi.data))) csi.data$Formula <- NA_character_
  if (!("IUPAC" %in% names(csi.data))) csi.data$IUPAC <- NA_character_
  if (!("Monoisotopic.Mass" %in% names(csi.data))) csi.data$Monoisotopic.Mass <- NA_real_

  # --- Filtering and Appending ---
  unique_in_csi <- setdiff(csi.data$feature.ID, existing_annotations$feature.ID)
  csi.data <- csi.data %>% dplyr::filter(feature.ID %in% unique_in_csi)

  # Set annotation metadata
  csi.data$annotation.type <- "CSI:FingerID"
  csi.data$confidence.level <- "3"

  csi.data$feature.ID <- as.numeric(csi.data$feature.ID)
  csi.data$confidence.score <- as.numeric(csi.data$confidence.score)

  # Append the filtered and processed data to the existing annotations
  updated_annotations <- existing_annotations %>%
    dplyr::bind_rows(csi.data)

  return(list(annotations = updated_annotations, cache = updated_cid_cache_df))
}
