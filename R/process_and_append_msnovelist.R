#' @title Process and Append MSNovelist Annotations
#' @description Cleans, standardizes, and appends MSNovelist de novo annotations to an existing data frame.
#'
#' @param msn.data Path to the MSNovelist tsv file.
#' @param existing_annotations A data frame of existing annotations.
#' @param msn.threshold A numeric threshold for MSNovelist ModelScore.
#' @param cid_cache_df A data frame used for caching CID lookups.
#' @param lipids.file Path to a lipids file for standardization.
#' @param cid_database_path Path to the PubChem CID database.
#' @param compute_id_prob A function to compute identification probability.
#' @param deduplicate_data A function to remove duplicate annotations.
#' @param standardise_annotation A function to standardize annotations.
#'
#' @return A list containing the updated annotations data frame and the updated CID cache.
#' @export
process_and_append_msnovelist <- function(
    msn.data,
    existing_annotations,
    msn.threshold,
    cid_cache_df,
    lipids.file,
    cid_database_path,
    compute_id_prob,
    deduplicate_data,
    standardise_annotation,
    standardisation
) {
  # Data Cleaning and Initial Processing
  msn.data <- read_checked_tsv(msn.data)

  # Use dplyr::select for robust column mapping instead of hardcoded indices
  msn.data <- msn.data %>%
    dplyr::select(
      confidence.score = ModelScore, # Adjust if you prefer SiriusScore here
      compound.name = name,
      smiles = smiles,
      feature.ID = mappingFeatureId
    )

  # Remove features already annotated in previous levels
  msn.data <- msn.data[!(msn.data$feature.ID %in% existing_annotations$feature.ID), ]

  # Filter based on ModelScore threshold
  msn.data$confidence.score <- as.numeric(msn.data$confidence.score)
  msn.data <- dplyr::filter(msn.data, confidence.score >= msn.threshold)

  # Add missing columns before processing to prevent errors
  missing_cols <- colnames(existing_annotations)[!colnames(existing_annotations) %in% colnames(msn.data)]
  for (col in missing_cols) {
    msn.data[[col]] <- NA
  }

  # Compute ID probability (Ensure this helper handles negative scores if using ModelScore!)
  msn.data <- compute_id_prob(msn.data, "confidence.score", msn.threshold)

  # Deduplicate and standardize
  msn.data <- deduplicate_data(msn.data, compound.name, confidence.score)
  result <- standardise_annotation(
    msn.data,
    name_col = "compound.name",
    smiles_col = "smiles",
    cid_cache_df = cid_cache_df,
    lipids.file = lipids.file,
    cid_database_path = cid_database_path,
    standardisation = standardisation
  )

  msn.data <- result$data
  updated_cid_cache_df <- result$cache

  # Add missing columns after standardization
  if (!("Formula" %in% names(msn.data))) msn.data$Formula <- NA_character_
  if (!("IUPAC" %in% names(msn.data))) msn.data$IUPAC <- NA_character_
  if (!("Monoisotopic.Mass" %in% names(msn.data))) msn.data$Monoisotopic.Mass <- NA_real_

  # --- Filtering and Appending ---
  unique_in_msn <- setdiff(msn.data$feature.ID, existing_annotations$feature.ID)
  msn.data <- msn.data %>% dplyr::filter(feature.ID %in% unique_in_msn)

  # Set annotation metadata
  msn.data$annotation.type <- "MSNovelist"
  msn.data$confidence.level <- "3"

  # Standardize types before binding
  existing_annotations$feature.ID <- as.numeric(existing_annotations$feature.ID)
  msn.data$feature.ID <- as.numeric(msn.data$feature.ID)

  existing_annotations$confidence.score <- as.numeric(existing_annotations$confidence.score)
  msn.data$confidence.score <- as.numeric(msn.data$confidence.score)

  msn.data$Formula <- as.character(msn.data$Formula)
  msn.data$Monoisotopic.Mass <- as.numeric(msn.data$Monoisotopic.Mass)
  msn.data$IUPAC <- as.character(msn.data$IUPAC)

  existing_annotations$gnps.shared.peaks <- as.numeric(existing_annotations$gnps.shared.peaks)
  msn.data$gnps.shared.peaks <- as.numeric(msn.data$gnps.shared.peaks)

  msn.data$library.name <- as.character(msn.data$library.name)
  msn.data$library.quality <- as.character(msn.data$library.quality)
  msn.data$NPC.pathway <- as.character(msn.data$NPC.pathway)
  msn.data$NPC.superclass <- as.character(msn.data$NPC.superclass)
  msn.data$gnps.library.usi <- as.character(msn.data$gnps.library.usi)
  msn.data$gnps.in.silico.bile.acid.info <- as.character(msn.data$gnps.in.silico.bile.acid.info)

  # Append the filtered and processed data to the existing annotations
  updated_annotations <- existing_annotations %>%
    dplyr::bind_rows(msn.data)

  return(list(annotations = updated_annotations, cache = updated_cid_cache_df))
}
