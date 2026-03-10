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
    cid_cache_df,
    lipids.file,
    cid_database_path,
    compute_id_prob,
) {
  # 1. Load and Initial Clean
  msn.df <- read_checked_tsv(msn.data)

  msn.df <- msn.df %>%
    dplyr::select(
      confidence.score = ModelScore,
      compound.name = name,
      smiles = smiles,
      feature.ID = mappingFeatureId
    ) %>%
    dplyr::mutate(confidence.score = as.numeric(confidence.score))

  # 2. Filter out features already successfully annotated in Levels 1-3
  msn.df <- msn.df[!(msn.df$feature.ID %in% existing_annotations$feature.ID), ]

  # 3. Sort by Score (Closest to Zero is Best)
  # Since ModelScores are negative, "descending" puts -1.5 above -10.0
  msn.df <- msn.df %>%
    dplyr::arrange(feature.ID, desc(confidence.score))

  # 4. Preparation for MAPS Pipeline
  missing_cols <- setdiff(colnames(existing_annotations), colnames(msn.df))
  for (col in missing_cols) msn.df[[col]] <- NA

  # 5. ID Probability (We pass a dummy threshold of -Inf since you want all hits)
  msn.df <- compute_id_prob(msn.df, "confidence.score", threshold = -Inf)

  final_msn <- msn.df

  # 8. Metadata and Type Correction
  final_msn$annotation.type <- "MSNovelist"
  final_msn$confidence.level <- "3" # Per your previous setting

  # Cast columns to ensure bind_rows doesn't crash on type mismatch
  final_msn <- final_msn %>%
    dplyr::mutate(
      feature.ID = as.numeric(feature.ID),
      confidence.score = as.numeric(confidence.score),
      Formula = as.character(Formula),
      IUPAC = as.character(IUPAC),
      Monoisotopic.Mass = as.numeric(Monoisotopic.Mass)
    )

  # 9. Final Merge
  updated_annotations <- existing_annotations %>%
    dplyr::mutate(feature.ID = as.numeric(feature.ID)) %>%
    dplyr::bind_rows(final_msn)

  return(list(annotations = updated_annotations))
}
