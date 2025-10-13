#' @title Standardize and Compute ID Probability for all Level 2 Annotations
#' @description Binds different Level 2 annotation data frames and applies a series
#'     of standardization and deduplication steps.
#'
#' @param gnps.data.lv2.low.conf A data frame of low confidence GNPS annotations.
#' @param ms2query.data.lv2 A data frame of Level 2 MS2Query annotations.
#' @param cid_cache_df A data frame containing the CID cache.
#' @param lipids.file A data frame of lipid information.
#' @param cid_database.path A character string. Path to the PubChem SQLite database.
#' @param ms2query.prob A numeric value. The probability threshold for MS2Query.
#'
#' @return A list with two elements: `data` (the processed Level 2 data frame)
#'     and `cache` (the updated CID cache data frame).
#'
#' @importFrom dplyr %>%
#'
#' @export
standardize_and_compute_all_lv2 <- function(
    gnps.data.lv2.low.conf,
    ms2query.data.lv2,
    cid_cache_df,
    lipids.file,
    cid_database.path,
    ms2query.prob) {

  lv2.annotations <- rbind(gnps.data.lv2.low.conf, ms2query.data.lv2)

  lv2.annotations$smiles <- trimws(lv2.annotations$smiles)

  lv2.annotations$CID <- NA
  lv2.annotations <- deduplicate_data(lv2.annotations, compound.name, confidence.score)

  # Reuse the `standardise_annotation` function
  result <- standardise_annotation(
    lv2.annotations,
    name_col = "compound.name",
    smiles_col = "smiles",
    cid_cache_df = cid_cache_df,
    lipids.file = lipids.file,
    cid_database_path = cid_database.path
  )

  lv2.annotations <- result$data
  cid_cache_df <- result$cache

  # Add new columns if they don't exist
  if (!("Formula" %in% names(lv2.annotations))) lv2.annotations$Formula <- NA_character_
  if (!("IUPAC" %in% names(lv2.annotations))) lv2.annotations$IUPAC <- NA_character_
  if (!("Monoisotopic.Mass" %in% names(lv2.annotations))) lv2.annotations$Monoisotopic.Mass <- NA_real_

  lv2.annotations <- deduplicate_data(lv2.annotations, compound.name, confidence.score)

  # Reuse the `compute_id_prob` function
  lv2.annotations <- compute_id_prob(lv2.annotations, "confidence.score", ms2query.prob)

  return(list(data = lv2.annotations, cache = cid_cache_df))
}
