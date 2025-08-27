#' Process and Append Annotations to a Master Data Frame
#'
#' A comprehensive function to handle the workflow of processing a new
#' annotation data source and appending it to a master annotation table.
#' It performs deduplication, standardization, probability calculation, and column
#' alignment before merging.
#'
#' @param df The data frame with new annotations to be processed.
#' @param master_df The master data frame to which the new annotations will be appended.
#' @param name_col A string specifying the column with compound names.
#' @param smiles_col A string specifying the column with SMILES strings.
#' @param score_col A string specifying the column with confidence scores.
#' @param prob_threshold A numeric value for the identification probability cutoff.
#' @param annotation_type A string describing the annotation source (e.g., "gnps", "csi").
#' @param confidence_level A string or number for the confidence level (e.g., "2", "3").
#' @param cid_cache_df A data frame used as a cache for CID lookups.
#' @param lipids.file A data frame containing LipidMaps data for standardization.
#' @param cid_database_path A string specifying the path to the PubChem SQLite database.
#'
#' @return A list containing the updated master data frame ('master_df') and the
#'   updated CID cache data frame ('cache').
#' @importFrom dplyr bind_rows
#' @importFrom rlang `[[`
#' @export
process_and_append_annotations <- function(
    df,
    master_df,
    name_col,
    smiles_col,
    score_col,
    prob_threshold,
    annotation_type,
    confidence_level,
    cid_cache_df,
    lipids.file,
    cid_database_path) {

  # Standardize and deduplicate pre- and post-standardization
  df <- deduplicate_by_score(df, id_col = feature.ID, score_col = score_col)

  result <- standardise_annotation(
    df,
    name_col = name_col,
    smiles_col = smiles_col,
    cid_cache_df = cid_cache_df,
    lipids.file = lipids.file,
    cid_database_path = cid_database_path
  )
  df <- result$data
  cid_cache_df <- result$cache

  df <- deduplicate_by_score(df, id_col = feature.ID, score_col = score_col)

  # Calculate ID probability
  df <- calculate_id_prob(df, id_col = feature.ID, score_col = score_col, prob_threshold = prob_threshold)

  # Add annotation metadata
  df$annotation.type <- annotation_type
  df$confidence.level <- confidence_level

  # Check and add missing columns before binding
  missing_cols <- setdiff(names(master_df), names(df))
  df[missing_cols] <- NA

  # Append to master dataframe
  master_df <- bind_rows(master_df, df)

  return(list(master_df = master_df, cache = cid_cache_df))
}
