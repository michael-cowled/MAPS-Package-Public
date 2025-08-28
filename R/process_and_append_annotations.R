#' Process Raw Annotations
#'
#' A comprehensive function to process a raw annotation data frame. It performs
#' pre-standardization deduplication, calls the standardization function,
#' performs post-standardization deduplication, and computes the identification
#' probability.
#'
#' @param raw_df The raw data frame of annotations.
#' @param score_col A string for the score column name.
#' @param threshold A numeric value for the confidence threshold.
#' @param name_col A string for the compound name column.
#' @param smiles_col A string for the SMILES column.
#' @param cid_cache_df A data frame used for CID caching.
#' @param lipids.file A data frame with lipid data for standardization.
#' @param db_path A string for the path to the PubChem database.
#'
#' @return A list with the processed data frame (`data`) and the updated
#'   CID cache (`cache`).
#' @export
#' @importFrom rlang sym !!
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#'   # Example usage assuming a dataset and all required files exist
#'   # processed_data <- process_annotations(
#'   #   raw_df = my_raw_data,
#'   #   score_col = "confidence.score",
#'   #   threshold = 0.7,
#'   #   name_col = "compound.name",
#'   #   smiles_col = "smiles",
#'   #   cid_cache_df = my_cache,
#'   #   lipids.file = my_lipids,
#'   #   db_path = "path/to/db"
#'   # )
#' }
process_annotations <- function(df, name_col, smiles_col, score_col, threshold,
                                cid_cache_df, lipids.file, db_path,
                                confidence.level=NULL, annotation.type=NULL) {
  # existing standardisation/deduplication/ID prob code...
  result <- standardise_with_cache(df, name_col, smiles_col, cid_cache_df, lipids.file, db_path)$data
  result <- deduplicate_data(result, name_col, score_col)
  result <- compute_id_prob(result, score_col, threshold)

  if(!is.null(confidence.level)) result$confidence.level <- confidence.level
  if(!is.null(annotation.type)) result$annotation.type <- annotation.type

  return(result)
}
