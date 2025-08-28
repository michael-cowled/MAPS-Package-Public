#' Standardise Compound Annotations with Caching
#'
#' This function acts as a wrapper for `standardise_annotation`, providing a
#' convenient interface to standardise compound names and SMILES strings
#' in a data frame. It resolves CIDs using a local cache and PubChem,
#' then updates the data with properties from a local SQLite database.
#'
#' @param df A data frame containing compound annotations.
#' @param name_col Name of the column in `df` containing compound names (as a string).
#' @param smiles_col Name of the column in `df` containing SMILES strings (as a string).
#' @param cid_cache_df A data frame to use as a cache for PubChem CID lookups. This cache will be updated.
#' @param lipids.file A data frame to use to lookup lipid names for PubChem CID lookups.
#' @param db_path Path to the "SQLite database" file containing PubChem data.
#'
#' @return A list with two elements: `data` (the updated data frame) and `cache` (the updated CID cache).
#' @export
#'
standardise_with_cache <- function(df, name_col, smiles_col, cid_cache_df, lipids.file, db_path) {
  # Call the primary standardisation function
  result <- standardise_annotation(
    data = df,
    name_col = name_col,
    smiles_col = smiles_col,
    cid_cache_df = cid_cache_df,
    lipids.file = lipids.file,
    cid_database_path = db_path
  )

  # Extract the updated data frame and cache from the result
  df <- result$data
  cid_cache_df <- result$cache

  # Ensure the expected output columns exist, even if standardise_annotation already handles this.
  for (col in c("Formula", "IUPAC", "Monoisotopic.Mass")) {
    if (!col %in% names(df)) {
      df[[col]] <- NA
    }
  }

  # Return the list with the updated data frame and cache
  return(list(data = df, cache = cid_cache_df))
}
