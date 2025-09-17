process_and_standardise <- function(data_frame, name_col, smiles_col, cid_cache, lipids_file, cid_db_path, confidence_col = "score") {

  # Ensure the data frame has the confidence column
  if (!confidence_col %in% names(data_frame)) {
    stop(paste("Data frame does not contain the required confidence column:", confidence_col))
  }

  # Pre-standardisation deduplication
  deduplicated_data <- deduplicate_data(data_frame, !!as.symbol(name_col), !!as.symbol(confidence_col))

  # Call the existing standardise_annotation function
  result <- standardise_annotation(
    deduplicated_data,
    name_col = name_col,
    smiles_col = smiles_col,
    cid_cache_df = cid_cache,
    lipids.file = lipids_file,
    cid_database_path = cid_db_path
  )

  # Update data and cache
  processed_data <- result$data
  updated_cache <- result$cache

  # Add any missing columns after standardization
  if (!("Formula" %in% names(processed_data))) processed_data$Formula <- NA_character_
  if (!("IUPAC" %in% names(processed_data))) processed_data$IUPAC <- NA_character_
  if (!("Monoisotopic.Mass" %in% names(processed_data))) processed_data$Monoisotopic.Mass <- NA_real_
  if (!("CID" %in% names(processed_data))) processed_data$CID <- NA_integer_

  # Post-standardisation deduplication
  processed_data <- deduplicate_data(processed_data, !!as.symbol(name_col), !!as.symbol(confidence_col))

  return(list(data = processed_data, cache = updated_cache))
}
