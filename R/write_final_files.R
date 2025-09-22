#' @title Write Final Processed Files
#' @description This function takes the final processed data frames and writes them to disk in specified locations. It also performs a final data verification check.
#'
#' @param final.annotation.df The main processed annotation data frame.
#' @param samples.df The processed samples data frame.
#' @param folder The base folder path for file output.
#' @param dataset.id The unique identifier for the dataset.
#' @param mzmine.data The initial MzMine data frame for verification.
#' @param cid_cache_df The CID cache data frame to be saved.
#' @param write_large_csv A function to write large CSV files.
#' @return A message indicating success or an error.
#' @importFrom dplyr %>% filter select full_join
#' @importFrom readr write_csv
#' @export
write_final_files <- function(
    final.annotation.df,
    samples.df,
    folder,
    dataset.id,
    mzmine.data,
    cid_cache_df,
    write_large_csv
) {
  #Step 1: Ensure dataset and other.data are correctly separated and combined
  dataset <- final.annotation.df %>% dplyr::filter(!is.na(smiles) & smiles != "N/A")
  other_data <- final.annotation.df %>% dplyr::filter(is.na(smiles) | smiles == "N/A")

  final.annotation.df2 <- dataset %>%
    dplyr::filter(!redundant) %>%
    dplyr::select(-redundant) %>%
    dplyr::bind_rows(other_data)

  write_large_csv(final.annotation.df2, paste0(folder, "/final-annotation-df.csv"))
  write_large_csv(final.annotation.df2, paste0("Y:/MA_BPA_Microbiome/Dataset-Annotations/", dataset.id, ".csv"))

  # Step 2: Write the samples data frame
  output_file <- paste0(dataset.id, "-samples-df.csv")
  write_large_csv(samples.df, output_file)

  # Step 3: Extract and write the top 10 features per sample
  final.annotation.df3 <- final.annotation.df2 %>%
    dplyr::select(feature.ID, mz, compound.name, smiles, confidence.level, Formula, IUPAC, Monoisotopic.Mass)

  samples.df.with.annotations <- samples.df %>%
    dplyr::full_join(final.annotation.df3, by = "feature.ID")

  top_10_features <- samples.df.with.annotations %>%
    dplyr::group_by(samples) %>%
    dplyr::top_n(10, area) %>%
    dplyr::ungroup()

  write_large_csv(top_10_features, paste0(dataset.id, "-top-10-features.csv"))

  # Step 4: Close connections and save cache
  closeAllConnections()
  tryCatch({
    readr::write_csv(cid_cache_df, "~/cid_cache.csv")
    message("[CACHE WRITE] Saved cache to: ", "~/cid_cache.csv")
  }, error = function(e) {
    warning("Failed to save cache: ", e$message)
  })

  # Step 5: Perform data verification check
  if (nrow(final.annotation.df) > nrow(mzmine.data)) {
    error.message <- paste0("!ATTENTION! Possible data mix-up: Check all data has been processed correctly for dataset ", dataset.id, "!")
    readr::write_csv(tibble::tibble(message = error.message), paste0("Y:/MA_BPA_Microbiome/Dataset-Annotations/", Sys.Date(), "_DATA_ERROR_", dataset.id, ".csv"))
    return(error.message)
  } else {
    message <- paste0("Script has finished: No issues detected in processing ", dataset.id, "!")
    return(message)
  }

  if (all(!duplicated(final.annotation.df$feature.ID))) {
    message("All feature.ID are unique ✅")
  } else {
    warning("Duplicates found ⚠️")
  }
}
