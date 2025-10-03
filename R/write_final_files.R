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
    dplyr::bind_rows(other_data) %>%
    dplyr::select(feature.ID, rt, mz, compound.name, smiles, annotation.type,
      confidence.level, confidence.score, id.prob, CID, HMDB.ID, Formula, IUPAC,
      Monoisotopic.Mass, mz.diff.ppm, feature.usi, gnps.library.usi,
      gnps.cluster.ID, gnps.in.silico.bile.acid.info,
      canopus.NPC.pathway, canopus.NPC.pathway.probability,
      canopus.NPC.superclass, canopus.NPC.superclass.probability,
      zodiac.formula, zodiac.confidence.score,
      Propagated.Feature.ID, Propagated.Annotation.Type,
      Propagated.Annotation.Class, Samples
    )

  # Step 2: Extract and write the top 10 features per sample
  top_10_features <- samples.df %>%
    dplyr::group_by(samples) %>%
    dplyr::top_n(10, area) %>%
    dplyr::ungroup()

  # Step 3: Writing all files
  if (Sys.getenv("USER_DOMAIN") == "unimelb") {
    write_large_csv(final.annotation.df2, paste0(folder, "/final-annotation-df.csv"))
    write_large_csv(final.annotation.df2, paste0("Y:/MA_BPA_Microbiome/Dataset-Annotations/", dataset.id, ".csv"))
    write_large_csv(samples.df, paste0("Y:/MA_BPA_Microbiome/Dataset-Abundances/", dataset.id, "-samples-df.csv"))
    write_large_csv(top_10_features, paste0("Y:/MA_BPA_Microbiome/Dataset-Abundances/", dataset.id, "-top-10-features.csv"))
  } else {
    write_large_csv(final.annotation.df2, paste0(folder, "/", dataset.id, ".csv"))
    write_large_csv(samples.df, paste0(folder, "/",dataset.id, "-samples-df.csv"))
    write_large_csv(top_10_features, paste0(folder, "/", dataset.id, "-top-10-features.csv"))
  }

  # Step 4: Close connections and save cache
  tryCatch({
    readr::write_csv(cid_cache_df, "~/MAPS/cid_cache.csv")
    message("[CACHE WRITE] Saved cache to: ", "~/MAPS/cid_cache.csv")
  }, error = function(e) {
    warning("Failed to save cache: ", e$message)
  })

  closeAllConnections()

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
