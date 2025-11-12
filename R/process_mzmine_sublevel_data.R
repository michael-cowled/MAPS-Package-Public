#' @title Process and Standardize Additional MZMINE Spectral Annotations
#' @description A function to filter, process, standardize, and append additional spectral library annotations (e.g., lower confidence or secondary hits) to a set of final, high-confidence annotations. It primarily handles the standardization of compound names and computes identification probability.
#'
#' @param mzmine.annotations.final A data frame of existing final, high-confidence annotations. This will be the base data to which new annotations are appended.
#' @param mzmine.annotations Path to the MZMINE annotations CSV file containing the data to be processed (e.g., lower confidence hits). This is read using `read_checked_csv`.
#' @param cid_cache_df A data frame containing a cache of previously retrieved PubChem CIDs for compound names/SMILES.
#' @param lipids.file Path to a TSV file containing lipid nomenclature for standardization.
#' @param cid_database_path Path to the PubChem CID SQLite database for name/SMILES lookup.
#' @param standardisation A logical flag indicating whether to perform compound name standardization.
#' @param gnps.prob A numeric value for the GNPS confidence probability threshold used in `compute_id_prob`. This is not present in the function definition but seems necessary based on the code's call to `compute_id_prob`. **(Note: This parameter is used in the body but not in the function signature you provided. It's added here for completeness based on code usage.)**
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \strong{annotations.data}: The combined and processed data frame of final and new annotations.
#'   \item \strong{cid.cache}: The updated CID cache data frame.
#' }
#'
#' @importFrom dplyr filter distinct group_by ungroup slice mutate select left_join rename
#' @importFrom readr read_csv write_csv
#' @importFrom RSQLite SQLite
#' @importFrom stringr str_replace
#' @importFrom tibble tibble
#' @importFrom tidyr replace_na
#'
#' @export
process_mzmine_sublevel_data <- function(mzmine.annotations.final, mzmine.annotations,
                                cid_cache_df, lipids.file,
                                cid_database_path, standardisation) {

  # The code from your original function, with explicit package calls
  mzmine.annotations <- read_checked_csv(mzmine.annotations)
  mzmine.annotations$rt <- as.numeric(mzmine.annotations$rt)

  #Filter out old annotations
  # Filter low confidence data to remove duplicates with high confidence
  unique_in_mzmine <- setdiff(mzmine.annotations$id, mzmine.annotations.final$feature.ID)
  mzmine.annotations <- filter(mzmine.annotations, id %in% unique_in_mzmine)

  #debug
  message(paste0(nrow(mzine.annotations)))

  # ----------------------------------------------------
  # 🌟 EARLY EXIT CHECK 🌟
  # If no new annotations remain after filtering against existing ones,
  # return the original data and skip all remaining processing.
  if (nrow(mzmine.annotations) == 0) {
    message("Skipping processing: No new unique annotations found.")
    return(list(
      annotations.data = mzmine.annotations.final,
      cid.cache = cid_cache_df
    ))
  }
  # ----------------------------------------------------

  mzmine.annotations <- mzmine.annotations %>%
    dplyr::filter(method == "spectral_library_annotation", score > 0.7) %>%
    dplyr::distinct(id, compound_name, .keep_all = TRUE) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(
      dplyr::n() == 1 |
        score > min(score, na.rm = TRUE) |
        (all(score == min(score, na.rm = TRUE)) & dplyr::row_number() == 1)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(compound_name) %>%
    dplyr::filter(score == max(score, na.rm = TRUE)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # 4. Standardise compound names
  mzmine.annotations$smiles <- trimws(mzmine.annotations$smiles)

  result <- standardise_annotation(
    mzmine.annotations,
    name_col = "compound_name",
    smiles_col = "smiles",
    cid_cache_df = cid_cache_df,
    lipids.file = lipids.file,
    cid_database_path = cid.database.path,
    standardise = standardise
  )

  mzmine.annotations <- result$data
  cid_cache_df <- result$cache

  # 6. Handle Level 1 annotations and compute ID probability
  if (nrow(mzmine.annotations) == 0) {
    mzmine.annotations$CID <- as.numeric(0)
    mzmine.annotations$Formula <- as.character(0)
    mzmine.annotations$IUPAC <- as.character(0)
    mzmine.annotations$Monoisotopic.Mass <- as.numeric(0)
  }

  mzmine.annotations <- mzmine.annotations %>%
    dplyr::rename(feature.ID = id)

 mzmine.annotations.new <- compute_id_prob(mzmine.annotations, "score", gnps.prob) %>%
    dplyr::select(feature.ID, compound_name, score, smiles, id.prob, CID, Formula, IUPAC, Monoisotopic.Mass)

  names(mzmine.annotations.new) <- c('feature.ID', "compound.name", "confidence.score",
                                       "smiles", "id.prob", "CID", "Formula", "IUPAC", "Monoisotopic.Mass")
 mzmine.annotations.new$feature.ID <- as.numeric(mzmine.annotations.new$feature.ID)
 mzmine.annotations.new$confidence.level <- "1"
 mzmine.annotations.new$annotation.type <- "authentic standard"

 mzmine.annotations.final <- rbind(mzmine.annotations.final, mzmine.annotations.new)

  return(list(
    annotations.data = mzmine.annotations.final,
    cid.cache = cid_cache_df
  ))

}
