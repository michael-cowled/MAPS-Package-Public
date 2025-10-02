#' @title Process MZMINE Data
#' @description A function to load, process, and enrich MZMINE annotation data.
#'
#' @param mzmine.annotations Path to the MZMINE annotations CSV file.
#' @param mzmine.data Path to the MZMINE data CSV file.
#' @param gnps.prob A numeric value for the GNPS confidence probability threshold.
#' @param cid.cache.path Path to the CID cache CSV file. Defaults to `~/cid_cache.csv`.
#' @param lipids.file.path Path to the lipids TSV file. Defaults to a specific network path.
#' @param cid.database.path Path to the PubChem CID SQLite database. Defaults to a specific network path.
#'
#' @return A list containing four data frames: `processed.data`, `lipids.data`, `annotations.data`, and `initial.mzmine.data`.
#'
#' @importFrom dplyr filter distinct group_by ungroup slice mutate select left_join
#' @importFrom readr read_csv write_csv
#' @importFrom readxl read_excel
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
#' @importFrom stringr str_replace
#' @importFrom tibble tibble
#' @importFrom tidyr replace_na
#'
#' @export
#'
process_mzmine_data <- function(mzmine.annotations, mzmine.data, gnps.prob,
                                cid.cache.path = "~/MAPS/cid_cache.csv",
                                lipids.file.path = "~/MAPS/lipids_expanded.tsv",
                                cid.database.path = "~/MAPS/PubChem_Indexed.sqlite") {

  # The code from your original function, with explicit package calls
  mzmine.annotations <- read_checked_csv(mzmine.annotations)
  mzmine.annotations$rt <- as.numeric(mzmine.annotations$rt)

  lipids <- mzmine.annotations %>%
    dplyr::filter(method == "matchedlipid", score > 0.7)

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

  # 3. Load database connections and cache
  # Added backup path check for the CID database
  if (!file.exists(cid.database.path)) {
    warning("CID database not found at local path. Using network path as a fallback.")
    cid.database.path <- "Y:/MA_BPA_Microbiome/Databases/PubChem/PubChem_Indexed.sqlite"
  }
  cid_db_con <- DBI::dbConnect(RSQLite::SQLite(), cid.database.path, flags = RSQLite::SQLITE_RO)

  cid_cache_df <- tryCatch({
    read_checked_csv(cid.cache.path)
  }, error = function(e) {
    message("[CACHE INIT] No cache file found. Initializing empty cache.")
    data.frame(LookupName = character(), SMILES = character(), CID = numeric(), stringsAsFactors = FALSE)
  })

  # Implement the file existence check here
  if (!file.exists(lipids.file.path)) {
    warning("Lipids file not found at local path. Using network path as a fallback.")
    lipids.file.path <- "Y:/MA_BPA_Microbiome/Databases/LipidMaps/lipids_expanded.tsv"
  }

  lipids.file <- read_checked_tsv(lipids.file.path)
  lipids.file <- lipids.file %>%
    dplyr::mutate(CID = as.numeric(CID))

  # 4. Standardise compound names
  mzmine.annotations$smiles <- trimws(mzmine.annotations$smiles)
  result <- standardise_annotation(
    mzmine.annotations,
    name_col = "compound_name",
    smiles_col = "smiles",
    cid_cache_df = cid_cache_df,
    lipids.file = lipids.file,
    cid_database_path = cid.database.path
  )

  mzmine.annotations <- result$data
  cid_cache_df <- result$cache

  # 5. Save cache
  tryCatch({
    readr::write_csv(cid_cache_df, cid.cache.path)
    message("[CACHE WRITE] Saved cache to: ", cid.cache.path)
  }, error = function(e) {
    warning("Failed to save cache: ", e$message)
  })

  # 6. Handle Level 1 annotations and compute ID probability
  if (nrow(mzmine.annotations) == 0) {
    mzmine.annotations$CID <- NA
  }

  mzmine.annotations <- mzmine.annotations %>%
    dplyr::rename(feature.ID = id)

  mzmine.annotations <- compute_id_prob(mzmine.annotations, "score", gnps.prob) %>%
    dplyr::select(feature.ID, compound_name, score, smiles, id.prob, CID, Formula, IUPAC, Monoisotopic.Mass)

  names(mzmine.annotations.final) <- c('feature.ID', "compound.name", "confidence.score",
                                       "smiles", "id.prob", "CID", "Formula", "IUPAC", "Monoisotopic.Mass")
  mzmine.annotations.final$feature.ID <- as.numeric(mzmine.annotations.final$feature.ID)
  mzmine.annotations.final$confidence.level <- "1"
  mzmine.annotations.final$annotation.type <- "authentic standard"

  # 7. Append to all mzmine features
  mzmine.data <- read_checked_csv(mzmine.data)
  sample.data <- mzmine.data

  if (!("spectral_db_matches.compound_name" %in% names(mzmine.data))) {
    mzmine.data$spectral_db_matches.compound_name <- NA
  }

  initial.mzmine.data <- mzmine.data

  mzmine.data <- mzmine.data %>%
    dplyr::select("id", "rt", "mz", "ion_identities.iin_id") %>%
    dplyr::mutate(dplyr::across(everything(), as.character)) %>%
    dplyr::mutate(dplyr::across(everything(), ~tidyr::replace_na(., "")))

  names(mzmine.data) <- c('feature.ID', "rt", "mz", "ion.identity.ID")
  mzmine.data$feature.ID <- as.numeric(mzmine.data$feature.ID)

  mzmine.data <- mzmine.data %>%
    dplyr::left_join(mzmine.annotations.final, by = "feature.ID")

  return(list(
    processed.data = mzmine.data,
    lipids.data = lipids,
    annotations.data = mzmine.annotations,
    initial.mzmine.data = initial.mzmine.data,
    cid.cache = cid_cache_df,
    lipids.file = lipids.file
  ))
}
