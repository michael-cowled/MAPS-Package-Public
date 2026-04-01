#' Standardise Compound Annotations
#'
#' Prioritises a fully vectorized cache lookup before falling back to individual
#' PubChem API lookups and local DB property retrieval.
#'
#' @export
standardise_annotation <- function(data,
                                   name_col = "compound_name",
                                   smiles_col = "smiles",
                                   cid_cache_df = NULL,
                                   lipids.file = NULL,
                                   cid_database_path = NULL,
                                   standardisation = TRUE,
                                   cache.location,
                                   enable_local_db = TRUE,
                                   enable_api = TRUE) {

  # --- Checks ---
  if (is.null(cid_cache_df)) stop("cid_cache_df must be provided.")
  if (!(name_col %in% names(data))) stop("Missing column: ", name_col)
  if (!(smiles_col %in% names(data))) stop("Missing column: ", smiles_col)
  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  if (enable_local_db && (is.null(cid_database_path) || !file.exists(cid_database_path))) {
    warning("Local PubChem database is missing or path is invalid. Forcing enable_local_db = FALSE.")
    enable_local_db <- FALSE
  }

  # --- Initialise ---
  data <- data[!grepl("candidate", data[[name_col]], ignore.case = TRUE), ]
  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  data[[name_col]] <- as.character(data[[name_col]])
  data[[smiles_col]] <- as.character(data[[smiles_col]])

  if (!"CID" %in% names(data)) data$CID <- NA_real_
  data$Formula <- NA_character_
  data$IUPAC <- NA_character_
  data$Monoisotopic.Mass <- NA_real_

  if (standardisation) {
    message("--- PASS 1: Resolving CIDs ---")

    # ======================================================================
    # BULK CACHE LOOKUP (Lightning Fast Vectorized Joins)
    # ======================================================================
    message("Performing bulk cache lookups...")
    initial_nas <- sum(is.na(data$CID))

    # 1. Match by Name
    cache_by_name <- cid_cache_df[!is.na(cid_cache_df$LookupName), c("LookupName", "CID")]
    # Ensure distinct names in cache to avoid row duplication in join
    cache_by_name <- cache_by_name[!duplicated(cache_by_name$LookupName), ]

    data <- data %>%
      dplyr::left_join(cache_by_name, by = stats::setNames("LookupName", name_col)) %>%
      dplyr::mutate(CID = dplyr::coalesce(CID.x, CID.y)) %>%
      dplyr::select(-CID.x, -CID.y)

    # 2. Match by SMILES (for remaining NAs)
    if ("SMILES" %in% names(cid_cache_df)) {
      cache_by_smiles <- cid_cache_df[!is.na(cid_cache_df$SMILES) & cid_cache_df$SMILES != "", c("SMILES", "CID")]
      cache_by_smiles <- cache_by_smiles[!duplicated(cache_by_smiles$SMILES), ]

      data <- data %>%
        dplyr::left_join(cache_by_smiles, by = stats::setNames("SMILES", smiles_col)) %>%
        dplyr::mutate(CID = dplyr::coalesce(CID.x, CID.y)) %>%
        dplyr::select(-CID.x, -CID.y)
    }

    resolved_from_cache <- initial_nas - sum(is.na(data$CID))
    message(sprintf("[CACHE HIT] Bulk resolved %d CIDs directly from cache.", resolved_from_cache))

    # ======================================================================
    # LOOP FOR MISSING ITEMS ONLY (API & Fallbacks)
    # ======================================================================
    missing_idx <- which(is.na(data$CID) | data$CID == "")

    if (length(missing_idx) > 0) {
      message(sprintf("Processing %d unresolved compounds...", length(missing_idx)))
      pb <- utils::txtProgressBar(min = 0, max = length(missing_idx), style = 3)

      for (idx in seq_along(missing_idx)) {
        i <- missing_idx[idx]
        name <- data[[name_col]][i]
        smiles <- data[[smiles_col]][i]

        if (is.na(name) || !nzchar(name)) {
          utils::setTxtProgressBar(pb, idx)
          next
        }

        # Call helper function (it handles lipids.file and API checks)
        pubchem_result <- get_cid_only_with_fallbacks(name, smiles, cid_cache_df, lipids.file, offline = !enable_api)
        data$CID[i] <- pubchem_result$CID
        cid_cache_df <- pubchem_result$cache

        utils::setTxtProgressBar(pb, idx)
      }
      close(pb)

      # Save cache ONCE at the end
      tryCatch({
        readr::write_csv(cid_cache_df, cache.location)
        message("Cache updated and saved successfully.")
      }, error = function(e) {
        warning("Failed to save cache: ", e$message, call. = FALSE)
      })
    } else {
      message("All compounds resolved from cache. Skipping API loop.")
    }

    # ======================================================================
    # --- PASS 2: Local DB (Optional) ---
    # ======================================================================
    if (enable_local_db) {
      message("\n--- PASS 2: Retrieving Properties from Local DB ---")
      db_con <- DBI::dbConnect(RSQLite::SQLite(), cid_database_path)
      cids_to_lookup <- unique(data$CID[!is.na(data$CID) & data$CID > 0])

      if (length(cids_to_lookup) > 0) {
        cid_str <- paste(cids_to_lookup, collapse = ", ")
        query <- sprintf("SELECT CID, Title, SMILES, Formula AS Formula_db, IUPAC AS IUPAC_db, `Monoisotopic.Mass` AS Monoisotopic_Mass_db FROM pubchem_data WHERE CID IN (%s) GROUP BY CID", cid_str)

        db_props <- tryCatch(DBI::dbGetQuery(db_con, query), error = function(e) stop(e$message, call. = FALSE))

        if (!is.null(db_props) && nrow(db_props) > 0) {
          if ("CID" %in% colnames(db_props)) db_props$CID <- as.numeric(db_props$CID)
          data <- data %>%
            dplyr::left_join(db_props, by = "CID") %>%
            dplyr::mutate(
              !!rlang::sym(name_col) := dplyr::coalesce(Title, IUPAC_db, !!rlang::sym(name_col)),
              !!rlang::sym(smiles_col) := dplyr::coalesce(SMILES, !!rlang::sym(smiles_col)),
              Formula = dplyr::coalesce(Formula_db, Formula),
              IUPAC = dplyr::coalesce(IUPAC_db, IUPAC),
              Monoisotopic.Mass = dplyr::coalesce(Monoisotopic_Mass_db, Monoisotopic.Mass)
            ) %>%
            dplyr::select(-Title, -SMILES, -Formula_db, -IUPAC_db, -Monoisotopic_Mass_db)
        }
      }
      DBI::dbDisconnect(db_con)
    }

    # ======================================================================
    # --- PASS 3: API Enrichment (Optional) ---
    # ======================================================================
    if (enable_api) {
      message("\n--- PASS 3: PubChem API Lookup ---")
      if (requireNamespace("jsonlite", quietly = TRUE)) {
        data <- update_compound_names(data, name_col = name_col, cid_col = "CID")
      }
    }
  }
  return(list(data = data, cache = cid_cache_df))
}
