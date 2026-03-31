#' Standardise Compound Annotations (Multi-Scenario)
#'
#' For each row in a data frame, attempts to resolve compound CIDs.
#' Features explicit toggles for local SQLite database lookups (Pass 2)
#' and live PubChem API enrichment (Pass 3) to support various network scenarios.
#'
#' @param data A data frame containing compound annotations.
#' @param name_col Name of the column in `data` containing compound names.
#' @param smiles_col Name of the column in `data` containing SMILES strings.
#' @param cid_cache_df A data frame to use as a cache for CID lookups.
#' @param lipids.file A data frame to use to lookup lipid names for CID lookups.
#' @param cid_database_path Path to the local "SQLite database" file.
#' @param standardisation Logical; whether to run the standardisation passes.
#' @param cache.location Path to save the updated cache.
#' @param enable_local_db Logical; if TRUE, executes Pass 2 (SQLite DB lookup).
#' @param enable_api Logical; if TRUE, allows API calls in Pass 1 and executes Pass 3.
#'
#' @return A list with two elements: `data` (the updated data frame) and `cache` (the updated CID cache).
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

  # Verify local DB status if requested
  if (enable_local_db) {
    if (is.null(cid_database_path) || !file.exists(cid_database_path)) {
      warning("Local PubChem database is missing or path is invalid. Forcing enable_local_db = FALSE.")
      enable_local_db <- FALSE
    }
  }

  # --- Filter and Initialise ---
  data <- data[!grepl("candidate", data[[name_col]], ignore.case = TRUE), ]
  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  data[[name_col]] <- as.character(data[[name_col]])
  data[[smiles_col]] <- as.character(data[[smiles_col]])

  if (!"CID" %in% names(data)) data$CID <- NA_real_
  data$Formula <- NA_character_
  data$IUPAC <- NA_character_
  data$Monoisotopic.Mass <- NA_real_

  # --- DB Connection (Conditional) ---
  if (enable_local_db) {
    message("[DB CONNECT] Connecting to local CID SQLite DB...")
    db_con <- DBI::dbConnect(RSQLite::SQLite(), cid_database_path)
    on.exit({
      if (DBI::dbIsValid(db_con)) {
        DBI::dbDisconnect(db_con)
        message("[DB DISCONNECT] Closed DB connection.")
      }
    }, add = TRUE)
  }

  if (standardisation) {
    # ======================================================================
    # --- PASS 1: Resolve CIDs ---
    # ======================================================================
    message("--- PASS 1: Resolving CIDs ---")
    pb <- utils::txtProgressBar(min = 0, max = nrow(data), style = 3)

    for (i in seq_len(nrow(data))) {
      name <- data[[name_col]][i]
      smiles <- data[[smiles_col]][i]

      if (!is.na(data$CID[i]) && data$CID[i] != "") {
        utils::setTxtProgressBar(pb, i)
        next
      }
      if (is.na(name) || !nzchar(name)) {
        utils::setTxtProgressBar(pb, i)
        next
      }

      # If enable_api is FALSE, we enforce offline mode in the helper function
      pubchem_result <- get_cid_only_with_fallbacks(name, smiles, cid_cache_df, lipids.file, offline = !enable_api)
      data$CID[i] <- pubchem_result$CID
      cid_cache_df <- pubchem_result$cache

      tryCatch({
        readr::write_csv(cid_cache_df, cache.location)
      }, error = function(e) {
        warning("Failed to save cache: ", e$message, call. = FALSE)
      })

      utils::setTxtProgressBar(pb, i)
    }
    close(pb)

    # ======================================================================
    # --- PASS 2: Retrieve Properties from Local DB ---
    # ======================================================================
    if (enable_local_db) {
      message("\n--- PASS 2: Retrieving Properties from Local DB ---")
      cids_to_lookup <- unique(data$CID[!is.na(data$CID) & data$CID > 0])

      if (length(cids_to_lookup) > 0) {
        cid_str <- paste(cids_to_lookup, collapse = ", ")
        query <- sprintf("SELECT CID, Title, SMILES,
                                 Formula AS Formula_db, IUPAC AS IUPAC_db,
                                 `Monoisotopic.Mass` AS Monoisotopic_Mass_db
                          FROM pubchem_data WHERE CID IN (%s) GROUP BY CID", cid_str)

        db_props <- tryCatch(
          DBI::dbGetQuery(db_con, query),
          error = function(e) {
            stop(sprintf("\n[FATAL DB ERROR] %s\nPlease check your connection to: %s", e$message, cid_database_path), call. = FALSE)
          }
        )

        if (!is.null(db_props) && nrow(db_props) > 0) {
          if ("CID" %in% colnames(db_props)) db_props$CID <- as.numeric(db_props$CID)

          message("[DB LOOKUP] Retrieved ", nrow(db_props), " rows")

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
        } else {
          message("[DB LOOKUP] No rows returned.")
        }
      }
    } else {
      message("\n--- PASS 2 SKIPPED: Local DB disabled. ---")
    }

    # ======================================================================
    # --- PASS 3: PubChem Lookup Integration (API) ---
    # ======================================================================
    if (enable_api) {
      message("\n--- PASS 3: PubChem API Lookup ---")
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        warning("PubChem lookup is enabled, but 'jsonlite' is not installed. Skipping.")
      } else {
        message("Starting live PubChem lookup...")
        data <- update_compound_names(data, name_col = name_col, cid_col = "CID")
        message("PubChem name enrichment complete.")
      }
    } else {
      message("\n--- PASS 3 SKIPPED: API disabled (Offline Mode). ---")
    }

  } else {
    message("Skipping standardisation entirely.")
  }

  return(list(data = data, cache = cid_cache_df))
}
