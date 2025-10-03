#' Standardise Compound Annotations Using PubChem
#'
#' For each row in a data frame, attempts to resolve compound CIDs using PubChem.
#' Then, uses the resolved CIDs to retrieve and update compound name (with PubChem Title),
#' SMILES, Formula, IUPAC, and Monoisotopic.Mass from a local SQLite database.
#' Rows containing "candidate" (case-insensitive) in the name column are excluded.
#' NOTE: The function is modified to proceed even if the local SQLite database is missing,
#' by skipping the property retrieval step (Pass 2).
#'
#' @param data A data frame containing compound annotations.
#' @param name_col Name of the column in `data` containing compound names (as a string).
#' @param smiles_col Name of the column in `data` containing SMILES strings (as a string).
#' @param cid_cache_df A data frame to use as a cache for PubChem CID lookups.
#' @param lipids.file A data frame to use to lookup lipid names for PubChem CID lookups.
#' @param cid_database_path Path to the "SQLite database" file.
#'
#' @return A list with two elements: `data` (the updated data frame) and `cache` (the updated CID cache).
#' @export
#'
standardise_annotation <- function(data,
                                   name_col = "compound_name",
                                   smiles_col = "smiles",
                                   cid_cache_df = NULL,
                                   lipids.file = NULL,
                                   cid_database_path = NULL) {

  # --- Libraries ---
  # NOTE: In R packages, use @importFrom in documentation instead of library() inside the function.
  # Leaving library calls here for standalone script execution, but it's bad package practice.
  # library(dplyr)
  # library(DBI)
  # library(RSQLite)
  # library(httr)

  # --- Checks ---
  if (is.null(cid_cache_df)) stop("cid_cache_df must be provided.")
  if (!(name_col %in% names(data))) stop("Missing column: ", name_col)
  if (!(smiles_col %in% names(data))) stop("Missing column: ", smiles_col)
  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  # Check for database path and existence to enable/disable DB lookup
  db_lookup_enabled <- TRUE
  if (is.null(cid_database_path) || !file.exists(cid_database_path)) {
    warning("Local PubChem database is missing or path is invalid. Properties (Formula, IUPAC, etc.) will NOT be retrieved from the local database (Pass 2 skipped).")
    db_lookup_enabled <- FALSE
    pubchem_lookup_enabled <- TRUE
  }

  # --- Filter and Initialise ---
  # Exclude 'candidate' features
  data <- data %>% dplyr::filter(!grepl("candidate", .data[[name_col]], ignore.case = TRUE))
  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  data[[name_col]] <- as.character(data[[name_col]])
  data[[smiles_col]] <- as.character(data[[smiles_col]])

  # Initialize the new columns (must be done regardless of DB status)
  data$CID <- NA_real_
  data$Formula <- NA_character_
  data$IUPAC <- NA_character_
  data$Monoisotopic.Mass <- NA_real_

  # --- DB Connection (Conditional) ---
  if (db_lookup_enabled) {
    message("[DB CONNECT] Connecting to CID SQLite DB...")
    db_con <- DBI::dbConnect(RSQLite::SQLite(), cid_database_path)
    # Ensure connection is closed upon function exit
    on.exit({
      if (DBI::dbIsValid(db_con)) {
        DBI::dbDisconnect(db_con)
        message("[DB DISCONNECT] Closed DB connection.")
      }
    }, add = TRUE)
  }

  # --- PASS 1: Resolve CIDs ---
  message("--- PASS 1: Resolving CIDs from PubChem & LipidMaps ---")
  # Assuming 'get_cid_only_with_fallbacks' is a separate function available in the package
  pb <- utils::txtProgressBar(min = 0, max = nrow(data), style = 3)

  for (i in seq_len(nrow(data))) {
    name <- data[[name_col]][i]
    smiles <- data[[smiles_col]][i]

    if (is.na(name) || !nzchar(name)) {
      utils::setTxtProgressBar(pb, i)
      next
    }

    # This function uses the PubChem PUG-REST API and the local lipid file
    pubchem_result <- get_cid_only_with_fallbacks(name, smiles, cid_cache_df, lipids.file)
    data$CID[i] <- pubchem_result$CID
    cid_cache_df <- pubchem_result$cache

    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  # --- PASS 2: Retrieve Properties from Local DB (Conditional) ---
  message("\n--- PASS 2: Retrieving Properties from Local DB ---")

  if (db_lookup_enabled) {
    cids_to_lookup <- unique(data$CID[!is.na(data$CID) & data$CID > 0])

    if (length(cids_to_lookup) > 0) {
      cid_str <- paste(cids_to_lookup, collapse = ", ")

      # Query the local database for properties
      # but we rename it to 'HMDB_db' for easy, non-conflicting use in R.
      query <- sprintf("SELECT CID, Title, SMILES,
                       Formula AS Formula_db, IUPAC AS IUPAC_db,
                       `Monoisotopic.Mass` AS Monoisotopic_Mass_db FROM pubchem_data WHERE CID IN (%s) GROUP BY CID", cid_str)

      message("[SQL QUERY]")
      message(query)

      db_props <- tryCatch(
        DBI::dbGetQuery(db_con, query),
        error = function(e) {
          message("[DB ERROR] ", e$message)
          return(NULL)
        }
      )

      if (!is.null(db_props) && nrow(db_props) > 0) {
        if ("CID" %in% colnames(db_props)) {
          db_props$CID <- as.numeric(db_props$CID)
        }

        message("[DB LOOKUP] Retrieved ", nrow(db_props), " rows")

        # Join retrieved properties and update/coalesce columns
        data <- data %>%
          dplyr::left_join(db_props, by = "CID") %>%
          dplyr::mutate(
            # Use Title/IUPAC from DB, fallback to original name
            !!rlang::sym(name_col) := dplyr::coalesce(Title, IUPAC_db, !!rlang::sym(name_col)),
            # Use SMILES from DB, fallback to original smiles
            !!rlang::sym(smiles_col) := dplyr::coalesce(SMILES, !!rlang::sym(smiles_col)),
            Formula = dplyr::coalesce(Formula_db, Formula),
            IUPAC = dplyr::coalesce(IUPAC_db, IUPAC),
            Monoisotopic.Mass = dplyr::coalesce(Monoisotopic_Mass_db, Monoisotopic.Mass)
          ) %>%
          # Clean up columns used for joining/coalescing
          dplyr::select(-Title, -SMILES, -Formula_db, -IUPAC_db, -Monoisotopic_Mass_db) # <--- CORRECTED LINE
      } else {
        message("[DB LOOKUP] No rows returned.")
      }
    }

# ======================================================================
# --- PASS 3: PubChem Lookup Integration ---
# ======================================================================

  # Check if the PubChem lookup feature is enabled
  }  else if (pubchem_lookup_enabled) {
    # Check if the required 'jsonlite' package is installed and load it
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      warning("PubChem lookup is enabled, but the 'jsonlite' package is not installed. Skipping lookup.")
    } else {
      library(jsonlite)
      message("Starting PubChem lookup (PASS 3)...")

      # Call the function to update the compound names in the data frame.
      # Ensure your data frame object is passed in and reassigned back.
      # Adjust 'name_col' and 'cid_col' if your columns are named differently.
      data <- update_compound_names(
        data, name_col, cid_col = "CID"            # Column containing the PubChem CID
      )

      message("PubChem name enrichment complete.")
    }
  }
  else {
    message("[DB AND PUBCHEM LOOKUP SKIPPED] Local database not available.")
  }
  return(list(data = data, cache = cid_cache_df))
}
