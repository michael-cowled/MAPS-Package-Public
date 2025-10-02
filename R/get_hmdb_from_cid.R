#' Retrieve HMDB IDs based on CID
#'
#' Takes a data frame containing a 'CID' column and queries the local PubChem
#' SQLite database to retrieve the corresponding HMDB.ID for each CID.
#' This is designed as a fallback or post-processing step if needed.
#'
#' @param final.annotation.df A data frame containing a column named 'CID'.
#' @param cid_database_path Path to the PubChem CID SQLite database (same as used in standardise_annotation).
#'
#' @return The input data frame with a new column 'HMDB.ID' added and populated.
#'
#' @importFrom dplyr select left_join mutate distinct any_of rename coalesce
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom RSQLite SQLite
#' @export
get_hmdb_from_cid <- function(final.annotation.df, cid_database_path) {

  # --- Checks ---
  if (!("CID" %in% names(final.annotation.df))) {
    stop("Input data frame must contain a 'CID' column.")
  }

  # --- 1. DB Connection and HMDB ID Lookup ---
  if (is.null(cid_database_path) || !file.exists(cid_database_path)) {
    warning("Local PubChem database is missing or path is invalid. Returning input data frame without HMDB IDs.")
    # Ensure HMDB.ID column exists as NA_character_ if no lookup was done
    if (!("HMDB.ID" %in% names(final.annotation.df))) {
      final.annotation.df$HMDB.ID <- NA_character_
    }
    return(final.annotation.df)
  }

  message("[DB CONNECT] Connecting to CID SQLite DB for HMDB lookup...")
  db_con <- DBI::dbConnect(RSQLite::SQLite(), cid_database_path, flags = RSQLite::SQLITE_RO)
  # Ensure connection is closed upon function exit
  on.exit(DBI::dbDisconnect(db_con), add = TRUE)

  # Filter for valid CIDs to query
  cids_to_lookup <- unique(final.annotation.df$CID[!is.na(final.annotation.df$CID) & final.annotation.df$CID > 0])

  if (length(cids_to_lookup) == 0) {
    message("No valid CIDs found to look up HMDB IDs.")
    if (!("HMDB.ID" %in% names(final.annotation.df))) {
      final.annotation.df$HMDB.ID <- NA_character_
    }
    return(final.annotation.df)
  }

  cid_str <- paste(cids_to_lookup, collapse = ", ")

  # Use an alias without a period (HMDB_ID) to prevent the SQL syntax error.
  # The column name in the DB is correctly backticked (`Primary.HMDB.ID`).
  query <- sprintf("SELECT CID, `Primary.HMDB.ID` AS HMDB_ID FROM pubchem_data WHERE CID IN (%s) GROUP BY CID", cid_str)

  message("[SQL QUERY] Retrieving HMDB IDs...")

  hmdb_map <- tryCatch(
    DBI::dbGetQuery(db_con, query) %>%
      # Ensure CID is numeric for the join
      dplyr::mutate(CID = as.numeric(CID)),
    error = function(e) {
      message("[DB ERROR] ", e$message)
      return(NULL)
    }
  )

  if (is.null(hmdb_map) || nrow(hmdb_map) == 0) {
    warning("HMDB lookup failed or returned no results. Returning input data frame.")
    if (!("HMDB.ID" %in% names(final.annotation.df))) {
      final.annotation.df$HMDB.ID <- NA_character_
    }
    return(final.annotation.df)
  }

  # Ensure the HMDB_ID column (the alias) is character
  hmdb_map$HMDB_ID <- as.character(hmdb_map$HMDB_ID)

  # --- 2. Perform the Lookup/Join ---

  # 1. Prepare the input data frame for joining and coalescing.
  # Rename the existing HMDB.ID column to a temporary name (HMDB_ID_existing)
  # This avoids collision with the new HMDB_ID column coming from the join.
  if ("HMDB.ID" %in% names(final.annotation.df)) {
    final.annotation.df <- final.annotation.df %>%
      dplyr::rename(HMDB_ID_existing = HMDB.ID)
  } else {
    # If it doesn't exist, create a placeholder NA column with the same temporary name
    final.annotation.df$HMDB_ID_existing <- NA_character_
  }

  # 2. Join the map (which contains CID and HMDB_ID) to the annotations
  final.annotation.df <- final.annotation.df %>%
    dplyr::left_join(hmdb_map, by = "CID") %>%

    # 3. Coalesce: Create the final HMDB.ID column (with period)
    # Coalesce: new value (HMDB_ID) | existing value (HMDB_ID_existing)
    dplyr::mutate(
      HMDB.ID = dplyr::coalesce(HMDB_ID, HMDB_ID_existing)
    ) %>%

    # 4. Clean up the temporary columns used for joining/coalescing
    dplyr::select(-dplyr::any_of(c("HMDB_ID", "HMDB_ID_existing")))

  message(paste("Successfully added HMDB IDs. Total unique HMDB IDs found:",
                sum(!is.na(final.annotation.df$HMDB.ID))))

  return(final.annotation.df)
}
