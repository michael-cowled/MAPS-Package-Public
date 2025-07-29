#' Standardise Compound Annotations Using PubChem
#'
#' For each row in a data frame, attempts to resolve compound annotations using PubChem.
#' Updates the compound name, SMILES, and adds a column for CID.
#' Rows containing "candidate" (case-insensitive) in the name column are excluded.
#'
#' @param data A data frame containing compound annotations.
#' @param name_col Name of the column in `data` containing compound names (as a string).
#' @param smiles_col Name of the column in `data` containing SMILES strings (as a string).
#' @param cid_database_path Path to the "SQLite database" file.
#' @param cid_cache_path Path to the "cid_cache.csv" file for loading/saving the cache.
#'
#' @return A data frame with standardised annotations and an additional PubChem info column:
#'   \itemize{
#'     \item \code{CID}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Example data
#' # data_to_process <- data.frame(
#' #   compound.name = c("glucose", "aspirin", "unknown compound", "caffeine", "candidate drug"),
#' #   smiles = c("C(C1C(C(C(C(O1)CO)O)O)O)O", "CC(=O)Oc1ccccc1C(=O)O", NA, "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", NA),
#' #   stringsAsFactors = FALSE
#' # )
#'
#' # df <- standardise_annotation(data_to_process, "compound.name", "smiles",
#' #                              "path/to/SQLite database",
#' #                              "path/to/cid_cache.csv")
#' }
standardise_annotation <- function(data,
                                   name_col = "compound_name",
                                   smiles_col = "smiles",
                                   cid_cache_df = NULL,
                                   cid_database_path = NULL) {

  library(dplyr)
  library(DBI)
  library(RSQLite)

  # --- Checks ---
  if (is.null(cid_cache_df)) stop("cid_cache_df must be provided.")
  if (!(name_col %in% names(data))) stop("Missing column: ", name_col)
  if (!(smiles_col %in% names(data))) stop("Missing column: ", smiles_col)
  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  # --- Filter out "candidate" names ---
  data <- data %>%
    filter(!grepl("candidate", .data[[name_col]], ignore.case = TRUE))

  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  data[[name_col]] <- as.character(data[[name_col]])
  data[[smiles_col]] <- as.character(data[[smiles_col]])
  data$CID <- NA_real_

  # --- Optional DB Connection ---
  db_con <- NULL
  if (!is.null(cid_database_path) && file.exists(cid_database_path)) {
    message("[DB CONNECT] Connecting to CID SQLite DB...")
    db_con <- dbConnect(SQLite(), cid_database_path)
    on.exit(if (!is.null(db_con)) dbDisconnect(db_con), add = TRUE)
  } else {
    message("[DB CONNECT] Database not found or path is NULL.")
  }

  # --- Internal function: safely query PubChem ---
  safe_get_pubchem <- function(name, smiles) {
    out <- tryCatch({
      get_cid_with_fallbacks(name, smiles)  # your own function
    }, error = function(e) {
      message("  [PUBCHEM ERROR] ", e$message)
      NULL
    })
    if (is.null(out)) {
      list(CID = -1, ResolvedName = NA_character_, SMILES = NA_character_)
    } else {
      out
    }
  }

  # --- Main Loop ---
  pb <- txtProgressBar(min = 0, max = nrow(data), style = 3)

  for (i in seq_len(nrow(data))) {
    name <- data[[name_col]][i]
    smiles <- data[[smiles_col]][i]

    if (is.na(name) || !nzchar(name)) {
      setTxtProgressBar(pb, i)
      next
    }

    resolved <- list(CID = NA_real_, ResolvedName = NA_character_, SMILES = NA_character_)

    # --- 1. Check Cache ---
    cached <- cid_cache_df %>% filter(LookupName == name & !is.na(ResolvedName) & !is.na(SMILES))
    if (nrow(cached) > 0 && !is.na(cached$CID[1])) {
      resolved <- list(
        CID = cached$CID[1],
        ResolvedName = cached$ResolvedName[1],
        SMILES = cached$SMILES[1]
      )
      message("  [CACHE HIT] CID: ", resolved$CID)
    }

    # --- 2. Check SQLite DB ---
    if (is.na(resolved$CID) && !is.null(db_con)) {
      query <- sprintf("SELECT CID, Title, SMILES FROM cid_data WHERE Title = '%s' COLLATE NOCASE",
                       gsub("'", "''", name))
      db_result <- tryCatch(dbGetQuery(db_con, query), error = function(e) NULL)

      if (!is.null(db_result) && nrow(db_result) > 0) {
        resolved <- list(
          CID = db_result$CID[1],
          ResolvedName = db_result$Title[1],
          SMILES = db_result$SMILES[1]
        )
        message("  [DB LOOKUP] Found CID: ", resolved$CID)
      }
    }

    # --- 3. Try PubChem ---
    if (is.na(resolved$CID)) {
      resolved <- safe_get_pubchem(name, smiles)

      if (!is.na(resolved$CID) && resolved$CID != -1) {
        message("  [PUBCHEM] Found CID: ", resolved$CID)
      } else {
        message("  [PUBCHEM] No CID found.")
      }
    }

    # --- 4. Fallback: fill missing from DB via CID ---
    if (!is.null(db_con) && !is.na(resolved$CID) && resolved$CID != -1 &&
        (is.na(resolved$ResolvedName) || is.na(resolved$SMILES))) {
      fallback_query <- sprintf("SELECT Title, SMILES FROM cid_data WHERE CID = %d", resolved$CID)
      fallback_result <- dbGetQuery(db_con, fallback_query)
      if (nrow(fallback_result) > 0) {
        if (is.na(resolved$ResolvedName)) resolved$ResolvedName <- fallback_result$Title[1]
        if (is.na(resolved$SMILES)) resolved$SMILES <- fallback_result$SMILES[1]
        message("  [DB FILL] Filled missing info from DB.")
      }
    }

    # --- Update Data ---
    data$CID[i] <- resolved$CID
    if (!is.na(resolved$ResolvedName)) data[[name_col]][i] <- resolved$ResolvedName
    if (!is.na(resolved$SMILES)) data[[smiles_col]][i] <- resolved$SMILES

    # --- Update Cache ---
    if (!name %in% cid_cache_df$LookupName) {
      cid_cache_df <- bind_rows(
        cid_cache_df,
        data.frame(
          LookupName = name,
          ResolvedName = resolved$ResolvedName,
          SMILES = resolved$SMILES,
          CID = resolved$CID,
          stringsAsFactors = FALSE
        )
      )
      message("  [CACHE ADD] ", name)
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)

  # --- Return ---
  return(list(data = data, cache = cid_cache_df))
}
