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

  if (is.null(cid_cache_df)) {
    stop("cid_cache_df must be provided")
  }

  # Helper function for fallback PubChem lookup (simplified placeholder)
  get_cid_with_fallbacks <- function(name, smiles) {
    # Example placeholder function:
    # Ideally, here you query PubChem API to get CID, ResolvedName, SMILES
    # Return a list with CID, ResolvedName, SMILES or NULL if not found
    # For now, simulate failure to find CID:
    return(NULL)
  }

  # Connect to SQLite DB if path provided
  db_con <- NULL
  if (!is.null(cid_database_path) && file.exists(cid_database_path)) {
    message("[DB CONNECT] Connecting to CID SQLite DB...")
    db_con <- dbConnect(SQLite(), cid_database_path)
  } else {
    message("[DB CONNECT] Database not found or path is NULL.")
  }

  # Return early if no rows
  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  # Check required columns
  if (!(name_col %in% names(data))) stop("Missing column: ", name_col)
  if (!(smiles_col %in% names(data))) stop("Missing column: ", smiles_col)

  # Coerce to character
  data[[name_col]] <- as.character(data[[name_col]])
  data[[smiles_col]] <- as.character(data[[smiles_col]])

  # Filter out 'candidate' rows
  data <- data[!grepl("candidate", data[[name_col]], ignore.case = TRUE), ]
  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  # Initialize CID column
  data$CID <- NA_real_

  pb <- txtProgressBar(min = 0, max = nrow(data), style = 3)

  for (i in seq_len(nrow(data))) {
    name <- data[[name_col]][i]
    smiles <- data[[smiles_col]][i]
    if (is.na(name) || !nzchar(name)) {
      setTxtProgressBar(pb, i)
      next
    }

    resolved_cid <- NA_real_
    resolved_name <- NA_character_
    resolved_smiles <- NA_character_

    # 1. Try cache first
    cached <- cid_cache_df %>% filter(LookupName == name & !is.na(ResolvedName) & !is.na(SMILES))
    if (nrow(cached) > 0 && !is.na(cached$CID[1])) {
      resolved_cid <- cached$CID[1]
      resolved_name <- cached$ResolvedName[1]
      resolved_smiles <- cached$SMILES[1]
      message(paste("  [CACHE HIT] CID:", resolved_cid))
    }

    # 2. Try local DB lookup if not in cache
    if (is.na(resolved_cid) && !is.null(db_con)) {
      query <- sprintf("SELECT CID, Title, SMILES FROM cid_data WHERE Title = '%s' COLLATE NOCASE", gsub("'", "''", name))
      db_result <- tryCatch(dbGetQuery(db_con, query), error = function(e) NULL)
      if (!is.null(db_result) && nrow(db_result) > 0) {
        resolved_cid <- db_result$CID[1]
        resolved_name <- db_result$Title[1]
        resolved_smiles <- db_result$SMILES[1]
        message(paste("  [DB LOOKUP] Found CID in DB:", resolved_cid))
      }
    }

    # 3. Fallback to PubChem (or any other external lookup)
    if (is.na(resolved_cid)) {
      info <- get_cid_with_fallbacks(name, smiles)
      if (!is.null(info) && !is.na(info$CID) && info$CID != -1) {
        resolved_cid <- info$CID
        resolved_name <- info$ResolvedName
        resolved_smiles <- info$SMILES
        message(paste("  [PUBCHEM] Found CID:", resolved_cid))
      } else {
        resolved_cid <- -1
        message("  [PUBCHEM] No CID found.")
      }
    }

    # 4. If CID found but missing name or smiles, try to fill from DB by CID
    if (!is.null(db_con) && !is.na(resolved_cid) && resolved_cid != -1 && (is.na(resolved_name) || is.na(resolved_smiles))) {
      query <- sprintf("SELECT Title, SMILES FROM cid_data WHERE CID = %d", resolved_cid)
      db_fallback <- dbGetQuery(db_con, query)
      if (nrow(db_fallback) > 0) {
        resolved_name <- ifelse(is.na(resolved_name), db_fallback$Title[1], resolved_name)
        resolved_smiles <- ifelse(is.na(resolved_smiles), db_fallback$SMILES[1], resolved_smiles)
        message("  [DB FILL] Filled missing info from DB.")
      }
    }

    # Update the row
    data$CID[i] <- resolved_cid
    if (!is.na(resolved_name)) data[[name_col]][i] <- resolved_name
    if (!is.na(resolved_smiles)) data[[smiles_col]][i] <- resolved_smiles

    # Update cache if new entry
    existing <- which(cid_cache_df$LookupName == name)
    if (length(existing) == 0) {
      cid_cache_df <- bind_rows(cid_cache_df,
                                data.frame(LookupName = name,
                                           ResolvedName = resolved_name,
                                           SMILES = resolved_smiles,
                                           CID = resolved_cid,
                                           stringsAsFactors = FALSE))
      message(paste("  [CACHE ADD]", name))
    }

    setTxtProgressBar(pb, i)
  }

  if (!is.null(db_con)) dbDisconnect(db_con)
  close(pb)

  return(list(data = data, cache = cid_cache_df))
}
