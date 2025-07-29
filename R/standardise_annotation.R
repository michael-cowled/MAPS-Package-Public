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
standardise_annotation <- function(df,
                                   name_col = "compound_name",
                                   smiles_col = "smiles",
                                   cid_cache_path = NULL,
                                   db_con = NULL,
                                   cid_database_path = NULL) {

  library(DBI)
  library(RSQLite)
  library(dplyr)

  # --- Handle DB connection ---
  internal_con <- FALSE
  if (!is.null(db_con) && dbIsValid(db_con)) {
    message("[DB CONNECT] Using externally provided DB connection.")
  } else if (!is.null(cid_database_path) && file.exists(cid_database_path)) {
    message("[DB CONNECT] Opening internal connection to: ", cid_database_path)
    db_con <- dbConnect(SQLite(), cid_database_path)
    internal_con <- TRUE
    on.exit({
      if (dbIsValid(db_con)) {
        dbDisconnect(db_con)
        message("[DB DISCONNECT] Closed internal DB connection.")
      }
    }, add = TRUE)
  } else {
    stop("No valid database connection or database path provided.")
  }

  # --- Load existing cache ---
  cache <- if (!is.null(cid_cache_path) && file.exists(cid_cache_path)) {
    read.csv(cid_cache_path, stringsAsFactors = FALSE)
  } else {
    data.frame(Query = character(), CID = numeric(), stringsAsFactors = FALSE)
  }

  # --- Annotate each row ---
  for (i in seq_len(nrow(df))) {
    compound_name <- trimws(df[[name_col]][i])
    smiles <- if (!is.null(smiles_col)) df[[smiles_col]][i] else NA

    if (is.na(compound_name) || grepl("candidate", compound_name, ignore.case = TRUE)) {
      next
    }

    message("Processing: ", compound_name)

    cached_row <- cache %>% filter(Query == compound_name)
    if (nrow(cached_row) > 0) {
      cid <- cached_row$CID[1]
      message("  [CACHE HIT] CID ", cid, " retrieved from cache.")
    } else {
      message("  [CACHE MISS] Querying CID for: ", compound_name)
      cid <- tryCatch({
        query <- paste0("SELECT CID FROM cid_lookup WHERE LookupName = ", dbQuoteString(db_con, compound_name))
        result <- dbGetQuery(db_con, query)
        if (nrow(result) > 0) result$CID[1] else -1
      }, error = function(e) {
        warning("  [DB ERROR] Failed to query CID for ", compound_name, ": ", e$message)
        -1
      })

      cache <- bind_rows(cache, data.frame(Query = compound_name, CID = cid, stringsAsFactors = FALSE))
    }

    df$CID[i] <- cid

    # --- Only do get_pubchem() if CID is valid ---
    if (!is.na(cid) && cid > 0) {
      pubchem_data <- tryCatch({
        get_pubchem(cid = cid, type = "cid", property = "properties", db_con = db_con)
      }, error = function(e) {
        warning("  [PUBCHEM ERROR] Failed for CID ", cid, ": ", e$message)
        NULL
      })

      if (!is.null(pubchem_data)) {
        df$SMILES[i] <- pubchem_data$smiles
        df$Title[i]  <- pubchem_data$title
      }
    } else {
      message("  [SKIP] CID = -1, skipping PubChem query.")
    }
  }

  # --- Save updated cache ---
  if (!is.null(cid_cache_path)) {
    write.csv(cache, cid_cache_path, row.names = FALSE)
    message("[CACHE] Written to: ", cid_cache_path)
  }

  return(df)
}
