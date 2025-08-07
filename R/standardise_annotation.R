standardise_annotation <- function(data,
                                   name_col = "compound_name",
                                   smiles_col = "smiles",
                                   cid_cache_df = NULL,
                                   cid_database_path = NULL) {

  # --- Libraries ---
  library(dplyr)
  library(DBI)
  library(RSQLite)

  # --- Checks ---
  if (is.null(cid_cache_df)) stop("cid_cache_df must be provided.")
  if (!(name_col %in% names(data))) stop("Missing column: ", name_col)
  if (!(smiles_col %in% names(data))) stop("Missing column: ", smiles_col)
  if (is.null(cid_database_path) || !file.exists(cid_database_path)) {
    stop("Database path must be provided and file must exist.")
  }
  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  # --- Filter and Initialise ---
  data <- data %>% filter(!grepl("candidate", .data[[name_col]], ignore.case = TRUE))
  if (nrow(data) == 0) return(list(data = data, cache = cid_cache_df))

  data[[name_col]] <- as.character(data[[name_col]])
  data[[smiles_col]] <- as.character(data[[smiles_col]])
  data$CID <- NA_real_

  # --- DB Connection ---
  message("[DB CONNECT] Connecting to CID SQLite DB...")
  db_con <- dbConnect(SQLite(), cid_database_path)
  on.exit({
    if (DBI::dbIsValid(db_con)) {
      dbDisconnect(db_con)
      message("[DB DISCONNECT] Closed DB connection.")
    }
  }, add = TRUE)

  # --- PASS 1: Resolve CIDs ---
  message("--- PASS 1: Resolving CIDs from PubChem ---")
  pb <- txtProgressBar(min = 0, max = nrow(data), style = 3)

  for (i in seq_len(nrow(data))) {
    name <- data[[name_col]][i]
    smiles <- data[[smiles_col]][i]

    if (is.na(name) || !nzchar(name)) {
      setTxtProgressBar(pb, i)
      next
    }

    pubchem_result <- get_cid_only_with_fallbacks(name, smiles, cid_cache_df)
    data$CID[i] <- pubchem_result$CID
    cid_cache_df <- pubchem_result$cache

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # --- PASS 2: Retrieve Properties from Local DB ---
  message("\n--- PASS 2: Retrieving Properties from Local DB ---")

  cids_to_lookup <- unique(data$CID[!is.na(data$CID) & data$CID > 0])

  if (length(cids_to_lookup) > 0) {
    cid_str <- paste(cids_to_lookup, collapse = ", ")

    query <- sprintf("
      SELECT CID, Title, SMILES, Formula AS Formula_db,
             IUPAC AS IUPAC_db, `Monoisotopic.Mass` AS Monoisotopic_Mass_db
      FROM pubchem_data
      WHERE CID IN (%s)
      GROUP BY CID", cid_str)

    message("[SQL QUERY]")
    message(query)

    db_props <- tryCatch(
      dbGetQuery(db_con, query),
      error = function(e) {
        message("[DB ERROR] ", e$message)
        return(NULL)
      }
    )

    if (!is.null(db_props) && nrow(db_props) > 0) {
      message("[DB LOOKUP] Retrieved ", nrow(db_props), " rows")

      # Merge properties into original data
      data <- data %>%
        left_join(db_props, by = "CID") %>%
        mutate(
          !!sym(name_col) := coalesce(Title, IUPAC_db, !!sym(name_col)),
          !!sym(smiles_col) := coalesce(SMILES, !!sym(smiles_col)),
          Formula = coalesce(Formula_db, Formula),
          IUPAC = coalesce(IUPAC_db, IUPAC),
          Monoisotopic.Mass = coalesce(Monoisotopic_Mass_db, Monoisotopic.Mass)
        ) %>%
        select(-Title, -SMILES, -Formula_db, -IUPAC_db, -Monoisotopic_Mass_db)
    } else {
      message("[DB LOOKUP] No rows returned.")
    }
  } else {
    message("[DB LOOKUP] No valid CIDs found.")
  }

  return(list(data = data, cache = cid_cache_df))
}
