#' Standardise Compound Annotations Using PubChem
#'
#' For each row in a data frame, attempts to resolve compound annotations using PubChem.
#' Updates the compound name, SMILES, and adds a column for CID.
#' Rows containing "candidate" (case-insensitive) in the name column are excluded.
#'
#' @param data A data frame containing compound annotations.
#' @param name_col Name of the column in `data` containing compound names (as a string).
#' @param smiles_col Name of the column in `data` containing SMILES strings (as a string).
#' @param filtered_merged_cid_data_path Path to the "Filtered_Merged_CID_Data.tsv" file.
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
#' #                              "path/to/Filtered_Merged_CID_Data.tsv",
#' #                              "path/to/cid_cache.csv")
#' }
standardise_annotation <- function(data, name_col, smiles_col, filtered_merged_cid_data_path, cid_cache_path = "cid_cache.csv") {
  # Load dplyr silently for bind_rows
  suppressMessages(library(dplyr))

  # --- Initialize or Load cid_cache_df ---
  # Define the structure for an empty cache (MolecularFormula and MonoisotopicMass removed)
  empty_cache_template <- data.frame(
    LookupName = character(), ResolvedName = character(), SMILES = character(),
    CID = numeric(),
    stringsAsFactors = FALSE
  )

  # Check if cid_cache_df exists in the global environment and is a data frame
  if (!exists("cid_cache_df", envir = .GlobalEnv) || !is.data.frame(get("cid_cache_df", envir = .GlobalEnv))) {
    message(paste("Initializing cid_cache_df. Attempting to load from:", cid_cache_path))
    if (file.exists(cid_cache_path)) {
      cid_cache_df_loaded <- tryCatch({
        read.csv(cid_cache_path, stringsAsFactors = FALSE)
      }, error = function(e) {
        message(paste("  Error loading cid_cache.csv:", e$message, "Initializing empty cache."))
        return(empty_cache_template) # Return template on error
      })
      # Ensure loaded data frame has all required columns from the new template, adding them as NA if missing
      for (col in names(empty_cache_template)) {
        if (!(col %in% names(cid_cache_df_loaded))) {
          # Add missing columns with appropriate empty vectors for consistency
          if (is.character(empty_cache_template[[col]])) {
            cid_cache_df_loaded[[col]] <- character(nrow(cid_cache_df_loaded))
          } else if (is.numeric(empty_cache_template[[col]])) {
            cid_cache_df_loaded[[col]] <- numeric(nrow(cid_cache_df_loaded))
          } else {
            cid_cache_df_loaded[[col]] <- NA # Fallback for other types
          }
        }
      }
      # Remove any columns from loaded cache that are no longer in the template
      cols_to_keep <- intersect(names(cid_cache_df_loaded), names(empty_cache_template))
      cid_cache_df <<- cid_cache_df_loaded[, cols_to_keep, drop = FALSE]
      # Reorder columns to match template
      cid_cache_df <<- cid_cache_df[, names(empty_cache_template), drop = FALSE]

    } else {
      message("  cid_cache.csv not found. Initializing empty cache.")
      cid_cache_df <<- empty_cache_template
    }
  } else {
    message("cid_cache_df already exists in global environment. Ensuring correct structure.")
    current_cid_cache_df <- get("cid_cache_df", envir = .GlobalEnv)

    # Use bind_rows to ensure all columns are present and correctly typed.
    # This handles both empty and non-empty data frames gracefully.
    # We bind the current cache with an empty version of the template to ensure all columns exist.
    # This will also drop columns not in the template.
    cid_cache_df <<- dplyr::bind_rows(current_cid_cache_df, empty_cache_template[0, , drop = FALSE])

    # Reorder columns to match template
    cid_cache_df <<- cid_cache_df[, names(empty_cache_template), drop = FALSE]
  }

  # --- New: Check if the input data frame is empty ---
  if (nrow(data) == 0) {
    message("Input data frame is empty. Returning as is.")
    return(data)
  }

  # --- New: Validate name_col and smiles_col existence and type ---
  if (!(name_col %in% names(data))) {
    stop(paste0("Error: Column '", name_col, "' not found in the input data frame."))
  }
  if (!(smiles_col %in% names(data))) {
    stop(paste0("Error: Column '", smiles_col, "' not found in the input data frame."))
  }

  # Ensure the columns are character type for grepl and other string operations
  if (!is.character(data[[name_col]])) {
    warning(paste0("Warning: Column '", name_col, "' is not character type. Attempting coercion."))
    data[[name_col]] <- as.character(data[[name_col]])
  }
  if (!is.character(data[[smiles_col]])) {
    warning(paste0("Warning: Column '", smiles_col, "' is not character type. Attempting coercion."))
    data[[smiles_col]] <- as.character(data[[smiles_col]])
  }


  # Filter out rows with "candidate"
  data <- data[!grepl("candidate", data[[name_col]], ignore.case = TRUE), ]

  n <- nrow(data)
  # If filtering results in an empty data frame, return it.
  if (n == 0) {
    message("After filtering 'candidate' rows, the data frame is empty. Returning as is.")
    return(data)
  }

  pb <- txtProgressBar(min = 0, max = n, style = 3)

  # Ensure output columns exist
  data[["CID"]] <- NA_real_ # Initialize as numeric NA

  # Load Filtered_Merged_CID_Data.tsv once outside the loop for efficiency
  library(DBI)
  library(RSQLite)

  # Connect to the SQLite database instead of reading TSV
  filtered_merged_db_path <- filtered_merged_cid_data_path  # Rename for clarity
  db_con <- NULL
  if (!is.null(filtered_merged_db_path) && file.exists(filtered_merged_db_path)) {
    message(paste("Connecting to CID database at:", filtered_merged_db_path))
    db_con <- dbConnect(SQLite(), filtered_merged_db_path)
  } else {
    message("  CID database path not provided or file does not exist. This fallback will not be used.")
  }


  for (i in seq_len(n)) {
    current_name <- data[[name_col]][i]
    current_smiles <- data[[smiles_col]][i]
    resolved_cid <- NA_real_ # To store the CID found for the current row
    resolved_name_from_lookup <- NA_character_ # To store the best resolved name from lookup
    resolved_smiles_from_lookup <- NA_character_ # To store the best resolved smiles from lookup

    # Added check for NA or empty string in current_name to prevent errors in subsequent steps
    if (is.na(current_name) || !nzchar(current_name)) {
      message(paste("Skipping row", i, ": current_name is NA or empty."))
      setTxtProgressBar(pb, i)
      next # Skip to next iteration if name is not valid
    }

    message(paste("Processing:", current_name))

    # Step 2: Check if CID is already in cid_cache_df for a specific LookupName.
    # Modified condition: only consider a cache hit if ResolvedName and SMILES are NOT NA
    cached_by_name <- cid_cache_df[!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == current_name &
                                     !is.na(cid_cache_df$ResolvedName) & !is.na(cid_cache_df$SMILES), ]

    if (nrow(cached_by_name) > 0 && !is.na(cached_by_name$CID[1])) {
      resolved_cid <- cached_by_name$CID[1]
      resolved_name_from_lookup <- cached_by_name$ResolvedName[1]
      resolved_smiles_from_lookup <- cached_by_name$SMILES[1]
      message(paste("  [CACHE HIT - LookupName] Resolved from cache for:", current_name, "(CID:", resolved_cid, ")"))
    }

    # Step 3: If no CID found from initial cache check, then check pubchem using get_cid_with_fallbacks
    info <- NULL
    if (is.na(resolved_cid)) {
      info <- get_cid_with_fallbacks(current_name, current_smiles)
      if (!is.null(info) && !is.na(info$CID)) {
        resolved_cid <- info$CID
        resolved_name_from_lookup <- info$ResolvedName
        resolved_smiles_from_lookup <- info$SMILES
        message(paste("  [PUBCHEM LOOKUP] Found CID:", resolved_cid, "for", current_name))
      } else {
        message(paste("  No CID found via PubChem for", current_name))
      }
    }

    # If a CID was resolved (either from initial cache or PubChem lookup)
    if (!is.na(resolved_cid) && resolved_cid != -1) {
      # Now, check the TSV file as a final fallback for name/SMILES if they are still NA from lookup
      # This is Step 5 from your original request, simplified and integrated.
      if ((is.na(resolved_name_from_lookup) || !nzchar(resolved_name_from_lookup)) &&
          (is.na(resolved_smiles_from_lookup) || !nzchar(resolved_smiles_from_lookup)) && # Check SMILES too
          !is.null(filtered_merged_data)) {
        tsv_entry <- NULL
        if (!is.null(db_con)) {
          query <- sprintf("SELECT Title, SMILES FROM cid_data WHERE CID = %d", resolved_cid)
          tsv_entry <- dbGetQuery(db_con, query)
        }
        if (nrow(tsv_entry) > 0) {
          # Use TSV Title and SMILES if lookup didn't provide them
          resolved_name_from_lookup <- tsv_entry$Title[1]
          resolved_smiles_from_lookup <- tsv_entry$SMILES[1]
          message(paste("  [TSV LOOKUP] Found CID:", resolved_cid, "in TSV. Using TSV data."))

          # Update cid_cache_df with TSV info if it's new or better
          new_cache_entry <- data.frame(
            LookupName = current_name, # Original name used for lookup
            ResolvedName = resolved_name_from_lookup,
            SMILES = resolved_smiles_from_lookup,
            CID = resolved_cid,
            stringsAsFactors = FALSE
          )

          # Add/update cache entry only if it's new or has more complete data
          existing_cache_idx <- which(cid_cache_df$LookupName == current_name & cid_cache_df$CID == resolved_cid)
          if (length(existing_cache_idx) > 0) {
            # Update existing entry if new data is more complete
            current_cache_row <- cid_cache_df[existing_cache_idx[1], ]
            if (is.na(current_cache_row$ResolvedName) && !is.na(new_cache_entry$ResolvedName[1])) {
              cid_cache_df$ResolvedName[existing_cache_idx[1]] <<- new_cache_entry$ResolvedName[1]
            }
            if (is.na(current_cache_row$SMILES) && !is.na(new_cache_entry$SMILES[1])) {
              cid_cache_df$SMILES[existing_cache_idx[1]] <<- new_cache_entry$SMILES[1]
            }
            message(paste("  [CACHE UPDATE] Updated existing CID cache entry from TSV for CID:", resolved_cid))
          } else {
            # Add as a new row
            cid_cache_df <<- dplyr::bind_rows(cid_cache_df, new_cache_entry)
            message(paste("  [CACHE ADD] Added new CID cache entry from TSV for CID:", resolved_cid))
          }

        } else {
          message(paste("  CID", resolved_cid, "not found in Filtered_Merged_CID_Data.tsv."))
        }
      } else if (is.null(filtered_merged_data)) {
        message("  Filtered_Merged_CID_Data.tsv not loaded, skipping TSV lookup.")
      }

      # Update the main data frame with the best resolved information
      data[["CID"]][i] <- resolved_cid
      if (!is.na(resolved_name_from_lookup) && nzchar(resolved_name_from_lookup)) {
        data[[name_col]][i] <- resolved_name_from_lookup
      }
      if (!is.na(resolved_smiles_from_lookup) && nzchar(resolved_smiles_from_lookup)) {
        data[[smiles_col]][i] <- resolved_smiles_from_lookup
      }

    } else {
      data[["CID"]][i] <- -1
      message(paste("  No CID could be resolved for", current_name))
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)
  return(data)

  if (!is.null(db_con)) {
    dbDisconnect(db_con)
  }
}
