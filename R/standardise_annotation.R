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
  # Check if cid_cache_df exists in the global environment and is a data frame
  if (!exists("cid_cache_df", envir = .GlobalEnv) || !is.data.frame(get("cid_cache_df", envir = .GlobalEnv))) {
    message(paste("Initializing cid_cache_df. Attempting to load from:", cid_cache_path))
    if (file.exists(cid_cache_path)) {
      cid_cache_df <<- tryCatch({
        read.csv(cid_cache_path, stringsAsFactors = FALSE)
      }, error = function(e) {
        message(paste("  Error loading cid_cache.csv:", e$message, "Initializing empty cache."))
        # Initialize with all expected columns, including those no longer used in new entries
        data.frame(
          LookupName = character(), ResolvedName = character(), SMILES = character(),
          CID = numeric(), MolecularFormula = character(), MonoisotopicMass = numeric(),
          stringsAsFactors = FALSE
        )
      })
    } else {
      message("  cid_cache.csv not found. Initializing empty cache.")
      # Initialize with all expected columns, including those no longer used in new entries
      cid_cache_df <<- data.frame(
        LookupName = character(), ResolvedName = character(), SMILES = character(),
        CID = numeric(), MolecularFormula = character(), MonoisotopicMass = numeric(),
        stringsAsFactors = FALSE
      )
    }
  } else {
    message("cid_cache_df already exists in global environment. Using existing cache.")
    # Ensure existing cid_cache_df has the necessary columns, add if missing
    required_cols <- c("LookupName", "ResolvedName", "SMILES", "CID", "MolecularFormula", "MonoisotopicMass")
    for (col in required_cols) {
      if (!(col %in% names(cid_cache_df))) {
        cid_cache_df[[col]] <<- NA # Add missing columns as NA
        message(paste("  Added missing column to cid_cache_df:", col))
      }
    }
  }


  # Filter out rows with "candidate"
  data <- data[!grepl("candidate", data[[name_col]], ignore.case = TRUE), ]

  n <- nrow(data)
  pb <- txtProgressBar(min = 0, max = n, style = 3)

  # Ensure output columns exist
  data[["CID"]] <- NA_real_ # Initialize as numeric NA

  # Load Filtered_Merged_CID_Data.tsv once outside the loop for efficiency
  filtered_merged_data <- NULL
  if (!is.null(filtered_merged_cid_data_path) && file.exists(filtered_merged_cid_data_path)) {
    message(paste("Loading Filtered_Merged_CID_Data.tsv from:", filtered_merged_cid_data_path))
    filtered_merged_data <- tryCatch({
      read.delim(filtered_merged_cid_data_path, sep = "\t", stringsAsFactors = FALSE)
    }, error = function(e) {
      message(paste("  Error loading Filtered_Merged_CID_Data.tsv:", e$message))
      return(NULL)
    })
    # Ensure expected columns are present
    if (!is.null(filtered_merged_data) && !all(c("CID", "Title", "SMILES") %in% names(filtered_merged_data))) {
      message("  Filtered_Merged_CID_Data.tsv must contain 'CID', 'Title', and 'SMILES' columns.")
      filtered_merged_data <- NULL
    }
  } else {
    message("  Filtered_Merged_CID_Data.tsv path not provided or file does not exist. This fallback will not be used.")
  }


  for (i in seq_len(n)) {
    current_name <- data[[name_col]][i]
    current_smiles <- data[[smiles_col]][i]
    resolved_cid <- NA_real_ # To store the CID found for the current row

    if (!is.na(current_name) && nzchar(current_name)) {
      message(paste("Processing:", current_name))

      # Step 2: Check if CID is already in cid_cache_df for a specific LookupName.
      cached_by_name <- cid_cache_df[!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == current_name, ]

      if (nrow(cached_by_name) > 0 && !is.na(cached_by_name$CID[1])) {
        resolved_cid <- cached_by_name$CID[1]
        data[[name_col]][i] <- cached_by_name$ResolvedName[1]
        data[[smiles_col]][i] <- cached_by_name$SMILES[1]
        message(paste("  [CACHE HIT - LookupName] Updated", name_col, ":", current_name, "→", data[[name_col]][i], "(CID:", resolved_cid, ")"))
      }

      # Step 3: If no CID found from initial cache check, then check pubchem using get_cid_with_fallbacks
      info <- NULL
      if (is.na(resolved_cid)) {
        info <- get_cid_with_fallbacks(current_name, current_smiles)
        if (!is.null(info) && !is.na(info$CID)) {
          resolved_cid <- info$CID
          data[[name_col]][i] <- info$ResolvedName
          data[[smiles_col]][i] <- info$SMILES
          message(paste("  [PUBCHEM LOOKUP] Found CID:", resolved_cid, "for", current_name))
        } else {
          message(paste("  No CID found via PubChem for", current_name))
        }
      }

      # If a CID was resolved (either from initial cache or PubChem lookup)
      if (!is.na(resolved_cid)) {
        # Step 4: Then for this new CID, check the cid_cache_df to make sure it is not already there,
        # and apply the logic of step 2, if it is.
        # This covers cases where get_cid_with_fallbacks might have returned a CID that wasn't fully cached
        # or was only partially cached.
        cached_by_cid <- cid_cache_df[!is.na(cid_cache_df$CID) & cid_cache_df$CID == resolved_cid, ]

        if (nrow(cached_by_cid) > 0 && !is.na(cached_by_cid$ResolvedName[1])) {
          data[[name_col]][i] <- cached_by_cid$ResolvedName[1]
          data[[smiles_col]][i] <- cached_by_cid$SMILES[1]
          message(paste("  [CACHE HIT - CID] Updated", name_col, ":", current_name, "→", data[[name_col]][i], "(CID:", resolved_cid, ")"))
        } else {
          # Step 5: If it isn't (i.e., CID not fully cached), then load in ONLY the row
          # corresponding to CID = row for Filtered_Merged_CID_Data.tsv.
          if (!is.null(filtered_merged_data)) {
            tsv_entry <- filtered_merged_data[filtered_merged_data$CID == resolved_cid, ]
            if (nrow(tsv_entry) > 0) {
              # Return the Title to compound.name and SMILES to smiles in the df of interest.
              data[[name_col]][i] <- tsv_entry$Title[1]
              data[[smiles_col]][i] <- tsv_entry$SMILES[1]
              message(paste("  [TSV LOOKUP] Found CID:", resolved_cid, "in TSV. Updated", name_col, ":", current_name, "→", data[[name_col]][i]))

              # Then update the cid_cache_df with the information pertaining to this new CID.
              # i.e. CID = CID, SMILES = smiles, ResolvedName = Title, and LookupName is the original name being queried in step 2.
              new_cache_entry <- data.frame(
                LookupName = current_name,
                ResolvedName = tsv_entry$Title[1],
                SMILES = tsv_entry$SMILES[1],
                CID = resolved_cid,
                stringsAsFactors = FALSE
              )
              # Ensure new columns for MolecularFormula and MonoisotopicMass are added as NA if they exist in cid_cache_df structure
              if ("MolecularFormula" %in% names(cid_cache_df) && !("MolecularFormula" %in% names(new_cache_entry))) {
                new_cache_entry$MolecularFormula <- NA_character_
              }
              if ("MonoisotopicMass" %in% names(cid_cache_df) && !("MonoisotopicMass" %in% names(new_cache_entry))) {
                new_cache_entry$MonoisotopicMass <- NA_real_
              }

              cid_cache_df <<- dplyr::bind_rows(cid_cache_df, new_cache_entry)
              message(paste("  [CACHE ADD] Added new CID cache entry from TSV for CID:", resolved_cid))
            } else {
              message(paste("  CID", resolved_cid, "not found in Filtered_Merged_CID_Data.tsv."))
            }
          } else {
            message("  Filtered_Merged_CID_Data.tsv not loaded, skipping TSV lookup.")
          }
        }
        data[["CID"]][i] <- resolved_cid
      } else {
        message(paste("  No CID could be resolved for", current_name))
      }
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)
  return(data)
}
