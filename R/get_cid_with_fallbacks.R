# get_cid_with_fallbacks.R
#' Get CID from PubChem with Fallbacks and Caching
#'
#' Attempts to resolve a PubChem CID using compound name, SMILES, and synonyms.
#' If successful, returns a list of annotation details, and updates a global cache.
#'
#' @param name A character string representing the compound name.
#' @param smiles A character string representing the compound SMILES (optional).
#'
#' @return A named list with keys: \code{CID}, \code{ResolvedName}, \code{SMILES},
#'   or \code{NULL} if all lookups fail.
#' @export
#'
#' @details This function expects a global object named \code{cid_cache_df} to exist,
#' and it will be updated via \code{<<-} with new entries when successful lookups are made.
#'
#' @examples
#' \dontrun{
#' # Assuming cid_cache_df is initialized globally as per your setup
#' get_cid_with_fallbacks("glucose")
#' get_cid_with_fallbacks("some compound", smiles = "CCO")
#' }
get_cid_with_fallbacks <- function(name, smiles = NA) {

  # --- Cache Lookup and Prioritization ---
  cached_entry <- NULL

  # 1. Check cache by LookupName
  name_matches <- cid_cache_df[!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name, ]
  if (nrow(name_matches) > 0) {
    # If multiple entries for the same LookupName, pick the most complete one (e.g., with CID)
    if (any(!is.na(name_matches$CID))) {
      cached_entry <- name_matches[!is.na(name_matches$CID), ][1, , drop = FALSE]
    } else {
      cached_entry <- name_matches[1, , drop = FALSE] # Take the first if no CID found
    }
  }

  # 2. If no strong match by name, try by SMILES (if provided and valid)
  if (is.null(cached_entry) && !is.na(smiles) && smiles != "") {
    smiles_matches <- cid_cache_df[!is.na(cid_cache_df$SMILES) & cid_cache_df$SMILES == smiles, ]
    if (nrow(smiles_matches) > 0) {
      if (any(!is.na(smiles_matches$CID))) {
        cached_entry <- smiles_matches[!is.na(smiles_matches$CID), ][1, , drop = FALSE]
      } else {
        cached_entry <- smiles_matches[1, , drop = FALSE]
      }
    }
  }

  # If a cached entry is found and it has a CID and all properties, return it.
  # This prevents unnecessary API calls for fully resolved entries.
  if (!is.null(cached_entry) &&
      !is.na(cached_entry$CID[1]) &&
      !is.na(cached_entry$ResolvedName[1]) &&
      !is.na(cached_entry$SMILES[1])) {
    message(paste("  [CACHE HIT] Fully resolved CID found for:", name, " (CID:", cached_entry$CID[1], ")"))
    return(as.list(cached_entry[1, ]))
  }

  # --- PubChem Lookup with Fallbacks ---
  resolved_cid <- NA_real_ # Initialize as numeric NA
  resolved_title <- NA_character_ # Initialize as character NA
  resolved_props <- list(SMILES = NA_character_) # Initialize with NAs

  # Start with cached CID if available, even if incomplete
  if (!is.null(cached_entry) && !is.na(cached_entry$CID[1])) {
    resolved_cid <- cached_entry$CID[1]
    message(paste("  [CACHE INCOMPLETE] Re-attempting PubChem for CID:", resolved_cid))
  }

  # 1. Try name
  if (is.na(resolved_cid)) {
    resolved_cid <- get_pubchem(name, "name", "cids")
  }

  # 2. Try SMILES
  if (is.na(resolved_cid) && !is.na(smiles) && smiles != "") {
    message(paste("  Name lookup failed. Trying SMILES:", smiles))
    resolved_cid <- get_pubchem(smiles, "smiles", "cids")
  }

  # 3. Try synonym
  if (is.na(resolved_cid)) {
    message(paste("  Name and SMILES failed. Trying synonym search for:", name))
    resolved_cid <- get_pubchem(name, "synonym", "cids")
  }

  # --- Retrieve Properties and Title based on the best CID found ---
  if (!is.na(resolved_cid)) {
    # Attempt to get properties
    temp_props <- get_pubchem(resolved_cid, "cid", "properties")

    # Ensure temp_props is a list before accessing elements, and initialize missing with NA
    if (is.list(temp_props)) {
      resolved_props$SMILES <- if (!is.null(temp_props$SMILES) && !is.na(temp_props$SMILES)) temp_props$SMILES else NA_character_
    } else {
      # If temp_props is NULL (due to get_pubchem error), ensure all are NA
      resolved_props$SMILES <- NA_character_
    }

    # Attempt to get title
    temp_title <- get_pubchem(resolved_cid, "cid", "title")
    if (is.character(temp_title) && !is.na(temp_title)) {
      resolved_title <- temp_title
    } else {
      resolved_title <- name # Fallback to original name if title retrieval fails or returns NULL/NA
    }

    # Merge with cached values if they are more complete
    if (!is.null(cached_entry)) {
      if (is.na(resolved_title) && !is.na(cached_entry$ResolvedName[1])) resolved_title <- cached_entry$ResolvedName[1]
      if (is.na(resolved_props$SMILES) && !is.na(cached_entry$SMILES[1])) resolved_props$SMILES <- cached_entry$SMILES[1]
    }

    # Create the new/updated entry
    final_entry <- data.frame(
      LookupName = name,
      ResolvedName = resolved_title,
      SMILES = resolved_props$SMILES,
      CID = resolved_cid,
      stringsAsFactors = FALSE
    )

    # Update or add to cache
    existing_row_idx <- which(
      (!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name) |
        (!is.na(cid_cache_df$SMILES) & cid_cache_df$SMILES == final_entry$SMILES[1] & !is.na(final_entry$SMILES[1])) |
        (!is.na(cid_cache_df$CID) & cid_cache_df$CID == final_entry$CID[1] & !is.na(final_entry$CID[1]))
    )

    if (length(existing_row_idx) > 0) {
      # Initialize updated_row with the current cached row's structure and values
      updated_row <- cid_cache_df[existing_row_idx[1], , drop = FALSE]

      # Now, update the values, preferentially filling NA values
      updated_row$ResolvedName[1] <- ifelse(is.na(updated_row$ResolvedName[1]), final_entry$ResolvedName[1], updated_row$ResolvedName[1])
      updated_row$SMILES[1] <- ifelse(is.na(updated_row$SMILES[1]), final_entry$SMILES[1], updated_row$SMILES[1])
      updated_row$CID[1] <- ifelse(is.na(updated_row$CID[1]), final_entry$CID[1], updated_row$CID[1])

      cid_cache_df[existing_row_idx[1], ] <<- updated_row # Update the global cache
      message(paste("  [CACHE UPDATE] Updated CID cache for:", name))
      return(as.list(updated_row[1, ]))
    } else {
      # Add as a new row if no matching entry found
      cid_cache_df <<- dplyr::bind_rows(cid_cache_df, final_entry)
      message(paste("  [CACHE ADD] Added new CID cache entry for:", name))
      return(as.list(final_entry[1, ]))
    }
  }

  message(paste("  Could not resolve CID for:", name))
  return(NULL)
}
