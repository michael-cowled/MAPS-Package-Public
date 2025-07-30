#' Get CID from PubChem with Fallbacks and Caching
#'
#' Attempts to resolve a PubChem CID using compound name, SMILES, and synonyms.
#' If successful, returns a list of annotation details, and updates a global cache.
#'
#' @param name A character string representing the compound name.
#' @param smiles A character string representing the compound SMILES (optional).
#'
#' @return A named list with keys: \code{CID}, \code{ResolvedName}, \code{SMILES},
#'   or \code{NULL} if all lookups fail. MolecularFormula and MonoisotopicMass are no longer returned.
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

  # NEW: If a cached entry is found and its CID is -1, return it immediately as unresolvable.
  if (!is.null(cached_entry) && !is.na(cached_entry$CID[1]) && cached_entry$CID[1] == -1) {
    message(paste("  [CACHE HIT] Unresolvable CID (-1) found for:", name, ". Skipping PubChem lookup."))
    return(as.list(cached_entry[1, c("CID", "ResolvedName", "SMILES")]))
  }

  # If a cached entry is found and it has a CID and all properties (excluding removed ones), return it.
  # This prevents unnecessary API calls for fully resolved entries.
  if (!is.null(cached_entry) &&
      !is.na(cached_entry$CID[1]) &&
      !is.na(cached_entry$ResolvedName[1]) &&
      !is.na(cached_entry$SMILES[1])) {
    message(paste("  [CACHE HIT] Fully resolved CID found for:", name, " (CID:", cached_entry$CID[1], ")"))
    return(as.list(cached_entry[1, c("CID", "ResolvedName", "SMILES")])) # Only return these columns
  }

  # --- PubChem Lookup with Fallbacks ---
  resolved_cid <- NA_real_ # Initialize as numeric NA
  # resolved_title will be determined more carefully below
  resolved_props <- list(SMILES = NA_character_) # Initialize with NAs

  # Start with cached CID if available, even if incomplete (but not -1)
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
  # Only proceed if a valid CID (not NA and not -1) was found
  if (!is.na(resolved_cid) && resolved_cid != -1) {
    # Attempt to get properties (now only SMILES)
    temp_props <- get_pubchem(resolved_cid, "cid", "properties")

    # Ensure temp_props is a list before accessing elements, and initialize missing with NA
    if (is.list(temp_props)) {
      resolved_props$SMILES <- if (!is.null(temp_props$SMILES) && !is.na(temp_props$SMILES)) temp_props$SMILES else NA_character_
    } else {
      # If temp_props is NULL (due to get_pubchem error), ensure SMILES is NA
      resolved_props$SMILES <- NA_character_
    }

    # Attempt to get title from PubChem
    temp_title <- get_pubchem(resolved_cid, "cid", "title")
    new_pubchem_title <- if (is.character(temp_title) && !is.na(temp_title)) temp_title else NA_character_

    # Determine the best ResolvedName (cleaned to remove any "| ..." suffix)
    resolved_title <- NA_character_

    if (!is.null(cached_entry) && !is.na(cached_entry$ResolvedName[1])) {
      if (!is.na(new_pubchem_title)) {
        resolved_title <- sub("\\s*\\|.*$", "", new_pubchem_title)  # Prefer new title, cleaned
      } else {
        resolved_title <- sub("\\s*\\|.*$", "", cached_entry$ResolvedName[1])  # Clean cached ResolvedName
      }
    } else if (!is.na(new_pubchem_title)) {
      resolved_title <- sub("\\s*\\|.*$", "", new_pubchem_title)  # Clean new title
    } else {
      resolved_title <- sub("\\s*\\|.*$", "", name)  # Fallback to input name, cleaned
    }


    # Ensure SMILES is also merged if cached_entry has a better one
    if (!is.null(cached_entry) && is.na(resolved_props$SMILES) && !is.na(cached_entry$SMILES[1])) {
      resolved_props$SMILES <- cached_entry$SMILES[1]
    }

    # Create the new/updated entry - Note: MolecularFormula and MonoisotopicMass are excluded
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

      # Directly assign final_entry's ResolvedName to updated_row's ResolvedName.
      # This ensures the best ResolvedName determined above is always applied to the cache.
      updated_row$ResolvedName[1] <- final_entry$ResolvedName[1]
      updated_row$SMILES[1] <- ifelse(is.na(updated_row$SMILES[1]), final_entry$SMILES[1], updated_row$SMILES[1])
      updated_row$CID[1] <- ifelse(is.na(updated_row$CID[1]), final_entry$CID[1], updated_row$CID[1])

      cid_cache_df[existing_row_idx[1], ] <<- updated_row # Update the global cache
      return(as.list(updated_row[1, c("CID", "ResolvedName", "SMILES")])) # Only return these columns
    } else {
      # Add as a new row if no matching entry found
      # No need to add MolecularFormula and MonoisotopicMass here, as they are removed from template.
      cid_cache_df <<- dplyr::bind_rows(cid_cache_df, final_entry)
      return(as.list(final_entry[1, c("CID", "ResolvedName", "SMILES")])) # Only return these columns
    }
  }

  message(paste("  Could not resolve CID for:", name))
  # If no valid CID was found or properties couldn't be retrieved, return NULL
  return(NULL)
}
