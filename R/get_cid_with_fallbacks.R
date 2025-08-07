#' Get CID from PubChem with Fallbacks and Caching
#'
#' Attempts to resolve a PubChem CID using compound name, SMILES, and synonyms.
#' If successful, returns a list of annotation details, and updates a provided cache.
#'
#' @param name A character string representing the compound name.
#' @param smiles A character string representing the compound SMILES (optional).
#' @param db_con An active DBI connection object to your SQLite database.
#' @param cid_cache_df The data frame to be used as a cache.
#'
#' @return A named list with keys: \code{CID}, \code{ResolvedName}, \code{SMILES},
#'   and the updated \code{cid_cache_df}, or \code{NULL} if all lookups fail.
#' @export
#'
get_cid_with_fallbacks <- function(name, smiles = NA, db_con = NULL, cid_cache_df) {

  # --- Cache Lookup and Prioritization ---
  cached_entry <- NULL

  # 1. Check cache by LookupName
  name_matches <- cid_cache_df[!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name, ]
  if (nrow(name_matches) > 0) {
    if (any(!is.na(name_matches$CID))) {
      cached_entry <- name_matches[!is.na(name_matches$CID), ][1, , drop = FALSE]
    } else {
      cached_entry <- name_matches[1, , drop = FALSE]
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

  # If a cached entry is found and its CID is -1, return it immediately as unresolvable.
  if (!is.null(cached_entry) && !is.na(cached_entry$CID[1]) && cached_entry$CID[1] == -1) {
    message(paste("  [CACHE HIT] Unresolvable CID (-1) found for:", name, ". Skipping PubChem lookup."))
    # Return a list containing both the resolved info and the original cache
    return(list(
      CID = cached_entry$CID[1],
      ResolvedName = cached_entry$ResolvedName[1],
      SMILES = cached_entry$SMILES[1],
      cache = cid_cache_df
    ))
  }

  # If a cached entry is found and it has a CID and all properties, return it.
  if (!is.null(cached_entry) &&
      !is.na(cached_entry$CID[1]) &&
      !is.na(cached_entry$ResolvedName[1]) &&
      !is.na(cached_entry$SMILES[1])) {
    message(paste("  [CACHE HIT] Fully resolved CID found for:", name, " (CID:", cached_entry$CID[1], ")"))
    # Return a list containing both the resolved info and the original cache
    return(list(
      CID = cached_entry$CID[1],
      ResolvedName = cached_entry$ResolvedName[1],
      SMILES = cached_entry$SMILES[1],
      cache = cid_cache_df
    ))
  }

  # --- PubChem Lookup with Fallbacks ---
  resolved_cid <- NA_real_
  resolved_props <- list(SMILES = NA_character_)

  # Start with cached CID if available, even if incomplete (but not -1)
  if (!is.null(cached_entry) && !is.na(cached_entry$CID[1])) {
    resolved_cid <- cached_entry$CID[1]
    message(paste("  [CACHE INCOMPLETE] Re-attempting PubChem for CID:", resolved_cid))
  }

  # 1. Try name
  if (is.na(resolved_cid)) {
    resolved_cid <- get_pubchem(name, "name", "cids", db_con)
  }

  # 2. Try SMILES
  if (is.na(resolved_cid) && !is.na(smiles) && smiles != "") {
    message(paste("  Name lookup failed. Trying SMILES:", smiles))
    resolved_cid <- get_pubchem(smiles, "smiles", "cids", db_con)
  }

  # 3. Try synonym
  if (is.na(resolved_cid)) {
    message(paste("  Name and SMILES failed. Trying synonym search for:", name))
    resolved_cid <- get_pubchem(name, "synonym", "cids", db_con)
  }

  # --- Retrieve Properties and Title based on the best CID found ---
  if (!is.na(resolved_cid) && resolved_cid != -1) {
    temp_props <- get_pubchem(resolved_cid, "cid", "properties", db_con)

    if (is.list(temp_props)) {
      resolved_props$SMILES <- if (!is.null(temp_props$SMILES) && !is.na(temp_props$SMILES)) temp_props$SMILES else NA_character_
    } else {
      resolved_props$SMILES <- NA_character_
    }

    temp_title <- get_pubchem(resolved_cid, "cid", "title", db_con)
    new_pubchem_title <- if (is.character(temp_title) && !is.na(temp_title)) temp_title else NA_character_

    # Determine the best ResolvedName
    resolved_title <- NA_character_

    if (!is.null(cached_entry) && !is.na(cached_entry$ResolvedName[1])) {
      if (!is.na(new_pubchem_title)) {
        resolved_title <- sub("\\s*\\|.*$", "", new_pubchem_title)
      } else {
        resolved_title <- sub("\\s*\\|.*$", "", cached_entry$ResolvedName[1])
      }
    } else if (!is.na(new_pubchem_title)) {
      resolved_title <- sub("\\s*\\|.*$", "", new_pubchem_title)
    } else {
      resolved_title <- sub("\\s*\\|.*$", "", name)
    }

    # Ensure SMILES is also merged if cached_entry has a better one
    if (!is.null(cached_entry) && is.na(resolved_props$SMILES) && !is.na(cached_entry$SMILES[1])) {
      resolved_props$SMILES <- cached_entry$SMILES[1]
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
      updated_row <- cid_cache_df[existing_row_idx[1], , drop = FALSE]
      updated_row$ResolvedName[1] <- final_entry$ResolvedName[1]
      updated_row$SMILES[1] <- ifelse(is.na(updated_row$SMILES[1]), final_entry$SMILES[1], updated_row$SMILES[1])
      updated_row$CID[1] <- ifelse(is.na(updated_row$CID[1]), final_entry$CID[1], updated_row$CID[1])
      cid_cache_df[existing_row_idx[1], ] = updated_row

      return(list(
        CID = updated_row$CID[1],
        ResolvedName = updated_row$ResolvedName[1],
        SMILES = updated_row$SMILES[1],
        cache = cid_cache_df
      ))
    } else {
      cid_cache_df <- dplyr::bind_rows(cid_cache_df, final_entry)

      return(list(
        CID = final_entry$CID[1],
        ResolvedName = final_entry$ResolvedName[1],
        SMILES = final_entry$SMILES[1],
        cache = cid_cache_df
      ))
    }
  }

  # If no valid CID was found or properties couldn't be retrieved, update cache with a failed entry
  message(paste("  Could not resolve CID for:", name))

  # Create a failed entry to add to cache
  failed_entry <- data.frame(
    LookupName = name,
    ResolvedName = NA_character_,
    SMILES = NA_character_,
    CID = -1,
    stringsAsFactors = FALSE
  )

  # Check if a failed entry for this LookupName already exists in cache
  existing_failed_idx <- which(!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name)

  if (length(existing_failed_idx) == 0) {
    cid_cache_df <- dplyr::bind_rows(cid_cache_df, failed_entry)
  }

  return(list(
    CID = -1,
    ResolvedName = NA_character_,
    SMILES = NA_character_,
    cache = cid_cache_df
  ))
}
