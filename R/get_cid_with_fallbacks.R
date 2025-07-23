#' Get CID from PubChem with Fallbacks and Caching
#'
#' Attempts to resolve a PubChem CID using compound name, SMILES, and synonyms.
#' If successful, returns a list of annotation details, and updates a global cache.
#'
#' @param name A character string representing the compound name.
#' @param smiles A character string representing the compound SMILES (optional).
#'
#' @return A named list with keys: \code{CID}, \code{ResolvedName}, \code{SMILES},
#'   \code{MolecularFormula}, and \code{MonoisotopicMass}, or \code{NULL} if all lookups fail.
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
  # Prioritize lookup by CID if available in cache for a given name/SMILES
  cached_entry <- NULL

  # 1. Check cache by LookupName
  name_matches <- cid_cache_df[!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name, ]
  if (nrow(name_matches) > 0) {
    # If multiple entries for the same LookupName, pick the most complete one (e.g., with CID)
    if (any(!is.na(name_matches$CID))) {
      cached_entry <- name_matches[!is.na(name_matches$CID), ][1, ]
    } else {
      cached_entry <- name_matches[1, ] # Take the first if no CID found
    }
  }

  # 2. If no strong match by name, try by SMILES (if provided and valid)
  if (is.null(cached_entry) && !is.na(smiles) && smiles != "") {
    smiles_matches <- cid_cache_df[!is.na(cid_cache_df$SMILES) & cid_cache_df$SMILES == smiles, ]
    if (nrow(smiles_matches) > 0) {
      if (any(!is.na(smiles_matches$CID))) {
        cached_entry <- smiles_matches[!is.na(smiles_matches$CID), ][1, ]
      } else {
        cached_entry <- smiles_matches[1, ]
      }
    }
  }

  # If a cached entry is found and it has a CID and all properties, return it.
  # This prevents unnecessary API calls for fully resolved entries.
  if (!is.null(cached_entry) &&
      !is.na(cached_entry$CID) &&
      !is.na(cached_entry$ResolvedName) &&
      !is.na(cached_entry$SMILES) &&
      !is.na(cached_entry$MolecularFormula) &&
      !is.na(cached_entry$MonoisotopicMass)) {
    message(paste("  [CACHE HIT] Fully resolved CID found for:", name, " (CID:", cached_entry$CID, ")"))
    return(as.list(cached_entry))
  }

  # --- PubChem Lookup with Fallbacks ---
  resolved_cid <- NA
  resolved_title <- NA
  resolved_props <- list()

  # Start with cached CID if available, even if incomplete
  if (!is.null(cached_entry) && !is.na(cached_entry$CID)) {
    resolved_cid <- cached_entry$CID
    message(paste("  [CACHE INCOMPLETE] Re-attempting PubChem for CID:", resolved_cid))
  }

  # 1. Try name if CID is still NA
  if (is.na(resolved_cid)) {
    resolved_cid <- get_pubchem(name, "name", "cids")
  }

  # 2. Try SMILES if CID is still NA and SMILES is provided
  if (is.na(resolved_cid) && !is.na(smiles) && smiles != "") {
    message(paste("  Name lookup failed. Trying SMILES:", smiles))
    resolved_cid <- get_pubchem(smiles, "smiles", "cids")
  }

  # 3. Try synonym if CID is still NA
  if (is.na(resolved_cid)) {
    message(paste("  Name and SMILES failed. Trying synonym search for:", name))
    resolved_cid <- get_pubchem(name, "synonym", "cids")
  }

  # --- Retrieve Properties and Title based on the best CID found ---
  if (!is.na(resolved_cid)) {
    # Attempt to get properties if not already fully cached or if properties are missing
    if (is.null(cached_entry) || is.na(cached_entry$SMILES) || is.na(cached_entry$MolecularFormula) || is.na(cached_entry$MonoisotopicMass)) {
      temp_props <- get_pubchem(resolved_cid, "cid", "properties")
      if (is.list(temp_props)) {
        resolved_props <- temp_props
      }
    } else {
      resolved_props <- list(
        SMILES = cached_entry$SMILES,
        MolecularFormula = cached_entry$MolecularFormula,
        MonoisotopicMass = cached_entry$MonoisotopicMass
      )
    }

    # Attempt to get title if not already fully cached or if title is missing
    if (is.null(cached_entry) || is.na(cached_entry$ResolvedName)) {
      temp_title <- get_pubchem(resolved_cid, "cid", "title")
      if (is.character(temp_title) && !is.na(temp_title)) {
        resolved_title <- temp_title
      } else {
        resolved_title <- name # Fallback to original name if title retrieval fails
      }
    } else {
      resolved_title <- cached_entry$ResolvedName
    }

    # Use existing cached values if they are more complete than newly fetched ones
    if (!is.null(cached_entry)) {
      resolved_title <- ifelse(is.na(resolved_title) && !is.na(cached_entry$ResolvedName), cached_entry$ResolvedName, resolved_title)
      resolved_props$SMILES <- ifelse(is.null(resolved_props$SMILES) && !is.na(cached_entry$SMILES), cached_entry$SMILES, resolved_props$SMILES)
      resolved_props$MolecularFormula <- ifelse(is.null(resolved_props$MolecularFormula) && !is.na(cached_entry$MolecularFormula), cached_entry$MolecularFormula, resolved_props$MolecularFormula)
      resolved_props$MonoisotopicMass <- ifelse(is.null(resolved_props$MonoisotopicMass) && !is.na(cached_entry$MonoisotopicMass), cached_entry$MonoisotopicMass, resolved_props$MonoisotopicMass)
    }

    # Create the new/updated entry
    final_entry <- data.frame(
      LookupName = name,
      ResolvedName = if (!is.null(resolved_title) && !is.na(resolved_title)) resolved_title else name,
      SMILES = if (!is.null(resolved_props$SMILES)) resolved_props$SMILES else if (!is.na(smiles)) smiles else NA,
      CID = resolved_cid,
      MolecularFormula = if (!is.null(resolved_props$MolecularFormula)) resolved_props$MolecularFormula else NA,
      MonoisotopicMass = if (!is.null(resolved_props$MonoisotopicMass)) resolved_props$MonoisotopicMass else NA,
      stringsAsFactors = FALSE
    )

    # Update or add to cache
    # First, try to find an existing entry to update
    existing_row_idx <- which(
      (!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name) |
        (!is.na(cid_cache_df$SMILES) & cid_cache_df$SMILES == final_entry$SMILES & !is.na(final_entry$SMILES)) |
        (!is.na(cid_cache_df$CID) & cid_cache_df$CID == final_entry$CID & !is.na(final_entry$CID))
    )

    if (length(existing_row_idx) > 0) {
      # Update the existing row with more complete information
      # Use `pmax` for numeric and `ifelse` for character to pick non-NA or current value
      current_cached_row <- cid_cache_df[existing_row_idx[1], ]

      updated_row <- data.frame(stringsAsFactors = FALSE)
      updated_row$LookupName <- current_cached_row$LookupName # Keep original lookup name
      updated_row$ResolvedName <- ifelse(is.na(current_cached_row$ResolvedName), final_entry$ResolvedName, current_cached_row$ResolvedName)
      updated_row$SMILES <- ifelse(is.na(current_cached_row$SMILES), final_entry$SMILES, current_cached_row$SMILES)
      updated_row$CID <- ifelse(is.na(current_cached_row$CID), final_entry$CID, current_cached_row$CID)
      updated_row$MolecularFormula <- ifelse(is.na(current_cached_row$MolecularFormula), final_entry$MolecularFormula, current_cached_row$MolecularFormula)
      updated_row$MonoisotopicMass <- ifelse(is.na(current_cached_row$MonoisotopicMass), final_entry$MonoisotopicMass, current_cached_row$MonoisotopicMass)

      cid_cache_df[existing_row_idx[1], ] <<- updated_row # Update the global cache
      message(paste("  [CACHE UPDATE] Updated CID cache for:", name))
      return(as.list(updated_row))
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
