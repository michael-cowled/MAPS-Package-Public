#' Get CID from PubChem with Fallbacks and Simplified Caching
#'
#' Attempts to resolve a PubChem CID using compound name and SMILES.
#' Updates a provided cache which only stores the CID, not other properties.
#'
#' @param name A character string representing the compound name.
#' @param smiles A character string representing the compound SMILES (optional).
#' @param cid_cache_df The data frame to be used as a cache, with columns LookupName and CID.
#'
#' @return A named list with keys: \code{CID} and the updated \code{cid_cache_df}.
#' @export
get_cid_only_with_fallbacks <- function(name, smiles = NA, cid_cache_df) {

  # Check cache for a pre-existing CID
  cached_entry <- cid_cache_df %>%
    filter(!is.na(LookupName) & LookupName == name) %>%
    slice(1)

  # Return cached CID if available
  if (nrow(cached_entry) > 0 && !is.na(cached_entry$CID[1])) {
    message(paste("  [CACHE HIT] CID found for:", name, " (CID:", cached_entry$CID[1], ")"))
    return(list(CID = cached_entry$CID[1], cache = cid_cache_df))
  }

  # --- PubChem Lookup with Fallbacks ---
  resolved_cid <- NA_real_

  # 1. Try name lookup
  resolved_cid <- get_pubchem_lite(name, "name", "cids")

  # 2. Try SMILES lookup
  if (is.na(resolved_cid) && !is.na(smiles) && smiles != "") {
    message(paste("  Name lookup failed. Trying SMILES:", smiles))
    resolved_cid <- get_pubchem_lite(smiles, "smiles", "cids")
  }

  # 3. Try synonym lookup
  if (is.na(resolved_cid)) {
    message(paste("  Name and SMILES failed. Trying synonym search for:", name))
    resolved_cid <- get_pubchem_lite(name, "synonym", "cids")
  }

  # --- Update Cache ---
  if (is.na(resolved_cid)) {
    message(paste("  [PUBCHEM] No CID found for:", name))
    # Create a failed entry to add to cache
    new_entry <- data.frame(LookupName = name, CID = -1, stringsAsFactors = FALSE)
  } else {
    new_entry <- data.frame(LookupName = name, CID = resolved_cid, stringsAsFactors = FALSE)
  }

  # Add the new entry to the cache
  cid_cache_df <- bind_rows(cid_cache_df, new_entry)

  return(list(CID = new_entry$CID[1], cache = cid_cache_df))
}
