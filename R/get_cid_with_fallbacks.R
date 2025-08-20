#' Get CID from PubChem with Fallbacks and Simplified Caching (Corrected)
#'
#' Attempts to resolve a PubChem CID using compound name and SMILES.
#' Updates a provided cache which only stores the CID, not other properties.
#'
#' @param name A character string representing the compound name.
#' @param smiles A character string representing the compound SMILES (optional).
#' @param cid_cache_df The data frame to be used as a cache, with columns LookupName and CID.
#' @param lipids.file A data frame to use to lookup lipid names for PubChem CID lookups.
#'
#' @return A named list with keys: \code{CID} and the updated \code{cid_cache_df}.
#' @export
get_cid_only_with_fallbacks <- function(name, smiles = NA, cid_cache_df, lipids.file) {

  # Check cache for a pre-existing CID
  cached_entry <- cid_cache_df %>%
    filter(!is.na(LookupName) & LookupName == name | !is.na(smiles) & SMILES == smiles) %>%
    slice(1)

  # --- 1. If a cached CID is found, return it immediately ---
  if (nrow(cached_entry) > 0 && !is.na(cached_entry$CID[1])) {
    message(paste("  [CACHE HIT] CID found for '", name, "' (CID:", cached_entry$CID[1], ")"))
    return(list(CID = cached_entry$CID[1], cache = cid_cache_df))
  }

  # --- 1b. Clean lipids.file to ensure unique CIDs ---
  lipids.file.clean <- lipids.file %>%
    group_by(CID) %>%
    summarise(
      Name = first(na.omit(Name)),
      Systematic.Name = first(na.omit(Systematic.Name)),
      Abbreviation = first(na.omit(Abbreviation)),
      HMDB.ID = first(na.omit(HMDB.ID)),
      Synonyms = paste(unique(na.omit(Synonyms)), collapse = "; "),
      smiles = first(na.omit(smiles)),
      .groups = "drop"
    )

  # --- 1c. Check lipids.file before moving on to PubChem ---
  lipid_match <- lipids.file.clean %>%
    rowwise() %>%
    filter(
      tolower(Name) == tolower(name) |
        tolower(Systematic.Name) == tolower(name) |
        tolower(Abbreviation) == tolower(name) |
        (!is.na(Synonyms) && tolower(name) %in% tolower(trimws(str_split(Synonyms, ";\\s*")[[1]]))) |
        (!is.na(smiles) && !is.na(.data$smiles) && .data$smiles == smiles)
    ) %>%
    ungroup()

  # If multiple matches, warn and pick the first
  if (nrow(lipid_match) > 1) {
    message(paste("  [LIPID DB WARNING] Multiple matches found for '", name,
                  "'. Using first match (CID:", lipid_match$CID[1], ")"))
    lipid_match <- lipid_match %>% slice(1)
  }

  if (nrow(lipid_match) > 0 && !is.na(lipid_match$CID[1])) {
    message(paste("  [LIPID DB] Found CID for '", name, "' in lipids.file (CID:", lipid_match$CID[1], ")"))
    # Add to cache
    new_entry <- data.frame(LookupName = name, CID = lipid_match$CID[1], stringsAsFactors = FALSE)
    cid_cache_df <- bind_rows(cid_cache_df, new_entry)
    return(list(CID = lipid_match$CID[1], cache = cid_cache_df))
  }

  # --- 2. If not in cache or lipids.file, perform PubChem lookup with fallbacks ---
  resolved_cid <- NA_real_

  # Try name lookup
  resolved_cid <- get_pubchem_lite(name, "name", "cids")

  # Try SMILES lookup
  if (is.na(resolved_cid) && !is.na(smiles) && smiles != "" && smiles != "N/A" && smiles != "NA") {
    message(paste("  Name lookup failed for '", name, "'. Trying SMILES:", smiles))
    resolved_cid <- get_pubchem_lite(smiles, "smiles", "cids")
  }

  # Try synonym lookup
  if (is.na(resolved_cid)) {
    message(paste("  Name and SMILES failed for '", name, "'. Trying synonym search..."))
    resolved_cid <- get_pubchem_lite(name, "synonym", "cids")
  }

  # --- 3. Update Cache based on lookup result ---
  if (!is.na(resolved_cid)) {
    message(paste("  [PUBCHEM] Found CID for '", name, "': ", resolved_cid))
    new_entry <- data.frame(LookupName = name, CID = resolved_cid, stringsAsFactors = FALSE)
  } else {
    message(paste("  [PUBCHEM] No CID found for '", name, "'."))
    new_entry <- data.frame(LookupName = name, CID = -1, stringsAsFactors = FALSE)
  }

  cid_cache_df <- bind_rows(cid_cache_df, new_entry)

  return(list(CID = new_entry$CID[1], cache = cid_cache_df))
}
