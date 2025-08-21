#' Get CID from PubChem with Fallbacks and LipidMaps Lookup
#'
#' Attempts to resolve a PubChem CID using compound name and SMILES, checking LipidMaps first,
#' and updates a provided cache which only stores the CID.
#'
#' @param name A character string representing the compound name.
#' @param smiles A character string representing the compound SMILES (optional).
#' @param cid_cache_df The data frame to be used as a cache, with columns LookupName and CID.
#' @param lipids.file A data frame to use to lookup lipid names for PubChem CID lookups.
#'
#' @return A named list with keys: \code{CID} and the updated \code{cid_cache_df}.
#' @export
get_cid_only_with_fallbacks <- function(name, smiles = NA, cid_cache_df, lipids.file) {

  # --- 0. Clean input for comparison ---
  name_clean <- trimws(tolower(name))
  smiles_clean <- if(!is.na(smiles)) trimws(tolower(smiles)) else NA

  # --- 1. Check cache first ---
  cached_entry <- cid_cache_df %>%
    filter(!is.na(LookupName) & LookupName == name) %>%
    slice(1)

  if (nrow(cached_entry) > 0 && !is.na(cached_entry$CID[1])) {
    message(paste0("  [CACHE HIT] CID found for '", name, "' (CID: ", cached_entry$CID[1], ")"))
    return(list(CID = cached_entry$CID[1], cache = cid_cache_df))
  }

  # --- 1b. Clean LipidMaps file ---
  lipids.file.clean <- lipids.file %>%
    mutate(CID_numeric = suppressWarnings(as.numeric(CID))) %>%
    filter(!is.na(CID_numeric) & CID_numeric > 0) %>%
    group_by(CID_numeric) %>%
    summarise(
      Name = first(na.omit(Name)),
      Systematic.Name = first(na.omit(Systematic.Name)),
      Abbreviation = first(na.omit(Abbreviation)),
      HMDB.ID = first(na.omit(HMDB.ID)),
      Synonyms = paste(unique(na.omit(Synonyms)), collapse = "; "),
      smiles = first(na.omit(smiles)),
      .groups = "drop"
    )

  # --- 1c. Check LipidMaps ---
  lipid_match <- lipids.file.clean %>%
    rowwise() %>%
    filter(
      tolower(trimws(Name)) == name_clean |
        tolower(trimws(Systematic.Name)) == name_clean |
        tolower(trimws(Abbreviation)) == name_clean |
        (!is.na(Synonyms) && name_clean %in% tolower(trimws(str_split(Synonyms, ";\\s*")[[1]]))) |
        (!is.na(smiles) && !is.na(.data$smiles) && tolower(trimws(.data$smiles)) == smiles_clean)
    ) %>%
    ungroup()

  if (nrow(lipid_match) > 1) {
    message(paste0("  [LIPID DB WARNING] Multiple matches found for '", name,
                   "'. Using first match (CID: ", lipid_match$CID_numeric[1], ")"))
    lipid_match <- lipid_match %>% slice(1)
  }

  if (nrow(lipid_match) > 0 && !is.na(lipid_match$CID_numeric[1])) {
    message(paste0("  [LIPID DB] Found CID for '", name, "' in lipids.file (CID: ", lipid_match$CID_numeric[1], ")"))
    new_entry <- data.frame(LookupName = name, CID = lipid_match$CID_numeric[1], stringsAsFactors = FALSE)
    cid_cache_df <- bind_rows(cid_cache_df, new_entry)
    return(list(CID = lipid_match$CID_numeric[1], cache = cid_cache_df))
  }

  # --- 2. PubChem lookup ---
  resolved_cid <- NA_real_

  resolved_cid <- get_pubchem_lite(name, "name", "cids")

  if (is.na(resolved_cid) && !is.na(smiles_clean) && smiles_clean != "" && smiles_clean != "N/A") {
    message(paste0("  Name lookup failed for '", name, "'. Trying SMILES: ", smiles))
    resolved_cid <- get_pubchem_lite(smiles, "smiles", "cids")
  }

  if (is.na(resolved_cid)) {
    message(paste0("  Name and SMILES failed for '", name, "'. Trying synonym search..."))
    resolved_cid <- get_pubchem_lite(name, "synonym", "cids")
  }

  # --- 3. Update cache ---
  if (!is.na(resolved_cid)) {
    message(paste0("  [PUBCHEM] Found CID for '", name, "': ", resolved_cid))
    new_entry <- data.frame(LookupName = name, CID = resolved_cid, stringsAsFactors = FALSE)
  } else {
    message(paste0("  [PUBCHEM] No CID found for '", name, "'."))
    new_entry <- data.frame(LookupName = name, CID = -1, stringsAsFactors = FALSE)
  }

  cid_cache_df <- bind_rows(cid_cache_df, new_entry)
  return(list(CID = new_entry$CID[1], cache = cid_cache_df))
}
