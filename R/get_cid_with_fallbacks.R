#' Get CID from PubChem with Fallbacks and LipidMaps Lookup
#'
#' Attempts to resolve a PubChem CID using compound name and SMILES, checking LipidMaps first,
#' and updates a provided cache which only stores the CID.
#'
#' @param name A character string representing the compound name.
#' @param smiles A character string representing the compound SMILES (optional).
#' @param cid_cache_df The data frame to be used as a cache, with columns LookupName and CID.
#' @param lipids.file A data frame to use to lookup lipid names for PubChem CID lookups.
#' @param offline Logical; if TRUE, strictly avoids live PubChem API lookups.
#'
#' @return A named list with keys: \code{CID} and the updated \code{cid_cache_df}.
#' @export
get_cid_only_with_fallbacks <- function(name, smiles = NA, cid_cache_df, lipids.file, offline = FALSE) {

  # --- 0. Clean input for comparison ---
  name_clean <- trimws(tolower(name))
  smiles_clean <- if(!is.na(smiles)) trimws(tolower(smiles)) else NA

  # --- 1. Check cache first ---
  cached_entry <- cid_cache_df %>%
    dplyr::filter((!is.na(LookupName) & LookupName == name) |
                    (!is.na(SMILES) & SMILES == smiles)) %>%
    dplyr::slice(1)

  if (nrow(cached_entry) > 0 && !is.na(cached_entry$CID[1])) {
    message(paste0("  [CACHE HIT] CID found for '", name, "' (CID: ", cached_entry$CID[1], ")"))
    return(list(CID = cached_entry$CID[1], cache = cid_cache_df))
  }

  # --- 1b. Clean LipidMaps file ---
  lipids.file.clean <- lipids.file %>%
    dplyr::mutate(CID_numeric = suppressWarnings(as.numeric(CID))) %>%
    dplyr::filter(!is.na(CID_numeric) & CID_numeric > 0) %>%
    dplyr::group_by(CID_numeric) %>%
    dplyr::summarise(
      Name = dplyr::first(stats::na.omit(Name)),
      Systematic.Name = dplyr::first(stats::na.omit(Systematic.Name)),
      Abbreviation = dplyr::first(stats::na.omit(Abbreviation)),
      HMDB.ID = dplyr::first(stats::na.omit(HMDB.ID)),
      Synonyms = paste(unique(stats::na.omit(Synonyms)), collapse = "; "),
      smiles = dplyr::first(stats::na.omit(smiles)),
      .groups = "drop"
    )

  # --- 1c. Check LipidMaps ---
  lipid_match <- lipids.file.clean %>%
    dplyr::rowwise() %>%
    dplyr::filter(
      tolower(trimws(Name)) == name_clean |
        tolower(trimws(Systematic.Name)) == name_clean |
        tolower(trimws(Abbreviation)) == name_clean |
        (!is.na(Synonyms) && name_clean %in% tolower(trimws(stringr::str_split(Synonyms, ";\\s*")[[1]]))) |
        (!is.na(smiles) && !is.na(.data$smiles) && tolower(trimws(.data$smiles)) == smiles_clean)
    ) %>%
    dplyr::ungroup()

  if (nrow(lipid_match) > 1) {
    message(paste0("  [LIPID DB WARNING] Multiple matches found for '", name,
                   "'. Using first match (CID: ", lipid_match$CID_numeric[1], ")"))
    lipid_match <- lipid_match %>% dplyr::slice(1)
  }

  if (nrow(lipid_match) > 0 && !is.na(lipid_match$CID_numeric[1])) {
    message(paste0("  [LIPID DB] Found CID for '", name, "' in lipids.file (CID: ", lipid_match$CID_numeric[1], ")"))
    new_entry <- data.frame(LookupName = name, CID = lipid_match$CID_numeric[1], stringsAsFactors = FALSE)
    cid_cache_df <- dplyr::bind_rows(cid_cache_df, new_entry)
    return(list(CID = lipid_match$CID_numeric[1], cache = cid_cache_df))
  }

  # ==========================================
  # --- NEW: OFFLINE CHECK BLOCK ---
  # ==========================================
  if (offline) {
    message(paste0("  [OFFLINE MODE] Skipping PubChem API for '", name, "'."))
    # Return NA_real_ to indicate it wasn't resolved, but do NOT add a -1 to the cache
    # because we haven't actually confirmed it doesn't exist on PubChem.
    return(list(CID = NA_real_, cache = cid_cache_df))
  }

  # --- 2. PubChem lookup (API) ---
  resolved_cid <- NA_real_

  # Send the raw, unencoded strings directly to the downstream function
  resolved_cid <- get_pubchem_lite(name, "name")

  if (is.na(resolved_cid) && !is.na(smiles_clean) && smiles_clean != "" && smiles_clean != "N/A") {
    message(paste0("  Name lookup failed for '", name, "'. Trying SMILES: ", smiles))
    resolved_cid <- get_pubchem_lite(smiles, "smiles")
  }

  if (is.na(resolved_cid)) {
    message(paste0("  Name and SMILES failed for '", name, "'. Trying synonym search..."))
    resolved_cid <- get_pubchem_lite(name, "synonym")
  }

  # --- 3. Update cache ---
  if (!is.na(resolved_cid)) {
    message(paste0("  [PUBCHEM] Found CID for '", name, "': ", resolved_cid))
    new_entry <- data.frame(LookupName = name, CID = resolved_cid, stringsAsFactors = FALSE)
  } else {
    message(paste0("  [PUBCHEM] No CID found for '", name, "'."))
    new_entry <- data.frame(LookupName = name, CID = -1, stringsAsFactors = FALSE)
  }

  cid_cache_df <- dplyr::bind_rows(cid_cache_df, new_entry)
  return(list(CID = new_entry$CID[1], cache = cid_cache_df))
}
