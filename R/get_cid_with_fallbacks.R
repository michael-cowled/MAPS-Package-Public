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

  name_clean <- trimws(tolower(name))
  smiles_clean <- if(!is.na(smiles)) trimws(tolower(smiles)) else NA

  # --- 1. Fast Cache Check (Base R) ---
  # Only grab the first matching row to mimic slice(1)
  cache_idx <- which((!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name) |
                       (!is.na(cid_cache_df$SMILES) & !is.na(smiles) & cid_cache_df$SMILES == smiles))

  if (length(cache_idx) > 0) {
    match_row <- cid_cache_df[cache_idx[1], ]
    if (!is.na(match_row$CID)) {
      # message(paste0("  [CACHE HIT] CID found for '", name, "' (CID: ", match_row$CID, ")"))
      return(list(CID = match_row$CID, cache = cid_cache_df))
    }
  }

  # --- 1b. Fast LipidMaps Check ---
  if (!is.null(lipids.file) && nrow(lipids.file) > 0) {
    # Using vectorized grepl instead of rowwise() + str_split()
    synonym_pattern <- paste0("(^|;\\s*)", stringr::str_escape(name_clean), "(\\s*;|$)")

    lipid_idx <- which(
      tolower(trimws(lipids.file$Name)) == name_clean |
        tolower(trimws(lipids.file$Systematic.Name)) == name_clean |
        tolower(trimws(lipids.file$Abbreviation)) == name_clean |
        (!is.na(lipids.file$Synonyms) & grepl(synonym_pattern, tolower(lipids.file$Synonyms))) |
        (!is.na(smiles_clean) & !is.na(lipids.file$smiles) & tolower(trimws(lipids.file$smiles)) == smiles_clean)
    )

    if (length(lipid_idx) > 0) {
      lipid_match <- lipids.file[lipid_idx[1], ]
      # message(paste0("  [LIPID DB] Found CID for '", name, "' (CID: ", lipid_match$CID_numeric, ")"))

      new_entry <- data.frame(LookupName = name, CID = lipid_match$CID_numeric, stringsAsFactors = FALSE)
      if ("SMILES" %in% names(cid_cache_df)) new_entry$SMILES <- smiles # keep schema aligned

      cid_cache_df <- rbind(cid_cache_df, new_entry)
      return(list(CID = lipid_match$CID_numeric, cache = cid_cache_df))
    }
  }

  # --- 2. Offline Check ---
  if (offline) {
    return(list(CID = NA_real_, cache = cid_cache_df))
  }

  # --- 3. PubChem API ---
  resolved_cid <- get_pubchem_lite(name, "name")
  if (is.na(resolved_cid) && !is.na(smiles_clean) && smiles_clean != "" && smiles_clean != "N/A") {
    resolved_cid <- get_pubchem_lite(smiles, "smiles")
  }
  if (is.na(resolved_cid)) {
    resolved_cid <- get_pubchem_lite(name, "synonym")
  }

  # --- 4. Update cache ---
  new_cid <- if (!is.na(resolved_cid)) resolved_cid else -1
  new_entry <- data.frame(LookupName = name, CID = new_cid, stringsAsFactors = FALSE)
  if ("SMILES" %in% names(cid_cache_df)) new_entry$SMILES <- smiles

  cid_cache_df <- rbind(cid_cache_df, new_entry)
  return(list(CID = new_entry$CID, cache = cid_cache_df))
}
