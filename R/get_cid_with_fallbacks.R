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
  cache_idx <- which((!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name) |
                       (!is.na(cid_cache_df$SMILES) & !is.na(smiles) & cid_cache_df$SMILES == smiles))

  if (length(cache_idx) > 0) {
    match_row <- cid_cache_df[cache_idx[1], ]
    if (!is.na(match_row$CID)) {
      message(paste0("  [CACHE HIT] CID found for '", name, "' (CID: ", match_row$CID, ")"))
      return(list(CID = match_row$CID, cache = cid_cache_df))
    }
  }

  # --- 1b. Fast LipidMaps Check ---
  if (!is.null(lipids.file) && nrow(lipids.file) > 0) {
    synonym_pattern <- paste0("(^|;\\s*)", stringr::str_escape(name_clean), "(\\s*;|$)")

    lipid_names <- if ("Name" %in% names(lipids.file)) tolower(trimws(lipids.file$Name)) else NA
    lipid_sys <- if ("Systematic.Name" %in% names(lipids.file)) tolower(trimws(lipids.file$Systematic.Name)) else NA
    lipid_abb <- if ("Abbreviation" %in% names(lipids.file)) tolower(trimws(lipids.file$Abbreviation)) else NA
    lipid_syn <- if ("Synonyms" %in% names(lipids.file)) lipids.file$Synonyms else NA
    lipid_smiles <- if ("smiles" %in% names(lipids.file)) tolower(trimws(lipids.file$smiles)) else NA

    lipid_idx <- which(
      (!is.na(lipid_names) & lipid_names == name_clean) |
        (!is.na(lipid_sys) & lipid_sys == name_clean) |
        (!is.na(lipid_abb) & lipid_abb == name_clean) |
        (!is.na(lipid_syn) & grepl(synonym_pattern, tolower(lipid_syn))) |
        (!is.na(smiles_clean) & !is.na(lipid_smiles) & lipid_smiles == smiles_clean)
    )

    if (length(lipid_idx) > 0) {
      if (length(lipid_idx) > 1) {
        message(paste0("  [LIPID DB WARNING] Multiple matches found for '", name, "'. Using first match."))
      }

      lipid_match <- lipids.file[lipid_idx[1], ]

      match_cid <- if ("CID_numeric" %in% names(lipid_match)) {
        lipid_match$CID_numeric
      } else if ("CID" %in% names(lipid_match)) {
        suppressWarnings(as.numeric(lipid_match$CID))
      } else {
        NA_real_
      }

      if (!is.na(match_cid) && match_cid > 0) {
        message(paste0("  [LIPID DB] Found CID for '", name, "' in lipids.file (CID: ", match_cid, ")"))
        new_entry <- data.frame(LookupName = name, CID = match_cid, stringsAsFactors = FALSE)
        if ("SMILES" %in% names(cid_cache_df)) new_entry$SMILES <- smiles

        cid_cache_df <- rbind(cid_cache_df, new_entry)
        return(list(CID = match_cid, cache = cid_cache_df))
      }
    }
  }

  # --- 2. Offline Check ---
  if (offline) {
    message(paste0("  [OFFLINE] Skipping PubChem API for '", name, "'."))
    return(list(CID = NA_real_, cache = cid_cache_df))
  }

  # --- 3. PubChem API ---
  resolved_cid <- get_pubchem_lite(name, "name")

  if (is.na(resolved_cid) && !is.na(smiles_clean) && smiles_clean != "" && smiles_clean != "N/A") {
    message(paste0("  Name lookup failed for '", name, "'. Trying SMILES: ", smiles))
    resolved_cid <- get_pubchem_lite(smiles, "smiles")
  }

  if (is.na(resolved_cid)) {
    message(paste0("  Name and SMILES failed for '", name, "'. Trying synonym search..."))
    resolved_cid <- get_pubchem_lite(name, "synonym")
  }

  # --- 4. Update cache ---
  if (!is.na(resolved_cid)) {
    message(paste0("  [PUBCHEM] Found CID for '", name, "': ", resolved_cid))
    new_cid <- resolved_cid
  } else {
    message(paste0("  [PUBCHEM] No CID found for '", name, "'."))
    new_cid <- -1
  }

  new_entry <- data.frame(LookupName = name, CID = new_cid, stringsAsFactors = FALSE)
  if ("SMILES" %in% names(cid_cache_df)) new_entry$SMILES <- smiles

  cid_cache_df <- rbind(cid_cache_df, new_entry)
  return(list(CID = new_entry$CID, cache = cid_cache_df))
}
