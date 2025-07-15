## 4. get_cid_with_fallbacks
get_cid_with_fallbacks <- function(name, smiles = NA) {
  # Prefer to match by exact name first
  name_match <- cid_cache_df[!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name, ]
  
  # If no match on name, try matching on SMILES (if available)
  if (nrow(name_match) == 0 && !is.na(smiles) && smiles != "") {
    name_match <- cid_cache_df[!is.na(cid_cache_df$SMILES) & cid_cache_df$SMILES == smiles, ]
  }
  
  if (nrow(name_match) > 0 && !is.na(name_match$CID[1])) {
    message(paste("  [CACHE] CID found for:", name))
    return(list(
      CID = name_match$CID[1],
      ResolvedName = name_match$ResolvedName[1],
      SMILES = name_match$SMILES[1],
      MolecularFormula = name_match$MolecularFormula[1],
      MonoisotopicMass = name_match$MonoisotopicMass[1]
    ))
  }
  
  # Step 1: Name → CID
  cid <- get_pubchem(name, "name", "cids")
  
  # Step 2: Try SMILES
  if (is.na(cid) && !is.na(smiles) && smiles != "") {
    message(paste("  Name lookup failed. Trying SMILES:", smiles))
    cid <- get_pubchem(smiles, "smiles", "cids")
  }
  
  # Step 3: Try synonym
  if (is.na(cid)) {
    message(paste("  Name and SMILES failed. Trying synonym search for:", name))
    cid <- get_pubchem(name, "synonym", "cids")
  }
  
  # If CID found, get details and cache
  if (!is.na(cid)) {
    props <- get_pubchem(cid, "cid", "properties")
    title <- get_pubchem(cid, "cid", "title")
    
    # Ensure props is a list before accessing
    if (!is.list(props)) props <- list()
    if (!is.character(title)) title <- NA
    
    new_row <- data.frame(
      LookupName = name,
      ResolvedName = if (!is.null(title) && !is.na(title)) title else name,
      SMILES = if (!is.null(props$SMILES)) props$SMILES else if (!is.na(smiles)) smiles else NA,
      CID = cid,
      MolecularFormula = if (!is.null(props$MolecularFormula)) props$MolecularFormula else NA,
      MonoisotopicMass = if (!is.null(props$MonoisotopicMass)) props$MonoisotopicMass else NA,
      stringsAsFactors = FALSE
    )
    
    cid_cache_df <<- bind_rows(cid_cache_df, new_row)
    return(as.list(new_row[1, ]))
  }
  
  return(NULL)
}