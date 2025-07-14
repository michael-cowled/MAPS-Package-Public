## 5. standardise_annotation
standardise_annotation <- function(data, name_col, smiles_col) {
  # Filter out rows where name_col contains "candidate" (case-insensitive)
  data <- data[!grepl("candidate", data[[name_col]], ignore.case = TRUE), ]
  
  n <- nrow(data)
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  # Ensure required output columns exist
  data[[paste0("CID")]] <- NA
  data[[paste0("MolecularFormula")]] <- NA
  data[[paste0("MonoisotopicMass")]] <- NA
  
  for (i in seq_len(n)) {
    current_name <- data[[name_col]][i]
    current_smiles <- data[[smiles_col]][i]
    
    if (!is.na(current_name) && nzchar(current_name)) {
      message(paste("Processing:", current_name))
      
      info <- get_cid_with_fallbacks(current_name, current_smiles)
      
      if (!is.null(info)) {
        data[[paste0("CID")]][i] <- info$CID
        if (!is.na(info$ResolvedName) && info$ResolvedName != current_name) {
          info$ResolvedName <- sub(" \\|.*", "", info$ResolvedName)
          data[[name_col]][i] <- info$ResolvedName
          message(paste("  Updated", name_col, ":", current_name, "→", info$ResolvedName))
        }
        data[[smiles_col]][i] <- info$SMILES
        data[[paste0("MolecularFormula")]][i] <- info$MolecularFormula
        data[[paste0("MonoisotopicMass")]][i] <- info$MonoisotopicMass
      } else {
        message(paste("  No CID found for", current_name))
      }
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(data)
}