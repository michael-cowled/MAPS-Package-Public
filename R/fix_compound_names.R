## 2. fix_compound_names 
fix_compound_names <- function(df, column) {
  # Convert to data.table for faster operations
  df <- data.table::as.data.table(df)
  
  # Define patterns and replacements for initial cleaning
  patterns <- c(
    "\\.(alpha|beta|gamma|omega)\\.-", "\\.(DELTA)\\.", "[[:space:]]*\\[IIN-based(?:[[:space:]]+on)?[[:space:]]*:[[:space:]]*[^)]*\\]",
    "NCG[[:alnum:]-]+_[[:alnum:]]+_", "NCG[^\\s!]+[\\s!]?",
    "Spectral Match to ", " from NIST14",
    "40.0 eV", "70.0 eV", "50.0 eV", ".00 eV", '"', "CollisionEnergy", "\\|.*",
    "_", "\\[", "\\]", "\\{", "\\}", " :102040", " :205060", " - ",
    "30.0 eV", "Massbank:[[:alnum:]]+\\s? ", "MassbankEU:[[:alnum:]]+\\s? ",
    "MoNA:[[:alnum:]]+\\s? ", " (Chimeric precursor selection)", "; \\(M\\+.*"
  )
  replacements <- c(
    "\\1-", "delta-", "", "", "", "", "", "", "", "", "", "", "", "",
    " ", "(", ")", "(", ")", "", "", "", "", "", "", "", "", "", ""
  )
  
  # Loop through columns and apply initial fixes
  for (col in column) {
    for (i in seq_along(patterns)) {
      df[grepl(patterns[i], get(col), ignore.case = TRUE), (col) := gsub(patterns[i], replacements[i], get(col))]
    }
    
    # --- New section for capitalization rules ---
    # Apply capitalization rules after initial cleaning
    
    # 1. If ALL CAPS, change to first letter capital, rest lowercase
    # Check if the entire string is uppercase (excluding spaces and non-alphabetic characters for the check)
    df[, (col) := ifelse(
      grepl("^[[:upper:]\\s\\d\\W]*$", get(col)) & grepl("[[:alpha:]]", get(col)), # Ensure it contains at least one letter
      paste0(toupper(substring(get(col), 1, 1)), tolower(substring(get(col), 2))),
      get(col)
    )]
    
    # 2. If ALL LOWERCASE, change first letter to uppercase
    # Check if the entire string is lowercase (excluding spaces and non-alphabetic characters for the check)
    df[, (col) := ifelse(
      grepl("^[[:lower:]\\s\\d\\W]*$", get(col)) & grepl("[[:alpha:]]", get(col)), # Ensure it contains at least one letter
      paste0(toupper(substring(get(col), 1, 1)), substring(get(col), 2)),
      get(col)
    )]
    
    # Mixed case names are left unchanged by these rules, as they won't match the ALL CAPS or ALL LOWERCASE patterns.
    # --- End of new section ---
  }
  
  return(as.data.frame(df))
}