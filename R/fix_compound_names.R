#' Clean and Standardize Compound Names
#'
#' Applies regex-based cleanup and casing corrections to specified column(s) of a data frame.
#'
#' @import data.table
#' @param df A data.frame or data.table with compound names.
#' @param columns A character vector of column names to clean.
#' @return A cleaned data frame.
#' @export
fix_compound_names <- function(df, columns) {
  requireNamespace("data.table", quietly = TRUE)
  dt <- data.table::as.data.table(df)

  # Added the "::" pattern at the beginning to strip metadata early
  patterns <- c(
    "[[:space:]]*::.*", "\\.(alpha|beta|gamma|omega)\\.-", "\\.(DELTA)\\.",
    "[[:space:]]*\\[IIN-based(?:[[:space:]]+on)?[[:space:]]*:[[:space:]]*[^)]*\\]",
    "NCG[[:alnum:]-]+_[[:alnum:]]+_", "NCG[^\\s!]+[\\s!]?",
    "Spectral Match to ", " from NIST14",
    "40.0 eV", "70.0 eV", "50.0 eV", ".00 eV", "\"", "CollisionEnergy", "\\|.*",
    "_", "\\[", "\\]", "\\{", "\\}", " :102040", " :205060",
    " - ", "30.0 eV", "Massbank:[[:alnum:]]+\\s? ", "MassbankEU:[[:alnum:]]+\\s? ",
    "MoNA:[[:alnum:]]+\\s? ", " \\(Chimeric precursor selection\\)",
    "; \\(M\\+.*", " M+formate"
  )

  replacements <- c(
    "", "\\1-", "delta-", "", "", "", "", "", "",
    "", "", "", "", "", " ", "(", ")", "(", ")", "",
    "", "", "", "", "", "", "", "", ""
  )

  for (col in columns) {
    # Vectorized cleanup is significantly faster than row-wise loops
    for (i in seq_along(patterns)) {
      dt[[col]] <- gsub(patterns[i], replacements[i], dt[[col]], ignore.case = TRUE)
    }

    # Trim whitespace
    dt[[col]] <- trimws(dt[[col]])

    # 1. ALL CAPS -> Capitalize first letter only
    all_caps_idx <- grepl("^[[:upper:]\\s\\d\\W]*$", dt[[col]]) & grepl("[[:alpha:]]", dt[[col]])
    if (any(all_caps_idx)) {
      vals <- dt[[col]][all_caps_idx]
      dt[all_caps_idx, (col) := paste0(toupper(substr(vals, 1, 1)), tolower(substr(vals, 2, nchar(vals))))]
    }

    # 2. all lowercase -> capitalize first letter
    all_lower_idx <- grepl("^[[:lower:]\\s\\d\\W]*$", dt[[col]]) & grepl("[[:alpha:]]", dt[[col]])
    if (any(all_lower_idx)) {
      vals <- dt[[col]][all_lower_idx]
      dt[all_lower_idx, (col) := paste0(toupper(substr(vals, 1, 1)), substr(vals, 2, nchar(vals)))]
    }
  }

  return(as.data.frame(dt))
}
