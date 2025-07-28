#' Clean and Standardize Compound Names
#'
#' Applies regex-based cleanup and casing corrections to specified column(s) of a data frame.
#'
#' @import data.table
#' @param df A data.frame or data.table with compound names.
#' @param column A character vector of column names to clean.
#' @return A cleaned data frame.
#' @export
fix_compound_names <- function(df, columns) {
  # Ensure data.table is loaded
  requireNamespace("data.table", quietly = TRUE)
  dt <- data.table::as.data.table(df)

  patterns <- c(
    "\\.(alpha|beta|gamma|omega)\\.-", "\\.(DELTA)\\.",
    "[[:space:]]*\\[IIN-based(?:[[:space:]]+on)?[[:space:]]*:[[:space:]]*[^)]*\\]",
    "NCG[[:alnum:]-]+_[[:alnum:]]+_", "NCG[^\\s!]+[\\s!]?",
    "Spectral Match to ", " from NIST14",
    "40.0 eV", "70.0 eV", "50.0 eV", ".00 eV", "\"", "CollisionEnergy", "\\|.*",
    "_", "\\[", "\\]", "\\{", "\\}", " :102040", " :205060",
    " - ", "30.0 eV", "Massbank:[[:alnum:]]+\\s? ", "MassbankEU:[[:alnum:]]+\\s? ",
    "MoNA:[[:alnum:]]+\\s? ", " \\(Chimeric precursor selection\\)",
    "; \\(M\\+.*"
  )
  replacements <- c(
    "\\1-", "delta-", "", "", "", "", "", "", "",
    "", "", "", "", "", " ", "(", ")", "(", ")", "",
    "", "", "", "", "", "", "", "", ""
  )

  for (col in columns) {
    for (i in seq_along(patterns)) {
      rows <- which(grepl(patterns[i], dt[[col]], ignore.case = TRUE))
      if (length(rows) > 0) {
        for (r in rows) {
          new_val <- gsub(patterns[i], replacements[i], dt[[col]][r], ignore.case = TRUE)
          data.table::set(dt, i = r, j = col, value = new_val)
        }
      }
    }

    # Trim whitespace before capitalization
    dt[[col]] <- trimws(dt[[col]])

    # Capitalization rules

    # 1. ALL CAPS -> Capitalize first letter only
    rows <- which(grepl("^[[:upper:]\\s\\d\\W]*$", dt[[col]]) & grepl("[[:alpha:]]", dt[[col]]))
    if (length(rows) > 0) {
      for (r in rows) {
        val <- dt[[col]][r]
        new_val <- paste0(toupper(substr(val, 1, 1)), tolower(substr(val, 2, nchar(val))))
        data.table::set(dt, i = r, j = col, value = new_val)
      }
    }

    # 2. all lowercase -> capitalize first letter
    rows <- which(grepl("^[[:lower:]\\s\\d\\W]*$", dt[[col]]) & grepl("[[:alpha:]]", dt[[col]]))
    if (length(rows) > 0) {
      for (r in rows) {
        val <- dt[[col]][r]
        new_val <- paste0(toupper(substr(val, 1, 1)), substr(val, 2, nchar(val)))
        data.table::set(dt, i = r, j = col, value = new_val)
      }
    }
  }

  # Return as data.frame for general compatibility
  return(as.data.frame(dt))
}
