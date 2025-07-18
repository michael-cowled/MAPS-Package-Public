#' Clean and Standardize Compound Names
#'
#' Applies regex-based cleanup and casing corrections to specified column(s) of a data frame.
#'
#' @import data.table
#' @param df A data.frame or data.table with compound names.
#' @param column A character vector of column names to clean.
#' @return A cleaned data frame.
#' @export
fix_compound_names <- function(df, column) {
  # Ensure we're working in a data.table-aware environment
  local({
    df <- data.table::as.data.table(df)

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

    for (col in column) {
      for (i in seq_along(patterns)) {
        df[grepl(patterns[i], get(col), ignore.case = TRUE), (col) := gsub(patterns[i], replacements[i], get(col))]
      }

      df[, (col) := ifelse(
        grepl("^[[:upper:]\\s\\d\\W]*$", get(col)) & grepl("[[:alpha:]]", get(col)),
        paste0(toupper(substring(get(col), 1, 1)), tolower(substring(get(col), 2))),
        get(col)
      )]

      df[, (col) := ifelse(
        grepl("^[[:lower:]\\s\\d\\W]*$", get(col)) & grepl("[[:alpha:]]", get(col)),
        paste0(toupper(substring(get(col), 1, 1)), substring(get(col), 2)),
        get(col)
      )]
    }

    return(as.data.frame(df))
  })
}
