## 14. deduplicate_data
# Default values for gnps.data (but should work for others)
#' Deduplicate Data by Highest Score per Name
#'
#' Selects the highest scoring entry per unique compound name.
#'
#' @param data Data frame containing the data.
#' @param name_col Column representing the compound name (unquoted).
#' @param score_col Column representing the score used to rank entries (unquoted).
#'
#' @return A data frame with duplicates removed, keeping only the highest score per name.
#' @export
#'
#' @examples
#' \dontrun{
#' deduplicated <- deduplicate_data(gnps.data, compound.name, confidence.score)
#' }
deduplicate_data <- function(data, name_col, score_col) {
  data %>%
    arrange(desc({{score_col}})) %>%
    group_by({{name_col}}) %>%
    slice_head(n = 1) %>%
    ungroup()
}
