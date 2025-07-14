## 14. deduplicate_data
# Default values for gnps.data (but should work for others)
deduplicate_data <- function(data, name_col, score_col) {
  data %>%
    arrange(desc({{score_col}})) %>%
    group_by({{name_col}}) %>%
    slice_head(n = 1) %>%
    ungroup()
}