#' Create a Sunburst Plot for Chemical Classifications
#'
#' This function takes a dataframe (typically from CANOPUS/SIRIUS) and generates
#' an interactive plotly sunburst chart based on NPC pathway and superclass.
#'
#' @param data A data frame containing chemical classification columns.
#' @param level1_col The column name for the top-level hierarchy (e.g., pathway). Default is "canopus.NPC.pathway".
#' @param level2_col The column name for the second-level hierarchy (e.g., superclass). Default is "canopus.NPC.superclass".
#' @param min_count Integer. Minimum number of occurrences for a Level 2 category to be included. Default is 10.
#' @param title Character. The title of the plot.
#'
#' @return A plotly object.
#' @import dplyr
#' @import plotly
#' @importFrom stringr str_replace_all str_to_title
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats setNames
#' @export
make_starburst <- function(data,
                                   level1_col = "canopus.NPC.pathway",
                                   level2_col = "canopus.NPC.superclass",
                                   min_count = 10,
                                   title = "Chemical Classification by Major Pathways") {

  # 1. Internal Helper for Label Cleaning
  clean_labels <- function(x) {
    x <- stringr::str_replace_all(x, "_", " ")
    x <- stringr::str_to_title(x)
    return(x)
  }

  # 2. Data Processing
  # We use .data[[col]] to handle string arguments for column names safely
  processed_data <- data %>%
    dplyr::select(Level1 = .data[[level1_col]],
                  Level2 = .data[[level2_col]]) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                ~ ifelse(is.na(.) | . == "" | . == "None", "Unclassified", .))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), clean_labels))

  # Filter based on user threshold
  filtered_data <- processed_data %>%
    dplyr::group_by(.data$Level1, .data$Level2) %>%
    dplyr::mutate(category_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$category_count >= min_count) %>%
    dplyr::filter(.data$Level1 != "Unclassified")

  # 3. Build Hierarchy Dataframe
  l1 <- filtered_data %>%
    dplyr::group_by(labels = .data$Level1) %>%
    dplyr::summarise(value = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(ids = .data$labels, parents = "")

  l2 <- filtered_data %>%
    dplyr::group_by(parents = .data$Level1, labels = .data$Level2) %>%
    dplyr::summarise(value = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(ids = paste0(.data$parents, " - ", .data$labels))

  hierarchy_df <- dplyr::bind_rows(l1, l2)

  # 4. Color Mapping
  unique_pathways <- unique(l1$labels)
  # Ensure we have enough colors; RColorBrewer Max for Set3 is 12
  n_colors <- length(unique_pathways)
  pal_size <- max(3, min(n_colors, 12))
  pathway_colors <- stats::setNames(
    RColorBrewer::brewer.pal(pal_size, "Set3")[1:n_colors],
    unique_pathways
  )

  # 5. Generate Plot
  p <- plotly::plot_ly(
    data = hierarchy_df,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~value,
    type = 'sunburst',
    branchvalues = 'total',
    marker = list(colors = pathway_colors)
  ) %>%
    plotly::layout(
      title = list(text = title, font = list(size = 18)),
      margin = list(l = 10, r = 10, b = 10, t = 60),
      paper_bgcolor = "white",
      font = list(family = "Arial", size = 14)
    )

  return(p)
}
