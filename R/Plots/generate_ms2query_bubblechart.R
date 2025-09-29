#' Generate a Bubble Chart of Metabolite Classes by Sample (MS2Query Data)
#'
#' Processes MS2Query classification and quantification data to generate a bubble chart
#' showing the sum of peak areas for each NPC superclass across samples.
#'
#' @param ms2query.path Path to the MS2Query output CSV file.
#' @param quant.data.path Path to the MZMINE quantification CSV file.
#' @param metadata.path Path to the metadata CSV file ("HGM/A - Analysis.csv").
#' @param folder Output directory path for saving the plot.
#' @param threshold MS2Query model prediction threshold (default 0.63).
#' @param noise_area_threshold Peak area threshold to filter noise (default 10000).
#'
#' @return The long format data frame used for plotting (invisibly).
#' @importFrom dplyr filter select mutate rename group_by summarise
#' @importFrom tidyr pivot_longer
#' @importFrom readr read_csv
#' @importFrom ggplot2 ggplot geom_point scale_size scale_color_manual labs guides theme_minimal element_rect element_text scale_y_discrete ggsave
generate_ms2query_bubblechart <- function(ms2query.path, quant.data.path, folder,
                                          metadata.path = "HGM/A - Analysis.csv",
                                          threshold = 0.63, noise_area_threshold = 10000) {

  # --- Helper Function for Summarization ---
  summarise_superclasss <- function(j, merged_data) {
    superclasss <- unique(merged_data$npc_superclass_results)
    temp_data_frame <- data.frame()

    for (i in superclasss) {
      # Filter merged_data based on npc_superclass_results == i
      filt_data <- merged_data[merged_data$npc_superclass_results == i, c(1, j)]

      # Calculate sum of the filtered data
      total_sum <- sum(filt_data[, 2], na.rm = TRUE)

      temp_df <- data.frame(npc_superclass_results = i, total_sum = total_sum)
      colnames(temp_df) <- c("npc_superclass_results", colnames(merged_data)[j])
      temp_data_frame <- rbind(temp_data_frame, temp_df)
    }
    return(temp_data_frame)
  }

  # 1. Load and Prepare Classification Data
  class.data <- read.csv(ms2query.path) %>%
    dplyr::filter(ms2query_model_prediction >= threshold)
  names(class.data)[10] <- "row.ID"
  class.data <- class.data %>%
    dplyr::select("row.ID", "npc_superclass_results") %>%
    dplyr::filter(npc_superclass_results != "None")

  # Tidy class data
  class.data <- class.data %>%
    dplyr::mutate(npc_superclass_results = ifelse(npc_superclass_results == "", "Unclassified", npc_superclass_results)) %>%
    dplyr::filter(npc_superclass_results != "Unclassified")

  # 2. Load and Prepare Quantification Data
  quant.data <- read.csv(quant.data.path)
  quant.data <- quant.data[, -c(2:13)] # Removes unneccessary columns

  # Apply noise threshold
  quant.data[, 2:ncol(quant.data)][quant.data[, 2:ncol(quant.data)] < noise_area_threshold] <- 0

  # Tidy filenames (remove .mzML.Peak.area suffix)
  colnames_quant <- colnames(quant.data)
  colnames_quant <- sub("\\.mzML\\.Peak\\.area$", "", colnames_quant)
  colnames(quant.data) <- colnames_quant

  # 3. Merge and Tidy
  QUERY_TABLE <- merge(class.data, quant.data, by = "row.ID")

  # Remove the ID column (and blanks if needed, though they are kept here as per original script comments)
  merged_data <- dplyr::select(QUERY_TABLE, -row.ID)

  # 4. Summarize by Superclass and Sample

  # Initialize the merged dataframe with the result from the first iteration (col 2)
  merged_df <- summarise_superclasss(j = 2, merged_data)

  # Loop over remaining sample columns (3 to end)
  for (j in 3:ncol(merged_data)) {
    tdf <- summarise_superclasss(j, merged_data)
    merged_df <- merge(merged_df, tdf, by = "npc_superclass_results", all = TRUE)
  }

  # Change from wide to long df format
  long <- merged_df %>%
    tidyr::pivot_longer(
      cols = -npc_superclass_results,
      names_to = "Sample",
      values_to = "Perc"
    )

  # 5. Generate Bubble Chart

  # Filter the data to remove rows where Perc is not more than 0
  filtered_long <- long %>%
    dplyr::filter(Perc > 0)

  category_colors <- rainbow(length(unique(filtered_long$npc_superclass_results)))

  plot_obj <- ggplot2::ggplot(filtered_long, ggplot2::aes(x = Sample, y = npc_superclass_results, size = Perc, color = npc_superclass_results)) +
    ggplot2::geom_point() +
    ggplot2::scale_size(range = c(1, 8)) +
    ggplot2::scale_color_manual(values = category_colors) +
    ggplot2::labs(x = "Fraction", y = "Metabolite Class", size = "Sum of Peak Areas") +
    ggplot2::guides(color = "none") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 16)
    ) +
    # Ensure y-axis order is reversed (highest classes at top)
    ggplot2::scale_y_discrete(limits = rev(unique(filtered_long$npc_superclass_results)))

  print(plot_obj)

  # Save the bubble chart
  ggplot2::ggsave(
    filename = paste0(folder, "/bubblechart.svg"),
    plot = plot_obj,
    width = 50,
    height = 40,
    units = "cm"
  )

  invisible(long)
}
