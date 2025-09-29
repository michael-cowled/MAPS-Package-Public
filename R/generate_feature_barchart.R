#' Generate a Barchart of Annotations per Sample/Fraction
#'
#' Counts the frequency of annotations (including propagated Ms/MS) across fractions
#' and creates a barchart. This function works with or without a metadata file.
#'
#' @param final.annotation.df Data frame containing feature annotations (from Annotation-table.R).
#' @param folder Output directory path for saving the CSV and plot.
#' @param metadata.path Optional path to the metadata CSV file (e.g., "HGM/A - Analysis.csv"). Defaults to NULL.
#'
#' @return The final merged data frame used for plotting (invisibly).
#' @importFrom dplyr select mutate group_by summarise full_join filter coalesce
#' @importFrom tidyr separate_rows
#' @importFrom readr write_csv
#' @importFrom ggplot2 ggplot geom_col aes xlab ylab scale_fill_manual theme element_text ylim ggsave
generate_feature_barchart <- function(final.annotation.df, folder, metadata.path = NULL) {

  # 1. Preparing Annotations.with.samples
  Samples <- final.annotation.df %>%
    dplyr::select(feature.ID, Samples, compound.name)
  names(Samples)[names(Samples) == "compound.name"] <- "annotation"

  # Convert 'annotation' column to presence/absence (0/1)
  Annotations.with.samples <- Samples %>%
    dplyr::mutate(annotation = ifelse(is.na(annotation), 0, 1))

  # 2. Summarising annotations per sample

  # Sum of all annotations per sample
  sum_annotations_per_sample <- Annotations.with.samples %>%
    tidyr::separate_rows(Samples, sep = "; ") %>%
    dplyr::group_by(Samples) %>%
    dplyr::summarise(sum_annotations = sum(annotation, na.rm = TRUE))

  # Sum of unique annotations per sample
  sum_unique_annotations_per_sample <- Annotations.with.samples %>%
    dplyr::mutate(Samples_count = stringr::str_count(Samples, ";") + 1) %>%
    dplyr::filter(Samples_count == 1) %>%
    tidyr::separate_rows(Samples, sep = "; ") %>%
    dplyr::group_by(Samples) %>%
    dplyr::summarise(sum_unique_annotations = sum(annotation, na.rm = TRUE))

  # 3. Merging annotations and calculating total features
  merged_df <- sum_annotations_per_sample %>%
    dplyr::full_join(sum_unique_annotations_per_sample, by = "Samples")

  # Calculate total number of features per sample
  total_features_per_sample <- Annotations.with.samples %>%
    tidyr::separate_rows(Samples, sep = "; ") %>%
    dplyr::group_by(Samples) %>%
    dplyr::summarise(total_features = n())

  # Merge the total features column into the final dataframe
  merged_df <- merged_df %>%
    dplyr::full_join(total_features_per_sample, by = "Samples")

  # 4. CONDITIONAL Mapping to metadata

  # Initialize the Sample column name for merging
  names(merged_df)[names(merged_df) == "Samples"] <- "Sample"

  if (!is.null(metadata.path) && file.exists(metadata.path)) {
    message("Metadata path found. Joining metadata for sample renaming.")

    # Use tryCatch to gracefully handle read errors
    metadata <- tryCatch({
      read.csv(metadata.path)
    }, error = function(e) {
      warning("Failed to read metadata file: ", e$message)
      return(NULL)
    })

    if (!is.null(metadata)) {
      # Proceed with merging only if metadata was read successfully
      names(metadata)[1] <- "Sample"

      # Check if Book.Code exists before selecting/coalescing
      if ("Book.Code" %in% names(metadata)) {
        metadata <- dplyr::select(metadata, Sample, Book.Code)

        merged_df <- merged_df %>%
          dplyr::full_join(metadata, by = "Sample") %>%
          dplyr::mutate(Sample = coalesce(Book.Code, Sample)) %>%
          dplyr::select(-Book.Code)
      } else {
        warning("Metadata file is missing 'Book.Code' column. Skipping sample renaming.")
      }
    }
  } else {
    message("No valid metadata path provided. Using raw sample names for plotting.")
  }

  # Ensure the final DF is clean and only contains rows with annotation data
  merged_df <- merged_df %>%
    dplyr::filter(!is.na(sum_annotations))

  # Save the resulting table
  readr::write_csv(merged_df, paste0(folder, "/features-per-sample.csv"))

  # 5. Barchart Generation
  category_colors <- rainbow(length(unique(merged_df$Sample)), alpha = 0.7)

  plot_obj2 <- ggplot2::ggplot(merged_df, ggplot2::aes(x = Sample, y = sum_annotations, fill = Sample)) +
    ggplot2::geom_col() +
    ggplot2::xlab("Sample") +
    ggplot2::ylab("Sum of Annotations") +
    ggplot2::scale_fill_manual(values = category_colors) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      text = ggplot2::element_text(size = 24),
      legend.position = "none"
    ) +
    ggplot2::ylim(0, 1000)

  print(plot_obj2)

  ggplot2::ggsave(
    filename = paste0(folder, "/feature-barchart.svg"),
    plot = plot_obj2,
    width = 50,
    height = 40,
    units = "cm"
  )

  invisible(merged_df)
}
