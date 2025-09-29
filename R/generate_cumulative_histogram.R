#' Generate a Cumulative Histogram of Unique Annotations Across Fractions
#'
#' Calculates the cumulative number of unique annotations discovered across ordered
#' fractions and creates a histogram.
#'
#' @param Annotations.with.samples Data frame from the barchart function (feature.ID, Samples, annotation).
#' @param folder Output directory path for saving the plot.
#' @param remove_sample Sample name to optionally exclude from the calculation (e.g., crude).
#'
#' @return The data frame containing the cumulative annotation count (invisibly).
#' @importFrom dplyr filter summarise n_distinct last add_row full_join
#' @importFrom tidyr fill
#' @importFrom stringr str_extract str_split str_detect
#' @importFrom ggplot2 ggplot geom_col labs theme element_text ggsave
generate_cumulative_histogram <- function(Annotations.with.samples, folder, remove_sample = "HGMA_1598") {

  # Optional Filtering (as suggested in original script comments)
  # Annotations.with.samples <- Annotations.with.samples %>%
  #   dplyr::filter(!stringr::str_detect(Samples, remove_sample))

  # Initialize an empty dataframe to store the results
  cumulative_annotations <- data.frame()

  # Get a sorted list of unique sample numbers
  sample_numbers <- sort(unique(as.numeric(stringr::str_extract(unlist(stringr::str_split(Annotations.with.samples$Samples, "; ")), "\\d+"))))

  # Loop through each sample number
  for (i in seq_along(sample_numbers)) {

    # Filter data for the current sample number
    # NOTE: This uses the sample number appended with "_" (e.g., _01, _02)
    current_sample_data <- Annotations.with.samples %>%
      dplyr::filter(stringr::str_detect(Samples, paste0("_", sample_numbers[i])))

    # Calculate the number of NEW annotations in the current sample
    # The original script's logic here was flawed as it recalculated all annotations
    # based on the remaining features. We rely on the filtering done in the loop.
    new_annotations <- sum(current_sample_data$annotation)

    # Add the new annotations to the cumulative count
    if (nrow(cumulative_annotations) == 0) {
      cumulative_annotations <- data.frame(
        Sample.Number = i,
        Cumulative.Unique.Annotations = new_annotations
      )
    } else {
      cumulative_annotations <- cumulative_annotations %>%
        tibble::add_row(
          Sample.Number = i,
          Cumulative.Unique.Annotations = new_annotations + dplyr::last(cumulative_annotations$Cumulative.Unique.Annotations)
        )
    }

    # Remove the processed features from the main dataframe for the next iteration
    # This simulates finding the features and removing them from the pool of 'available' features.
    Annotations.with.samples <- Annotations.with.samples %>%
      dplyr::filter(!feature.ID %in% current_sample_data$feature.ID)
  }

  # Ensure all chronological sample numbers are present
  all_sample_numbers <- data.frame(Sample.Number = seq_along(sample_numbers))
  cumulative_annotations <- dplyr::full_join(all_sample_numbers, cumulative_annotations, by = "Sample.Number")

  # Fill missing Cumulative.Unique.Annotations with the previous value
  cumulative_annotations <- cumulative_annotations %>%
    tidyr::fill(Cumulative.Unique.Annotations, .direction = "down")

  # Print the result
  print(cumulative_annotations)

  # Create the cumulative histogram
  plot_obj3 <- ggplot2::ggplot(cumulative_annotations, ggplot2::aes(x = Sample.Number, y = Cumulative.Unique.Annotations)) +
    ggplot2::geom_col() +
    ggplot2::labs(
      x = "Fraction Number",
      y = "Cumulative Number of Unique Annotations"
    ) +
    ggplot2::theme(text = ggplot2::element_text(size = 24))

  print(plot_obj3)

  # Save the histogram
  ggplot2::ggsave(
    filename = paste0(folder, "/cumulative-histogram.svg"),
    plot = plot_obj3,
    width = 50,
    height = 40,
    units = "cm"
  )

  invisible(cumulative_annotations)
}
