#' Append annotation data frames with automatic column alignment
#'
#' This function appends a new annotation data frame to an existing
#' annotation table. If the new data contains columns that are not
#' present in the existing annotations, those columns are added to
#' the existing annotations and filled with \code{NA}.
#'
#' @param data.to.append A data frame of new annotations to be appended.
#' @param existing.annotations A data frame of existing annotations.
#'
#' @return A combined data frame containing both the new and existing
#' annotations with aligned columns.
#'
#' @examples
#' \dontrun{
#' new <- data.frame(feature.ID = 1:2, compound.name = c("A", "B"))
#' old <- data.frame(feature.ID = 3, compound.name = "C", confidence.score = 0.95)
#' combined <- append_annotations(new, old)
#' }
#'
#' @export
append_annotations <- function(data.to.append, existing.annotations) {
  # Identify columns present in new data but not in existing annotations
  missing_cols <- colnames(data.to.append)[
    !colnames(data.to.append) %in% colnames(existing.annotations)
  ]

  # Add missing columns to existing annotations and fill with NA
  for (col in missing_cols) {
    existing.annotations[[col]] <- NA
  }

  # Bind rows while preserving column order
  combined_annotations <- dplyr::bind_rows(data.to.append, existing.annotations)

  return(combined_annotations)
}
