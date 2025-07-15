#' Filter Data Based on Authentic Standards Using Retention Time
#'
#' Filters a dataset to exclude features that do not closely match authentic standards'
#' retention times within a specified tolerance.
#'
#' @param data Data frame to filter, e.g. GNPS or CSI data.
#' @param mzmine.data Data frame of authentic standards, with columns including `compound.name`, `CID`, `feature.ID`, and `rt`.
#' @param compound_col Character string. Column name in \code{data} containing compound IDs (e.g. CID).
#' @param feature_id_col Character string. Column name identifying features, defaults to `"feature.ID"`.
#' @param rt_tolerance Numeric. Maximum allowed retention time difference to keep features, default 0.1.
#'
#' @return Filtered data frame with features excluded that do not match authentic standards within the RT tolerance.
#' @export
#'
#' @examples
#' \dontrun{
#' filtered_data <- filter_data_by_standards(
#'   data = gnps.data,
#'   mzmine.data = mzmine_standards,
#'   compound_col = "CID"
#' )
#' }
filter_data_by_standards <- function(
    data,
    mzmine.data,
    compound_col,
    feature_id_col = "feature.ID",
    rt_tolerance = 0.1
) {
  # Step 1: Filter authentic standards from mzmine.data
  mzmine.standards <- dplyr::filter(mzmine.data, !is.na(compound.name))

  # Step 2: Unique CIDs to check (from standards)
  mzmine.to.check <- unique(mzmine.standards$CID)

  # Step 3: Join RT from mzmine.data, ensure numeric
  data.cleaned <- dplyr::left_join(
    data,
    dplyr::select(mzmine.data, dplyr::all_of(feature_id_col), rt),
    by = feature_id_col
  ) %>%
    dplyr::mutate(rt = suppressWarnings(as.numeric(as.character(rt))))

  # Step 4: Loop through each CID in standards to filter data.cleaned
  for (i in mzmine.to.check) {
    temp <- dplyr::filter(data.cleaned, .data[[compound_col]] == i)
    if (nrow(temp) == 0) next

    temp2 <- dplyr::filter(mzmine.standards, CID == i)
    temp3 <- dplyr::filter(temp, .data[[feature_id_col]] != temp2$feature.ID[1])

    rt.std.raw <- temp2$rt[1]
    rt.std <- suppressWarnings(as.numeric(as.character(rt.std.raw)))

    if (length(rt.std) == 0 || is.na(rt.std)) {
      warning(paste("Skipping CID", i, "- invalid standard rt:", rt.std.raw))
      next
    }

    # Remove entries with RT difference > rt_tolerance
    for (j in temp3[[feature_id_col]]) {
      rt.j.raw <- data.cleaned$rt[data.cleaned[[feature_id_col]] == j]
      rt.j <- suppressWarnings(as.numeric(as.character(rt.j.raw)))
      if (length(rt.j) == 0 || all(is.na(rt.j))) next

      if (abs(rt.std - rt.j[1]) > rt_tolerance) {
        data.cleaned <- dplyr::filter(data.cleaned, .data[[feature_id_col]] != j)
      }
    }

    # Remove remaining entries within tolerance but not the standard itself
    temp <- dplyr::filter(data.cleaned, .data[[compound_col]] == i)
    temp3 <- dplyr::filter(temp, .data[[feature_id_col]] != temp2$feature.ID[1])

    rt.others <- vapply(temp3[[feature_id_col]], function(fid) {
      rt.val <- data.cleaned$rt[data.cleaned[[feature_id_col]] == fid]
      if (length(rt.val) == 0 || all(is.na(rt.val))) return(NA_real_)
      suppressWarnings(as.numeric(as.character(rt.val[1])))
    }, numeric(1))

    valid_idx <- !is.na(rt.others)

    close_match_ids <- temp3[[feature_id_col]][valid_idx & abs(rt.std - rt.others[valid_idx]) <= rt_tolerance]

    data.cleaned <- dplyr::filter(data.cleaned, !(.data[[feature_id_col]] %in% close_match_ids))
  }

  return(data.cleaned)
}
