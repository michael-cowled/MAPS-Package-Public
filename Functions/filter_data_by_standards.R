## 13. filter_data_by_standards
filter_data_by_standards <- function(
    data,              # either gnps.data or csi.data
    mzmine.data,       # standards data frame (same as before)
    compound_col,      # column name for compound ID or CID in 'data'
    feature_id_col = "feature.ID",
    rt_tolerance = 0.1
) {
  library(dplyr)
  
  # Step 1: Filter authentic standards from mzmine.data
  mzmine.standards <- filter(mzmine.data, !is.na(compound.name))
  
  # Step 2: Unique CIDs to check (from standards)
  mzmine.to.check <- unique(mzmine.standards$CID)
  
  # Step 3: Join RT from mzmine.data, ensure numeric
  data.cleaned <- data %>%
    left_join(mzmine.data %>% select(all_of(feature_id_col), rt), by = feature_id_col) %>%
    mutate(rt = suppressWarnings(as.numeric(as.character(rt))))
  
  # Step 4: Loop through each CID in standards to filter data.cleaned
  for (i in mzmine.to.check) {
    temp <- filter(data.cleaned, .data[[compound_col]] == i)
    if (nrow(temp) == 0) next
    
    temp2 <- filter(mzmine.standards, CID == i)
    temp3 <- filter(temp, .data[[feature_id_col]] != temp2$feature.ID[1])  # standard feature.ID
    
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
        data.cleaned <- filter(data.cleaned, .data[[feature_id_col]] != j)
      }
    }
    
    # Remove remaining entries within tolerance but not standard
    temp <- filter(data.cleaned, .data[[compound_col]] == i)
    temp3 <- filter(temp, .data[[feature_id_col]] != temp2$feature.ID[1])
    
    rt.others <- vapply(temp3[[feature_id_col]], function(fid) {
      rt.val <- data.cleaned$rt[data.cleaned[[feature_id_col]] == fid]
      if (length(rt.val) == 0 || all(is.na(rt.val))) return(NA_real_)
      suppressWarnings(as.numeric(as.character(rt.val[1])))
    }, numeric(1))
    
    valid_idx <- !is.na(rt.others)
    
    close_match_ids <- temp3[[feature_id_col]][valid_idx & abs(rt.std - rt.others[valid_idx]) <= rt_tolerance]
    
    data.cleaned <- filter(data.cleaned, !(.data[[feature_id_col]] %in% close_match_ids))
  }
  
  return(data.cleaned)
}