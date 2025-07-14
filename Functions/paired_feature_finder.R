## 7. paired_feature_finder
# Function: Finds the features linked by an edge in a GNPS cluster
paired_feature_finder <- function(ID) {
  filtered_pairs <- gnps.cluster.pairs %>%
    filter(CLUSTERID1 == ID | CLUSTERID2 == ID) %>%
    arrange(desc(Cosine)) # Sort by decreasing Cosine
  
  # Extract paired values
  paired_values <- filtered_pairs %>%
    mutate(paired_value = ifelse(CLUSTERID1 == ID, CLUSTERID2, CLUSTERID1)) %>%
    pull(paired_value)
  
  return(paired_values)
}