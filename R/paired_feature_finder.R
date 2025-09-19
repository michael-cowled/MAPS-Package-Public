#' Find Paired Features in GNPS Cluster
#'
#' Given a GNPS cluster ID, returns a vector of paired feature IDs (i.e., directly connected nodes)
#' sorted by cosine similarity (descending).
#'
#' @param ID A numeric or character GNPS cluster ID to look up.
#' @param gnps.cluster.pairs A data frame containing GNPS cluster pairs, with
#'   columns \code{CLUSTERID1}, \code{CLUSTERID2}, and \code{Cosine}.
#'
#' @return A vector of paired GNPS cluster IDs connected to the input \code{ID}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming gnps.cluster.pairs is a data frame in your environment
#' paired <- paired_feature_finder(332, gnps.cluster.pairs)
#' }
paired_feature_finder <- function(ID, gnps.cluster.pairs) {
  filtered_pairs <- gnps.cluster.pairs %>%
    dplyr::filter(CLUSTERID1 == ID | CLUSTERID2 == ID) %>%
    dplyr::arrange(desc(Cosine)) # Sort by decreasing cosine similarity

  paired_values <- filtered_pairs %>%
    dplyr::mutate(paired_value = ifelse(CLUSTERID1 == ID, CLUSTERID2, CLUSTERID1)) %>%
    dplyr::pull(paired_value)

  return(paired_values)
}
