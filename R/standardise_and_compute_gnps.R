#' Standardise and compute GNPS annotations
#'
#' This function standardises GNPS Level 2 high-confidence annotations
#' by resolving SMILES to PubChem CIDs, filtering duplicates, and computing
#' identification probabilities.
#'
#' @param gnps.data.lv2.high.conf A data frame of GNPS Level 2 high-confidence annotations.
#' @param cid_cache_df A data frame cache for storing SMILES-to-CID lookups.
#' @param lipids.file Path to lipid reference file used for annotation standardisation.
#' @param cid_database_path Path to the SQLite database containing CID lookups.
#' @param ms2query.prob Minimum probability threshold for ID assignment.
#'
#' @return A list with:
#' \describe{
#'   \item{data}{A data frame of standardised GNPS annotations with computed ID probabilities.}
#'   \item{cache}{The updated CID cache data frame.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- standardise_and_compute_gnps(
#'   gnps.data.lv2.high.conf,
#'   cid_cache_df,
#'   "lipids.csv",
#'   "cid_lookup.sqlite",
#'   ms2query.prob = 0.7
#' )
#' }
#'
#' @export
standardise_and_compute_gnps <- function(
    gnps.data.lv2.high.conf,
    cid_cache_df,
    lipids.file,
    cid_database_path,
    ms2query.prob
) {

  # Pre-standardisation filtering for duplicates
  gnps.data.lv2.high.conf <- deduplicate_data(
    gnps.data.lv2.high.conf,
    compound.name,
    confidence.score
  )

  # Standardize annotations
  result <- standardise_annotation(
    gnps.data.lv2.high.conf,
    name_col = "compound.name",
    smiles_col = "smiles",
    cid_cache_df = cid_cache_df,
    lipids.file = lipids.file,
    cid_database_path = cid_database_path
  )

  gnps.data.lv2.high.conf <- result$data
  cid_cache_df <- result$cache

  # Add new columns if missing
  if (!("Formula" %in% names(gnps.data.lv2.high.conf)))
    gnps.data.lv2.high.conf$Formula <- NA_character_
  if (!("IUPAC" %in% names(gnps.data.lv2.high.conf)))
    gnps.data.lv2.high.conf$IUPAC <- NA_character_
  if (!("Monoisotopic.Mass" %in% names(gnps.data.lv2.high.conf)))
    gnps.data.lv2.high.conf$Monoisotopic.Mass <- NA_real_

  # Post-standardisation filtering of duplicates
  gnps.data.lv2.high.conf <- deduplicate_data(
    gnps.data.lv2.high.conf,
    compound.name,
    confidence.score
  )

  # Compute ID probability
  gnps.data.lv2.high.conf <- compute_id_prob(
    gnps.data.lv2.high.conf,
    "confidence.score",
    ms2query.prob
  )

  return(list(
    data = gnps.data.lv2.high.conf,
    cache = cid_cache_df
  ))
}
