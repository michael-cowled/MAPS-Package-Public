#' Get CID from PubChem with Fallbacks and Caching
#'
#' Attempts to resolve a PubChem CID using compound name, SMILES, and synonyms.
#' If successful, returns a list of annotation details, and updates a global cache.
#'
#' @param name A character string representing the compound name.
#' @param smiles A character string representing the compound SMILES (optional).
#'
#' @return A named list with keys: \code{CID}, \code{ResolvedName}, \code{SMILES},
#'   \code{MolecularFormula}, and \code{MonoisotopicMass}, or \code{NULL} if all lookups fail.
#' @export
#'
#' @details This function expects a global object named \code{cid_cache_df} to exist,
#' and it will be updated via \code{<<-} with new entries when successful lookups are made.
#'
#' @examples
#' \dontrun{
#' get_cid_with_fallbacks("glucose")
#' get_cid_with_fallbacks("some compound", smiles = "CCO")
#' }
get_cid_with_fallbacks <- function(name, smiles = NA) {
  # Check cache: by name
  name_match <- cid_cache_df[!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name, ]

  # If not found, try matching on SMILES
  if (nrow(name_match) == 0 && !is.na(smiles) && smiles != "") {
    name_match <- cid_cache_df[!is.na(cid_cache_df$SMILES) & cid_cache_df$SMILES == smiles, ]
  }

  # Return cached result if available
  if (nrow(name_match) > 0 && !is.na(name_match$CID[1])) {
    message(paste("  [CACHE] CID found for:", name))
    return(list(
      CID = name_match$CID[1],
      ResolvedName = name_match$ResolvedName[1],
      SMILES = name_match$SMILES[1],
      MolecularFormula = name_match$MolecularFormula[1],
      MonoisotopicMass = name_match$MonoisotopicMass[1]
    ))
  }

  # Step 1: Try name
  cid <- get_pubchem(name, "name", "cids")

  # Step 2: Try SMILES
  if (is.na(cid) && !is.na(smiles) && smiles != "") {
    message(paste("  Name lookup failed. Trying SMILES:", smiles))
    cid <- get_pubchem(smiles, "smiles", "cids")
  }

  # Step 3: Try synonym
  if (is.na(cid)) {
    message(paste("  Name and SMILES failed. Trying synonym search for:", name))
    cid <- get_pubchem(name, "synonym", "cids")
  }

  # Step 4: If CID found, get properties and title
  if (!is.na(cid)) {
    props <- get_pubchem(cid, "cid", "properties")
    title <- get_pubchem(cid, "cid", "title")

    if (!is.list(props)) props <- list()
    if (!is.character(title)) title <- NA

    new_row <- data.frame(
      LookupName = name,
      ResolvedName = if (!is.null(title) && !is.na(title)) title else name,
      SMILES = if (!is.null(props$SMILES)) props$SMILES else if (!is.na(smiles)) smiles else NA,
      CID = cid,
      MolecularFormula = if (!is.null(props$MolecularFormula)) props$MolecularFormula else NA,
      MonoisotopicMass = if (!is.null(props$MonoisotopicMass)) props$MonoisotopicMass else NA,
      stringsAsFactors = FALSE
    )

    cid_cache_df <<- dplyr::bind_rows(cid_cache_df, new_row)
    return(as.list(new_row[1, ]))
  }

  return(NULL)
}
