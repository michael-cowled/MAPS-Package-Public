#' @title Check and Install Packages
#' @description An internal helper function to check for, install, and load R packages from CRAN or GitHub.
#'
#' @param packages A character vector of package names to check and install from CRAN.
#' @param github_packages A named list where keys are package names and values are their GitHub repositories.
#'
#' @keywords internal
#' @noRd
#'
check_and_install <- function(packages, github_packages = list()) {
  # Install 'remotes' if needed for GitHub installs
  if (!requireNamespace("remotes", quietly = TRUE)) {
    utils::install.packages("remotes")
  }

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% names(github_packages)) {
        message(paste("Installing", pkg, "from GitHub:", github_packages[[pkg]]))
        remotes::install_github(github_packages[[pkg]])
      } else {
        message(paste("Installing", pkg, "from CRAN"))
        utils::install.packages(pkg, dependencies = TRUE)
      }
    }
    suppressPackageStartupMessages(base::library(pkg, character.only = TRUE))
  }
}
