## 1. check_and_install
# Function to check and install missing packages, then load them
check_and_install <- function(packages) {
  to_install <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(to_install) > 0) {
    message(paste("Installing:", paste(to_install, collapse = ", ")))
    install.packages(to_install, dependencies = TRUE)
    library(packages)
  }
  # Always load the packages, regardless of whether they were just installed
  message(paste("Loading:", paste(packages, collapse = ", ")))
  invisible(lapply(packages, library, character.only = TRUE))
}