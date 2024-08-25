# Define the list of packages to install
packages <- c("dplyr", "readxl", "writexl", "foreach", "doParallel", 
              "optparse", "openxlsx", "stringr", "readr")

# Function to check and install packages
install_if_needed <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, repos = 'https://cloud.r-project.org')
    }
  }
}

# Call the function with the list of packages
install_if_needed(packages)
