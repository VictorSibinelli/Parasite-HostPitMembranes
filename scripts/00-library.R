######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Scrpt 00 - Instaling and loading packages needed in the analysis + packages version control
#################################################################

# Package names
packages <- c("tidyverse", "here", "plotly", "data.table",
              "readr", "groundhog", "car", "htmlwidgets",
              "lattice","lme4","nlme")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

#########################################################################
# If you desire to run the analysis with the same package versions, remove the comments in the lines below and run it

# groundhog.library(packages, "2024-07-01")