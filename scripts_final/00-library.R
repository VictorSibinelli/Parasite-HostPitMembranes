######################################################################
#
# Author: Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# Date: 13/07/2024
# Script 00 - Installing Required Packages for Analysis
#           + Package Version Control (Optional)
######################################################################

# ----------------------------
# Package Names (Alphabetically Ordered)
# ----------------------------
# List of required packages for the analysis in alphabetical order
packages <- c(
  "car", "data.table", "devtools", "DHARMa", "emmeans", "FactoMiner", "factoextra", 
  "ggbiplot", "ggnewscale", "ggpubr", "ggrepel", "ggstatsplot", "GGally", "groundhog", 
  "here", "htmlwidgets", "MASS", "matrixStats", "nlme", "performance", "plotly", "predictmeans", 
  "sjPlot", "stringr", "tidyverse", "viridis", "patchwork"
)

# ----------------------------
# Install Packages if Not Already Installed
# ----------------------------
# Check which packages are not installed and install them
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  # Installing missing packages
  install.packages(packages[!installed_packages])
}

# Confirmation message for package installation
print("Packages Installed Successfully")

# ----------------------------
# Optional: Control Package Versions with 'groundhog'
# ----------------------------
# Uncomment the following lines if you want to ensure specific versions
# for package consistency across sessions.

# groundhog.library(packages, "2024-11-01")

# ----------------------------
# Clean Environment
# ----------------------------
# Clear all objects from the environment to avoid conflicts
rm(list = ls())

# End of Script
######################################################################
