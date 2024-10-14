######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Scrpt 02 - Test assumptions test
#################################################################
# load packages
library(here)
source(here("scripts", "01-DataWrangling.R"))
rm(list=ls())
# Load data
pitO <- read.csv(here("data", "processed", "PitOdata.csv"))
pitdata <- read.csv(here("data", "processed", "pitdata.csv"))
source(here("scripts","Functions.R"))
# Relevel the factors for each data frame
relevel_factors(ls())

# Variables of interest for the first dataframe
variables1 <- c("pitavg", "pcavg", "peavg", "pcd")
species <- unique(pitdata$ssp)

# Initialize data frame to store outlier information for the first dataframe
outlier_info1 <- data.frame(
  species = species,
  pitavg = numeric(length(species)),
  pcavg = numeric(length(species)),
  peavg = numeric(length(species)),
  pcd = numeric(length(species)),
  stringsAsFactors = FALSE
)

# Perform Shapiro-Wilk test for each species and variable in the first dataframe
pitresults1 <- perform_normality_test(pitdata, species, variables1)
print(pitresults1)

# Directory for saving plots and table
output_dir_figs <- here("outputs", "figs", "assumptions")
output_dir_tables <- here("outputs", "tables", "assumptions")


## Iterate over variables to replace outliers with NA and check for normality in the first dataframe
for (var in variables1) {
  # Set up the PNG file for output
  png(filename = file.path(output_dir_figs, paste0(var, "_density.png")),
      width = 8000, height = 6000, res = 600)
  
  # Set up the plotting area
  par(mfrow = c(2, 4), oma = c(4, 4, 2, 2), mar = c(4, 4, 2, 1))
  
  # Iterate over each species
  for (ssp in species) {
    result <- replace_outliers_test_normality(pitdata, var, ssp)
    pitdata <- result$data
    
    # Plot density comparison
    p <- plot_density_comparison(result$original_density, result$updated_density, ssp, var)
    p
  
    # Update outlier_info1 with the count of replaced data points
    outlier_info1[outlier_info1$species == ssp, var] <- result$outlier_count
  }
  
  # Add variable name as a title for the entire page
  mtext(var, side = 3, outer = TRUE, cex = 1.5, line = 1)
  
  # Close the PNG device
  dev.off()
}

# Print outlier info for first dataframe
print(outlier_info1)

# Perform analysis for the second dataframe
species <- unique(pitO$ssp)
variables2 <- c("PitOpening", "PitDiameter")

# Initialize data frame to store outlier information for the second dataframe
outlier_info2 <- data.frame(
  species = species,
  PitOpening = numeric(length(species)),
  PitDiameter = numeric(length(species)),
  stringsAsFactors = FALSE
)

# Perform Shapiro-Wilk test for each species and variable in the second dataframe
pitresults2 <- perform_normality_test(pitO, species, variables2)
print(pitresults2)

# Iterate over variables to replace outliers with NA and check for normality in the second dataframe
for (var in variables2) {
  png(filename = file.path(output_dir_figs, paste0(var, "_density.png")),
      width = 8000, height = 6000, res=600)
  par(mfrow = c(2, 4), oma = c(4, 4, 2, 2), mar = c(4, 4, 2, 1))
  
  for (ssp in species) {
    result <- replace_outliers_test_normality(pitO, var, ssp)
    pitO <- result$data
    
    plot_density_comparison(result$original_density, result$updated_density, ssp,var)
    
    # Update outlier_info2 with the count of replaced data points
    outlier_info2[outlier_info2$species == ssp, var] <- result$outlier_count
  }
  
  mtext(var, side = 3, outer = TRUE, cex = 1.5, line = 1)
  dev.off()
}

# Print outlier info for second dataframe
print(outlier_info2)

# Merge and save outlier information
outlier_info_combined <- merge(outlier_info1, outlier_info2, by = "species")
print(outlier_info_combined)
write.csv(outlier_info_combined, file = file.path(output_dir_tables, "pit_outlier.csv"), row.names = FALSE)

# Save cleaned data
fwrite(pitdata, file = here("data", "processed", "pitdata_clean.csv"))
fwrite(pitO, file = here("data", "processed", "pitO_clean.csv"))

rm(list=ls())

# homogeneity of variance test
###############################################


##Re-upload original data , can be done with cleanes data
pitO <- read.csv(here("data", "processed", "PitOdata.csv"))
pitdata <- read.csv(here("data", "processed", "pitdata_clean.csv"))
source(here("scripts","Functions.R"))
# Relevel the factors for each data frame
relevel_factors(ls())

variables1 <- c("pitavg", "pcavg", "peavg", "pcd")
variables2 <- c("PitOpening", "PitDiameter")

for (var in variables1) {
  p <- ggplot(data = pitdata, aes_string(x = "ssp", y = var)) +
    geom_boxplot(na.rm = TRUE) +
    labs(
      title = paste("Boxplot of", var, "variance across Species"),
      x = "Species",
      y = var
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  print(p)
  plot_name <- paste(var,"_Variance.png",sep="")
  ggsave(here("outputs","figs","assumptions",plot_name), plot = p, dpi = 600, width = 10, height = 7)
}

for (var in variables2) {
  p <- ggplot(data = pitO, aes_string(x = "ssp", y = var)) +
    geom_boxplot(na.rm = TRUE) +
    labs(
      title = paste("Boxplot of", var, "variance across Species"),
      x = "Species",
      y = var
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  print(p)
  plot_name <- paste(var,"_Variance.png",sep="")
  ggsave(here("outputs","figs","assumptions",plot_name), plot = p, dpi = 600, width = 10, height = 7)
}
cat("Wall Thickness Test assumptions Complete
    \nSummary:
    \n1. Small deviation from normality for most of traits
    \n2. Aprox. homocesdasticity observed in pitAVG and pitOpening
    \n3. Heteroscedasticity observed in pcd and pit DIameter")

