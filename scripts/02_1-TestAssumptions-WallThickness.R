######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 02.2 - Test assumptions test
######################################################################

# Load the data wrangling script
library(here)
source(here("scripts", "01-DataWrangling.R"))

# Remove all objects from the environment
rm(list = ls())

# Load data
wdata <- read.csv(here("data", "processed", "wdata.csv"))
source(here("scripts", "Functions.R"))
dataframes <- ls()

# Relevel the factors for each data frame
for (df_name in dataframes) {
  df <- get(df_name) # Get the data frame by name
  
  if ("ssp" %in% colnames(df)) { # Check if 'ssp' column exists
    df$ssp <- factor(df$ssp, levels = c(
      "Psittacanthus robustus", "Vochysia thyrsoidea",
      "Phoradendron perrotettii", "Tapirira guianensis",
      "Struthanthus rhynchophyllus", "Tipuana tipu",
      "Viscum album", "Populus nigra"
    ))
    assign(df_name, df) # Assign the modified data frame back to its name
  }
  rm(df, df_name) # Remove duplicated dataframe
}

# Variables of interest
variables <- c("wthickness")

# Perform Shapiro-Wilk test for each species and variable
species <- unique(wdata$ssp)

# Directory for saving plots and tables
output_dir_figs <- here("outputs", "figs", "assumptions")
output_dir_tables <- here("outputs", "tables", "assumptions")

# Initialize data frame to store outlier information
outlier_info <- data.frame(
  species = character(),
  wthickness = numeric(),
  stringsAsFactors = FALSE
)

# Set up the plotting area with 2 rows and 4 columns
png(filename = file.path(output_dir_figs, "WthicknessDensity.png"), 
    width = 6400, height = 4800, res = 600)  # Increase image size for bigger graphs

# Set the layout for multiple plots per page
par(mfrow = c(2, 4), oma = c(4, 4, 4, 2), mar = c(5, 5, 3, 2))  # Adjusted margins

# Iterate over species and variables to replace outliers with NA and check for normality
for (ssp in species) {
  outlier_counts <- numeric(length(variables))
  names(outlier_counts) <- variables
  
  for (var in variables) {
    # Use the revised version of `replace_outliers_test_normality`
    result <- replace_outliers_test_normality(wdata, var, ssp)
    wdata <- result$data
    
    # Plot density comparison using the original and updated density from the result
    plot_density_comparison(result$original_density, result$updated_density, ssp, var)
    
    # Store the number of points replaced with NA in the outlier count for this variable
    outlier_counts[var] <- result$outlier_count
  }
  
  # Add the counts to the outlier_info dataframe for each species
  outlier_info <- rbind(outlier_info, data.frame(
    species = ssp,
    wthickness = outlier_counts["wthickness"],
    stringsAsFactors = FALSE
  ))
}

# Close the PNG device after plotting
dev.off()

# Save the outlier information as a CSV file
write.csv(outlier_info, file = file.path(output_dir_tables, "w_outlier.csv"), row.names = FALSE)
print(outlier_info)

# Save the cleaned data after outlier replacement
fwrite(wdata, file = here("data", "processed", "wdata_clean.csv"))

# Homogeneity of variance test
h <- ggplot(wdata, aes(x = ssp, y = wthickness)) +
  geom_boxplot(na.rm = TRUE) +
  labs(
    title = "Boxplot of Wall Thickness Across Species",
    x = "Species", 
    y = "Wall Thickness"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
h

# Save graph
ggsave(here(output_dir_figs, "WthicknessVariance.png"), plot = h, dpi = 600, width = 10, height = 7)

# Homogeneity of variance found inside pairs
