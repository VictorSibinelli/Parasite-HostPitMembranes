######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Scrpt 02 - Test assumptions test
#################################################################
# load packages
library(here)
source(here("scripts", "01-DataWrangling.R"))

# Load data
pitO <- read.csv(here("data", "processed", "PitOdata.csv"))
pitdata <- read.csv(here("data", "processed", "pitdata.csv"))

# Function to plot density before and after outlier removal on the same graph
plot_density_comparison <- function(original_density, updated_density, variable) {
  plot(original_density, main = variable, 
       xlab = variable, ylab = "Density", lwd = 2, col = "red", xlim = range(c(original_density$x, updated_density$x), na.rm = TRUE))
  lines(updated_density, lwd = 2, col = "blue")
  legend("topright", legend = c("Before", "After"), col = c("red", "blue"), lwd = 2)
}

# Function to replace outliers with NA and test for normality
replace_outliers_test_normality <- function(data, variable, species, max_iter = 10) {
  outlier_count <- 0
  removed_outliers <- numeric()
  
  original_data <- data[[variable]][data$ssp == species]
  shapiro_p <- NA
  
  for (i in 1:max_iter) {
    outliers <- boxplot.stats(data[[variable]][data$ssp == species])$out
    if (length(outliers) == 0) break
    x1 <- outliers[which.max(abs(outliers - median(data[[variable]][data$ssp == species], na.rm = TRUE)))]
    data[[variable]][data[[variable]] == x1 & data$ssp == species] <- NA
    outlier_count <- outlier_count + 1
    removed_outliers <- c(removed_outliers, x1)
    shapiro_test <- shapiro.test(data[[variable]][data$ssp == species])
    shapiro_p <- shapiro_test$p.value
    if (shapiro_p > 0.05) break
  }
  
  original_density <- density(original_data, na.rm = TRUE)
  updated_density <- density(data[[variable]][data$ssp == species], na.rm = TRUE)
  
  return(list(data = data, outlier_count = outlier_count, shapiro_p = shapiro_p, removed_outliers = removed_outliers,
              original_density = original_density, updated_density = updated_density))
}

# Function to perform normality test across species and variables
perform_normality_test <- function(data, species_list, variables) {
  results <- data.frame(species = character(), stringsAsFactors = FALSE)
  for (ssp in species_list) {
    shapiro_results <- sapply(variables, function(var) {
      shapiro.test(data[[var]][data$ssp == ssp])$p.value
    })
    results <- rbind(results, c(ssp, shapiro_results))
  }
  colnames(results) <- c("species", variables)
  return(results)
}

# Variables of interest
variables <- c("pitavg", "pcavg", "peavg", "pcd")

# Perform Shapiro-Wilk test for each species and variable
species <- unique(pitdata$ssp)
pitresults <- perform_normality_test(pitdata, species, variables)

# Print initial normality test results
print(pitresults)

# Directory for saving plots and table
output_dir_figs <- here("outputs", "figs", "assumptions")
output_dir_tables <- here("outputs", "tables", "assumptions")
dir.create(output_dir_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_tables, recursive = TRUE, showWarnings = FALSE)

# Initialize data frame to store outlier information
outlier_info <- data.frame(
  species = character(),
  pitavg = numeric(),
  pcavg = numeric(),
  peavg = numeric(),
  pcd = numeric(),
  stringsAsFactors = FALSE
)

# Iterate over species and variables to replace outliers with NA and check for normality
for (ssp in species) {
  # Set up a 2x2 plotting layout
  
 # png(filename = file.path(output_dir_figs, paste0("pit_density_", ssp, ".png")), width = 800, height = 800)
  par(mfrow = c(2, 2), oma = c(4, 4, 2, 2), mar = c(4, 4, 2, 1))
  
  outlier_counts <- numeric(length(variables))
  names(outlier_counts) <- variables
  
  for (var in variables) {
    result <- replace_outliers_test_normality(pitdata, var, ssp)
    pitdata <- result$data
    
    # Plot density comparison
    plot_density_comparison(result$original_density, result$updated_density, var)
    
    # Store the number of points substituted for NA
    outlier_counts[var] <- result$outlier_count
  }
  
  # Add the counts to the outlier_info dataframe
  outlier_info <- rbind(outlier_info, data.frame(
    species = ssp,
    pitavg = outlier_counts["pitavg"],
    pcavg = outlier_counts["pcavg"],
    peavg = outlier_counts["peavg"],
    pcd = outlier_counts["pcd"],
    stringsAsFactors = FALSE
  ))
  
  # Add a main title for each species
  mtext(ssp, side = 3, outer = TRUE, cex = 1.5, line = 1)
  
  # Close the PNG device
  #dev.off()
}

# Save the outlier information as a CSV file
write.csv(outlier_info, file = file.path(output_dir_tables, "pit_outlier.csv"), row.names = FALSE)

# Print final data (optionally save it to a file for further analysis)
print(pitdata)
# save clean data
fwrite(pitdata, file = here("data", "processed", "pitdata_clean.csv"))

# homogeneity of variance test
###############################################

ggplot(pitdata, aes(x = ssp, y = pitavg, fill = ssp)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of pitavg Across Species",
    x = "Species", y = "pitavg"
  ) +
  theme_minimal()


### Graph shows heterocedasticity between groups

ggplot(pitdata, aes(x = ssp, y = peavg, fill = ssp)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of peavg Across Species",
    x = "Species", y = "peavg"
  ) +
  theme_minimal()

ggplot(pitdata, aes(x = ssp, y = pcavg, fill = ssp)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of pcavg Across Species",
    x = "Species", y = "pcavg"
  ) +
  theme_minimal()


ggplot(pitdata, aes(x = ssp, y = pcd, fill = ssp)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of pcd Across Species",
    x = "Species", y = "pcd"
  ) +
  theme_minimal()

### heterocedasticity present in all pit traits


#################################################
#Pit Opening data
outliers_removed_df <- data.frame(
  Species = character(),
  Outliers_Removed = numeric(),
  stringsAsFactors = FALSE
)

# Function to plot density before and after outlier removal
plot_density_comparison <- function(original_data, updated_data, species) {
  # Original density before outlier removal
  original_density <- density(original_data$Length, na.rm = TRUE)
  
  # Updated density after outlier removal
  updated_density <- density(updated_data$Length, na.rm = TRUE)
  
  # Plot density
  plot(original_density, main = paste("Density Plot for", species), 
       xlab = "Length", ylab = "Density", lwd = 2, col = "red", xlim = range(c(original_density$x, updated_density$x)))
  lines(updated_density, lwd = 2, col = "blue")
  legend("topright", legend = c("Before", "After"), col = c("red", "blue"), lwd = 2)
}

# Function to replace outliers and check for normality
process_species_data <- function(data, species) {
  # Subset data for the current species
  subset_data <- data[data$ssp == species, ]
  
  # Initialize variables
  original_data <- subset_data
  outliers_removed <- 0
  shapiro_p <- shapiro.test(subset_data$Length)$p.value
  
  # Iterate until data becomes normally distributed or only one data point is left
  while (shapiro_p < 0.05 && nrow(subset_data) > 1) {
    # Find and remove the outlier
    outliers <- boxplot.stats(subset_data$Length)$out
    if (length(outliers) == 0) break
    x1 <- outliers[which.max(abs(outliers - median(subset_data$Length, na.rm = TRUE)))]
    subset_data$Length[subset_data$Length == x1] <- NA
    outliers_removed <- outliers_removed + 1
    
    # Recalculate Shapiro-Wilk test
    shapiro_p <- shapiro.test(subset_data$Length)$p.value
  }
  
  # Plot densities before and after removal
  plot_density_comparison(original_data, subset_data, species)
  
  # Add to the DataFrame
  outliers_removed_df <<- rbind(outliers_removed_df, data.frame(
    Species = species,
    Outliers_Removed = outliers_removed,
    stringsAsFactors = FALSE
  ))
  
  return(subset_data)
}

# Process each species and plot
species_list <- unique(pitOdata$ssp)
for (species in species_list) {
  # Plot and process data
  processed_data <- process_species_data(pitOdata, species)
}

# Print the DataFrame showing the number of outliers removed
print(outliers_removed_df)

