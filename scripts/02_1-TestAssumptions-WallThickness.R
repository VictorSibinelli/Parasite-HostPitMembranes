######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 02.2 - Test assumptions test
######################################################################


# Load the data wrangling script
source(here("scripts", "01-DataWrangling.R"))

# Load data
wdata <- read.csv(here("data", "processed", "wdata.csv"))

# Function to plot density before and after outlier removal on the same graph
plot_density_comparison <- function(original_density, updated_density, species, variable) {
  plot(original_density, main = paste(species, "-", variable), 
       xlab = variable, ylab = "Density", lwd = 2, col = "red", 
       xlim = range(c(original_density$x, updated_density$x), na.rm = TRUE))
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

# Variables of interest
variables <- c("wthickness")

# Perform Shapiro-Wilk test for each species and variable
species <- unique(wdata$ssp)

# Directory for saving plots and table
output_dir_figs <- here("outputs", "figs", "assumptions")
output_dir_tables <- here("outputs", "tables", "assumptions")
dir.create(output_dir_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_tables, recursive = TRUE, showWarnings = FALSE)

# Initialize data frame to store outlier information
outlier_info <- data.frame(
  species = character(),
  wthickness = numeric(),
  stringsAsFactors = FALSE
)

# Set up the plotting area with 2 rows and 4 columns
png(filename = file.path(output_dir_figs, "WthicknessDensity.png"), width = 1600, height = 800)
par(mfrow = c(2, 4), oma = c(4, 4, 2, 2), mar = c(4, 4, 2, 1))

# Iterate over species and variables to replace outliers with NA and check for normality
for (ssp in species) {
  outlier_counts <- numeric(length(variables))
  names(outlier_counts) <- variables
  
  for (var in variables) {
    result <- replace_outliers_test_normality(wdata, var, ssp)
    wdata <- result$data
    
    # Plot density comparison
    plot_density_comparison(result$original_density, result$updated_density, ssp, var)
    
    # Store the number of points substituted for NA
    outlier_counts[var] <- result$outlier_count
  }
  
  # Add the counts to the outlier_info dataframe
  outlier_info <- rbind(outlier_info, data.frame(
    species = ssp,
    wthickness = outlier_counts["wthickness"],
    stringsAsFactors = FALSE
  ))
}

# Close the PNG device
dev.off()

# Save the outlier information as a CSV file
write.csv(outlier_info, file = file.path(output_dir_tables, "w_outlier.csv"), row.names = FALSE)
print(outlier_info)
# Save the cleaned data
fwrite(wdata, file = here("data", "processed", "wdata_clean.csv"))

# Print final data (optionally save it to a file for further analysis)
print(wdata)


# Homogeneity of variance test
ggplot(wdata, aes(x = ssp, y = wthickness, fill = ssp)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of Wall Thickness Across Species",
    x = "Species", y = "Wall Thickness"
  ) +
  theme_minimal()
#hogeneity of variance found inside pairs
