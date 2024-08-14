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

# Relevel the factors for each data frame
dataframes <- ls()
for (df_name in dataframes) {
  df <- get(df_name) # Get the data frame by name
  
  # Check if 'ssp' column exists and is a factor or character
  if ("ssp" %in% colnames(df)) {
    if (!is.factor(df$ssp)) {
      df$ssp <- as.factor(df$ssp) # Convert to factor if not already
    }
    
    # Relevel the factor with the specified levels
    df$ssp <- factor(df$ssp, levels = c(
      "Psittacanthus robustus", "Vochysia thyrsoidea",
      "Phoradendron perrotettii", "Tapirira guianensis",
      "Struthanthus rhynchophyllus", "Tipuana tipu",
      "Viscum album", "Populus nigra"
    ))
    
    # Assign the modified data frame back to its name
    assign(df_name, df, envir = .GlobalEnv)
  }
}



# Function to plot density before and after outlier removal on the same graph
plot_density_comparison <- function(original_density, updated_density, title) {
  plot(original_density, main = title, 
       xlab = "", ylab = "Density", lwd = 2, col = "red", xlim = range(c(original_density$x, updated_density$x), na.rm = TRUE))
  lines(updated_density, lwd = 2, col = "blue")
  legend("topright", legend = c("Before", "After"), col = c("red", "blue"), lwd = 2)
}

# Function to replace outliers with NA and test for normality
replace_outliers_test_normality <- function(data, variable, species, max_iter = 10) {
  outlier_count <- 0
  
  original_data <- data[[variable]][data$ssp == species]
  shapiro_p <- NA
  
  for (i in 1:max_iter) {
    outliers <- boxplot.stats(data[[variable]][data$ssp == species])$out
    if (length(outliers) == 0) break
    # Update outliers count
    outlier_count <- length(outliers)
    data[[variable]][data[[variable]] %in% outliers & data$ssp == species] <- NA
    shapiro_test <- shapiro.test(data[[variable]][data$ssp == species])
    shapiro_p <- shapiro_test$p.value
    if (shapiro_p > 0.05) break
  }
  
  original_density <- density(original_data, na.rm = TRUE)
  updated_density <- density(data[[variable]][data$ssp == species], na.rm = TRUE)
  
  return(list(data = data, outlier_count = outlier_count, shapiro_p = shapiro_p, 
              original_density = original_density, updated_density = updated_density))
}

# Function to perform normality test across species and variables
perform_normality_test <- function(data, species_list, variables) {
  results <- data.frame(species = character(), stringsAsFactors = FALSE)
  for (ssp in species_list) {
    shapiro_results <- sapply(variables, function(var) {
      if (sum(!is.na(data[[var]][data$ssp == ssp])) > 2) {
        return(shapiro.test(data[[var]][data$ssp == ssp])$p.value)
      } else {
        return(NA)
      }
    })
    results <- rbind(results, c(ssp, shapiro_results))
  }
  colnames(results) <- c("species", variables)
  return(results)
}

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
dir.create(output_dir_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_tables, recursive = TRUE, showWarnings = FALSE)

# Iterate over variables to replace outliers with NA and check for normality in the first dataframe
for (var in variables1) {
  png(filename = file.path(output_dir_figs, paste0(var, "_density.png")), width = 1200, height = 800)
  par(mfrow = c(2, 4), oma = c(4, 4, 2, 2), mar = c(4, 4, 2, 1))
  
  for (ssp in species) {
    result <- replace_outliers_test_normality(pitdata, var, ssp)
    pitdata <- result$data
    
    plot_density_comparison(result$original_density, result$updated_density, ssp)
    
    # Update outlier_info1 with the count of replaced data points
    outlier_info1[outlier_info1$species == ssp, var] <- result$outlier_count
  }
  
  mtext(var, side = 3, outer = TRUE, cex = 1.5, line = 1)
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
  png(filename = file.path(output_dir_figs, paste0(var, "_density.png")), width = 1200, height = 800)
  par(mfrow = c(2, 4), oma = c(4, 4, 2, 2), mar = c(4, 4, 2, 1))
  
  for (ssp in species) {
    result <- replace_outliers_test_normality(pitO, var, ssp)
    pitO <- result$data
    
    plot_density_comparison(result$original_density, result$updated_density, ssp)
    
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
fwrite(pitO[,-5], file = here("data", "processed", "pitO_clean.csv"))

# Print final data
print(pitdata)
print(pitO)



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

ggplot(pitO, aes(x = ssp, y = PitDiameter, fill = ssp)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of PitDiameter Across Species",
    x = "Species", y = "Pit Diameter"
  ) +
  theme_minimal()

ggplot(pitO, aes(x = ssp, y = PitOpening, fill = ssp)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of Pit Opening Across Species",
    x = "Species", y = "Pit Opening"
  ) +
  theme_minimal()
##Only pit opening has homogeneity of variance

