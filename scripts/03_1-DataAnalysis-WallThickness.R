######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.1 - Data Analysis - Vessel walls
#################################################################
library(here)
source(here("scripts", "02_1-TestAssumptions-WallThickness.R"))
rm(list=ls())

wdata <- read.csv(here("data", "processed", "wdata.csv"))
wdata_clean <- read.csv(here("data", "processed", "wdata_clean.csv"))



# List of data frame names
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
  rm(df, df_name) # remove duplicated dataframe
}

# Data subsets by species
species_data <- list(
  "Psittacanthus robustus" = wdata_clean[wdata_clean$ssp == "Psittacanthus robustus", ],
  "Vochysia thyrsoidea" = wdata_clean[wdata_clean$ssp == "Vochysia thyrsoidea", ],
  "Phoradendron perrotettii" = wdata_clean[wdata_clean$ssp == "Phoradendron perrotettii", ],
  "Tapirira guianensis" = wdata_clean[wdata_clean$ssp == "Tapirira guianensis", ],
  "Struthanthus rhynchophyllus" = wdata_clean[wdata_clean$ssp == "Struthanthus rhynchophyllus", ],
  "Tipuana tipu" = wdata_clean[wdata_clean$ssp == "Tipuana tipu", ],
  "Viscum album" = wdata_clean[wdata_clean$ssp == "Viscum album", ],
  "Populus nigra" = wdata_clean[wdata_clean$ssp == "Populus nigra", ]
)

# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

# Initialize results data frame
VWall_results <- data.frame(
  Parasite = character(),
  Host = character(),
  ParasiteMean = numeric(),
  HostMean = numeric(),
  MeanDifference = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

# Perform t-tests for vessel wall thickness
for (pair in species_pairs) {
  # Extract data for each species in the pair
  parasite_data <- species_data[[pair[1]]]
  host_data <- species_data[[pair[2]]]
  
  # Check if there is data for both species
  if (nrow(parasite_data) == 0 || nrow(host_data) == 0) {
    cat("No data available for species pair:", pair[1], "and", pair[2], "\n")
    next # Skip this pair if data is missing
  }
  
  # Use tryCatch to handle errors during the t-test
  tryCatch({
    # Perform t-test for vessel wall thickness
    ttest <- t.test(parasite_data$wthickness, host_data$wthickness, var.equal = TRUE)
    
    # Append t-test results to the data frame
    VWall_results <- rbind(VWall_results, data.frame(
      Parasite = pair[1],
      Host = pair[2],
      ParasiteMean = ttest$estimate[1],
      HostMean = ttest$estimate[2],
      MeanDifference = diff(ttest$estimate),
      pvalue = ttest$p.value,
      stringsAsFactors = FALSE
    ))
    row.names(VWall_results) <- NULL
  }, error = function(e) {
    # Handle the error, print a message and continue the loop
    cat("Error in t-test for species pair:", pair[1], "and", pair[2], "\n")
    print(e)
  })
}

# Display the results
print(VWall_results)

# Write results to a CSV file
fwrite(VWall_results, file = here("outputs", "tables", "VWall_results.csv"))


###LME model

# Initialize a dataframe to store results
VWall_AIC <- data.frame(
  PairTested = character(),
  ParasiteMean = numeric(),
  HostMean = numeric(),
  EstimatedDifference = numeric(),
  REVariance = numeric(),
  PValue = numeric(),
  DeltaAIC = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each species pair and fit the model
for (pair in species_pairs) {
  
  # Subset data for the current species pair
  subset_data <- wdata_clean[wdata_clean$ssp %in% pair, ]
  
  # Check if subset has at least two levels of 'ssp'
  if (length(unique(subset_data$ssp)) < 2) {
    cat("\nSkipping pair", paste(pair, collapse = " vs "), "due to insufficient levels in 'ssp'.\n")
    next
  }
  
  # Try to fit the models using tryCatch
  tryCatch({
    # Fit the linear mixed-effects model
    full_model <- lmer(wthickness ~ ssp + (1 | label), data = subset_data)
    reduced_model <- lmer(wthickness ~ (1 | label), data = subset_data)
    
    # Print summaries of the models
    cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
    print(summary(full_model))
    
    # Perform the Likelihood Ratio Test
    lrt <- anova(reduced_model, full_model)
    
    # Extract information from the full model
    fixed_effects <- summary(full_model)$coefficients
    random_effects <- as.data.frame(VarCorr(full_model))
    residual_variance <- sigma(full_model)^2
    
    # Extract fixed effects and variances
    estimated_difference <- fixed_effects[2,1] # Adjust if necessary
    variance_explained_by_label <- random_effects$vcov[1] / sum(random_effects$vcov)
    lrt_p_value <- lrt$`Pr(>Chisq)`[2]
    
    # Calculate Delta AIC
    full_model_aic <- AIC(full_model)
    reduced_model_aic <- AIC(reduced_model)
    delta_aic <- full_model_aic - reduced_model_aic
    
    # Append results to the dataframe
    VWall_AIC <- rbind(VWall_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = fixed_effects[1,1], # Assumes the first level is the reference
      HostMean = fixed_effects[1,1] + fixed_effects[2,1], # Adjust if necessary
      EstimatedDifference = estimated_difference,
      REVariance = variance_explained_by_label,
      PValue = lrt_p_value,
      DeltaAIC = delta_aic,
      stringsAsFactors = FALSE
    ))
    
    # Standardized residuals
    residuals <- resid(full_model)
    fitted_values <- fitted(full_model)
    standardized_residuals <- residuals / sqrt(residual_variance)
    
    # Create a data frame for plotting
    residuals_df <- data.frame(
      FittedValues = fitted_values,
      StandardizedResiduals = standardized_residuals,
      Pair = paste(pair, collapse = " vs ")
    )
    # Create and save the plot
    plot <- ggplot(residuals_df, aes(x = residuals)) +
      geom_histogram(binwidth = 0.1, color = "black", fill = "lightblue") +
      ggtitle(paste("Residuals Distribution for", paste(pair, collapse = " vs "))) +
      xlab("Residuals") +
      ylab("Frequency")
    print(plot)
    
    # Create and save the plot
    plot <- ggplot(residuals_df, aes(x = FittedValues, y = StandardizedResiduals)) +
      geom_point() +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      ggtitle(paste("Standardized Residuals vs Fitted Values for", paste(pair, collapse = " vs "))) +
      xlab("Fitted Values") +
      ylab("Standardized Residuals") +
      theme_minimal()
    print(plot)
    
  }, error = function(e) {
    # Handle the error
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}


# Print the results dataframe
print(VWall_AIC) # No difference was found on wall thickness

# Save the results dataframe to a CSV file
fwrite(VWall_AIC, file = here("outputs", "tables", "Wdata_AIC.csv"))

