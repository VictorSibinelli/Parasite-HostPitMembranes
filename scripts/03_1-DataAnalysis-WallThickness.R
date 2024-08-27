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
# Data subsets
species_data <- list(
  "P. robustus" = wdata_clean[wdata_clean$ssp == "Psittacanthus robustus", ],
  "V. thyrsoidea" = wdata_clean[wdata_clean$ssp == "Vochysia thyrsoidea", ],
  "P. perrotettii" = wdata_clean[wdata_clean$ssp == "Phoradendron perrotettii", ],
  "T. guianensis" = wdata_clean[wdata_clean$ssp == "Tapirira guianensis", ],
  "S. rhynchophyllus" = wdata_clean[wdata_clean$ssp == "Struthanthus rhynchophyllus", ],
  "T. tipu" = wdata_clean[wdata_clean$ssp == "Tipuana tipu", ],
  "V. album" = wdata_clean[wdata_clean$ssp == "Viscum album", ],
  "P. nigra" = wdata_clean[wdata_clean$ssp == "Populus nigra", ])

VWall_results <- data.frame(Parasite = character(),
                            Host = character(),
                            ParasiteMean = numeric(),
                            HostMean = numeric(),
                            MeanDifference = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)

# Perform t-tests for Pit membrane thickness
for (pair in list(c("P. robustus", "V. thyrsoidea"), 
                  c("P. perrotettii", "T. guianensis"), 
                  c("S. rhynchophyllus", "T. tipu"), 
                  c("V. album", "P. nigra"))){
  
  parasite_data <- species_data[[pair[1]]]
  host_data <- species_data[[pair[2]]]
  
  # Perform t-test for vessel wall thickness
  ttest <- t.test(parasite_data$wthickness, host_data$wthickness, var.equal = T)
  VWall_results <- rbind(VWall_results, data.frame(
    Parasite = pair[1],
    Host = pair[2],
    ParasiteMean = ttest$estimate[1],
    HostMean = ttest$estimate[2],
    MeanDifference = diff(ttest$estimate),
    pvalue= ttest$p.value,
    stringsAsFactors = FALSE
  ))
}
VWall_results
fwrite(VWall_results, file = here("outputs", "tables", "VWall_results.csv"))



# List of species pairs
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

# Initialize a dataframe to store results
results_df <- data.frame(
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
  results_df <- rbind(results_df, data.frame(
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
}

# Print the results dataframe
print(results_df)

# Save the results dataframe to a CSV file
fwrite(results_df, file = here("data", "processed", "Wdata_AIC.csv"))
