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
source(here("scripts","Functions.R"))
# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu")#,
  #c("Viscum album", "Populus nigra")
)
relevel_factors(ls())


###LME model
# Initialize a dataframe to store results
VWall_AIC <- data.frame(
  PairTested = character(),
  ParasiteMean = numeric(),
  HostMean = numeric(),
  EstimatedDifference = numeric(),
  REVariance = numeric(),
  RelDiff = numeric(),
  DeltaAIC = numeric(),
  stringsAsFactors = FALSE
)

full_model <- lmer(wthickness ~ parasitism + (1 | label), data = wdata)
summary(full_model)
reduced_model <- lmer(wthickness ~ 1+ (1 | label), data = wdata)
summary(reduced_model)
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
full_model$levels
# Append results to the dataframe
VWall_AIC <- rbind(VWall_AIC, data.frame(
  PairTested = paste(levels(wdata$parasitism), collapse = " vs "),
  ParasiteMean = fixed_effects[1,1], # Assumes the first level is the reference
  HostMean = fixed_effects[1,1] + fixed_effects[2,1], # Adjust if necessary
  EstimatedDifference = estimated_difference,
  REVariance = variance_explained_by_label,
  RelDiff = estimated_difference/fixed_effects[1,1],
  DeltaAIC = delta_aic,
  stringsAsFactors = FALSE
))


# Loop through each species pair and fit the model
for (pair in species_pairs) {
  
  # Subset data for the current species pair
  subset_data <- wdata[wdata$ssp %in% pair, ]
  
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
      RelDiff = estimated_difference/fixed_effects[1,1],
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
    
    testDispersion(full_model)
    title(main=paste(pair, collapse = " vs "),line= 1,cex.main=1)
    simulationOutput <- simulateResiduals(fittedModel = full_model, plot = T)
    title(main = paste(pair, collapse = " vs "),line=1,cex.main=1)
    
  }, error = function(e) {
    # Handle the error
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}

# Print the results dataframe
print(VWall_AIC) # No difference was found on wall thickness

# Save the results dataframe to a CSV file
fwrite(VWall_AIC, file = here("outputs", "tables", "Wdata_AIC.csv"))


cat("Summary\n
1. Model results:\nSignificant differences found in grouped Parasite x Host and S. rhynchophyllus x T. tipu \n
Low effect size (relative difference) in both cases\n\n
2.Residual analysys:\n
    S. rhynchophyllus x T. tipu showed a significant, but very weak deviation of homegeneity
    in DHARMa's simulation of scaled residuals.")
