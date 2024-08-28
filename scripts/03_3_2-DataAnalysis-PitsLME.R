######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.3.2 - Data Analysis - Pit membranes
#################################################################
library(here)
source(here("scripts", "02_3-TestAssumptions-Pit.R"))
rm(list = ls())


# Load data
pitdata_clean <- read.csv(here("data", "processed", "pitdata_clean.csv"))
pitOdata_clean <- read.csv(here("data", "processed", "pitO_clean.csv"))
pitdata <- read.csv(here("data", "processed", "pitdata.csv"))
pitOdata <- read.csv(here("data", "processed", "pitOdata.csv"))

# Relevel factors
dataframes <- c("pitdata_clean", "pitOdata_clean", "pitdata", "pitOdata")
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
}

# Data subsets
species_data <- list(
  "P. robustus" = pitdata_clean[pitdata_clean$ssp == "Psittacanthus robustus", ],
  "V. thyrsoidea" = pitdata_clean[pitdata_clean$ssp == "Vochysia thyrsoidea", ],
  "P. perrotettii" = pitdata_clean[pitdata_clean$ssp == "Phoradendron perrotettii", ],
  "T. guianensis" = pitdata_clean[pitdata_clean$ssp == "Tapirira guianensis", ],
  "S. rhynchophyllus" = pitdata_clean[pitdata_clean$ssp == "Struthanthus rhynchophyllus", ],
  "T. tipu" = pitdata_clean[pitdata_clean$ssp == "Tipuana tipu", ],
  "V. album" = pitdata_clean[pitdata_clean$ssp == "Viscum album", ],
  "P. nigra" = pitdata_clean[pitdata_clean$ssp == "Populus nigra", ]
)

species_data2 <- list(
  "P. robustus" = pitOdata_clean[pitOdata_clean$ssp == "Psittacanthus robustus", ],
  "V. thyrsoidea" = pitOdata_clean[pitOdata_clean$ssp == "Vochysia thyrsoidea", ],
  "P. perrotettii" = pitOdata_clean[pitOdata_clean$ssp == "Phoradendron perrotettii", ],
  "T. guianensis" = pitOdata_clean[pitOdata_clean$ssp == "Tapirira guianensis", ],
  "S. rhynchophyllus" = pitOdata_clean[pitOdata_clean$ssp == "Struthanthus rhynchophyllus", ],
  "T. tipu" = pitOdata_clean[pitOdata_clean$ssp == "Tipuana tipu", ],
  "V. album" = pitOdata_clean[pitOdata_clean$ssp == "Viscum album", ],
  "P. nigra" = pitOdata_clean[pitOdata_clean$ssp == "Populus nigra", ]
)

# List of species pairs
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

# Initialize result data frames
PitDiameter_results <- data.frame(
  PairTested = character(),
  ParasiteMean = numeric(),
  HostMean = numeric(),
  EstimatedDifference = numeric(),
  REVariance = numeric(),
  PValue = numeric(),
  DeltaAIC = numeric(),
  stringsAsFactors = FALSE
)

PitOpening_results <- data.frame(
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
  subset_data <- pitOdata[pitOdata$ssp %in% pair, ]
  
  
  # Fit the linear mixed-effects model
  full_model <- lme(PitDiameter ~ ssp, 
                    random = ~ 1 | label, 
                    data = subset_data, 
                    weights = varIdent(form = ~ 1 | ssp))
  
  reduced_model <- lme(PitDiameter ~ 1, 
                       random = ~ 1 | label, 
                       data = subset_data, 
                       weights = varIdent(form = ~ 1 | ssp))
  
  # Print summaries of the models
  cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
  print(summary(full_model))
  
  # Perform the Likelihood Ratio Test
  lrt <- anova(reduced_model, full_model)
  
  # Extract information from the full model
  fixed_effects <- as.data.frame(summary(full_model)$tTable)
  random_effects <- as.numeric(VarCorr(full_model)[1])
  residual_variance <- sigma(full_model)^2
  
  # Extract fixed effects and variances
  estimated_difference <- fixed_effects$Value[2]
  variance_explained_by_label <- random_effects /sigma(full_model)^2
  lrt_p_value <- lrt$`p-value`
  
  # Calculate Delta AIC
  full_model_aic <- AIC(full_model)
  reduced_model_aic <- AIC(reduced_model)
  delta_aic <- full_model_aic - reduced_model_aic
  
  # Append results to the dataframe
  PitDiameter_results <- rbind(PitDiameter_results, data.frame(
    PairTested = paste(pair, collapse = " vs "),
    ParasiteMean = fixed_effects$Value[1], 
    HostMean = fixed_effects$Value[1] + fixed_effects$Value[2], 
    EstimatedDifference = estimated_difference,
    REVariance = variance_explained_by_label,
    PValue = lrt_p_value,
    DeltaAIC = delta_aic,
    stringsAsFactors = FALSE
  ))
  
  # Standardized residuals
  residuals <- resid(full_model, na.omit=T)
  fitted_values <- fitted(full_model)
  standardized_residuals <- residuals / sigma(full_model)^2
  
  # Create a data frame for plotting
  residuals_df <- data.frame(
    FittedValues = fitted_values,
    StandardizedResiduals = standardized_residuals,
    Pair = paste(pair, collapse = " vs ")
  )
  
  # Create and save the plot
  plot1 <- ggplot(residuals_df, aes(x = StandardizedResiduals)) +
    geom_histogram(binwidth = 0.1, color = "black", fill = "grey") +
    ggtitle(paste("Pit Diameter Residuals Distribution for", paste(pair, collapse = " vs "))) +
    xlab("Standardized Residuals") +
    ylab("Frequency")
  print(plot1)
  
  # Create and save the residuals vs. fitted values plot
  plot2 <- ggplot(residuals_df, aes(x = FittedValues, y = StandardizedResiduals)) +
    geom_point() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ggtitle(paste("Pit Diameter Standardized Residuals vs Fitted Values for", paste(pair, collapse = " vs "))) +
    xlab("Fitted Values") +
    ylab("Standardized Residuals") +
    theme_minimal()
  print(plot2)
}



############# repeat for Pit Opening
for (pair in species_pairs) {
  # Subset data for the current species pair
  subset_data <- pitOdata[pitOdata$ssp %in% pair, ]
  
  
  # Fit the linear mixed-effects model
  full_model <- lme(PitOpening ~ ssp, 
                    random = ~ 1 | label, 
                    data = subset_data, 
                    weights = varIdent(form = ~ 1 | ssp))
  
  reduced_model <- lme(PitOpening ~ 1, 
                       random = ~ 1 | label, 
                       data = subset_data, 
                       weights = varIdent(form = ~ 1 | ssp))
  
  # Print summaries of the models
  cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
  print(summary(full_model))
  
  # Perform the Likelihood Ratio Test
  lrt <- anova(reduced_model, full_model)
  
  # Extract information from the full model
  fixed_effects <- as.data.frame(summary(full_model)$tTable)
  random_effects <- as.numeric(VarCorr(full_model)[1])
  residual_variance <- sigma(full_model)^2
  
  # Extract fixed effects and variances
  estimated_difference <- fixed_effects$Value[2]
  variance_explained_by_label <- random_effects /sigma(full_model)^2
  lrt_p_value <- lrt$`p-value`
  
  # Calculate Delta AIC
  full_model_aic <- AIC(full_model)
  reduced_model_aic <- AIC(reduced_model)
  delta_aic <- full_model_aic - reduced_model_aic
  
  # Append results to the dataframe
  PitOpening_results <- rbind(PitOpening_results, data.frame(
    PairTested = paste(pair, collapse = " vs "),
    ParasiteMean = fixed_effects$Value[1], 
    HostMean = fixed_effects$Value[1] + fixed_effects$Value[2], 
    EstimatedDifference = estimated_difference,
    REVariance = variance_explained_by_label,
    PValue = lrt_p_value,
    DeltaAIC = delta_aic,
    stringsAsFactors = FALSE
  ))
  
  # Standardized residuals
  residuals <- resid(full_model, na.omit=T)
  fitted_values <- fitted(full_model)
  standardized_residuals <- residuals / sigma(full_model)^2
  
  # Create a data frame for plotting
  residuals_df <- data.frame(
    FittedValues = fitted_values,
    StandardizedResiduals = standardized_residuals,
    Pair = paste(pair, collapse = " vs ")
  )
  
  # Create and save the plot
  plot1 <- ggplot(residuals_df, aes(x = StandardizedResiduals)) +
    geom_histogram(binwidth = 0.1, color = "black", fill = "grey") +
    ggtitle(paste("Pit Opening Residuals Distribution for", paste(pair, collapse = " vs "))) +
    xlab("Standardized Residuals") +
    ylab("Frequency")
  print(plot1)
  
  # Create and save the residuals vs. fitted values plot
  plot2 <- ggplot(residuals_df, aes(x = FittedValues, y = StandardizedResiduals)) +
    geom_point() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ggtitle(paste("Pit Opening Standardized Residuals vs Fitted Values for", paste(pair, collapse = " vs "))) +
    xlab("Fitted Values") +
    ylab("Standardized Residuals") +
    theme_minimal()
  print(plot2)
}


PitDiameter_results[,2:7] <- round(PitDiameter_results[,2:7],3) #round values
PitDiameter_results <- PitDiameter_results %>%  drop_na() #remove redundant lines
fwrite(PitDiameter_results, file=here("outputs","tables","PitDiameter_AIC")) #save data frame


PitOpening_results[,2:7] <- round(PitOpening_results[,2:7],3) #round values
PitOpening_results <- PitOpening_results %>%  drop_na() #remove redundant lines
fwrite(PitOpening_results, file=here("outputs","tables","PitOpening_AIC")) #save data frame
