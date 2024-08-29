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
  c("Struthanthus rhynchophyllus", "Tipuana tipu")  #,
  #c("Viscum album", "Populus nigra")
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

for (pair in species_pairs) {
  # Subset data for the current species pair
  subset_data <- pitOdata[pitOdata$ssp %in% pair, ]
  
  if (nrow(subset_data) < 1) {
    cat("No data for species pair:", paste(pair, collapse = " vs "), "\n")
    next
  }
  
  # Fit the linear mixed-effects model
  full_model <- tryCatch({
    lme(PitDiameter ~ ssp, 
        random = ~ 1 | label, 
        data = subset_data, 
        weights = varIdent(form = ~ 1 | ssp))
  }, error = function(e) {
    cat("Error in fitting model for pair:", paste(pair, collapse = " vs "), "\n")
    return(NULL)
  })
  
  if (is.null(full_model)) {
    next
  }
  
  # Fit the reduced model
  reduced_model <- tryCatch({
    lme(PitDiameter ~ 1, 
        random = ~ 1 | label, 
        data = subset_data, 
        weights = varIdent(form = ~ 1 | ssp))
  }, error = function(e) {
    cat("Error in fitting reduced model for pair:", paste(pair, collapse = " vs "), "\n")
    return(NULL)
  })
  
  if (is.null(reduced_model)) {
    next
  }
  
  # Print summaries of the models
  cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
  print(summary(full_model))
  
  # Perform the Likelihood Ratio Test
  lrt <- anova(reduced_model, full_model)
  
  # Extract information from the full model
  fixed_effects <- as.data.frame(summary(full_model)$tTable)
  
  if (nrow(fixed_effects) < 2) {
    cat("Insufficient fixed effects for pair:", paste(pair, collapse = " vs "), "\n")
    next
  }
  
  random_effects <- as.numeric(VarCorr(full_model)[1])
  
  # Calculate metrics
  estimated_difference <- fixed_effects$Value[2]
  variance_explained_by_label <- random_effects / sigma(full_model)^2
  lrt_p_value <- lrt$`p-value`[2] # Use correct index
  delta_aic <- AIC(full_model) - AIC(reduced_model)
  
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
  
  # Residual plots for model analysis
  par(mar=c(5,5,5,5))
  resid_plot <- residplot(full_model,newwd = F)
  title(sub = paste("PitDiameter", paste(pair, collapse = " vs ")), adj=0.5, line = 4,cex.sub=0.9)
  
  cookd_plot <- CookD(full_model,newwd = F)
  title(main = paste("PitDiameter", paste(pair, collapse = " vs ")), adj=0.5, line = 0.5)
  
  
}


for (pair in species_pairs) {
  # Subset data for the current species pair
  subset_data <- pitOdata[pitOdata$ssp %in% pair, ]
  
  if (nrow(subset_data) < 1) {
    cat("No data for species pair:", paste(pair, collapse = " vs "), "\n")
    next
  }
  
  # Fit the linear mixed-effects model
  full_model <- tryCatch({
    lme(PitOpening ~ ssp, 
        random = ~ 1 | label, 
        data = subset_data, 
        weights = varIdent(form = ~ 1 | ssp))
  }, error = function(e) {
    cat("Error in fitting full model for pair:", paste(pair, collapse = " vs "), "\n")
    return(NULL)
  })
  
  if (is.null(full_model)) {
    next
  }
  
  # Fit the reduced model
  reduced_model <- tryCatch({
    lme(PitOpening ~ 1, 
        random = ~ 1 | label, 
        data = subset_data, 
        weights = varIdent(form = ~ 1 | ssp))
  }, error = function(e) {
    cat("Error in fitting reduced model for pair:", paste(pair, collapse = " vs "), "\n")
    return(NULL)
  })
  
  if (is.null(reduced_model)) {
    next
  }
  
  # Print summaries of the models
  cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
  print(summary(full_model))
  
  # Perform the Likelihood Ratio Test
  lrt <- anova(reduced_model, full_model)
  
  # Extract information from the full model
  fixed_effects <- as.data.frame(summary(full_model)$tTable)
  
  if (nrow(fixed_effects) < 2) {
    cat("Insufficient fixed effects for pair:", paste(pair, collapse = " vs "), "\n")
    next
  }
  
  random_effects <- as.numeric(VarCorr(full_model)[1])
  
  # Calculate metrics
  estimated_difference <- fixed_effects$Value[2]
  variance_explained_by_label <- random_effects / sigma(full_model)^2
  lrt_p_value <- lrt$`p-value`[2] # Use correct index
  delta_aic <- AIC(full_model) - AIC(reduced_model)
  
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
  
#residual plot analysis for the models
  par(mar=c(5,5,5,1))
  resid_plot <- residplot(full_model, newwd = F)
  title(sub = paste("PitOpening", paste(pair, collapse = " vs ")), adj=0.5, line = 4, cex.sub=0.9)
  
  cookd_plot <- CookD(full_model, newwd = F)
  title(main = paste("PitOpening", paste(pair, collapse = " vs ")), adj=0.5, line = 0.5)
}




#save result tables
fwrite(PitDiameter_results, file=here("outputs","tables","PitDiameter_AIC")) #save data frame

fwrite(PitOpening_results, file=here("outputs","tables","PitOpening_AIC")) #save data frame
