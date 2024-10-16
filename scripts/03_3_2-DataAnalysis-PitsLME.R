######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.3.2 - Data Analysis - Pit membranes
#####################################################################

library(here)
source(here("scripts", "02_3-TestAssumptions-Pit.R")) # Source test assumptions script
rm(list=ls())
source(here("scripts", "Functions.R")) # Source custom functions

# Load data
pitdata <- read.csv(here("data", "processed", "pitdata.csv"))
pitOdata <- read.csv(here("data", "processed", "pitOdata.csv"))

# Modify outlier values to the second-highest value
n <- length(pitOdata$PitDiameter)
pitOdata[209, "PitDiameter"] <- sort(pitOdata$PitDiameter, partial = n - 1)[n - 1]
pitOdata[465, "PitOpening"] <- sort(pitOdata$PitOpening, partial = n - 1)[n - 1]

# Replace the 5 largest PitOpening values for Tipuana tipu with the 6th largest
tipuana_indices <- pitOdata$ssp == "Tipuana tipu"
top_5_tipuanas <- tail(sort(pitOdata$PitOpening[tipuana_indices]), 5)
pitOdata$PitOpening[tipuana_indices & pitOdata$PitOpening %in% top_5_tipuanas] <- 
  sort(pitOdata$PitOpening[tipuana_indices])[6]

# Relevel factors
relevel_factors(ls())

# List of species pairs
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album","Populus nigra")
)

# Initialize result data frames
PitDiameter_results <- data.frame(PairTested = character(),
                                  ParasiteMean = numeric(),
                                  HostMean = numeric(),
                                  EstimatedDifference = numeric(),
                                  REVariance = numeric(),
                                  PValue = numeric(),
                                  DeltaAIC = numeric(),
                                  stringsAsFactors = FALSE) %>% as.tibble

PitOpening_results <- data.frame(PairTested = character(),
                                 ParasiteMean = numeric(),
                                 HostMean = numeric(),
                                 EstimatedDifference = numeric(),
                                 REVariance = numeric(),
                                 PValue = numeric(),
                                 DeltaAIC = numeric(),
                                 stringsAsFactors = FALSE) %>% as.tibble()


PitMembrane_results <- data.frame(Parasite = character(),
                                  Host = character(),
                                  ParasiteMean = numeric(),
                                  HostMean = numeric(),
                                  MeanDifference = numeric(),
                                  pvalue = numeric(),
                                  stringsAsFactors = FALSE) %>% as.tible

pcd_results <- data.frame(Parasite = character(),
                          Host = character(),
                          ParasiteMean = numeric(),
                          HostMean = numeric(),
                          MeanDifference = numeric(),
                          pvalue = numeric(),
                          stringsAsFactors = FALSE) %>% as.tibble()

# Fit models for PitDiameter and perform ANOVA
full_model_diameter <- lme(PitDiameter ~ parasitism, random = ~ 1 | ssp/label, data = pitOdata,
                           weights = varIdent(form = ~ 1 | ssp))
reduced_model_diameter <- lme(PitDiameter ~ 1, random = ~ 1 | ssp/label, data = pitOdata,
                              weights = varIdent(form = ~ 1 | ssp))
pd <- anova(full_model_diameter, reduced_model_diameter) # Group-level difference non-significant
pd
# Fit models for PitOpening and perform ANOVA
full_model_opening <- lme(PitOpening ~ parasitism, random = ~ 1 | ssp/label, data = pitOdata,
                          weights = varIdent(form = ~ 1 | ssp))
reduced_model_opening <- lme(PitOpening ~ 1, random = ~ 1 | ssp/label, data = pitOdata,
                             weights = varIdent(form = ~ 1 | ssp))
po <- anova(full_model_opening, reduced_model_opening)
po# Pit Opening difference significant
summary(full_model_opening) # Display summary of the full model

# Append results to the PitDiameter_results dataframe
PitDiameter_results <- rbind(PitDiameter_results, data.frame(
  PairTested = paste(levels(pitdata$parasitism), collapse = " vs "),
  ParasiteMean =summary(full_model_diameter)$tTable[1], 
  HostMean = summary(full_model_diameter)$tTable[1] + summary(full_model_diameter)$tTable[2], 
  EstimatedDifference = summary(full_model_diameter)$tTable[2],
  REVariance =(as.numeric(VarCorr(full_model_diameter)[[2]]) + as.numeric(VarCorr(full_model_diameter)[[4]])) /
    as.numeric(VarCorr(full_model_diameter)[[5]]),
  PValue = pd$`p-value`[2],
  DeltaAIC = AIC(full_model_diameter)-AIC(reduced_model_diameter),
  stringsAsFactors = FALSE
))

PitOpening_results <- rbind(PitOpening_results, data.frame(
  PairTested = paste(levels(pitdata$parasitism), collapse = " vs "),
  ParasiteMean =summary(full_model_opening)$tTable[1], 
  HostMean = summary(full_model_opening)$tTable[1] + summary(full_model_opening)$tTable[2], 
  EstimatedDifference = summary(full_model_opening)$tTable[2],
  REVariance =(as.numeric(VarCorr(full_model_opening)[[2]]) + as.numeric(VarCorr(full_model_opening)[[4]])) /
    as.numeric(VarCorr(full_model_opening)[[5]]),
  PValue = po$`p-value`[2],
  DeltaAIC = AIC(full_model_opening)-AIC(reduced_model_opening),
  stringsAsFactors = FALSE
))

# Loop through each species pair to analyze PitDiameter
for (pair in species_pairs) {
  # Subset data for the current species pair
  subset_data <- pitOdata[pitOdata$ssp %in% pair, ]
  
  # Check if there is data for the species pair
  if (nrow(subset_data) < 1) {
    cat("No data for species pair:", paste(pair, collapse = " vs "), "\n")
    next
  }
  
  # Fit the full model and handle potential errors
  full_model <- tryCatch({
    lme(PitDiameter ~ ssp, random = ~ 1 | label, data = subset_data, 
        weights = varIdent(form = ~ 1 | ssp))
  }, error = function(e) {
    cat("Error in fitting model for pair:", paste(pair, collapse = " vs "), "\n")
    return(NULL)
  })
  
  if (is.null(full_model)) next  # Skip to the next pair if model fitting failed
  
  # Fit the reduced model and handle potential errors
  reduced_model <- tryCatch({
    lme(PitDiameter ~ 1, random = ~ 1 | label, data = subset_data, 
        weights = varIdent(form = ~ 1 | ssp))
  }, error = function(e) {
    cat("Error in fitting reduced model for pair:", paste(pair, collapse = " vs "), "\n")
    return(NULL)
  })
  
  if (is.null(reduced_model)) next  # Skip to the next pair if model fitting failed
  
  # Print summaries of the models
  cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
  print(summary(full_model))
  
  # Perform the Likelihood Ratio Test
  lrt <- anova(reduced_model, full_model)
  
  # Extract fixed effects and check for sufficient data
  fixed_effects <- as.data.frame(summary(full_model)$tTable)
  if (nrow(fixed_effects) < 2) {
    cat("Insufficient fixed effects for pair:", paste(pair, collapse = " vs "), "\n")
    next
  }
  
  # Calculate metrics
  estimated_difference <- fixed_effects$Value[2]
  variance_explained_by_label <- as.numeric(VarCorr(full_model)[1]) / 
    as.numeric(VarCorr(full_model)[2])
  lrt_p_value <- lrt$`p-value`[2]  # Correct index for p-value
  delta_aic <- AIC(full_model) - AIC(reduced_model)
  
  # Append results to the PitDiameter_results dataframe
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
  
  # Calculate standardized residuals for plotting
  residuals <- resid(full_model, na.omit = TRUE)
  fitted_values <- fitted(full_model)
  standardized_residuals <- residuals / sigma(full_model)^2
  
  # Create a data frame for plotting residuals
  residuals_df <- data.frame(
    FittedValues = fitted_values,
    StandardizedResiduals = standardized_residuals,
    Pair = paste(pair, collapse = " vs ")
  )
  
  # Residual plots for model analysis
  par(mar = c(5, 5, 5, 5))  # Set margins for the plot
  resid_plot <- residplot(full_model, newwd = FALSE)
  title(sub = paste("PitDiameter", paste(pair, collapse = " vs ")), 
        adj = 0.5, line = 4, cex.sub = 0.9)
  print(check_model(full_model))
  
  # Cook's distance plot
  cookd_plot <- CookD(full_model, newwd = FALSE)
  title(main = paste("PitDiameter", paste(pair, collapse = " vs ")), 
        adj = 0.5, line = 0.5)
}

#########################################################
## Pit Opening Analysis
for (pair in species_pairs) {
  # Subset data for the current species pair
  subset_data <- pitOdata[pitOdata$ssp %in% pair, ]
  
  # Check if there is data for the species pair
  if (nrow(subset_data) < 1) {
    cat("No data for species pair:", paste(pair, collapse = " vs "), "\n")
    next
  }
  
  # Fit the full model and handle potential errors
  full_model <- tryCatch({
    lme(PitOpening ~ ssp, random = ~ 1 | label, data = subset_data, 
        weights = varIdent(form = ~ 1 | ssp))
  }, error = function(e) {
    cat("Error in fitting model for pair:", paste(pair, collapse = " vs "), "\n")
    return(NULL)
  })
  
  if (is.null(full_model)) next  # Skip to the next pair if model fitting failed
  
  # Fit the reduced model and handle potential errors
  reduced_model <- tryCatch({
    lme(PitOpening ~ 1, random = ~ 1 | label, data = subset_data, 
        weights = varIdent(form = ~ 1 | ssp))
  }, error = function(e) {
    cat("Error in fitting reduced model for pair:", paste(pair, collapse = " vs "), "\n")
    return(NULL)
  })
  
  if (is.null(reduced_model)) next  # Skip to the next pair if model fitting failed
  
  # Print summaries of the models
  cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
  print(summary(full_model))
  
  # Perform the Likelihood Ratio Test
  lrt <- anova(reduced_model, full_model)
  
  # Extract fixed effects and check for sufficient data
  fixed_effects <- as.data.frame(summary(full_model)$tTable)
  if (nrow(fixed_effects) < 2) {
    cat("Insufficient fixed effects for pair:", paste(pair, collapse = " vs "), "\n")
    next
  }
  
  # Calculate metrics
  estimated_difference <- fixed_effects$Value[2]
  variance_explained_by_label <- as.numeric(VarCorr(full_model)[1]) / 
    as.numeric(VarCorr(full_model)[2])
  lrt_p_value <- lrt$`p-value`[2]  # Correct index for p-value
  delta_aic <- AIC(full_model) - AIC(reduced_model)
  
  # Append results to the PitOpening_results dataframe
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
  
  # Calculate standardized residuals for plotting
  residuals <- resid(full_model, na.omit = TRUE)
  fitted_values <- fitted(full_model)
  standardized_residuals <- residuals / sigma(full_model)^2
  
  # Create a data frame for plotting residuals
  residuals_df <- data.frame(
    FittedValues = fitted_values,
    StandardizedResiduals = standardized_residuals,
    Pair = paste(pair, collapse = " vs ")
  )
  
  # Residual plots for model analysis
  par(mar = c(5, 5, 5, 5))  # Set margins for the plot
  resid_plot <- residplot(full_model, newwd = FALSE)
  title(sub = paste("PitOpening", paste(pair, collapse = " vs ")), 
        adj = 0.5, line = 4, cex.sub = 0.9)
  print(check_model(full_model, residual_type = "simulated"))
  
  # Cook's distance plot
  cookd_plot <- CookD(full_model, newwd = FALSE)
  title(main = paste("PitOpening", paste(pair, collapse = " vs ")), 
        adj = 0.5, line = 0.5)
}

######Pit chamber
ssps <- unlist(species_pairs)
# Initialize empty data frame to store results
PeXPc_df <- data.frame(Difference = numeric(), PValue = numeric(), stringsAsFactors = FALSE)

for (s in ssps) {
  subset_data <- pitdata[pitdata$ssp == s, ]
  
  if (nrow(subset_data) > 1) { # Ensure there's enough data for the t-test
    t <- t.test(subset_data$peavg, subset_data$pcavg, var.equal = TRUE)
    result <- data.frame(Difference = diff(t$estimate), PValue = t$p.value)
    PeXPc_df <- rbind(PeXPc_df, result)
    rownames(PeXPc_df)[nrow(PeXPc_df)] <- s
  } else {
    cat("Not enough data for species:", s, "\n")
  }
}
# Display the result
PeXPc_df


full_model <- lme(pcd ~ parasitism, random = ~ 1 | ssp, data = pitdata)
reduced_model <- lme(pcd ~ 1, random = ~ 1 | ssp, data = pitdata)
a <- anova(full_model,reduced_model)#non siginificant
full_model$coefficients$fixed

pcd_results <- rbind(pcd_results, data.frame(
  Parasite = "Parasites",
  Host = "Hosts",
  ParasiteMean = full_model$coefficients$fixed[1]*1000,
  HostMean = sum(full_model$coefficients$fixed) * 1000,
  MeanDifference = full_model$coefficients$fixed[2] * 1000,
  pvalue = a$`p-value`[2],
  stringsAsFactors = FALSE
))


full_model <- lme(pitavg ~ parasitism, random = ~ 1 | ssp,na.omit(pitdata))
reduced_model <- lme(pitavg ~ 1, random = ~ 1 | ssp,na.omit(pitdata))
anova(full_model,reduced_model)#Significant
PitMembrane_results <- rbind(PitMembrane_results, data.frame(
  Parasite = "Parasites",
  Host = "Hosts",
  ParasiteMean = full_model$coefficients$fixed[1]*1000,
  HostMean = sum(full_model$coefficients$fixed) * 1000,
  MeanDifference = full_model$coefficients$fixed[2] * 1000,
  pvalue = a$`p-value`[2],
  stringsAsFactors = FALSE
))

for (pair in species_pairs) {
  subset_data <- pitdata[pitdata$ssp %in% pair, ]
  t1 <- t.test(pcd~ssp,var.equal=F,data=subset_data)
  pcd_results <- rbind(pcd_results, data.frame(
    Parasite = pair[1],
    Host = pair[2],
    ParasiteMean = t1$estimate[1] * 1000,
    HostMean = t1$estimate[2] * 1000,
    MeanDifference = diff(t1$estimate) * 1000,
    pvalue = t1$p.value,
    stringsAsFactors = FALSE
  ))
  t2 <- t.test(pitavg~ssp,var.equal=T,data = subset_data)
  PitMembrane_results <- rbind(PitMembrane_results, data.frame(
    Parasite = pair[1],
    Host = pair[2],
    ParasiteMean = t2$estimate[1],
    HostMean = t2$estimate[2],
    MeanDifference = diff(t2$estimate),
    pvalue = t2$p.value,
    stringsAsFactors = FALSE
  ))
}


# Convert data frames to tibbles
PitDiameter_results <- as_tibble(PitDiameter_results)
PitOpening_results <- as_tibble(PitOpening_results)
PitMembrane_results <- as_tibble(PitMembrane_results)
pcd_results <- as_tibble(pcd_results)

##########################################################

#save result tables
fwrite(PitDiameter_results, file=here("outputs","tables","PitDiameter_AIC")) #save data frame
fwrite(PitOpening_results, file=here("outputs","tables","PitOpening_AIC")) #save data frame
fwrite(pcd_results, file=here("outputs","tables","PCD_AIC")) #save data frame
fwrite(PitMembrane_results, file=here("outputs","tables","PitMembrane_AIC")) #save data frame
