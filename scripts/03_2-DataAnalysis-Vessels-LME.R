######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.2 - Data Analysis - Vessels LME
#################################################################
# Load necessary libraries and custom functions
library(here)
source(here("scripts", "00-library.R"))  # Assumptions test script
# Clear the environment
rm(list = ls())
source(here("scripts", "Functions.R"))
# Load processed data
HydraulicData <- read.csv(here("data", "processed", "HydraulicData.csv"))
vdata <- read.csv(here("data", "processed", "vdata.csv"))
vadata <- read.csv(here("data", "processed", "vadata.csv"))

### Substitute the most influential values in HydraulicData



HydraulicData[c(177,178,180),"HydraulicDiameter"] <- sort(HydraulicData$HydraulicDiameter[HydraulicData$ssp=="Vochysia thyrsoidea"])[4]
HydraulicData[c(102,98,101),"vdensity"] <- rev(sort(HydraulicData$vdensity[HydraulicData$ssp=="Tapirira guianensis"]))[4]
HydraulicData[c(172),"vdensity"] <- rev(sort(HydraulicData$vdensity[HydraulicData$ssp=="Vochysia thyrsoidea"]))[2]
HydraulicData[c(74),"VesselFraction"] <- rev(sort(HydraulicData$VesselFraction[HydraulicData$ssp=="Struthanthus rhynchophyllus"]))[1]
HydraulicData[c(51),"VesselFraction"] <- rev(sort(HydraulicData$VesselFraction[HydraulicData$ssp=="Psittacanthus robustus"]))[3]
HydraulicData[c(51),"Kmax"] <- rev(sort(HydraulicData$Kmax[HydraulicData$ssp=="Psittacanthus robustus"]))[2]
HydraulicData[c(138,141,140,134),"Kmax"] <- rev(sort(HydraulicData$Kmax[HydraulicData$ssp=="Tipuana tipu"]))[5]

### Fit Linear Mixed-Effects Models (LME) for Parasite vs Host comparisons

# List of species pairs for analysis
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

# Variables to model
vars <- colnames(HydraulicData)[c(4:6, 8)]
relevel_factors(ls())
# Initialize empty lists to store models and results
model_list <- list()
resul_table_list <- list()

# Loop over each variable to fit the models and store results
for (v in vars) {
  model_name <- paste("full_model", v, sep = "_")
  table_name <- paste(v, "results", sep = "_")
  
  # Fit the full model (with parasitism effect) and reduced model (without parasitism)
  full_model <- lme(as.formula(paste(v, "~ parasitism")), 
                    random = ~ 1 | ssp/indiv, 
                    data = HydraulicData,
                    control = list(maxIter = 100, msMaxIter = 100),
                    weights = varIdent(form = ~ 1 | ssp))
  
  reduced_model <- lme(as.formula(paste(v, "~ 1")), 
                       random = ~ 1 | ssp/indiv, 
                       data = HydraulicData,
                       control = list(maxIter = 100, msMaxIter = 100),
                       weights = varIdent(form = ~ 1 | ssp))

  
  fixed_effects <- round(fixed.effects(full_model),digits = 2)
  re <- as.numeric(VarCorr(full_model)[, "Variance"])
  # Extract residuals
  residuals<- round(residuals(full_model),digits=2)
  
  CoV_parasite <-round(sd(residuals[HydraulicData$parasitism == "Parasite"])/fixed_effects[1],digits = 2)
  CoV_host <- round(sd(residuals[HydraulicData$parasitism == "Host"])/sum(fixed_effects),digits=2)
  # Calculate specific metrics
  estimated_difference <- fixed_effects[2]
  variance_explained_by_RE <-(re[2]+re[4])/ 
    (re[5]+re[2]+re[4])
  
  # Perform LRT and calculate metrics
  lrt <- AIC(reduced_model,full_model)
  
  delta_aic <- diff(lrt$AIC)
  
  # Confidence intervals for fixed effects
  conf_intervals <- intervals(full_model)
  

  # Store the full model
  
  ms <- summary(full_model)

  model_list[[model_name]] <- full_model

  # Create a results table summarizing the model output
  resul_table_list[[table_name]] <- data.frame(
    PairTested = paste(unique(HydraulicData$parasitism), collapse = " x "),
    ParasiteMean = paste(fixed_effects[1],"(",CoV_parasite, ")",sep = ""),
    HostMean = paste(sum(fixed_effects)," (",CoV_host, ")",sep=""),
    REVariance = variance_explained_by_RE*100,
    RelDiff = estimated_difference/fixed_effects[1],
    DeltaAIC = delta_aic,
    stringsAsFactors = FALSE
  ) %>% as_tibble()
  
  # Print model summaries
  cat(paste("----------------------------------------------------------------------------------------
  Summary for full_model", v, "Parasites x Host
  ---------------------------------------------------------------------------------------------------", sep = " "))
  print(summary(full_model))
  
  # Perform ANOVA on full and reduced models
  pd <- anova(full_model, reduced_model)
  print(pd)

  # Residual plots for model diagnostics
  par(mar = c(5, 5, 5, 5))  # Set plot margins
  resid_plot <- residplot(full_model, newwd = FALSE)
  title(sub = paste("Residuals of", v, "Parasite x Host"), adj = 0.5, line = 4, cex.sub = 0.9)

  
  print(check_model(full_model))

  # Calculate Cook's Distance
  cookd_plot <- CookD(full_model, newwd = FALSE, idn = 10)
  
  # Add title and thresholds
  title(main = paste("Cook's Distance for", v, "Parasite x Host"), adj = 0.5, line = 0.5)
  abline(h = 4 / nrow(HydraulicData), col = "red", lty = 2)  # Standard threshold
  abline(h = 0.05, col = "blue", lty = 2)  # 3x mean threshold
  
  }



# Display the formatted confidence intervals

resul_table_list
HydraulicData[174,"VesselFraction"] <- rev(sort(HydraulicData$VesselFraction[HydraulicData$ssp=="Vochysia thyrsoidea"]))[2]
HydraulicData[51,"vdensity"] <- rev(sort(HydraulicData$vdensity[HydraulicData$ssp=="Psittacanthus robustus"]))[2]
HydraulicData[172,"vdensity"] <- rev(sort(HydraulicData$vdensity[HydraulicData$ssp=="Vochysia thyrsoidea"]))[2]
HydraulicData[179,"HydraulicDiameter"] <- rev(sort(HydraulicData$HydraulicDiameter[HydraulicData$ssp=="Vochysia thyrsoidea"]))[2]
HydraulicData[80,"vdensity"] <- rev(sort(HydraulicData$vdensity[HydraulicData$ssp=="Struthanthus rhynchophyllus"]))[2]
HydraulicData[90,"Kmax"] <- rev(sort(HydraulicData$Kmax[HydraulicData$ssp=="Struthanthus rhynchophyllus"]))[2]
HydraulicData[c(127),"HydraulicDiameter"] <- sort(HydraulicData$HydraulicDiameter[HydraulicData$ssp=="Tipuana tipu"])[10]
HydraulicData[c(132,134),"HydraulicDiameter"] <- rev(sort(HydraulicData$HydraulicDiameter[HydraulicData$ssp=="Tipuana tipu"]))[7]
HydraulicData[112,"Kmax"] <- rev(sort(HydraulicData$Kmax[HydraulicData$ssp=="Tapirira guianensis"]))[2]
HydraulicData[112,"VesselFraction"] <- rev(sort(HydraulicData$VesselFraction[HydraulicData$ssp=="Tapirira guianensis"]))[2]
HydraulicData[c(93,97,98,99,101,102),"vdensity"] <- rev(sort(HydraulicData$vdensity[HydraulicData$ssp=="Tapirira guianensis"]))[8]
#########################################################################################################
for (pair in species_pairs) {
  
  # Filter data for the current species pair
  data <- HydraulicData %>% subset(ssp %in% pair)
  
  for (v in vars) {
    
    # Dynamic names for models and result tables
    model_name <- paste(paste(pair, collapse = " vs "), "full_model", v, sep = "_")
    table_name <- paste(v, "results", sep = "_")
    
    # Fit the full model (with species effect) and reduced model (without species effect)
    full_model <- lme(as.formula(paste(v, "~ ssp")), 
                      random = ~ 1 | indiv, 
                      data = data,  # Use filtered data for current species pair
                      control = list(maxIter = 100, msMaxIter = 100),
                      weights = varIdent(form = ~ 1 | ssp))
    
    reduced_model <- lme(as.formula(paste(v, "~ 1")), 
                         random = ~ 1 | indiv, 
                         data = data,  # Use filtered data for current species pair
                         control = list(maxIter = 100, msMaxIter = 100),
                         weights = varIdent(form = ~ 1 | ssp))
    
    
    fixed_effects <- round(fixed.effects(full_model),digits = 2)
    re <- as.numeric(VarCorr(full_model)[, "Variance"])
    # Extract residuals
    residuals<- round(residuals(full_model),digits=2)
    
    CoV_parasite <-round(sd(residuals[data$parasitism == pair[1]])/fixed_effects[1],digits = 2)
    CoV_host <- round(sd(residuals[data$parasitism == pair[2]])/sum(fixed_effects),digits=2)
    # Calculate specific metrics
    estimated_difference <- fixed_effects[2]
    variance_explained_by_RE <-(re[2]+re[4])/ 
      (re[5]+re[2]+re[4])
    
    # Perform LRT and calculate metrics
    lrt <- AIC(reduced_model,full_model)
    
    delta_aic <- diff(lrt$AIC)
    
    # Confidence intervals for fixed effects
    conf_intervals <- intervals(full_model)
    # Store the full model in the model_list
    model_list[[model_name]] <- full_model
    
    # Create a results table summarizing the model output for the current variable and pair
    result <- data.frame(
      PairTested = paste(unique(data$ssp), collapse = " x "),
      ParasiteMean = paste(fixed_effects[1],"(",CoV_parasite, ")",sep = ""),
      HostMean = paste(sum(fixed_effects)," (",CoV_host, ")",sep=""),
      REVariance = variance_explained_by_RE*100,
      RelDiff = estimated_difference/fixed_effects[1],
      DeltaAIC = delta_aic,
      stringsAsFactors = FALSE
    ) %>% as_tibble()
    
    # Append results to the table or initialize if NULL
    if (!is.null(resul_table_list[[table_name]])) {
      resul_table_list[[table_name]] <- rbind(resul_table_list[[table_name]], result)
    } else {
      resul_table_list[[table_name]] <- result
    }
    
    # Print model summaries for the current variable
    cat(paste("----------------------------------------------------------------------------------------
  Summary for full_model", v, "Parasites x Host
  ---------------------------------------------------------------------------------------------------", sep = " "))
    print(summary(full_model))
    
    # Perform ANOVA on full and reduced models, and print the results
    pd <- anova(full_model, reduced_model)
    print(pd)

    par(mar = c(5, 5, 5, 5))  # Set plot margins
    resid_plot <- residplot(full_model, newwd = FALSE)
    title(sub = paste("Residuals of", v, paste(pair,collapse = " x ")), adj = 0.5, line = 4, cex.sub = 0.9)
    
    print(check_model(full_model))

    # Cook's distance plot to assess influential points
    cookd_plot <- CookD(full_model, newwd = FALSE,idn=10)
    title(main = paste("Cook's Distance for", v, paste(pair,collapse = " x ")), adj = 0.5, line = 0.5)
    abline(h=4/nrow(data),col="red")
    
  }  
} 
resul_table_list

# Loop through each table and save it
for (t in seq_along(resul_table_list)) {
  # Create the file name for each table
  table_name <- paste(names(resul_table_list)[t], ".csv", sep = "")
  
  # Use fwrite to save each table to the specified directory
  fwrite(resul_table_list[[t]], file = here("outputs", "tables", table_name))
  
  # Print a success message
  cat(paste(table_name, "successfully saved\n"))
}

