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
##HydraulicDiameter

HydraulicData[c(177,178,179,180),"HydraulicDiameter"] <- sort(HydraulicData$HydraulicDiameter[HydraulicData$ssp=="Vochysia thyrsoidea"])[4]
HydraulicData[c(172),"vdensity"] <- rev(sort(HydraulicData$vdensity[HydraulicData$ssp=="Vochysia thyrsoidea"]))[2]

HydraulicData[c(51),"Kmax"] <- rev(sort(HydraulicData$Kmax[HydraulicData$ssp=="Psittacanthus robustus"]))[2]


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
  
  # Store the full model
  
  ms <- summary(full_model)

  model_list[[model_name]] <- full_model

  # Create a results table summarizing the model output
  resul_table_list[[table_name]] <- data.frame(
    PairTested = paste(unique(HydraulicData$parasitism), collapse = " x "),
    ParasiteMean = paste0(round(full_model$coefficients$fixed[1], 3), " ± ", round(ms$tTable[,"Std.Error"][1], 3)),
    HostMean = paste0(round(sum(full_model$coefficients$fixed), 3), " ± ", round(ms$tTable[,"Std.Error"][2], 3)),
    EstimatedDifference = round(full_model$coefficients$fixed[2], 3),
    REVariance = (as.numeric(VarCorr(full_model)[2]) + as.numeric(VarCorr(full_model)[4])) / 
      as.numeric(VarCorr(full_model)[5]),
    PValue = anova(full_model, reduced_model)$`p-value`[2],
    DeltaAIC = AIC(full_model) - AIC(reduced_model),
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

  cookd_plot <- CookD(full_model, newwd = FALSE,idn = 10)
  title(main = paste("Cook's Distance for",v, "Parasite x Host"), adj = 0.5, line = 0.5)
  abline(h=4/nrow(HydraulicData),col="red")
  }

# Display the formatted confidence intervals

resul_table_list

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
    
    # Store the full model in the model_list
    model_list[[model_name]] <- full_model
    
    # Create a results table summarizing the model output for the current variable and pair
    result <- data.frame(
      PairTested = paste(unique(data$ssp), collapse = " x "),  # Use filtered data
      ParasiteMean = paste0(round(full_model$coefficients$fixed[1], 3), " ± ", round(ms$tTable[,"Std.Error"][1], 3)),
      HostMean = paste0(round(sum(full_model$coefficients$fixed), 3), " ± ", round(ms$tTable[,"Std.Error"][2], 3)),
      EstimatedDifference = round(full_model$coefficients$fixed[2], 3),
      REVariance =as.numeric(VarCorr(full_model)[,"Variance"][1])/sum(as.numeric(VarCorr(full_model)[,"Variance"])),
      PValue = anova(full_model, reduced_model)$`p-value`[2],
      DeltaAIC = AIC(full_model) - AIC(reduced_model),
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
    cookd_plot <- CookD(full_model, newwd = FALSE)
    title(main = paste("Cook's Distance for", v, paste(pair,collapse = " x ")), adj = 0.5, line = 0.5)
    abline(h=4/nrow(data),col="red")
  }  
} 
resul_table_list

# Loop through each table and save it
for (t in seq_along(resul_table_list)) {
  # Create the file name for each table
  table_name <- paste(names(resul_table_list)[t], ".R", sep = "")
  
  # Use fwrite to save each table to the specified directory
  fwrite(resul_table_list[[t]], file = here("outputs", "tables", table_name))
  
  # Print a success message
  cat(paste(table_name, "successfully saved\n"))
}

