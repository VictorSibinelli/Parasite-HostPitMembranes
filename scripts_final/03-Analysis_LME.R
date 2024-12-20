### LME analysis



library(tidyverse)
library(here)
library(nlme)
library(predictmeans)
library(performance)
source(here("scripts","00-library.R"))

# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data<- read.csv(here("data", "processed", "VesselDiameter_data.csv"))
Hydraulic_data<- read.csv(here("data", "processed", "HydraulicData.csv"))
PitFraction_data<- read.csv(here("data", "processed", "PitFraction_data.csv"))
PitDiOp_data<- read.csv(here("data", "processed", "PitDiOp_data.csv"))
PitMembrane_data<- read.csv(here("data", "processed", "PitMembrane_data.csv"))
source(here("scripts", "Functions.R"))

# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

relevel_factors(ls())


####Wall thickness
# Initialize data frames for results
VWall_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
Wall_full <- lme(
  WallThickness ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Wall_data,                 # Data frame
  method = "ML",
  weights = varIdent(form = ~ 1 | ssp)
)

Wall_null <- lme(
  WallThickness ~ 1,                # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Wall_data,                 # Data frame
  method = "ML",
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(Wall_full), digits = 2)
re <- as.numeric(VarCorr(Wall_full)[, "Variance"])
residuals <- round(residuals(Wall_full), digits = 2)

# Calculate CoV for parasitism groups
CoV_parasite <- round(sd(residuals[Wall_data$parasitism == "Parasite"]) / fixed_effects[1], digits = 2)
CoV_host <- round(sd(residuals[Wall_data$parasitism == "Host"]) / sum(fixed_effects), digits = 2)

# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(Wall_null, Wall_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

# Append global model results
VWall_AIC <- rbind(VWall_AIC, data.frame(
  PairTested = paste(levels(Wall_data$parasitism), collapse = " vs "),
  ParasiteMean = paste(fixed_effects[1], "(", CoV_parasite, ")", sep = ""),
  HostMean = paste(sum(fixed_effects), "(", CoV_host, ")", sep = ""),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = fixed_effects[2] / fixed_effects[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(Wall_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      WallThickness ~ ssp,           # Fixed effects
      random = ~ 1 | ssp / indiv,    # Random effects
      data = subset_data,            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      WallThickness ~ 1,             # Fixed effects
      random = ~ 1 | ssp / indiv,    # Random effects
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    # Print model summary
    cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
    print(summary(full_model))
    print(anova(full_model, reduced_model))
    
    # Extract fixed effects, variance, and residuals
    fixed_effects <- round(fixed.effects(full_model), digits = 2)
    re <- as.numeric(VarCorr(full_model)[, "Variance"])
    residuals <- round(residuals(full_model), digits = 2)
    
    # Calculate CoV and other metrics
    CoV_parasite <- round(sd(residuals[subset_data$ssp == pair[1]]) / fixed_effects[1], digits = 2)
    CoV_host <- round(sd(residuals[subset_data$ssp == pair[2]]) / sum(fixed_effects), digits = 2)
    variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])
    
    # Perform likelihood ratio test for the pair
    lrt <- anova(full_model, reduced_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Append results for the species pair
    VWall_AIC <- rbind(VWall_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste(fixed_effects[1], "(", CoV_parasite, ")", sep = ""),
      HostMean = paste(sum(fixed_effects), "(", CoV_host, ")", sep = ""),
      REVariance = variance_explained_by_RE * 100,
      RelDiff = fixed_effects[2] / fixed_effects[1],
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}
rm(VWall_AIC,Wall_data,Wall_full,Wall_null,full_model,reduced_model,lrt)
#############################################################
#Vessel Diameter
VWall_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
Wall_full <- lme(
  WallThickness ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Wall_data,                 # Data frame
  method = "ML",
  weights = varIdent(form = ~ 1 | ssp)
)

Wall_null <- lme(
  WallThickness ~ 1,                # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Wall_data,                 # Data frame
  method = "ML",
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(Wall_full), digits = 2)
re <- as.numeric(VarCorr(Wall_full)[, "Variance"])
residuals <- round(residuals(Wall_full), digits = 2)

# Calculate CoV for parasitism groups
CoV_parasite <- round(sd(residuals[Wall_data$parasitism == "Parasite"]) / fixed_effects[1], digits = 2)
CoV_host <- round(sd(residuals[Wall_data$parasitism == "Host"]) / sum(fixed_effects), digits = 2)

# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(Wall_null, Wall_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

# Append global model results
VWall_AIC <- rbind(VWall_AIC, data.frame(
  PairTested = paste(levels(Wall_data$parasitism), collapse = " vs "),
  ParasiteMean = paste(fixed_effects[1], "(", CoV_parasite, ")", sep = ""),
  HostMean = paste(sum(fixed_effects), "(", CoV_host, ")", sep = ""),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = fixed_effects[2] / fixed_effects[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(Wall_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      WallThickness ~ ssp,           # Fixed effects
      random = ~ 1 | ssp / indiv,    # Random effects
      data = subset_data,            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      WallThickness ~ 1,             # Fixed effects
      random = ~ 1 | ssp / indiv,    # Random effects
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    # Print model summary
    cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
    print(summary(full_model))
    print(anova(full_model, reduced_model))
    
    # Extract fixed effects, variance, and residuals
    fixed_effects <- round(fixed.effects(full_model), digits = 2)
    re <- as.numeric(VarCorr(full_model)[, "Variance"])
    residuals <- round(residuals(full_model), digits = 2)
    
    # Calculate CoV and other metrics
    CoV_parasite <- round(sd(residuals[subset_data$ssp == pair[1]]) / fixed_effects[1], digits = 2)
    CoV_host <- round(sd(residuals[subset_data$ssp == pair[2]]) / sum(fixed_effects), digits = 2)
    variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])
    
    # Perform likelihood ratio test for the pair
    lrt <- anova(full_model, reduced_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Append results for the species pair
    VWall_AIC <- rbind(VWall_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste(fixed_effects[1], "(", CoV_parasite, ")", sep = ""),
      HostMean = paste(sum(fixed_effects), "(", CoV_host, ")", sep = ""),
      REVariance = variance_explained_by_RE * 100,
      RelDiff = fixed_effects[2] / fixed_effects[1],
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}
