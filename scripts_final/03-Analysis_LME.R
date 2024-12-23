### LME analysis



library(tidyverse)
library(here)
library(nlme)
library(predictmeans)
library(performance)
source(here("scripts","00-library.R"))

# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data<- read.csv(here("data", "processed", "VesselDiameter_data.csv")) %>% group_by(ssp) %>%
  filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
  ungroup()
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
    print(anova(reduced_model,full_model))
    
    # Extract fixed effects, variance, and residuals
    fixed_effects <- round(fixed.effects(full_model), digits = 2)
    re <- as.numeric(VarCorr(full_model)[, "Variance"])
    residuals <- round(residuals(full_model), digits = 2)
    
    # Calculate CoV and other metrics
    CoV_parasite <- round(sd(residuals[subset_data$ssp == pair[1]]) / fixed_effects[1], digits = 2)
    CoV_host <- round(sd(residuals[subset_data$ssp == pair[2]]) / sum(fixed_effects), digits = 2)
    variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])
    
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model,full_model)
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
residplot(full_model)
check_model(full_model)
CookD(full_model,newwd = F, idn = 10)
abline(h=4/nrow(subset_data))

rm(VWall_AIC,Wall_data,Wall_full,Wall_null,full_model,reduced_model,lrt)
#############################################################
####################################################################


#Vessel Diameter


VDiameter_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)




# Fit global models
VesselDiameter_full <- lme(
  VesselDiameter ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = VesselDiameter_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

VesselDiameter_null <- lme(
  VesselDiameter ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = VesselDiameter_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(VesselDiameter_full), digits = 2)
re <- as.numeric(VarCorr(VesselDiameter_full)[, "Variance"])
residuals <- round(residuals(VesselDiameter_full), digits = 2)

# Calculate CoV for parasitism groups
CoV_parasite <- round(sd(residuals[VesselDiameter_data$parasitism == "Parasite"]) / fixed_effects[1], digits = 2)
CoV_host <- round(sd(residuals[VesselDiameter_data$parasitism == "Host"]) / sum(fixed_effects), digits = 2)

# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(VesselDiameter_null,VesselDiameter_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

# Append global model results

residplot(VesselDiameter_full,newwd = F)
title(sub = "VDiameter PxH")
check_model(VesselDiameter_full,show_dots = F)
CookD(VesselDiameter_full, idn = 20,newwd = F)
abline(h=4/nrow(VesselDiameter_data),col="red")

#High deviation form normality of resisuals, but low leverage
# Set up parallel backend
library(parallel)
library(doParallel)
library(foreach)
num_cores <- detectCores() - 2 # Use one less core than available
cl <- makeCluster(num_cores)  # Create a cluster
registerDoParallel(cl)
# Initialize results dataframe
VDiameter_AIC <- foreach(pair = species_pairs, .combine = rbind, .packages = c("nlme", "performance")) %dopar% {
  subset_data <- subset(VesselDiameter_data, ssp %in% pair)
  
  tryCatch({
    # Fit models
    full_model <- lme(
      VesselDiameter ~ ssp,
      random = ~ 1 | ssp / indiv,
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      VesselDiameter ~ 1,
      random = ~ 1 | ssp / indiv,
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    # Extract metrics
    fixed_effects <- round(fixed.effects(full_model), digits = 2)
    re <- as.numeric(VarCorr(full_model)[, "Variance"])
    residuals <- round(residuals(full_model), digits = 2)
    CoV_parasite <- round(sd(residuals[subset_data$ssp == pair[1]]) / fixed_effects[1], digits = 2)
    CoV_host <- round(sd(residuals[subset_data$ssp == pair[2]]) / sum(fixed_effects), digits = 2)
    variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])
    
    # Likelihood ratio test
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Log progress
    cat("\nCompleted pair:", paste(pair, collapse = " vs "), "\n")
    
    # Return results
    data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste(fixed_effects[1], "(", CoV_parasite, ")", sep = ""),
      HostMean = paste(sum(fixed_effects), "(", CoV_host, ")", sep = ""),
      REVariance = variance_explained_by_RE * 100,
      RelDiff = fixed_effects[2] / fixed_effects[1],
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
    NULL
  })
}
VDiameter_AIC <- rbind(VDiameter_AIC, data.frame(
  PairTested = paste(levels(VesselDiameter_data$parasitism), collapse = " vs "),
  ParasiteMean = paste(fixed_effects[1], "(", CoV_parasite, ")", sep = ""),
  HostMean = paste(sum(fixed_effects), "(", CoV_host, ")", sep = ""),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = fixed_effects[2] / fixed_effects[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

##very big model assumptions deviations
###############################################
################################################

#Hydraulic diameter

HDiameter_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
Hydraulic_full <- lme(
  HydraulicDiameter ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

Hydraulic_null <- lme(
  HydraulicDiameter ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(Hydraulic_full), digits = 2)
re <- as.numeric(VarCorr(Hydraulic_full)[, "Variance"])
residuals <- round(residuals(Hydraulic_full), digits = 2)

# Calculate CoV for parasitism groups
CoV_parasite <- round(sd(residuals[Hydraulic_data$parasitism == "Parasite"]) / fixed_effects[1], digits = 2)
CoV_host <- round(sd(residuals[Hydraulic_data$parasitism == "Host"]) / sum(fixed_effects), digits = 2)

# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(Hydraulic_null,Hydraulic_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

residplot(Hydraulic_full,newwd = F)
title(sub = "HDiameter PxH")
check_model(Hydraulic_full,show_dots = F)
CookD(Hydraulic_full, idn = 20,newwd = F)
abline(h=4/nrow(Hydraulic_data),col="red")
# Append global model results
HDiameter_AIC <- rbind(HDiameter_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
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
  subset_data <- subset(Hydraulic_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      HydraulicDiameter ~ ssp,           # Fixed effects
      random = ~ 1 | ssp / indiv,    # Random effects
      data = subset_data,            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      HydraulicDiameter ~ 1,             # Fixed effects
      random = ~ 1 | ssp / indiv,    # Random effects
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    # Print model summary
    cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
    print(summary(full_model))
    print(anova( reduced_model,full_model))
    
    # Extract fixed effects, variance, and residuals
    fixed_effects <- round(fixed.effects(full_model), digits = 2)
    re <- as.numeric(VarCorr(full_model)[, "Variance"])
    residuals <- round(residuals(full_model), digits = 2)
    
    # Calculate CoV and other metrics
    CoV_parasite <- round(sd(residuals[subset_data$ssp == pair[1]]) / fixed_effects[1], digits = 2)
    CoV_host <- round(sd(residuals[subset_data$ssp == pair[2]]) / sum(fixed_effects), digits = 2)
    variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])
    
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    residplot(full_model,newwd = F)
    title(sub = pair)
    print(check_model(full_model,show_dots = F))
    print(CookD(full_model, idn = 20,newwd = F))
    abline(h=4/nrow(Hydraulic_data),col="red")
    title(sub=paste0(pair,collacpse=" x "))
    # Append results for the species pair
    HDiameter_AIC <- rbind(HDiameter_AIC, data.frame(
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
#large deviations from assumptions
#####





#####Vessel Dendity
VDensity_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
VDensity_full <- lme(
  VesselDensity ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

VDensity_null <- lme(
  VesselDensity ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(VDensity_full), digits = 2)
re <- as.numeric(VarCorr(VDensity_full)[, "Variance"])
residuals <- round(residuals(VDensity_full), digits = 2)

# Calculate CoV for parasitism groups
CoV_parasite <- round(sd(residuals[Hydraulic_data$parasitism == "Parasite"]) / fixed_effects[1], digits = 2)
CoV_host <- round(sd(residuals[Hydraulic_data$parasitism == "Host"]) / sum(fixed_effects), digits = 2)

# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(VDensity_null,VDensity_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

# Append global model results
VDensity_AIC <- rbind(VDensity_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
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
  subset_data <- subset(Hydraulic_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      VesselDensity ~ ssp,           # Fixed effects
      random = ~ 1 | ssp / indiv,    # Random effects
      data = subset_data,            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      VesselDensity ~ 1,             # Fixed effects
      random = ~ 1 | ssp / indiv,    # Random effects
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    # Print model summary
    cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
    print(summary(full_model))
    print(anova( reduced_model,full_model))
    
    # Extract fixed effects, variance, and residuals
    fixed_effects <- round(fixed.effects(full_model), digits = 2)
    re <- as.numeric(VarCorr(full_model)[, "Variance"])
    residuals <- round(residuals(full_model), digits = 2)
    
    # Calculate CoV and other metrics
    CoV_parasite <- round(sd(residuals[subset_data$ssp == pair[1]]) / fixed_effects[1], digits = 2)
    CoV_host <- round(sd(residuals[subset_data$ssp == pair[2]]) / sum(fixed_effects), digits = 2)
    variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])
    
    
    residplot(full_model,newwd = F)
    title(sub = pair)
    print(check_model(full_model,show_dots = F))
    print(CookD(full_model, idn = 20,newwd = F))
    abline(h=4/nrow(Hydraulic_data),col="red")
    title(sub=paste0(pair,collacpse=" x "))
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Append results for the species pair
    VDensity_AIC <- rbind(VDensity_AIC, data.frame(
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

