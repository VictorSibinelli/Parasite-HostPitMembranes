### LME analysis


library(glmmTMB)
library(tidyverse)
library(here)
library(nlme)
#library(predictmeans)
library(emmeans)
library(performance)
source(here("scripts","00-library.R"))

# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data<- read.csv(here("data", "processed", "VesselDiameter_data.csv")) %>% group_by(ssp,indiv) %>%
  filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
  ungroup()
Hydraulic_data<- read.csv(here("data", "processed", "HydraulicData.csv"))
PitFraction_data<- read.csv(here("data", "processed", "PitFraction_data.csv"))
PitDiOp_data<- read.csv(here("data", "processed", "PitDiOp_data.csv"))

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


# Calculate 95% Confidence Intervals using intervals() for the lme model
conf_int <- emmeans(Wall_full,~parasitism) %>% summary()

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect



# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(Wall_null, Wall_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

# Append global model results
VWall_AIC <- rbind(VWall_AIC, data.frame(
  PairTested = paste(levels(Wall_data$parasitism), collapse = " vs "),
  ParasiteMean = paste0(CI_Parasite,collapse= "-"),
  HostMean = paste0(CI_Host, collapse  = "-"),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = fixed_effects[2] / fixed_effects[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))


# residplot(Wall_full,newwd = F)
# title(sub = "Wall Thickness PxH")
# check_model(Wall_full,show_dots = F)
# CookD(Wall_full, idn = 20,newwd = F)
# abline(h=4/nrow(Wall_data),col="red")

Wall_data[510,"WallThickness"] <- rev(sort(Wall_data$WallThickness[Wall_data$indiv=="Tipuana tipu 2"]))[2]

# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(Wall_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      WallThickness ~ ssp,           # Fixed effects
      random = ~ 1 |indiv,    # Random effects
      data = subset_data,            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      WallThickness ~ 1,             # Fixed effects
      random = ~ 1 | indiv,    # Random effects
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
    re <- VarCorr(full_model)[1,"Variance"] %>% as.numeric()
    residuals <- VarCorr(full_model)[2,"Variance"] %>% as.numeric()
    variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))
    
    # Calculate 95% Confidence Intervals using intervals() for the lme model
    conf_int <- emmeans(full_model,~ssp) %>% summary()
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
    CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect
    
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model,full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    # 
    # residplot(full_model,newwd = F)
    # title(sub = paste0(pair,collapse = " x "))
    # print(check_model(full_model,show_dots = F))
    # print(CookD(full_model, idn = 20,newwd = F))
    # abline(h=4/nrow(subset_data),col="red")
    # title(sub=paste0(pair,collacpse=" x "))
    # 
    # Append results for the species pair
    VWall_AIC <- rbind(VWall_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste0(CI_Parasite,collapse = "-"),
      HostMean = paste0(CI_Host,collapse = "-"),
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


print(VWall_AIC)



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
  log(VesselDiameter) ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = VesselDiameter_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

VesselDiameter_null <- lme(
  log(VesselDiameter) ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = VesselDiameter_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(VesselDiameter_full), digits = 2)
re <- as.numeric(VarCorr(VesselDiameter_full)[, "Variance"])
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Calculate 95% Confidence Intervals using intervals() for the lme model
conf_int <- emmeans(VesselDiameter_full,~parasitism) %>% summary()

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect


# Perform likelihood ratio test (LRT)
lrt <- anova(VesselDiameter_null,VesselDiameter_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

residplot(VesselDiameter_full,newwd = F)
title(sub = "VDiameter PxH")
check_model(VesselDiameter_full,show_dots = F)
CookD(VesselDiameter_full, idn = 20,newwd = F)
abline(h=4/nrow(VesselDiameter_data),col="red")


VDiameter_AIC <- rbind(VDiameter_AIC, data.frame(
  PairTested = paste(levels(VesselDiameter_data$parasitism), collapse = " vs "),
  ParasiteMean = paste0(CI_Parasite,collapse = "-"),
  HostMean = paste0(CI_Host,collapse = "-"),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = fixed_effects[2] / fixed_effects[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))





for (pair in species_pairs) {
  subset_data <- subset(VesselDiameter_data, ssp %in% pair) 
  tryCatch({
    # Fit models
    full_model <- lme(
      log(VesselDiameter) ~ ssp,
      random = ~ 1 | indiv,
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      log(VesselDiameter) ~ 1,
      random = ~ 1 | indiv,
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    

    
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))
    print(CookD(full_model, idn = 20,newwd = F))
    abline(h=4/nrow(subset_data),col="red")
    title(sub=paste0(pair,collacpse=" x "))
    
    fixed_effects <- round(fixed.effects(full_model), digits = 2)
    re <- VarCorr(full_model)[1,"Variance"] %>% as.numeric()
    residuals <- VarCorr(full_model)[2,"Variance"] %>% as.numeric()
    variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))
    conf_int <- emmeans(full_model,~ssp) %>% summary()
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
    CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect
    
    # Likelihood ratio test
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Log progress
    cat("\nCompleted pair:", paste(pair, collapse = " vs "), "\n")
    
    # Return results
    VDiameter_AIC <-    rbind(VDiameter_AIC,data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste0(CI_Parasite,collapse = "-"),
      HostMean = paste0(CI_Host,collapse = "-"),
      REVariance = variance_explained_by_RE * 100,
      RelDiff = fixed_effects[2] / fixed_effects[1],
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
    NULL
  })
}

 



###############################################
################################################

#Hydraulic diameter
# Initialize HDiameter_AIC with appropriate column names and types
HDiameter_AIC <- data.frame(
  PairTested = character(),
  ParasiteMean = character(),
  HostMean = character(),
  REVariance = numeric(),
  RelDiff = numeric(),
  DeltaAIC = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
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
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

conf_int <- emmeans(Hydraulic_full,~parasitism) %>% summary()

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect

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
  ParasiteMean = paste0(CI_Parasite,collapse = "-"),
  HostMean = paste0(CI_Host,collapse = "-"),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = fixed_effects[2] / fixed_effects[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))
Hydraulic_data[61,"HydraulicDiameter"] <- rev(sort(Hydraulic_data$HydraulicDiameter[Hydraulic_data$indiv=="Populus nigra 1"]))[2]
# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(Hydraulic_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      HydraulicDiameter ~ ssp,           # Fixed effects
      random = ~ 1 |indiv,    # Random effects
      data = subset_data,            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      HydraulicDiameter ~ 1,             # Fixed effects
      random = ~ 1 | indiv,    # Random effects
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    # Extract fixed effects, variance, and residuals
    fixed_effects <- round(fixed.effects(full_model), digits = 2)
    re <- VarCorr(full_model)[1,"Variance"] %>% as.numeric()

    variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))
    conf_int <- emmeans(full_model,~ssp) %>% summary()
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
    CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect
    
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))
    print(CookD(full_model, idn = 20,newwd = F))
    abline(h=4/nrow(subset_data),col="red")
    title(sub=paste0(pair,collacpse=" x "))
    
    # Append results for the species pair
    new_row <- data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste0(CI_Parasite,collapse = "-"),
      HostMean = paste0(CI_Host,collapse = "-"),
      REVariance = variance_explained_by_RE * 100,
      RelDiff = fixed_effects[2] / fixed_effects[1],
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    )
    
    HDiameter_AIC <- rbind(HDiameter_AIC, new_row)
    
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}







#####Vessel Dendity
VDensity_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

VDensity_glm <- glmmTMB(
  VesselDensity ~ parasitism + (1 | ssp / indiv),  # Fixed and random effects
  data = Hydraulic_data,                          # Data frame
  family = nbinom2(),                             # Poisson family for count data
  control = glmmTMBControl(optCtrl = list(iter.max = 150, eval.max = 150))
)
VDensity_glm_null <- glmmTMB(
  VesselDensity ~ 1 + (1 | ssp / indiv),  # Fixed and random effects
  data = Hydraulic_data,                          # Data frame
  family = nbinom2(),                             # Poisson family for count data
  control = glmmTMBControl(optCtrl = list(iter.max = 150, eval.max = 150))
)
lrt <- anova(VDensity_glm_null,VDensity_glm)

# Extract variance components
var_components <- VarCorr(VDensity_glm)$cond

# Variance explained by "indiv:ssp"
indiv_var <- var_components$`indiv:ssp`[1, 1]
ssp_var <- var_components$`ssp`[1, 1]
v <- summary(VDensity_glm)
# Variance explained by "ssp"
resid_var <-v$sigma

# Residual variance (fixed at 1 for negative binomial on the log scale)

# Calculate total variance
total_variance <- indiv_var+resid_var+ssp_var

# Percentage of variance explained by each random effect
variance_explained_by_RE <- (indiv_var+ssp_var)/total_variance


fixed_effects <- fixef(VDensity_glm)$cond

# Confidence intervals for fixed effects
ci_fixed <- confint(VDensity_glm, parm = "beta_", level = 0.95)

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- ci_fixed[1,] %>% exp()%>% round(digits = 2)  # CI for Parasite (Intercept)
CI_Parasite <- c(CI_Parasite[1], CI_Parasite[3], CI_Parasite[2])
CI_Host <-(ci_fixed[1,]+ci_fixed[2,]) %>% exp() %>% round(digits=2)
CI_Host <- c(CI_Host[1], CI_Host[3], CI_Host[2])
# Variance explained by random effects


# Perform likelihood ratio test (LRT)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`Pr(>Chisq)`)

# Append global model results
VDensity_AIC <- rbind(VDensity_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
  ParasiteMean = paste0(CI_Parasite,collapse = "-"),
  HostMean = paste0(CI_Host,collapse = "-"),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = 100* (CI_Parasite["Estimate"]-CI_Host["Estimate"])/CI_Parasite["Estimate"] ,
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

residplot(VDensity_glm,newwd = F)
title(sub = "Vdens PxH")
check_model(VDensity_glm,show_dots = F)

simulateResiduals(fittedModel = VDensity_glm) %>% plot()
testDispersion(simulateResiduals(fittedModel = VDensity_glm))
# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(Hydraulic_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- glmmTMB(
      VesselDensity ~ ssp + (1 | indiv),  # Fixed and random effects
      data = subset_data,                          # Data frame
      family = nbinom2(),                             # Poisson family for count data
      control = glmmTMBControl(optCtrl = list(iter.max = 150, eval.max = 150))
    )
    
    reduced_model <-glmmTMB(
      VesselDensity ~ 1 + (1 |indiv),  # Fixed and random effects
      data = subset_data,                          # Data frame
      family = nbinom2(),                             # Poisson family for count data
      control = glmmTMBControl(optCtrl = list(iter.max = 150, eval.max = 150))
    )
    
    lrt <- anova(reduced_model,full_model)
    
    
    # Extract variance components
    var_components <- VarCorr(full_model)$cond
    
    # Variance explained by "indiv:ssp"
    indiv_var <- var_components$`indiv`[1, 1]

    v <- summary(full_model)
    # Variance explained by "ssp"
    resid_var <-v$sigma
    
    # Residual variance (fixed at 1 for negative binomial on the log scale)
    
    # Calculate total variance
    total_variance <- indiv_var+resid_var
    
    # Percentage of variance explained by each random effect
    variance_explained_by_RE <- (indiv_var)/total_variance
    

    
    
    fixed_effects <- fixef(full_model)$cond
    
    # Confidence intervals for fixed effects
    ci_fixed <- confint(full_model, parm = "beta_", level = 0.95)
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- ci_fixed[1,] %>% exp()%>% round(digits = 2)  # CI for Parasite (Intercept)
    CI_Parasite <- c(CI_Parasite[1], CI_Parasite[3], CI_Parasite[2])
    CI_Host <-(ci_fixed[1,]+ci_fixed[2,]) %>% exp() %>% round(digits=2)
    CI_Host <- c(CI_Host[1], CI_Host[3], CI_Host[2])
    # Variance explained by random effects
    
    
    # Perform likelihood ratio test (LRT)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`Pr(>Chisq)`)
    
    # Likelihood ratio test
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))

    simulateResiduals(fittedModel = full_model) %>% plot() %>% print()
    testDispersion(simulateResiduals(fittedModel = full_model)) %>% print()
    
    # Append global model results
    VDensity_AIC <- rbind(VDensity_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste0(CI_Parasite,collapse = "-"),
      HostMean = paste0(CI_Host,collapse = "-"),
      REVariance = variance_explained_by_RE*100,
      RelDiff = 100* (CI_Parasite["Estimate"]-CI_Host["Estimate"])/CI_Parasite["Estimate"] ,
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))

  })
}








###################################################
###Vessel Fraction
VFraction_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
VFraction_full <- lme(
  VesselFraction ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

VFraction_null <- lme(
  VesselFraction ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(VFraction_full), digits = 2)
re <- as.numeric(VarCorr(VFraction_full)[, "Variance"])
conf_int <- emmeans(VFraction_full,~parasitism) %>% summary()

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect

# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(VFraction_null,VFraction_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

residplot(VFraction_full,newwd = F)
title(sub = "VFraction PxH")
check_model(VFraction_full,show_dots = F)
CookD(VFraction_full, idn = 20,newwd = F)
abline(h=4/nrow(Hydraulic_data),col="red")

# Append global model results
VFraction_AIC <- rbind(VFraction_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
  ParasiteMean = paste0(CI_Parasite,collapse = "-"),
  HostMean = paste0(CI_Host,collapse = "-"),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = fixed_effects[2] / fixed_effects[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

Hydraulic_data[84,"VesselFraction"] <-
  rev(sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Psittacanthus robustus 2"]))[2]
Hydraulic_data[222,"VesselFraction"] <- 
  rev(sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Vochysia thyrsoidea 1"]))[2]
Hydraulic_data[c(237),"VesselFraction"] <- 
  rev(sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Vochysia thyrsoidea 2"]))[2]
Hydraulic_data[c(232,236,238),"VesselFraction"] <- 
  sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Vochysia thyrsoidea 2"])[4]
Hydraulic_data[c(239,240),"VesselFraction"] <- 
  rev(sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Vochysia thyrsoidea 3"]))[3]
Hydraulic_data[135,"VesselFraction"] <- 
  rev(sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Tapirira guianensis 2"]))[2]
Hydraulic_data[c(127,132),"VesselFraction"] <- 
sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Tapirira guianensis 2"])[3] 
# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(Hydraulic_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      VesselFraction ~ ssp,           # Fixed effects
      random = ~ 1 | indiv,    # Random effects
      data = subset_data,            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      VesselFraction ~ 1,             # Fixed effects
      random = ~ 1 |indiv,    # Random effects
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
    re <- VarCorr(full_model)[1,"Variance"] %>% as.numeric()

    variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))
    
    conf_int <- emmeans(full_model,~ssp) %>% summary()
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
    CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect
    
    
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))
    print(CookD(full_model, idn = 20,newwd = F))
    abline(h=4/nrow(Hydraulic_data),col="red")
    title(sub=paste0(pair,collacpse=" x "))
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Append results for the species pair
    VFraction_AIC <- rbind(VFraction_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste0(CI_Parasite,collapse = "-"),
      HostMean = paste0(CI_Host,collapse = "-"),
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






#######################################
#Kmax
Kmax_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
Kmax_full <- lme(
  log(Kmax) ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

Kmax_null <- lme(
 log(Kmax) ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(Kmax_full), digits = 2)
re <- as.numeric(VarCorr(Kmax_full)[, "Variance"])
residuals <- round(residuals(Kmax_full), digits = 2)

conf_int <- emmeans(Kmax_full,~parasitism) %>% summary()

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect

# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(Kmax_null,Kmax_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

# Append global model results
Kmax_AIC <- rbind(Kmax_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism %>% as.factor()), collapse = " vs "),
  ParasiteMean = paste0(CI_Parasite,collapse = "-"),
  HostMean = paste0(CI_Host,collapse = "-"),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = fixed_effects[2] / fixed_effects[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))
residplot(VFraction_full,newwd = F)
title(sub = "VFraction PxH")
check_model(VFraction_full,show_dots = F)
CookD(VFraction_full, idn = 20,newwd = F)
abline(h=4/nrow(Hydraulic_data),col="red")
# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(Hydraulic_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      log(Kmax) ~ ssp,           # Fixed effects
      random = ~ 1 |indiv,    # Random effects
      data = subset_data,            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      log(Kmax) ~ 1,             # Fixed effects
      random = ~ 1 | indiv,    # Random effects
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
    re <- VarCorr(full_model)[1,"Variance"] %>% as.numeric()

    variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))
    conf_int <- emmeans(full_model,~ssp) %>% summary()
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
    CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect
    
    # Likelihood ratio test
    
    
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))
    print(CookD(full_model, idn = 20,newwd = F))
    abline(h=4/nrow(Hydraulic_data),col="red")
    title(sub=paste0(pair,collacpse=" x "))
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Append results for the species pair
    Kmax_AIC <- rbind(Kmax_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste0(CI_Parasite,collapse = "-"),
      HostMean = paste0(CI_Host,collapse = "-"),
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




############################################
####Pit diameter
PitDiameter_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)
PitDiOp_data[59,]


# Fit global models
PitDiameter_full <- lme(
  PitDiameter ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = PitDiOp_data %>% na.omit(),                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

PitDiameter_null <- lme(
  PitDiameter ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = PitDiOp_data %>% na.omit,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(PitDiameter_full), digits = 2)
re <- as.numeric(VarCorr(PitDiameter_full)[, "Variance"])
conf_int <- emmeans(PitDiameter_full,~parasitism) %>% summary()

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect

# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(PitDiameter_null,PitDiameter_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

# Append global model results
PitDiameter_AIC <- rbind(PitDiameter_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
  ParasiteMean = paste0(CI_Parasite,collapse = "-"),
  HostMean = paste0(CI_Host,collapse = "-"),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = fixed_effects[2] / fixed_effects[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

residplot(PitDiameter_full,newwd = F)
title(sub = paste0("Pit diam PxH",collapse = " x "))
print(check_model(PitDiameter_full,show_dots = F))
print(CookD(PitDiameter_full, idn = 20,newwd = F))
abline(h=4/nrow(PitDiOp_data),col="red")
title(sub=paste0(pair,collacpse=" x "))

# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(PitDiOp_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      PitDiameter ~ ssp,           # Fixed effects
      random = ~ 1 | indiv,    # Random effects
      data = subset_data %>% na.omit(),            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
     
    )
    
    reduced_model <- lme(
      PitDiameter ~ 1,             # Fixed effects
      random = ~ 1 |indiv,    # Random effects
      data = subset_data %>% na.omit(),
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
    re <- VarCorr(full_model)[1,"Variance"] %>% as.numeric()

    variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))
    # Calculate CoV and other metrics

    conf_int <- emmeans(full_model,~ssp) %>% summary()
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
    CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect
    
    # Likelihood ratio test
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))
    print(CookD(full_model, idn = 20,newwd = F))
    abline(h=4/nrow(Hydraulic_data),col="red")
    title(sub=paste0(pair,collacpse=" x "))
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Append results for the species pair
    PitDiameter_AIC <- rbind(PitDiameter_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste0(CI_Parasite,collapse = "-"),
      HostMean = paste0(CI_Host,collapse = "-"),
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


###########################################
##PitOpening
PitOpening_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
PitOpening_full <- lme(
  PitOpening ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data =PitDiOp_data %>% na.omit(),                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

PitOpening_null <- lme(
  PitOpening ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = PitDiOp_data %>% na.omit(),                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(PitOpening_full), digits = 2)
re <- as.numeric(VarCorr(PitOpening_full)[, "Variance"])
residuals <- round(residuals(PitOpening_full), digits = 2)
conf_int <- emmeans(PitOpening_full,~parasitism) %>% summary()

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect

# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(PitOpening_null,PitOpening_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

# Append global model results
PitOpening_AIC <- rbind(PitOpening_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
  ParasiteMean = paste0(CI_Parasite,collapse = "-"),
  HostMean = paste0(CI_Host,collapse = "-"),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = fixed_effects[2] / fixed_effects[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(PitDiOp_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      PitOpening ~ ssp,           # Fixed effects
      random = ~ 1 | indiv,    # Random effects
      data = subset_data %>% na.omit,            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      PitOpening ~ 1,             # Fixed effects
      random = ~ 1 | indiv,    # Random effects
      data = subset_data %>% na.omit(),
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
    re <- VarCorr(full_model)[1,"Variance"] %>% as.numeric()
    conf_int <- emmeans(full_model,~ssp) %>% summary()
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
    CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect
    
    # Likelihood ratio test
    variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))

    
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))
    print(CookD(full_model, idn = 20,newwd = F))
    abline(h=4/nrow(Hydraulic_data),col="red")
    title(sub=paste0(pair,collacpse=" x "))
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Append results for the species pair
    PitOpening_AIC <- rbind(PitPitDiameter_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste0(CI_Parasite,collapse = "-"),
      HostMean = paste0(CI_Host,collapse = "-"),
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

###############################################
##pitFraction


PitFraction_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
PitFraction_full <- lme(
  PitFraction ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

PitFraction_null <- lme(
  PitFraction ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(PitFraction_full), digits = 2)
re <- VarCorr(full_model)[1,"Variance"] %>% as.numeric()
residuals <- VarCorr(full_model)[2,"Variance"] %>% as.numeric()
variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))
conf_int <- emmeans(PitFraction_full,~parasitism) %>% summary()

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect

# Variance explained by random effects
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Perform likelihood ratio test (LRT)
lrt <- anova(PitFraction_null,PitFraction_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

# Append global model results
PitFraction_AIC <- rbind(PitFraction_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
  ParasiteMean = paste0(CI_Parasite,collapse = "-"),
  HostMean = paste0(CI_Host,collapse = "-"),
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
      PitFraction ~ ssp,           # Fixed effects
      random = ~ 1 |indiv,    # Random effects
      data = subset_data,            # Subset data
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      PitFraction ~ 1,             # Fixed effects
      random = ~ 1 |indiv,    # Random effects
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
    re <- VarCorr(full_model)[1,"Variance"] %>% as.numeric()
  
    variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))
    conf_int <- emmeans(full_model,~ssp) %>% summary()
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
    CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect
    
    
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))
    print(CookD(full_model, idn = 20,newwd = F))
    abline(h=4/nrow(Hydraulic_data),col="red")
    title(sub=paste0(pair,collacpse=" x "))
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Append results for the species pair
    PitFraction_AIC <- rbind(PitFraction_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste0(CI_Parasite,collapse = "-"),
      HostMean = paste0(CI_Host,collapse = "-"),
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