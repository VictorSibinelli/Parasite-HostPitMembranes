
### LME analysis
library(glmmTMB)
library(tidyverse)
library(here)
library(nlme)
library(emmeans)
library(performance)
library(predictmeans)
library(DHARMa)
library(sjPlot)
# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data<- read.csv(here("data", "processed", "VesselDiameter_data.csv")) 
Hydraulic_data<- read.csv(here("data", "processed", "HydraulicData.csv"))
PitFraction_data<- read.csv(here("data", "processed", "PitFraction_data.csv"))
PitDiOp_data<- read.csv(here("data", "processed", "PitDiOp_data.csv"))

source(here("scripts", "Functions.R"))
output_dir <- here("outputs", "tables", "Models")

# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

relevel_factors(ls())
library(report)

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


residplot(Wall_full,newwd = F)
title(sub = "Wall Thickness PxH")
check_model(Wall_full,show_dots = F)
CookD(Wall_full, idn = 20,newwd = F)
abline(h=4/nrow(Wall_data),col="red")



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
    file_name <-paste("Tvw",
                      paste0(
                        paste(str_extract(pair, "^\\S+"), collapse = "_"), 
                        "-Tvw.doc"
                      ))
                      file_path <- here(output_dir, file_name)
                      
                      # Create and save the model table
                      tab_m <- sjPlot::tab_model(
                        reduced_model, 
                        full_model,
                        title = paste("Vessel Wall Thickness -", paste(pair, collapse = " vs ")),
                        file = file_path,
                        show.aic = TRUE,
                        show.r2
                        show.reflvl = TRUE,
                        pred.labels = c("Intercept", "Species Effect"),
                        dv.labels = c("Null Model", "Full Model"),
                        p.style = "numeric_stars"
                      )
                      
                      # Explicitly print the table to save it
                      print(tab_m)
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}


print(VWall_AIC)

sjPlot::tab_model(Wall_null,Wall_full,show.aic = T,show.reflvl = T,
                  title = "Vessel Wall Thickness - Parasites vs Hosts",
                  pred.labels = c("Parasitism","Host"),
                  dv.labels = c("Null Model","Full Model"),
                  p.style = "numeric_stars",
                  file = here("outputs","tables","Tvw_tab_model.doc")
                  )

#######################
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

# residplot(VesselDiameter_full,newwd = F)
# title(sub = "VDiameter PxH")
# # check_model(VesselDiameter_full,show_dots = F)


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
    VesselDiameter ~ ssp,
      random = ~ 1 | indiv,
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      VesselDiameter ~ 1,
      random = ~ 1 | indiv,
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    


    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
     # print(check_model(full_model,show_dots = F))

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
plot(VesselDiameter_data$VesselDiameter~VesselDiameter_data$parasitism)
print(VDiameter_AIC)
##############################################################################################
 

TopVDiameter_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)




# Fit global models
TopVesselDiameter_full <- lme(
  log(VesselDiameter) ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = VesselDiameter_data %>% group_by(ssp,indiv) %>%
    filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
    ungroup(),                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

TopVesselDiameter_null <- lme(
  log(VesselDiameter) ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = VesselDiameter_data%>% group_by(ssp,indiv) %>%
    filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
    ungroup(),                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(TopVesselDiameter_full), digits = 2)
re <- as.numeric(VarCorr(TopVesselDiameter_full)[, "Variance"])
variance_explained_by_RE <- (re[2] + re[4]) / (re[5] + re[2] + re[4])

# Calculate 95% Confidence Intervals using intervals() for the lme model
conf_int <- emmeans(TopVesselDiameter_full,~parasitism) %>% summary()

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% exp() %>% round(digits = 2)
CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>%exp() %>%  round(digits = 2) 


# Perform likelihood ratio test (LRT)
lrt <- anova(TopVesselDiameter_null,TopVesselDiameter_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`) 

# residplot(TopVesselDiameter_full,newwd = F)
# title(sub = "TopVDiameter PxH")
# check_model(TopVesselDiameter_full,show_dots = F)
# CookD(TopVesselDiameter_full, idn = 20,newwd = F)
# abline(h=4/nrow(TopVesselDiameter_data),col="red")
# 

TopVDiameter_AIC <- rbind(TopVDiameter_AIC, data.frame(
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
      data = subset_data%>% group_by(ssp,indiv) %>%
        filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
        ungroup(),
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      log(VesselDiameter) ~ 1,
      random = ~ 1 | indiv,
      data = subset_data%>% group_by(ssp,indiv) %>%
        filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
        ungroup(),
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    
    # 
    # residplot(full_model,newwd = F)
    # title(sub = paste0(pair,collapse = " x "))
    # print(check_model(full_model,show_dots = F))
    # print(CookD(full_model, idn = 20,newwd = F))
    # abline(h=4/nrow(subset_data),col="red")
    # title(sub=paste0(pair,collacpse=" x "))
    # 
    fixed_effects <- round(fixed.effects(full_model), digits = 2)
    re <- VarCorr(full_model)[1,"Variance"] %>% as.numeric()
    residuals <- VarCorr(full_model)[2,"Variance"] %>% as.numeric()
    variance_explained_by_RE <-re/sum(as.numeric(VarCorr(full_model)[,"Variance"]))
    conf_int <- emmeans(full_model,~ssp) %>% summary()
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% exp() %>% round(digits = 2)
    CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>%exp() %>%  round(digits = 2) # Adding Intercept to Host effect
    
    # Likelihood ratio test
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Log progress
    cat("\nCompleted pair:", paste(pair, collapse = " vs "), "\n")
    
    # Return results
    TopVDiameter_AIC <-rbind(TopVDiameter_AIC,data.frame(
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

# residplot(Hydraulic_full,newwd = F)
# title(sub = "HDiameter PxH")
# check_model(Hydraulic_full,show_dots = F)
# CookD(Hydraulic_full, idn = 20,newwd = F)
# abline(h=4/nrow(Hydraulic_data),col="red")
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
    # residplot(full_model,newwd = F)
    # title(sub = paste0(pair,collapse = " x "))
    # print(check_model(full_model,show_dots = F))
    # print(CookD(full_model, idn = 20,newwd = F))
    # abline(h=4/nrow(subset_data),col="red")
    # title(sub=paste0(pair,collacpse=" x "))
    # 
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
# 
# residplot(VDensity_glm,newwd = F)
# title(sub = "Vdens PxH")
# check_model(VDensity_glm,show_dots = F)
# 
# simulateResiduals(fittedModel = VDensity_glm) %>% plot()
# testDispersion(simulateResiduals(fittedModel = VDensity_glm))
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
    
    # # Likelihood ratio test
    # residplot(full_model,newwd = F)
    # title(sub = paste0(pair,collapse = " x "))
    # print(check_model(full_model,show_dots = F))
    # 

    
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
VFraction_full <-  glmmTMB(
  (VesselFraction/100) ~ parasitism + (1 | ssp / indiv),  # Fixed effect: species, Random effect: individual
  family = beta_family(link = "logit"),  # Beta regression with logit link
  data = Hydraulic_data
) 

VFraction_null <-  glmmTMB(
  (VesselFraction/100) ~ (1 | ssp / indiv),  # Fixed effect: species, Random effect: individual
  family = beta_family(link = "logit"),  # Beta regression with logit link
  data = Hydraulic_data
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


res <- simulateResiduals(VFraction_full)
plot(res)  # Check for residual patterns

# Model fit assessment
check_model(VFraction_full)
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
# 
# Hydraulic_data[84,"VesselFraction"] <-
#   rev(sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Psittacanthus robustus 2"]))[2]
# Hydraulic_data[222,"VesselFraction"] <- 
#   rev(sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Vochysia thyrsoidea 1"]))[2]
# Hydraulic_data[c(237),"VesselFraction"] <- 
#   rev(sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Vochysia thyrsoidea 2"]))[2]
# Hydraulic_data[c(232,236,238),"VesselFraction"] <- 
#   sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Vochysia thyrsoidea 2"])[4]
# Hydraulic_data[c(239,240),"VesselFraction"] <- 
#   rev(sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Vochysia thyrsoidea 3"]))[3]
# Hydraulic_data[135,"VesselFraction"] <- 
#   rev(sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Tapirira guianensis 2"]))[2]
# Hydraulic_data[c(127,132),"VesselFraction"] <- 
# sort(Hydraulic_data$VesselFraction[Hydraulic_data$indiv=="Tapirira guianensis 2"])[3] 
# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(Hydraulic_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- glmmTMB(
      (VesselFraction/100) ~ ssp + (1 | indiv),  # Fixed effect: species, Random effect: individual
      family = beta_family(link = "logit"),  # Beta regression with logit link
      data = Hydraulic_data
    ) 
    
    reduced_model <- glmmTMB(
      (VesselFraction/100) ~ (1 | indiv),  # Fixed effect: species, Random effect: individual
      family = beta_family(link = "logit"),  # Beta regression with logit link
      data = Hydraulic_data
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
    
    
    res <- simulateResiduals(full_model)
    plot(res)  # Check for residual patterns
    
    check_model(full_model)
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
# Fit the model using glmmTMB

Kmax_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)
Kmax_full <- glmmTMB(
  Kmax ~ parasitism + (1 | ssp/indiv), # Fixed and random effects
  data = Hydraulic_data,              # Data frame
  family = gaussian(link = "log"),    # Use a log link for the response
  control = glmmTMBControl(
    optCtrl = list(iter.max = 300, eval.max = 300) # Max iterations for fitting
  ),
  dispformula = ~ ssp                  # Variance structure for dispersion
)

Kmax_null <- glmmTMB(
  Kmax ~ 1 + (1 | ssp/indiv), # Fixed and random effects
  data = Hydraulic_data,              # Data frame
  family = gaussian(link = "log"),    # Use a log link for the response
  control = glmmTMBControl(
    optCtrl = list(iter.max = 300, eval.max = 300) # Max iterations for fitting
  ),
  dispformula = ~ ssp                  # Variance structure for dispersion
)



# Extract fixed effects, random effects variance, and residuals

fixed_effects <- fixef(Kmax_full)
fixed_effects$cond
# Extract variance components
var_components <- VarCorr(Kmax_full) 

re <- var_components$cond %>% as.data.frame() %>% sum()
residual_variance <- var(exp(residuals(Kmax_full, type = "pearson")), na.rm = TRUE)

total_variance <- re + residual_variance
variance_explained_by_RE <-re/total_variance 

# Extract confidence intervals for fixed effects only
ci_fixed <- confint(Kmax_full, parm = "beta_")
ci_fixed <- ci_fixed[grep("^cond", rownames(ci_fixed)), ]
# View the fixed effects confidence intervals
print(ci_fixed)

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- ci_fixed[1,] %>% exp()%>% round(digits = 2)  # CI for PaDevasite (Intercept)
CI_Host <-(ci_fixed[2, ] + ci_fixed[1, ] )%>% exp() %>% round(digits = 2)

CI_Parasite <- c(CI_Parasite[1], CI_Parasite[3], CI_Parasite[2])
CI_Host <- c(CI_Host[1], CI_Host[3], CI_Host[2])
# Variance explained by random effects


# Perform likelihood ratio test (LRT)
lrt <- anova(Kmax_null,Kmax_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`Pr(>Chisq)`)
residplot(Kmax_full,newwd = F)
check_model(Kmax_full,show_dots = F)

simulateResiduals(fittedModel = Kmax_full) %>% plot() %>% print()
Kmax_AIC <- rbind(Kmax_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism %>% as.factor()), collapse = " vs "),
  ParasiteMean = paste0(CI_Parasite,collapse = "-"),
  HostMean = paste0(CI_Host,collapse = "-"),
  REVariance = variance_explained_by_RE * 100,
  RelDiff = exp(fixed_effects$cond[2]) / exp(fixed_effects$cond[1]),
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))



for (pair in species_pairs) {
  subset_data <- subset(Hydraulic_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- glmmTMB(
      Kmax ~ ssp + (1 |indiv), # Fixed and random effects
      data = subset_data,              # Data frame
      family = gaussian(link = "log"),    # Use a log link for the response
      control = glmmTMBControl(
        optCtrl = list(iter.max = 300, eval.max = 300) # Max iterations for fitting
      ),
      dispformula = ~ indiv                  # Variance structure for dispersion
    )
    
    reduced_model <- glmmTMB(
      Kmax ~ 1+ (1 |indiv), # Fixed and random effects
      data = subset_data,              # Data frame
      family = gaussian(link = "log"),    # Use a log link for the response
      control = glmmTMBControl(
        optCtrl = list(iter.max = 300, eval.max = 300) # Max iterations for fitting
      ),
      dispformula = ~ indiv                  # Variance structure for dispersion
    )
    
    # Print model summary
    cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
    print(summary(full_model))
    print(anova( reduced_model,full_model))
    
    fixed_effects <- fixef(full_model)
    
    # Extract variance components
    var_components <- VarCorr(full_model) 
    
    re <- var_components$cond %>% as.data.frame() %>% sum()
    residual_variance <- var(exp(residuals(full_model, type = "pearson")), na.rm = TRUE)
    
    total_variance <- re + residual_variance
    variance_explained_by_RE <-re/total_variance 
    
    # Extract confidence intervals for fixed effects only
    ci_fixed <- confint(full_model, parm = "beta_")
    ci_fixed <- ci_fixed[grep("^cond", rownames(ci_fixed)), ]
    # View the fixed effects confidence intervals
    print(ci_fixed)
    
    # Extract the CI for the intercept (Parasite) and parasitismHost (Host)
    CI_Parasite <- ci_fixed[1,] %>% exp()%>% round(digits = 2)  # CI for PaDevasite (Intercept)
    CI_Host <-(ci_fixed[2, ] + ci_fixed[1, ] )%>% exp() %>% round(digits = 2)
    
    CI_Parasite <- c(CI_Parasite[1], CI_Parasite[3], CI_Parasite[2])
    CI_Host <- c(CI_Host[1], CI_Host[3], CI_Host[2])
    # Variance explained by random effects
    
    
    # Perform likelihood ratio test (LRT)
    lrt <- anova(reduced_model,full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`Pr(>Chisq)`)
    
    
    # residplot(full_model,newwd = F)
    # title(sub = paste0(pair,collapse = " x "))
    # print(check_model(full_model,show_dots = F))
    # simulateResiduals(fittedModel = full_model) %>% plot() %>% print()
    # title(sub = paste0(pair,collapse = " x "))
    # 
    # Append results for the species pair
    Kmax_AIC <- rbind(Kmax_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = paste0(CI_Parasite,collapse = "-"),
      HostMean = paste0(CI_Host,collapse = "-"),
      REVariance = variance_explained_by_RE * 100,
      RelDiff = exp(fixed_effects$cond[2]) / exp(fixed_effects$cond[1]),
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
    # residplot(full_model,newwd = F)
    # title(sub = paste0(pair,collapse = " x "))
    # print(check_model(full_model,show_dots = F))
    # print(CookD(full_model, idn = 20,newwd = F))
    # abline(h=4/nrow(Hydraulic_data),col="red")
    # title(sub=paste0(pair,collacpse=" x "))
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
residplot(PitOpening_full,newwd = F)
title(sub = paste0("Pit Opening PxH",collapse = " x "))
print(check_model(PitOpening_full,show_dots = F))
print(CookD(PitOpening_full, idn = 20,newwd = F))
abline(h=4/nrow(PitDiOp_data),col="red")
title(sub=paste0(pair,collacpse=" x "))
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


PitDiOp_data[615,"PitOpening"] <- rev(sort(PitDiOp_data$PitOpening[PitDiOp_data$indiv=="Tapirira guianensis 1"]))[2]
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

    
    # residplot(full_model,newwd = F)
    # title(sub = paste0(pair,collapse = " x "))
    # print(check_model(full_model,show_dots = F))
    # print(CookD(full_model, idn = 20,newwd = F))
    # abline(h=4/nrow(Hydraulic_data),col="red")
    # title(sub=paste0(pair,collacpse=" x "))
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    # Append results for the species pair
    PitOpening_AIC <- rbind(PitOpening_AIC, data.frame(
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
  data = PitFraction_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

PitFraction_null <- lme(
  PitFraction ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = PitFraction_data,                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Extract fixed effects, random effects variance, and residuals
fixed_effects <- round(fixed.effects(PitFraction_full), digits = 2)
# Extract variances from VarCorr and remove NA values
re <- VarCorr(PitFraction_full)[,"Variance"] %>%
  as.numeric() %>%
  .[!is.na(.)]

# Print the result


residuals <- re[3]
variance_explained_by_RE <-(re[1]+re[2])/sum(re)
conf_int <- emmeans(PitFraction_full,~parasitism) %>% summary()

# Extract the CI for the intercept (Parasite) and parasitismHost (Host)
CI_Parasite <- c(conf_int$lower.CL[1],fixed_effects[1],conf_int$upper.CL[1]) %>% round(digits = 2)
CI_Host <-c(conf_int$lower.CL[2],sum(fixed_effects),conf_int$upper.CL[2]) %>% round(digits = 2) # Adding Intercept to Host effect




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
residplot(PitFraction_full,newwd = F)
title(sub = paste0("Pit fraction PxH",collapse = " x "))
print(check_model(PitFraction_full,show_dots = F))
print(CookD(PitFraction_full, idn = 20,newwd = F))
abline(h=4/nrow(PitDiOp_data),col="red")
title(sub=paste0(pair,collacpse=" x "))


PitFraction_data[161,"PitFraction"] <- rev(sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Viscum album 1"]))[2]
PitFraction_data[81,"PitFraction"] <- rev(sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Struthanthus rhynchophyllus 2"]))[2]
PitFraction_data[85,"PitFraction"] <- sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Struthanthus rhynchophyllus 3"])[2]
PitFraction_data[c(146,151),"PitFraction"] <- sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Tipuana tipu 3"])[3]
PitFraction_data[19,"PitFraction"] <- sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Phoradendron perrotettii 2"])[2]
PitFraction_data[109,"PitFraction"] <- sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Tapirira guianensis 2"])[2]

PitFraction_data[34,"PitFraction"] <- sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Psittacanthus robustus 1"])[2]
PitFraction_data[56,"PitFraction"] <- rev(sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Psittacanthus robustus 3"]))[2]
PitFraction_data[c(44,47),"PitFraction"] <- sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Psittacanthus robustus 2"])[3]
PitFraction_data[c(185,195),"PitFraction"] <- sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Vochysia thyrsoidea 1"])[3]
PitFraction_data[191,"PitFraction"] <- rev(sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Vochysia thyrsoidea 1"]))[2]
PitFraction_data[c(202),"PitFraction"] <- sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Vochysia thyrsoidea 2"])[2]
PitFraction_data[c(206),"PitFraction"] <- rev(sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Vochysia thyrsoidea 2"]))[2]
PitFraction_data[c(208,209,210),"PitFraction"] <- rev(sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Vochysia thyrsoidea 3"]))[4]
PitFraction_data[c(214),"PitFraction"] <- sort(PitFraction_data$PitFraction[PitFraction_data$indiv=="Vochysia thyrsoidea 3"])[2]







# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(PitFraction_data, ssp %in% pair)
  
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



AIC_tables <- list(VWall_AIC,
                   VDiameter_AIC,
                   TopVDiameter_AIC,
                   HDiameter_AIC,
                   VDensity_AIC,
                   VFraction_AIC,
                   Kmax_AIC,
                   PitDiameter_AIC,
                   PitOpening_AIC,
                   PitFraction_AIC
)
# Assign names to each element in the list
names(AIC_tables) <- c(
  "VWall_AIC",
  "VDiameter_AIC",
  "TopVDiameter_AIC",
  "HDiameter_AIC",
  "VDensity_AIC",
  "VFraction_AIC",
  "Kmax_AIC",
  "PitDiameter_AIC",
  "PitOpening_AIC",
  "PitFraction_AIC"
)

# Save each element in the list as a CSV file
lapply(names(AIC_tables), function(name) {
  write.csv(
    AIC_tables[[name]],                 # The data to write
    file = here("outputs","tables",paste0(name, ".csv")),        # File name based on the element name
    row.names = FALSE                   # Exclude row names
  )
})
