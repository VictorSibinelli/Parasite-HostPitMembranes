
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
library(stringr)

# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data<- read.csv(here("data", "processed", "VesselDiameter_data.csv")) 
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

output_dir <- here("outputs", "tables", "Models")

# Initialize results dataframe
VWall_AIC <- data.frame(
  PairTested = character(),
  RelDiff = numeric(),
  DeltaAIC = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
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
summary(Wall_full)

lrt <- anova(Wall_null, Wall_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

residplot(Wall_full,newwd = F)
title(sub = "Wall Thickness PxH")
check_model(Wall_full,show_dots = F)
CookD(Wall_full, idn = 20,newwd = F)
abline(h=4/nrow(Wall_data),col="red")


# Append global model results
VWall_AIC <- rbind(VWall_AIC, data.frame(
  PairTested = paste(levels(Wall_data$parasitism), collapse = " vs "),
  RelDiff = Wall_full$coefficients$fixed[2] / Wall_full$coefficients$fixed[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

sjPlot::tab_model(Wall_null,Wall_full,show.aic = T,show.reflvl = T,
                  title = "Vessel Wall Thickness - Parasites vs Hosts",
                  pred.labels = c("Parasitism","Host"),
                  dv.labels = c("Null Model","Full Model"),
                  p.style = "numeric_stars",
                  file = here("outputs","tables","Models","Tvw_model.doc")
)


Wall_data[510,"WallThickness"] <-
  rev(sort(Wall_data$WallThickness[Wall_data$indiv=="Tipuana tipu 2"]))[2]


# Iterate through species pairs
for (pair in species_pairs) {
  cat("\n", strrep("-", 50), "\n")
  cat("Processing pair:", paste(pair, collapse = " vs "), "\n")
  
  # Subset data with validation
  subset_data <- subset(Wall_data, ssp %in% pair)
  
  tryCatch({
    # --------------------------------------------------------------------------
    # Model Fitting
    # --------------------------------------------------------------------------
    full_model <- lme(
      WallThickness ~ ssp,
      random = ~ 1 | indiv,
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      WallThickness ~ 1,
      random = ~ 1 | indiv,
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    
    lrt <- anova(reduced_model, full_model)
    full_model$coefficients$fixed[2]
    VWall_AIC <- rbind(VWall_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      RelDiff = full_model$coefficients$fixed[2] / full_model$coefficients$fixed[1],
      DeltaAIC = diff(lrt$AIC),
      p_value = na.omit(lrt$`p-value`),
      stringsAsFactors = FALSE
    ))
    
    # Inside the loop's tryCatch() after model fitting:
    
    # Generate filename and path
    file_name <- paste0(
      "Tvw_",
      paste(str_extract(pair, "^\\S+"), collapse = "_"), 
      ".doc"
    )
    file_path <- here(output_dir, file_name)
    
    # FORCE RENDER AND SAVE
    suppressWarnings({
      tab <- sjPlot::tab_model(
        reduced_model,
        full_model,
        title = paste("Vessel Wall Thickness -", paste(pair, collapse = " vs ")),
        file = file_path,
        show.aic = TRUE,
        show.reflvl = TRUE,
        pred.labels = c("Intercept", "Species Effect"),
        dv.labels = c("Null Model", "Full Model"),
        p.style = "numeric_stars"
      )
      
      print(tab)
    })
    
    # Double-check file existence
    Sys.sleep(0.5) # Short pause for file I/O
    
    
  }, error = function(e) {
    cat("\n❌ ERROR:", conditionMessage(e), "\n")
  })
}
#####################################################################################
########################################################
# Initialize results dataframe
VDiameter_AIC <- data.frame(
  PairTested = character(),
  RelDiff = numeric(),
  DeltaAIC = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Fit global models
VesselDiameter_full <- lme(
  VesselDiameter ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,        # Random effects
  data = VesselDiameter_data,        # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

VesselDiameter_null <- lme(
  VesselDiameter ~ 1,                # Fixed effects
  random = ~ 1 | ssp / indiv,        # Random effects
  data = VesselDiameter_data,        # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

lrt <- anova(VesselDiameter_null, VesselDiameter_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

VDiameter_AIC <- rbind(VDiameter_AIC, data.frame(
  PairTested = paste(levels(VesselDiameter_data$parasitism), collapse = " vs "),
  RelDiff = fixed.effects(VesselDiameter_full)[2] / fixed.effects(VesselDiameter_full)[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))
residplot(VesselDiameter_full,newwd = F)
title(sub = "VDiameter PxH")
 check_model(VesselDiameter_full,show_dots = F)

sjPlot::tab_model(VesselDiameter_null, VesselDiameter_full, show.aic = TRUE, show.reflvl = TRUE,
                  title = "Vessel Diameter - Parasites vs Hosts",
                  pred.labels = c("Parasitism", "Host"),
                  dv.labels = c("Null Model", "Full Model"),
                  p.style = "numeric_stars",
                  file = here("outputs","tables","Models","Dv_model.doc"))

# Iterate through species pairs
for (pair in species_pairs) {
  cat("\n", strrep("-", 50), "\n")
  cat("Processing pair:", paste(pair, collapse = " vs "), "\n")
  
  subset_data <- subset(VesselDiameter_data, ssp %in% pair)
  
  tryCatch({
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
    
    
    # residplot(full_model,newwd = F)
    # title(sub = paste0(pair,collapse = " x "))
    # print(check_model(full_model,show_dots = F))
    # print(CookD(full_model, idn = 20,newwd = F))
    # abline(h=4/nrow(subset_data),col="red")
    # title(sub=paste0(pair,collacpse=" x "))

    lrt <- anova(reduced_model, full_model)
    
    VDiameter_AIC <- rbind(VDiameter_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      RelDiff = full_model$coefficients$fixed[2] / full_model$coefficients$fixed[1],
      DeltaAIC = diff(lrt$AIC),
      p_value = na.omit(lrt$`p-value`),
      stringsAsFactors = FALSE
    ))
    
    file_name <- paste0("VDiameter_", paste(str_extract(pair, "^\\S+"), collapse = "_"), ".doc")
    file_path <- here(output_dir, file_name)
    
    suppressWarnings({
      tab <- sjPlot::tab_model(
        reduced_model, full_model,
        title = paste("Vessel Diameter -", paste(pair, collapse = " vs ")),
        file = file_path,
        show.aic = TRUE,
        show.reflvl = TRUE,
        pred.labels = c("Intercept", "Species Effect"),
        dv.labels = c("Null Model", "Full Model"),
        p.style = "numeric_stars"
      )
      print(tab)
    })
    
    Sys.sleep(0.5) # Short pause for file I/O
    
  }, error = function(e) {
    cat("\n❌ ERROR:", conditionMessage(e), "\n")
  })
}


################################################################################
###############################################################################
# Initialize results dataframe
TopVDiameter_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)


TopVesselDiameter_full <- glmmTMB(
  VesselDiameter ~ parasitism +( 1 | ssp / indiv),        # Random effects
  data = VesselDiameter_data %>% group_by(ssp,indiv) %>%
    filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
    ungroup(),        # Data frame
  family = lognormal(link = "log"),  # Log-normal family with log link
  control = glmmTMBControl(optCtrl = list(maxiter = 150))
)

TopVesselDiameter_null <- glmmTMB(
  VesselDiameter ~ ( 1 | ssp / indiv),        # Random effects
  data = VesselDiameter_data %>% group_by(ssp,indiv) %>%
    filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
    ungroup(),        # Data frame
  family = lognormal(link = "log"),  # Log-normal family with log link
  control = glmmTMBControl(optCtrl = list(maxiter = 150))
)

# Model checks
residplot(TopVesselDiameter_full, newwd = FALSE)
title(sub = "TopVDiameter PxH")
check_model(TopVesselDiameter_full, show_dots = FALSE)


# Perform likelihood ratio test (LRT)
lrt <- anova(TopVesselDiameter_null, TopVesselDiameter_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`Pr(>Chisq)`)[1]


TopVDiameter_AIC <- rbind(TopVDiameter_AIC, data.frame(
  PairTested = paste(levels(VesselDiameter_data$parasitism), collapse = " vs "),
  RelDiff =fixef(TopVesselDiameter_full)$cond[2]/fixef(TopVesselDiameter_full)$cond[1] ,
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

sjPlot::tab_model(TopVesselDiameter_null, TopVesselDiameter_full, 
                  show.aic = TRUE, show.reflvl = TRUE,
                  title = "log Top 10% Vessel Diameter - Parasites vs Hosts",
                  pred.labels = c("Parasitism", "Host"),
                  file= here(output_dir,"TopD-model.doc"),
                  dv.labels = c("Null Model", "Full Model"),
                  p.style = "numeric_stars")

# Iterate through species pairs
for (pair in species_pairs) {
  cat("\n", strrep("-", 50), "\n")
  cat("Processing pair:", paste(pair, collapse = " vs "), "\n")
  
  subset_data <- subset(VesselDiameter_data, ssp %in% pair) %>%
    group_by(ssp, indiv) %>%
    filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
    ungroup()
  
  tryCatch({
    full_model <- glmmTMB(
      VesselDiameter ~ ssp +( 1 |indiv),        # Random effects
      data = subset_data,        # Data frame
      family = lognormal(link = "log"),  # Log-normal family with log link
      control = glmmTMBControl(optCtrl = list(maxiter = 150))
    )
    
    reduced_model <- glmmTMB(
      VesselDiameter ~ ( 1 |indiv),        # Random effects
      data = subset_data,        # Data frame
      family = lognormal(link = "log"),  # Log-normal family with log link
      control = glmmTMBControl(optCtrl = list(maxiter = 150))
    )
    
    # Model checks
    residplot(full_model, newwd = FALSE)
    title(sub = paste0(pair, collapse = " x "))
    check_model(full_model, show_dots = FALSE)

    
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`Pr(>Chisq)`)[1]
    
    TopVDiameter_AIC <- rbind(TopVDiameter_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      RelDiff = fixef(full_model)$cond[2]/fixef(reduced_model)$cond[1],
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
    
    file_name <- paste0("TopVDiameter_", paste(str_extract(pair, "^\\S+"), collapse = "_"), ".doc")
    file_path <- here(output_dir, file_name)
    
    suppressWarnings({
      tab <- sjPlot::tab_model(
        reduced_model, full_model,
        title = paste("log Top Vessel Diameter -", paste(pair, collapse = " vs ")),
        file = file_path,
        show.aic = TRUE,
        show.reflvl = TRUE,
        pred.labels = c("Intercept", "Species Effect"),
        dv.labels = c("Null Model", "Full Model"),
        p.style = "numeric_stars"
      )
      print(tab)
    })
    
    Sys.sleep(0.5) # Short pause for file I/O
    
  }, error = function(e) {
    cat("\n❌ ERROR:", conditionMessage(e), "\n")
  })
}
#######################################################################
# Initialize results dataframe
HydraulicDiameter_AIC <- data.frame(
  PairTested = character(), ParasiteMean = numeric(), HostMean = numeric(),
  REVariance = numeric(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
HydraulicDiameter_full <- lme(
  HydraulicDiameter ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,              # Data frame without grouping
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

HydraulicDiameter_null <- lme(
  HydraulicDiameter ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = Hydraulic_data,              # Data frame without grouping
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Model checks
residplot(HydraulicDiameter_full, newwd = FALSE)
title(sub = "HydraulicDiameter PxH")
check_model(HydraulicDiameter_full, show_dots = FALSE)
CookD(HydraulicDiameter_full, idn = 20, newwd = FALSE)
abline(h = 4 / nrow(Hydraulic_data), col = "red")

# Perform likelihood ratio test (LRT)
lrt <- anova(HydraulicDiameter_null, HydraulicDiameter_full)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)

HydraulicDiameter_AIC <- rbind(HydraulicDiameter_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
  RelDiff = fixed.effects(HydraulicDiameter_full)[2] / fixed.effects(HydraulicDiameter_full)[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

sjPlot::tab_model(HydraulicDiameter_null, HydraulicDiameter_full, 
                  show.aic = TRUE, show.reflvl = TRUE,
                  title = "Hydraulic Diameter - Parasites vs Hosts",
                  pred.labels = c("Parasitism", "Host"),
                  file = here(output_dir, "HydraulicDiameter-model.doc"),
                  dv.labels = c("Null Model", "Full Model"),
                  p.style = "numeric_stars")

Hydraulic_data[61,"HydraulicDiameter"] <- 
rev(sort(Hydraulic_data$HydraulicDiameter[Hydraulic_data$indiv=="Populus nigra 1"]))[2]
# Iterate through species pairs
for (pair in species_pairs) {
  cat("\n", strrep("-", 50), "\n")
  cat("Processing pair:", paste(pair, collapse = " vs "), "\n")
  
  subset_data <- subset(Hydraulic_data, ssp %in% pair)  # No grouping
  
  tryCatch({
    full_model <- lme(
      HydraulicDiameter ~ ssp,
      random = ~ 1 | indiv,
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    reduced_model <- lme(
      HydraulicDiameter ~ 1,
      random = ~ 1 | indiv,
      data = subset_data,
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method = "ML"
    )
    
    # Model checks
    residplot(full_model, newwd = FALSE)
    title(sub = paste0(pair, collapse = " x "))
    check_model(full_model, show_dots = FALSE)
    CookD(full_model, idn = 20, newwd = FALSE)
    abline(h = 4 / nrow(subset_data), col = "red")
    
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)
    
    HydraulicDiameter_AIC <- rbind(HydraulicDiameter_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      RelDiff = fixed.effects(full_model)[2] / fixed.effects(full_model)[1],
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
    
    file_name <- paste0("HydraulicDiameter_", paste(str_extract(pair, "^\\S+"),
                                                    collapse = "_"), ".doc")
    file_path <- here(output_dir, file_name)
    
    suppressWarnings({
      tab <- sjPlot::tab_model(
        reduced_model, full_model,
        title = paste("Hydraulic Diameter -", paste(pair, collapse = " vs ")),
        file = file_path,
        show.aic = TRUE,
        show.reflvl = TRUE,
        pred.labels = c("Intercept", "Species Effect"),
        dv.labels = c("Null Model", "Full Model"),
        p.style = "numeric_stars"
      )
      print(tab)
    })
    
    Sys.sleep(0.5) # Short pause for file I/O
    
  }, error = function(e) {
    cat("\n❌ ERROR:", conditionMessage(e), "\n")
  })
}
############################################################################
################################################################################

# Initialize results dataframe
VDensity_AIC <- data.frame(
  PairTested = character(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
VDensity_glm <- glmmTMB(
  round(VesselDensity)  ~ parasitism + (1 | ssp / indiv),  # Fixed and random effects
  data = Hydraulic_data,                          # Data frame
  family = nbinom2(),                             # Negative binomial family for count data
  control = glmmTMBControl(optCtrl = list(iter.max = 150, eval.max = 150))
)

VDensity_glm_null <- glmmTMB(
  round(VesselDensity) ~ 1 + (1 | ssp / indiv),  # Null model with no fixed effects
  data = Hydraulic_data,                          # Data frame
  family = nbinom2(),                             # Negative binomial family for count data
  control = glmmTMBControl(optCtrl = list(iter.max = 150, eval.max = 150))
)

# Model checks
residplot(VDensity_glm, newwd = FALSE)
title(sub = "Vessel Density PxH")
check_model(VDensity_glm, show_dots = FALSE)

# Perform likelihood ratio test (LRT)
lrt <- anova(VDensity_glm_null, VDensity_glm)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`Pr(>Chisq)`)[1]

# Append global model results
VDensity_AIC <- rbind(VDensity_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
  RelDiff = fixef(VDensity_glm)$cond[2]/fixef(VDensity_glm)$cond[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

# Save model summary to a Word document
sjPlot::tab_model(VDensity_glm_null, VDensity_glm, 
                  show.aic = TRUE, show.reflvl = TRUE,
                  title = "Vessel Density - Parasites vs Hosts",
                  pred.labels = c("Parasitism", "Host"),
                  file = here(output_dir, "VesselDensity-model.doc"),
                  dv.labels = c("Null Model", "Full Model"),
                  p.style = "numeric_stars")

# Iterate through species pairs
for (pair in species_pairs) {
  cat("\n", strrep("-", 50), "\n")
  cat("Processing pair:", paste(pair, collapse = " vs "), "\n")
  
  subset_data <- subset(Hydraulic_data, ssp %in% pair)  # No grouping
  
  tryCatch({
    full_model <- glmmTMB(
      round(VesselDensity) ~ ssp + (1 | indiv),  # Fixed and random effects
      data = subset_data,                          # Data frame
      family = nbinom2(),                             # Negative binomial family for count data
      control = glmmTMBControl(optCtrl = list(iter.max = 150, eval.max = 150))
    )
    
    reduced_model <- glmmTMB(
      round(VesselDensity) ~ 1 + (1 | indiv),  # Null model with no fixed effects
      data = subset_data,                          # Data frame
      family = nbinom2(),                             # Negative binomial family for count data
      control = glmmTMBControl(optCtrl = list(iter.max = 150, eval.max = 150))
    )
    
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`Pr(>Chisq)`)[1]
    
    # Append results for species pair
    VDensity_AIC <- rbind(VDensity_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      RelDiff = fixef(full_model)$cond[2]/fixef(reduced_model)$cond[1],
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
    residplot(VDensity_glm, newwd = FALSE)
    title(sub = paste("VesselDensity",pair))
    print(check_model(VDensity_glm, show_dots = FALSE))
    res <- simulateResiduals(full_model)
    plot(res)
    file_name <- paste0("VesselDensity_", paste(str_extract(pair, "^\\S+"), collapse = "_"), ".doc")
    file_path <- here(output_dir, file_name)
    
    suppressWarnings({
      tab <- sjPlot::tab_model(
        reduced_model, full_model,
        title = paste("Vessel Density -", paste(pair, collapse = " vs ")),
        file = file_path,
        show.aic = TRUE,
        show.reflvl = TRUE,
        pred.labels = c("Intercept", "Species Effect"),
        dv.labels = c("Null Model", "Full Model"),
        p.style = "numeric_stars"
      )
      print(tab)
    })
    
    Sys.sleep(0.5) # Short pause for file I/O
    
  }, error = function(e) {
    cat("\n❌ ERROR:", conditionMessage(e), "\n")
  })
}
######################################################################################

##vessel Fraction
VFraction_AIC <- data.frame(
  PairTested = character(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)


VFraction_full <- glmmTMB(
  (VesselFraction/100) ~ parasitism + (1 | ssp / indiv),  # Fixed effect: species, Random effect: individual
  family = beta_family(link = "logit"),  # Beta regression with logit link
  data = Hydraulic_data
) 
summary(VFraction_full)
VFractio_null <- glmmTMB(
  (VesselFraction/100) ~  + (1 | ssp / indiv),  # Fixed effect: species, Random effect: individual
  family = beta_family(link = "logit"),  # Beta regression with logit link
  data = Hydraulic_data
) 

# Perform likelihood ratio test (LRT)
lrt <- anova(VFraction_full,VFractio_null)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`Pr(>Chisq)`)[1]


residplot(VFraction_full, newwd = FALSE)
title(sub = "VFraction PxH")
check_model(VFraction_full, show_dots = FALSE)

res <- simulateResiduals(VFraction_full)
plot(res) 

VFraction_AIC <- rbind(VFraction_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
  RelDiff = plogis(fixef(VFraction_full)$cond[1] + fixef(VFraction_full)$cond[2]) / 
    plogis(fixef(VFraction_full)$cond[1]),
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

sjPlot::tab_model(VFractio_null, VFraction_full, 
                  show.aic = TRUE, show.reflvl = TRUE,
                  title = "Vessel Fraction - Parasites vs Hosts",
                  pred.labels = c("Parasitism", "Host"),
                  file = here(output_dir, "VesselFraction-model.doc"),
                  dv.labels = c("Null Model", "Full Model"),
                  p.style = "numeric_stars")

for (pair in species_pairs) {
  subset_data <- subset(Hydraulic_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for the species pair
    full_model <- glmmTMB(
      (VesselFraction/100) ~ ssp + (1 | indiv),  # Fixed effect: species, Random effect: individual
      family = beta_family(link = "logit"),  # Beta regression with logit link
      data = subset_data
    ) 
    
    reduced_model <- glmmTMB(
      (VesselFraction/100) ~ (1 | indiv),  # Fixed effect: species, Random effect: individual
      family = beta_family(link = "logit"),  # Beta regression with logit link
      data = subset_data
    ) 
    # Print model summary
    cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
    print(summary(full_model))
    print(anova( reduced_model,full_model))

    res <- simulateResiduals(full_model)
    plot(res)  
    title(sub = paste0(pair,collapse = " x "))
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))
    title(sub=paste0(pair,collacpse=" x "))
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`Pr(>Chisq)`)[1]
    
    # Append results for the species pair
    VFraction_AIC <- rbind(VFraction_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      RelDiff =plogis(fixef(full_model)$cond[1] + fixef(full_model)$cond[2]) / 
        plogis(fixef(full_model)$cond[1])
      ,
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
    
    file_name <- paste0("VesselFraction_", paste(str_extract(pair, "^\\S+"), collapse = "_"), ".doc")
    file_path <- here(output_dir, file_name)
    suppressWarnings({
      tab <- sjPlot::tab_model(
        reduced_model, full_model,
        title = paste("Vessel Fraction -", paste(pair, collapse = " vs ")),
        file = file_path,
        show.aic = TRUE,
        show.reflvl = TRUE,
        pred.labels = c("Intercept", "Species Effect"),
        dv.labels = c("Null Model", "Full Model"),
        p.style = "numeric_stars"
      )
      print(tab)
    })
    
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}


#################################################################################
##################################################################################
# Initialize results dataframe

Kmax_AIC <- data.frame(
  PairTested = character(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit full and null models
Kmax_full <- glmmTMB(
  Kmax ~ parasitism + (1 | ssp/indiv),
  family = Gamma(link = "log"),
  data = Hydraulic_data,
  control = glmmTMBControl(
    optimizer = optim,
    optArgs = list(method = "BFGS")
  )
)

Kmax_null <- glmmTMB(
  Kmax ~ (1 | ssp/indiv),
  family = Gamma(link = "log"),
  data = Hydraulic_data,
  control = glmmTMBControl(
    optimizer = optim,
    optArgs = list(method = "BFGS")
  )
)

# Perform likelihood ratio test (LRT)
lrt <- anova(Kmax_full, Kmax_null)
delta_aic <- AIC(Kmax_full) - AIC(Kmax_null)
pv <- lrt$`Pr(>Chisq)`[2]

# Model checks
res <- simulateResiduals(Kmax_full)
plot(res) 

# Append results for the global model
Kmax_AIC <- rbind(Kmax_AIC, data.frame(
  PairTested = paste(levels(Hydraulic_data$parasitism), collapse = " vs "),
  RelDiff =exp(fixef(Kmax_full)$cond[2]) /exp(fixef(Kmax_full)$cond[1]),
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

# Output model summary
sjPlot::tab_model(Kmax_null, Kmax_full, 
                  show.aic = TRUE, show.reflvl = TRUE,
                  title = "Kmax - Parasites vs Hosts",
                  pred.labels = c("Parasitism", "Host"),
                  file = here(output_dir, "Kmax-model.doc"),
                  dv.labels = c("Null Model", "Full Model"),
                  p.style = "numeric_stars")

# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(Hydraulic_data, ssp %in% pair)
  
  tryCatch({
    # Fit models for species pair
    full_model <- glmmTMB(
      Kmax ~ ssp + (1 | indiv),
      family = Gamma(link = "log"),
      data = subset_data,
      control = glmmTMBControl(
        optimizer = optim,
        optArgs = list(method = "BFGS")
      )
    )
    
    reduced_model <- glmmTMB(
      Kmax ~ (1 | indiv),
      family = Gamma(link = "log"),
      data = subset_data,
      control = glmmTMBControl(
        optimizer = optim,
        optArgs = list(method = "BFGS")
      )
    )
    
    # Perform LRT for species pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- AIC(full_model) - AIC(reduced_model)
    pv <- lrt$`Pr(>Chisq)`[2]
    
    # Append results
    Kmax_AIC <- rbind(Kmax_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      RelDiff = exp(fixef(full_model)$cond[2]) /exp(fixef(reduced_model)$cond[1]),
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
    
    # Residual diagnostics
    res <- simulateResiduals(full_model)
    plot(res, sub = paste0(pair, collapse = " x "))
    
    # Output model results
    file_name <- paste0("Kmax_", paste(pair, collapse = "_"), ".doc")
    file_path <- here(output_dir, file_name)
    
    sjPlot::tab_model(
      reduced_model, full_model,
      title = paste("Kmax -", paste(pair, collapse = " vs ")),
      file = file_path,
      show.aic = TRUE,
      show.reflvl = TRUE,
      pred.labels = c("Intercept", "Species Effect"),
      dv.labels = c("Null Model", "Full Model"),
      p.style = "numeric_stars"
    )
    Sys.sleep(0.5)
  }, error = function(e) {
    cat("\nError in species pair:", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}

# Output final results
tab_model(Kmax_full)

#############################################################################
############################################################################
# Initialize results dataframe
PitDiameter_AIC <- data.frame(
  PairTested = character(), RelDiff = numeric(), DeltaAIC = numeric(),
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

PitDiameter_null <-  lme(
  PitDiameter ~ 1,       # Fixed effects
  random = ~ 1 | ssp / indiv,       # Random effects
  data = PitDiOp_data %>% na.omit(),                 # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Perform likelihood ratio test (LRT)
lrt <- anova(PitDiameter_full, PitDiameter_null)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)[1]

# Model checks
residplot(PitDiameter_full, newwd = FALSE)
title(sub = "PitDiameter PxH")
check_model(PitDiameter_full, show_dots = FALSE)
CookD(PitDiameter_full)
abline(h = 4 / nrow(PitDiOp_data), col = "red")
title(sub = "Pit DIameter PxH")
# Append results for the global model
PitDiameter_AIC <- rbind(PitDiameter_AIC, data.frame(
  PairTested = paste(levels(PitDiOp_data$parasitism), collapse = " vs "),
  RelDiff = fixed.effects(PitDiameter_full)[2] / fixed.effects(PitDiameter_full)[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

# Output using sjPlot for the global model
sjPlot::tab_model(PitDiameter_null, PitDiameter_full, 
                  show.aic = TRUE, show.reflvl = TRUE,
                  title = "Pit Diameter - Parasites vs Hosts",
                  pred.labels = c("Parasitism", "Host"),
                  file = here(output_dir, "PitDiameter-model.doc"),
                  dv.labels = c("Null Model", "Full Model"),
                  p.style = "numeric_stars")


PitDiOp_data[623,"PitDiameter"] <- 
rev(sort(PitDiOp_data$PitDiameter[PitDiOp_data$indiv=="Tapirira guianensis 1"]))[2]
PitDiOp_data[362,"PitDiameter"] <- 
  rev(sort(PitDiOp_data$PitDiameter[PitDiOp_data$indiv=="Psittacanthus robustus 2"]))[2]
# Iterate through species pairs

for (pair in species_pairs) {
  subset_data <- subset(PitDiOp_data, ssp %in% pair) %>% na.omit()
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      PitDiameter ~ ssp,       # Fixed effects
      random = ~ 1 |indiv,       # Random effects
      data = subset_data,                 # Data frame
      method = "ML",
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp)
    )
    
    reduced_model <- lme(
      PitDiameter ~ 1,       # Fixed effects
      random = ~ 1 | indiv,       # Random effects
      data = subset_data,                 # Data frame
      method = "ML",
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp)
    )
    
    # Perform likelihood ratio test (LRT) for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)[1]
    
    # Append results for the species pair
    PitDiameter_AIC <- rbind(PitDiameter_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      RelDiff = fixed.effects(full_model)[2] / fixed.effects(full_model)[1],
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
    
    # Model checks
    residplot(full_model, newwd = FALSE)
    title(sub = paste0(pair, collapse = " x "))
    check_model(full_model, show_dots = FALSE)
    CookD(full_model,newwd = F)
    abline(h = 4 / nrow(subset_data), col = "red")
    title(sub = paste0(pair, collapse = " x "))
    # Output using sjPlot for the species pair models
    file_name <- paste0("PitDiameter_", paste(pair, collapse = "_"), ".doc")
    file_path <- here(output_dir, file_name)
    suppressWarnings({
      tab <-  sjPlot::tab_model(
        reduced_model, full_model,
        title = paste("PitDiameter -", paste(pair, collapse = " vs ")),
        file = file_path,
        show.aic = TRUE,
        show.reflvl = TRUE,
        pred.labels = c("Intercept", "Species Effect"),
        dv.labels = c("Null Model", "Full Model"),
        p.style = "numeric_stars"
      )
      
      print(tab)
    })
    Sys.sleep(0.5)
    
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}

# Check the resulting data frame
head(PitDiameter_AIC)

########################################################################################################
###################################################################################

# Initialize results dataframe
PitOpening_AIC <- data.frame(
  PairTested = character(), RelDiff = numeric(), DeltaAIC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
PitOpening_full <- lme(
  PitOpening ~ parasitism,       # Fixed effects
  random = ~ 1 | ssp / indiv,    # Random effects
  data = PitDiOp_data %>% na.omit(),  # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

PitOpening_null <- lme(
  PitOpening ~ 1,                # Fixed effects
  random = ~ 1 | ssp / indiv,    # Random effects
  data = PitDiOp_data %>% na.omit(),  # Data frame
  method = "ML",
  control = list(maxIter = 150, msMaxIter = 150),
  weights = varIdent(form = ~ 1 | ssp)
)

# Perform likelihood ratio test (LRT)
lrt <- anova(PitOpening_full, PitOpening_null)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`p-value`)[1]

# Model checks
residplot(PitOpening_full, newwd = FALSE)
title(sub = "Pit Opening PxH")
check_model(PitOpening_full, show_dots = FALSE)
CookD(PitOpening_full,newwd = F)
abline(h = 4 / nrow(PitDiOp_data), col = "red")
title(sub = "Pit Opening PxH")

# Append results for the global model
PitOpening_AIC <- rbind(PitOpening_AIC, data.frame(
  PairTested = paste(levels(PitDiOp_data$parasitism), collapse = " vs "),
  RelDiff = fixed.effects(PitOpening_full)[2] /fixed.effects(PitOpening_full)[1],
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))

# Output using sjPlot for the global model
sjPlot::tab_model(PitOpening_null, PitOpening_full, 
                  show.aic = TRUE, show.reflvl = TRUE,
                  title = "Pit Opening - Parasites vs Hosts",
                  pred.labels = c("Parasitism", "Host"),
                  file = here(output_dir, "PitOpening-model.doc"),
                  dv.labels = c("Null Model", "Full Model"),
                  p.style = "numeric_stars")



PitDiOp_data[615,"PitOpening"] <- 
  rev(sort(PitDiOp_data$PitOpening[PitDiOp_data$indiv=="Tapirira guianensis 1"]))[2]
PitDiOp_data[797,"PitOpening"] <- 
  rev(sort(PitDiOp_data$PitOpening[PitDiOp_data$indiv=="Tipuana tipu 1"]))[2]
# Iterate through species pairs
for (pair in species_pairs) {
  subset_data <- subset(PitDiOp_data, ssp %in% pair) %>% na.omit()
  
  tryCatch({
    # Fit models for the species pair
    full_model <- lme(
      PitOpening ~ ssp,       # Fixed effects
      random = ~ 1 | indiv,   # Random effects
      data = subset_data,     # Data frame
      method = "ML",
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp)
    )
    
    reduced_model <- lme(
      PitOpening ~ 1,         # Fixed effects
      random = ~ 1 | indiv,   # Random effects
      data = subset_data,     # Data frame
      method = "ML",
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp)
    )
    
    # Perform likelihood ratio test (LRT) for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`p-value`)[1]
    
    # Append results for the species pair
    PitOpening_AIC <- rbind(PitOpening_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      RelDiff = fixed.effects(full_model)[2] / fixed.effects(full_model)[1],
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
    
    # Model checks
    residplot(full_model, newwd = FALSE)
    title(sub = paste0(pair, collapse = " x "))
    check_model(full_model, show_dots = FALSE)
    CookD(full_model,newwd = F)
    abline(h = 4 / nrow(subset_data), col = "red")
    title(sub = paste0(pair, collapse = " x "))
    
    # Output using sjPlot for the species pair models
    file_name <- paste0("PitOpening_", paste(pair, collapse = "_"), ".doc")
    file_path <- here(output_dir, file_name)
    suppressWarnings({
      tab <-  sjPlot::tab_model(
        reduced_model, full_model,
        title = paste("PitOpening -", paste(pair, collapse = " vs ")),
        file = file_path,
        show.aic = TRUE,
        show.reflvl = TRUE,
        pred.labels = c("Intercept", "Species Effect"),
        dv.labels = c("Null Model", "Full Model"),
        p.style = "numeric_stars"
      )
      
      print(tab)
    })
    Sys.sleep(0.5)
    
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}
#############################################################################################
#############################################################################################

PitFraction_AIC <- data.frame( RelDiff = numeric(), DeltaAIC = numeric(),
                               p_value = numeric(), stringsAsFactors = FALSE
)

# Fit global models
PitFraction_full <- glmmTMB(
  (PitFraction/100) ~ parasitism + (1 | ssp / indiv),  # Fixed effect: species, Random effect: individual
  family = beta_family(link = "logit"),  # Beta regression with logit link
  data = PitFraction_data
) 

PitFraction_null <- glmmTMB(
  (PitFraction/100) ~  (1 | ssp / indiv),  # Fixed effect: species, Random effect: individual
  family = beta_family(link = "logit"),  # Beta regression with logit link
  data = PitFraction_data
) 


lrt <- anova(PitFraction_full,PitFraction_null)
delta_aic <- diff(lrt$AIC)
pv <- na.omit(lrt$`Pr(>Chisq)`)[1]


residplot(PitFraction_full, newwd = FALSE)
title(sub = "PitFraction PxH")
check_model(PitFraction_full, show_dots = FALSE)

res <- simulateResiduals(PitFraction_full)
plot(res) 

PitFraction_AIC <- rbind(PitFraction_AIC, data.frame(
  PairTested = paste(levels(PitDiOp_data$parasitism), collapse = " vs "),
  RelDiff = plogis(fixef(PitFraction_full)$cond[1] + fixef(PitFraction_full)$cond[2]) / 
    plogis(fixef(PitFraction_full)$cond[1])
  ,
  DeltaAIC = delta_aic,
  p_value = pv,
  stringsAsFactors = FALSE
))




sjPlot::tab_model(PitFraction_null, PitFraction_full, 
                  show.aic = TRUE, show.reflvl = TRUE,
                  title = "Pit Fraction - Parasites vs Hosts",
                  pred.labels = c("Parasitism", "Host"),
                  file = here(output_dir, "PitFraction-model.doc"),
                  dv.labels = c("Null Model", "Full Model"),
                  p.style = "numeric_stars")

for (pair in species_pairs) {
  subset_data <- subset(PitFraction_data, ssp %in% pair) %>% na.omit()
  
  tryCatch({
    # Fit models for the species pair
    full_model <- glmmTMB(
      (PitFraction/100) ~ ssp + (1 | indiv),  # Fixed effect: species, Random effect: individual
      family = beta_family(link = "logit"),  # Beta regression with logit link
      data = subset_data
    ) 
    
    reduced_model <- glmmTMB(
      (PitFraction/100) ~ (1 | indiv),  # Fixed effect: species, Random effect: individual
      family = beta_family(link = "logit"),  # Beta regression with logit link
      data = subset_data
    ) 
    # Print model summary
    cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
    print(summary(full_model))
    print(anova( reduced_model,full_model))
    
    res <- simulateResiduals(full_model)
    plot(res)  
    title(sub = paste0(pair,collapse = " x "))
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))
    title(sub=paste0(pair,collacpse=" x "))
    # Perform likelihood ratio test for the pair
    lrt <- anova(reduced_model, full_model)
    delta_aic <- diff(lrt$AIC)
    pv <- na.omit(lrt$`Pr(>Chisq)`)[1]
    
    # Append results for the species pair
    PitFraction_AIC <- rbind(PitFraction_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      RelDiff =plogis(fixef(full_model)$cond[1] + fixef(full_model)$cond[2]) / 
        plogis(fixef(full_model)$cond[1])
      ,
      DeltaAIC = delta_aic,
      p_value = pv,
      stringsAsFactors = FALSE
    ))
    
    file_name <- paste0("PitFraction_", paste(str_extract(pair, "^\\S+"), collapse = "_"), ".doc")
    file_path <- here(output_dir, file_name)
    suppressWarnings({
      tab <- sjPlot::tab_model(
        reduced_model, full_model,
        title = paste("Pit Fraction -", paste(pair, collapse = " vs ")),
        file = file_path,
        show.aic = TRUE,
        show.reflvl = TRUE,
        pred.labels = c("Intercept", "Species Effect"),
        dv.labels = c("Null Model", "Full Model"),
        p.style = "numeric_stars"
      )
      print(tab)
    })
    
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}




AIC_tables <- list(VWall_AIC,
                   VDiameter_AIC,
                   TopVDiameter_AIC,
                   HydraulicDiameter_AIC,
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

