### LME analysis




library(here)
source(here("scripts","00-library.R"))

# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDensity_data <- read.csv(here("data", "processed", "VDensity_data.csv"))
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
Wall_full <- lme(
  WallThickness ~ parasitism,      # Fixed effects
  random = ~ 1 | ssp/indiv,     # Random effects
  data = Wall_data,                  # Data frame
  method="ML",
  weights = varIdent(form = ~ 1 | ssp)
)
Wall_null <- lme(
  WallThickness ~ 1,      # Fixed effects
  random = ~ 1 | ssp/indiv,     # Random effects
  data = Wall_data,                  # Data frame
  method="ML",
  weights = varIdent(form = ~ 1 | ssp)
)
anova(Wall_full,Wall_null)
emmeans(Wall_full,~parasitism,infer=T)
residplot(Wall_full,newwd = F,id=T)
CookD(Wall_full)
title("Wall Thickness",cex=0.5,line=3)
abline(h=4/length(Wall_data$WallThickness),col="red")



for (pair in species_pairs) {
  
  subset_data <- subset(Wall_data,ssp %in% pair)
  
  if (length(unique(subset_data$ssp)) < 2) {
    cat("\nSkipping pair", paste(pair, collapse = " vs "), "due to insufficient levels in 'ssp'.\n")
    next
  }
  
  tryCatch({
    full_model <- full_model <- lme(
      Wthickness ~ ssp,      # Fixed effects
      random = ~ 1 | ssp/label,     # Random effects
      data = subset_data,                  # Data frame
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method="ML"
    )
    reduced_model <- lme(
      Wthickness ~ 1,      # Fixed effects
      random = ~ 1 | ssp/label,     # Random effects
      data = subset_data,   
      control = list(maxIter = 150, msMaxIter = 150),
      weights = varIdent(form = ~ 1 | ssp),
      method="ML"
    )
    
    cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
    print(summary(full_model))
    
    
    # 
    #     print(check_model(full_model))
    #     print(CookD(full_model,neww=F,idn=10))
    #     abline(h=4/length(wdata$Wthickness),col="red")
    #     title(main=paste(pair,collapse = " x "),cex=0.5,line=3)
    
    
    
    # Extract fixed effects and random effects variance
    fixed_effects <- round(fixed.effects(full_model),digits = 2)
    re <- as.numeric(VarCorr(full_model)[, "Variance"])
    # Extract residuals
    residuals<- round(residuals(full_model),digits=2)
    
    CoV_parasite <-round(sd(residuals[subset_data$ssp==pair[1]])/fixed_effects[1],digits = 2)
    CoV_host <- round(sd(residuals[subset_data$ssp == pair[2]])/sum(fixed_effects),digits=2)
    
    # Calculate specific metrics
    estimated_difference <- fixed_effects[2]
    variance_explained_by_RE <-(re[2]+re[4])/ 
      (re[5]+re[2]+re[4])
    
    # Perform LRT and calculate metrics
    lrt <- AIC(reduced_model,full_model)
    
    delta_aic <- diff(lrt$AIC)
    
    # Confidence intervals for fixed effects
    conf_intervals <- intervals(full_model,which = "fixed")
    
    # Append results
    VWall_AIC <- rbind(VWall_AIC, data.frame(
      PairTested = paste(unique(subset_data$ssp), collapse = " vs "),
      ParasiteMean =paste(fixed_effects[1],"(",CoV_parasite, ")",sep = ""),
      HostMean = paste(sum(fixed_effects)," (",CoV_host, ")",sep=""), 
      REVariance = variance_explained_by_RE*100 ,
      RelDiff = estimated_difference/fixed_effects[1],
      DeltaAIC = delta_aic,
      stringsAsFactors = FALSE
    ))
    
    # Confidence intervals table for the current pair
    CI95_pair <-data.frame(
      Group = unique(subset_data$ssp),
      Lower = c(conf_intervals$fixed[1,"lower"],sum(conf_intervals$fixed[,"lower"])),
      Estimate=c(conf_intervals$fixed[1,"est."],sum(conf_intervals$fixed[,"est."])),
      Upper = c(conf_intervals$fixed[1,"upper"],sum(conf_intervals$fixed[,"upper"]))
    )
    
    # Plot for the current species pair
    
    CI95 <- rbind(CI95,CI95_pair)
    
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}



#############################################################
