######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.1 - Data Analysis - Vessel walls
#################################################################
library(here)
source(here("scripts","00-library.R"))

# Load data
wdata <- read.csv(here("data", "processed", "wdata.csv")) %>% drop_na()
source(here("scripts", "Functions.R"))

# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)



relevel_factors(ls())



### Initialize result data frames
VWall_AIC <- data.frame(PairTested = character(), ParasiteMean = numeric(),
                        HostMean = numeric(),
                        REVariance = numeric(), RelDiff = numeric(),
                        DeltaAIC = numeric(), stringsAsFactors = FALSE)

CI95 <- data.frame(Group = character(), Estimate = numeric(), 
                   Lower = numeric(), Upper = numeric(), stringsAsFactors = FALSE)


full_model <- lme(
  Wthickness ~ parasitism,      # Fixed effects
  random = ~ 1 | ssp/label,     # Random effects
  data = wdata,                  # Data frame
  method="ML"
)
reduced_model <- lme(
  Wthickness ~ 1,                # Fixed effect: Only the intercept
  random = ~ 1 | ssp/label,      # Random intercepts for ssp and label (nested)
  data = wdata,                  # Dataset
  method="ML"
)

 summary(full_model)
# summary(reduced_model)
# (lrt <- anova(reduced_model, full_model))
# check_model(full_model)
# CookD(full_model,neww=F,idn=10)
# abline(h=4/length(wdata$Wthickness),col="red")

# Extract fixed effects and random effects variance
fixed_effects <- round(fixed.effects(full_model),digits = 2)
re <- as.numeric(VarCorr(full_model)[, "Variance"])
# Extract residuals
residuals<- round(residuals(full_model),digits=2)

CoV_parasite <-round(sd(residuals[wdata$parasitism == "Parasite"])/fixed_effects[1],digits = 2)
CoV_host <- round(sd(residuals[wdata$parasitism == "Host"])/sum(fixed_effects),digits=2)
# Calculate specific metrics
estimated_difference <- fixed_effects[2]
variance_explained_by_RE <-(re[2]+re[4])/ 
  (re[5]+re[2]+re[4])

# Perform LRT and calculate metrics
lrt <- AIC(reduced_model,full_model)

delta_aic <- diff(lrt$AIC)

# Confidence intervals for fixed effects
conf_intervals <- intervals(full_model)



# Append results to the dataframe
VWall_AIC <- rbind(VWall_AIC, data.frame(
  PairTested = paste(levels(wdata$parasitism), collapse = " vs "),
  ParasiteMean =paste(fixed_effects[1],"(",CoV_parasite, ")",sep = ""),
  HostMean = paste(sum(fixed_effects)," (",CoV_host, ")",sep=""), 
  REVariance = variance_explained_by_RE*100 ,
  RelDiff = estimated_difference/fixed_effects[1],
  DeltaAIC = delta_aic,
  stringsAsFactors = FALSE
))

# Confidence intervals table
CI95 <- data.frame(
  Group = c("Parasite", "Host"),
  Lower = c(conf_intervals$fixed[1,"lower"],sum(conf_intervals$fixed[,"lower"])),
  Estimate=c(conf_intervals$fixed[1,"est."],sum(conf_intervals$fixed[,"est."])),
  Upper = c(conf_intervals$fixed[1,"upper"],sum(conf_intervals$fixed[,"upper"]))
)
print(VWall_AIC)
print(CI95)



wdata$Wthickness[200] <- rev(sort(wdata$Wthickness[wdata$ssp=="Struthanthus rhynchophyllus"]))[7]
 wdata$Wthickness[c(221,214,227,257,242)] <- sort(wdata$Wthickness[wdata$ssp=="Struthanthus rhynchophyllus"])[6]

 wdata$Wthickness[323] <- sort(wdata$Wthickness[wdata$ssp=="Tapirira guianensis"])[3]
 wdata$Wthickness[c(359,378,407,406,412,442,360,379,389,390,409,439)]<- rev(sort(wdata$Wthickness[wdata$ssp=="Tipuana tipu"]))[13]

 wdata$Wthickness[c(401,402,431,437,372,377)] <- sort(wdata$Wthickness[wdata$ssp=="Tipuana tipu"])[7]


# Loop through each species pair and fit the model
for (pair in species_pairs) {
  
  subset_data <- subset(wdata,ssp %in% pair)
  
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

row.names(VWall_AIC) <- NULL
print(VWall_AIC)

write.csv(VWall_AIC, file = here("outputs", "tables", "Wdata_AIC.csv"))
