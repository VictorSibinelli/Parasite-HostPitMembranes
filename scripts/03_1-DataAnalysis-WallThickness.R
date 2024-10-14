######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.1 - Data Analysis - Vessel walls
#################################################################

library(here)
source(here("scripts", "02_1-TestAssumptions-WallThickness.R"))
rm(list = ls())

# Load data
wdata <- read.csv(here("data", "processed", "wdata.csv"))
wdata_clean <- read.csv(here("data", "processed", "wdata_clean.csv"))
source(here("scripts", "Functions.R"))

# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu")
)

relevel_factors(ls())

### Initialize result data frames
VWall_AIC <- data.frame(PairTested = character(), ParasiteMean = numeric(),
                        HostMean = numeric(), EstimatedDifference = numeric(),
                        REVariance = numeric(), RelDiff = numeric(),
                        DeltaAIC = numeric(), stringsAsFactors = FALSE)

CI95 <- data.frame(Group = character(), Estimate = numeric(), 
                   Lower = numeric(), Upper = numeric(), stringsAsFactors = FALSE)


full_model <- lmer(wthickness ~ parasitism + (1 | ssp/label), data = wdata)
reduced_model <- lmer(wthickness ~ 1+ (1 | ssp/label), data = wdata)
summary(full_model)
summary(reduced_model)
lrt <- anova(reduced_model, full_model)

# Extract information from the full model
fixed_effects <- summary(full_model)$coefficients
random_effects <- as.data.frame(VarCorr(full_model))
residual_variance <- sigma(full_model)^2

# Extract fixed effects and variances
estimated_difference <- fixed_effects[2,1] # Adjust if necessary
variance_explained_by_RE <- random_effects$vcov[1]+random_effects$vcov[2] / sum(random_effects$vcov)
lrt_p_value <- lrt$`Pr(>Chisq)`[2]
delta_aic <- diff(lrt$AIC)

# Confidence intervals for fixed effects
conf_intervals <- confint(full_model, level = 0.95)
conf_intervals <- conf_intervals[rownames(conf_intervals) %in% rownames(summary(full_model)$coefficients), ]


# Append results to the dataframe
VWall_AIC <- rbind(VWall_AIC, data.frame(
  PairTested = paste(levels(wdata$parasitism), collapse = " vs "),
  ParasiteMean = fixed_effects[1,1], # Assumes the first level is the reference
  HostMean = fixed_effects[1,1] + fixed_effects[2,1], # Adjust if necessary
  EstimatedDifference = estimated_difference,
  REVariance = variance_explained_by_RE,
  RelDiff = estimated_difference/fixed_effects[1,1],
  DeltaAIC = delta_aic,
  stringsAsFactors = FALSE
))

# Confidence intervals table
CI95 <- data.frame(
  Group = c("Parasite", "Host"),
  Estimate = c(summary(full_model)$coef[1, 1], sum(summary(full_model)$coef[, 1])),
  Lower = c(conf_intervals[1, 1], conf_intervals[1, 1] + conf_intervals[2, 1]),
  Upper = c(conf_intervals[1, 2], conf_intervals[1, 2] + conf_intervals[2, 2])
)


print(check_model(full_model))
testDispersion(full_model)
title(main=paste(levels(wdata$parasitism), collapse = " vs "),line= 1,cex.main=1)
simulationOutput <- simulateResiduals(fittedModel = full_model, plot = T)
title(main = paste(levels(wdata$parasitism), collapse = " vs "),line=1,cex.main=1)

# Loop through each species pair and fit the model
for (pair in species_pairs) {
  
  subset_data <- wdata[wdata$ssp %in% pair, ]
  
  if (length(unique(subset_data$ssp)) < 2) {
    cat("\nSkipping pair", paste(pair, collapse = " vs "), "due to insufficient levels in 'ssp'.\n")
    next
  }
  
  tryCatch({
    full_model <- lmer(wthickness ~ ssp + (1 | label), data = subset_data)
    reduced_model <- lmer(wthickness ~ (1 | label), data = subset_data)
    
    cat("\nModel Summary for Full Model (", paste(pair, collapse = " vs "), "):\n")
    print(summary(full_model))
    
    # Likelihood Ratio Test
    anova(reduced_model, full_model)
    
    # Confidence intervals
    conf_intervals <- confint(full_model, level = 0.95)
    conf_intervals <- conf_intervals[rownames(conf_intervals) %in% rownames(summary(full_model)$coefficients), ]
    
    delta_aic <- AIC(full_model) - AIC(reduced_model)
    variance_explained_by_RE <- as.data.frame(VarCorr(full_model))$vcov[1] / sum(as.data.frame(VarCorr(full_model))$vcov)
    
    # Append results
    VWall_AIC <- rbind(VWall_AIC, data.frame(
      PairTested = paste(pair, collapse = " vs "),
      ParasiteMean = summary(full_model)$coef[1, 1],
      HostMean = summary(full_model)$coef[1, 1] + summary(full_model)$coef[2, 1],
      EstimatedDifference = summary(full_model)$coef[2, 1],
      REVariance = variance_explained_by_RE,
      RelDiff = summary(full_model)$coef[2, 1] / summary(full_model)$coef[1, 1],
      DeltaAIC = delta_aic,
      stringsAsFactors = FALSE
    ))
    
    # Confidence intervals table for the current pair
    CI95_pair <- data.frame(
      Group = unique(subset_data$ssp),
      Estimate = c(summary(full_model)$coef[1, 1], sum(summary(full_model)$coef[, 1])),
      Lower = c(conf_intervals[1, 1], conf_intervals[1, 1] + conf_intervals[2, 1]),
      Upper = c(conf_intervals[1, 2], conf_intervals[1, 2] + conf_intervals[2, 2])
    )
    
    # Plot for the current species pair
   
    CI95 <- rbind(CI95,CI95_pair)
    
    print(check_model(full_model))
    testDispersion(full_model)
    title(main=paste(pair, collapse = " vs "),line= 1,cex.main=1)
    simulationOutput <- simulateResiduals(fittedModel = full_model, plot = T)
    title(main = paste(pair, collapse = " vs "),line=1,cex.main=1)
    
    
  }, error = function(e) {
    cat("\nAn error occurred for pair", paste(pair, collapse = " vs "), ":", e$message, "\n")
  })
}

# Set the desired order for the groups
desired_order <- c("Parasite", "Host",
                   "Psittacanthus robustus", "Vochysia thyrsoidea",
                   "Phoradendron perrotettii", "Tapirira guianensis",
                   "Struthanthus rhynchophyllus", "Tipuana tipu")

# Convert Group to a factor with specified levels
CI95$Group <- factor(CI95$Group, levels = desired_order)

# Create the plot
CI95 %>%
  ggplot(aes(Group, Estimate)) +
  geom_point(size = 4, aes(color = Group)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  coord_flip() +
  labs(title = "Estimates and 95% Confidence Intervals",
       x = "Effect", y = "Estimate") +
  scale_color_manual(values = rep(c("firebrick", "black"), length.out = nlevels(CI95$Group))) +  # Alternate colors
  theme_minimal()


write.csv(VWall_AIC, file = here("outputs", "tables", "Wdata_AIC.csv"))
write.csv(CI95, file = here("outputs", "tables", "WTCI95.csv"))
cat("Summary:\n
1. Model results: Significant differences were found in Struthanthus rhynchophyllus vs Tipuana tipu and Parasite vs Host.
2. Residual analysis: All models were checked. Deviations were either non-existent or small enough to be disregarded.")
