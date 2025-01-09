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
    
    
    residplot(full_model,newwd = F)
    title(sub = paste0(pair,collapse = " x "))
    print(check_model(full_model,show_dots = F))
    simulateResiduals(fittedModel = full_model) %>% plot() %>% print()
    title(sub = paste0(pair,collapse = " x "))
    
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

