# Install and load necessary packages
install.packages("glmmTMB")
install.packages("DHARMa")  # For model diagnostics
library(broom.mixed)  # For checking model fit

library(glmmTMB)
library(DHARMa)
library(performance)


# Fit the GLMM using Beta regression
model <- glmmTMB(
  (VesselFraction/100) ~ parasitism + (1 | ssp / indiv),  # Fixed effect: species, Random effect: individual
  family = beta_family(link = "logit"),  # Beta regression with logit link
  data = Hydraulic_data
) 

# Model summary
summary(model)
# Model diagnostics
res <- simulateResiduals(model)
plot(res)  # Check for residual patterns

# Model fit assessment
check_model(model)

# Likelihood ratio test (Compare full model with null model)
null_model <- glmmTMB(
  (VesselFraction/100) ~ (1 | ssp / indiv),  # Fixed effect: species, Random effect: individual
  family = beta_family(link = "logit"),  # Beta regression with logit link
  data = Hydraulic_data
) 
anova(model, null_model, test = "Chisq")

# Estimated effects and confidence intervals
confint(model)

# Predicted values
gg_data<- Hydraulic_data[,c("VesselFraction","parasitism")]
gg_data$pedicted <- predict(model, type = "response")

# Plot observed vs. predicted values
library(ggplot2)
gg_data %>% ggplot(aes(x = VesselFraction, y = predicted, color = parasitism)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Observed vs Predicted Tissue Fraction", x = "Observed", y = "Predicted")


summary_df <- broom.mixed::tidy(model, effects = "fixed")
summary_df$estimate_fraction <- exp(summary_df$estimate) / (1 + exp(summary_df$estimate))
