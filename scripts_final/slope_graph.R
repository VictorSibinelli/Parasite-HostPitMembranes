library(viridis)
library(GGally)
library(ggstatsplot)
library(devtools)
library(ggbiplot)
library(here)
library(patchwork)
library(MASS)  # for robust regression (rlm)
library(tidyverse)
# Read data
Median_data <- read.csv(here("data", "processed", "Median_data.csv"))
#Median_data[17,"Kmax"] <- NA
# Rename the columns for traits
var_names <- c("Dpo", "Dpit", "Tvw", "D", "Fpit", "Dh", "VD", "Fv", "Kmax")
colnames(Median_data)[3:11] <- var_names

# Create the ggpairs plot
ggpairs(
  Median_data,
  columns = 3:11,
  aes(colour = parasitism, fill = parasitism),  # Map both colour and fill to parasitism
  lower = list(
    continuous = wrap("smooth", method = "rlm")  # Use robust regression in lower panels
  ),
  upper = list(continuous = wrap("cor", method = "spearman"))  # Use Spearman correlation in upper panels
) +
  scale_colour_manual(
    values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("grey", 0.7))  # Customize line and point colors
  ) +
  scale_fill_manual(
    values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("grey", 0.7))  # Customize fill colors
  ) +
  theme_minimal()  # Apply a minimal theme

# Follow-up test: Test if the slope of the relationship between traits differs for Parasites and Hosts
# Let's test the relationship for each pair of traits


trait_pairs <- combn(var_names, 2, simplify = TRUE)

# Create the slope_differences data frame with p_value column
slope_differences <- data.frame(trait_pair = character(0), p_value = numeric(0))

# Loop through each trait pair
for (pair in 1:ncol(trait_pairs)) {
  trait1 <- trait_pairs[1, pair]
  trait2 <- trait_pairs[2, pair]
  
  # Fit robust linear model with interaction between trait and parasitism
  model <- rlm(as.formula(paste(trait1, "~", trait2, "* parasitism")), data = Median_data)
  
  # Check the model summary to identify the correct interaction term name
  model_summary <- summary(model)
  
  # Extract p-value for the interaction term (it should be something like trait2:parasitismHost or trait2:parasitismParasite)
  interaction_term <- grep(paste(trait2, ":parasitism", sep = ""), rownames(model_summary$coefficients), value = TRUE)
  
  if (length(interaction_term) > 0) {
    interaction_p_value <- model_summary$coefficients[interaction_term, 4]
    # Add the result to the slope_differences dataframe
    slope_differences <- rbind(slope_differences, data.frame(trait_pair = paste(trait1, trait2, sep = "_"), p_value = interaction_p_value))
  }
}

library(MASS)  # For rlm function
library(car)   # For Anova function
library(ggplot2)
library(dplyr)
library(patchwork)

# Filter the slope_differences for significant p_values (less than 0.05)
significant_pairs <- slope_differences %>% filter(p_value < 0.05)

# Loop through each significant trait pair and create a regression plot
regression_plots <- list()

for (pair in significant_pairs$trait_pair) {
  # Split the trait_pair into trait1 and trait2
  traits <- strsplit(pair, "_")[[1]]
  trait1 <- traits[1]
  trait2 <- traits[2]
  
  # Fit robust linear model for plotting
  model <- rlm(as.formula(paste(trait1, "~", trait2, "* parasitism")), data = Median_data)
  
  # Use Anova from car package to get p-values and statistics
  anova_results <- Anova(model, type = 3)
  p_value <- anova_results$`Pr(>F)`[2]  # p-value for interaction term (trait2:parasitism)
  
  # Extract slope (interaction term) and calculate R²
  model_summary <- summary(model)
  interaction_term <- grep(paste(trait2, ":parasitism", sep = ""), rownames(model_summary$coefficients), value = TRUE)
  
  if (length(interaction_term) > 0) {
    slope <- model_summary$coefficients[interaction_term, 1]
    r_squared <- 1 - (model_summary$rss / sum((Median_data[[trait1]] - mean(Median_data[[trait1]]))^2))  # R² calculation
  } else {
    slope <- NA
    r_squared <- NA
  }
  
  # Create a regression plot for this pair
  plot <- ggplot(Median_data, aes_string(x = trait2, y = trait1, colour = "parasitism", fill = "parasitism")) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "rlm", se = TRUE, aes(colour = parasitism)) +  # Add robust regression line
    labs(x = trait2, y = trait1, title = paste("Regression between", trait1, "and", trait2)) +
    annotate("text", x = Inf, y = Inf, label = paste("Slope: ", round(slope, 2), "\nR²: ", round(r_squared, 2), "\np-value: ", round(p_value, 4)),
             hjust = 1.1, vjust = 1.1, size = 3, color = "black", fontface = "italic") +
    scale_colour_manual(values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("grey", 0.7))) +
    scale_fill_manual(values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("grey", 0.7))) +
    theme_minimal() +
    theme(legend.position = "none")  # Remove the legend
  
  # Store the plot in the list
  regression_plots[[pair]] <- plot
}

# Display all regression plots (you can use patchwork to arrange them)
patchwork::wrap_plots(regression_plots)
