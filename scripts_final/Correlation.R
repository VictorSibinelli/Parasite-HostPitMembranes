####covariance
library(viridis)
library(GGally)
library(ggstatsplot)
library(devtools)
library(ggbiplot)
library(here)
library(patchwork)
library(car)
library(tidyverse)
library(MASS)  # for rlm function

Median_data<- read.csv(here("data", "processed", "Median_data.csv"))
Median_data$Kmax <- log(Median_data$Kmax)

head(Median_data)
var_names <- c("Dpo","Dpit","Tvw","D","Dtop","Fpit","Dh","VD","Fv","Kmax")
colnames(Median_data)[3:11] <- var_names
# Create the ggpairs plot

# Scatter plot with separate regression lines for each group
grouped_ggscatterstats(
  data = Median_data,
  x = "Dtop",
  y = "VD",
  xlab = "Top Diameter",
  ylab = "Vessel Diameter",
  type = "robust",
  grouping.var = parasitism,
  show.ci = TRUE,         # Show confidence intervals for regression lines
  ggtheme = ggplot2::theme_minimal()
)


spearman_plot <- ggpairs(
  Median_data,
  columns = 3:11,
  aes(colour = parasitism, fill = parasitism),  # Map both colour and fill to parasitism
  lower = list(
    continuous = wrap("smooth", method = "rlm")  # Use robust regression in lower panels
  ),
  upper = list(continuous = wrap("cor", method = "spearman",size=6))  # Correlation in upper panels
) +
  scale_colour_manual(
    values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("grey", 0.8))  # Customize line and point colors
  ) +
  scale_fill_manual(
    values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("grey", 0.8))  # Customize fill colors
  ) +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 20),   # Axis labels
    strip.text = element_text(size = 20),  # Facet strip text
    legend.text = element_text(size = 20)  # Legend text
  )  # Apply a minimal theme

print(spearman_plot)


# # Create the correlation matrix for "Host"
# host_corr <- ggcorrmat(
#   data = subset(Median_data, parasitism == "Host"),
#   type = "robust",
#   label = TRUE,
#   p.adjust.method = "BH"
# ) +
#   scale_fill_viridis_c(option = "viridis", alpha = 0.8) +
#   labs(title = "Host")
# 
# # Create the correlation matrix for "Parasite"
# parasite_corr <- ggcorrmat(
#   data = subset(Median_data, parasitism == "Parasite"),
#   type = "np",
#   method = "spearman",
#   label = TRUE
# ) +
#   scale_fill_viridis_c(option = "viridis", alpha = 0.8) +
#   labs(title = "Parasite") +
#   theme(legend.position = "none") +  # Remove the legend for the second plot
#   guides(fill = "none")             # Remove the color scale for the second plot
# 
# # Combine the two plots
# combined_plot <- host_corr + parasite_corr + 
#   plot_annotation(title = "spearman's Correlation Matrix")&
#   theme(plot.title = element_text(size = 20, hjust = 0.5))  
# 
# # Display the combined plot
# print(combined_plot)


trait_pairs <- combn(var_names, 2, simplify = TRUE)

# Create the slope_differences data frame with p_value column
slope_differences <- data.frame(trait_pair = character(0), p_value = numeric(0))

# Loop through each trait pair
for (pair in 1:ncol(trait_pairs)) {
  trait1 <- trait_pairs[1, pair]
  trait2 <- trait_pairs[2, pair]
  
  # Fit robust linear model with interaction between trait and parasitism
  model <- rlm(as.formula(paste(trait1, "~", trait2, "* parasitism")), data = Median_data,maxit = 100)
  
  # Check the model summary to identify the correct interaction term name
  model_summary <- summary(model)
  
  # Extract p-value for the interaction term (it should be something like trait2:parasitismHost or trait2:parasitismParasite)
  interaction_term <- grep(paste(trait2, ":parasitism", sep = ""), rownames(model_summary$coefficients), value = TRUE)
  anova_results <- Anova(model, type = 3) 
  if (length(interaction_term) > 0) {
    interaction_p_value <- anova_results$`Pr(>F)`[4]
    # Add the result to the slope_differences dataframe
    slope_differences <- rbind(slope_differences, data.frame(trait_pair = paste(trait1, trait2, sep = "_"), p_value = interaction_p_value))
  }
}
sig_slopes <- slope_differences %>% filter(p_value<0.05)

# Loop through each significant trait pair and create a regression plot
regression_plots <- list()

for (pair in sig_slopes$trait_pair) {
  # Split the trait_pair into trait1 and trait2
  traits <- strsplit(pair, "_")[[1]]
  trait1 <- traits[1]
  trait2 <- traits[2]
  
  # Fit the robust linear model
  model <- rlm(as.formula(paste(trait1, "~", trait2, "* parasitism")),
               data = Median_data, maxit = 100)
  
  # Perform Type III ANOVA
  anova_results <- Anova(model, type = 3)
  model$coefficients
  
  # Extract the p-value for the interaction term
  p_value <- anova_results$`Pr(>F)`[4]
  print(paste("P-value for", trait1, "vs", trait2, "interaction:", p_value))
  
  parasite_slope <- coef(model)[trait2] + coef(model)[4]
  host_slope <- coef(model)[trait2]
  
  # Generate the plot
  plot <- Median_data %>%
    ggplot(aes(x = .data[[trait2]], y = .data[[trait1]])) +  # Map column names dynamically
    geom_point(aes(color = parasitism),size=3) +
    stat_smooth(method = function(formula, data, weights = weight) 
      rlm(formula, data, weights = weight, method = "MM"),
      aes(color = parasitism),size=3) +
    scale_colour_manual(values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("black", 0.7))) +
    scale_fill_manual(values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("black", 0.7))) +
    annotate("text", 
             x = max(Median_data[[trait2]])*0.75,  # Move annotation outside the plot area
             y = max(Median_data[[trait1]])*1.35 , 
             label = paste0("Parasite slope: ", round(parasite_slope, 3)),
             hjust = 0, vjust = 1, color = "red",size=6) +
    annotate("text", 
             x = max(Median_data[[trait2]]) * 0.75,  # Move annotation outside the plot area
             y = max(Median_data[[trait1]]) * 1.3, 
             label = paste0("Host slope: ", round(host_slope, 3)),
             hjust = 0, vjust = 1, color = "black",size=6) +
    annotate("text", 
             x = max(Median_data[[trait2]]) * 0.75,  # Move annotation outside the plot area
             y = max(Median_data[[trait1]]) * 1.25, 
             label = paste0("p-value: ", round(p_value, 4)),
             hjust = 0, vjust = 1, color = "black",size=6) +
    theme_minimal() +
    theme(legend.position = "none",
            axis.title = element_text(size = 20),  # Axis titles
            axis.text = element_text(size = 17),
           )+
    xlim(c(0, max(Median_data[[trait2]]) * 1.2))+
    ylim(c(0,max(Median_data[[trait1]])*1.35))# Adjust margin to fit annotations outside
  
  # Print the plot
  print(plot)
  regression_plots[[pair]] <- plot
}

# Combine all plots and save
regression_plots <- patchwork::wrap_plots(regression_plots)

ggsave(filename = here("outputs", "figs", "trait_regression_plot.png"),regression_plots,
       dpi = 600, units = "in", height = 30, width = 30)
ggsave(filename = here("outputs", "figs", "spearman_plot.png"),spearman_plot,
       dpi = 600, units = "in", height = 15, width = 18)
