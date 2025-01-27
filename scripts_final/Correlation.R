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
var_names <- c("Dpo","Dpit","Tvw","D","Dtop","Fpit","Dh","VD","Fv","logKmax")
colnames(Median_data)[3:12] <- var_names

source(here("scripts", "Functions.R"))

# Create the ggpairs plot with robust regression and Spearman correlation
spearman_plot <- ggpairs(
  Median_data,
  columns = 3:12,
  aes(colour = parasitism, fill = parasitism),  # Map both colour and fill to parasitism
  lower = list(
    continuous = wrap("smooth", method = "rlm", maxit = 50)  # Robust regression
  ),
  upper = list(
    continuous = wrap("cor", method = "spearman", size = 7)  # Spearman correlation
  )
) +
  ggplot2::scale_color_manual(
    labels = c("P", "H"),  # Custom legend labels
    values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("grey", 0.8))
  ) +
  ggplot2::scale_fill_manual(
    labels = c("P", "H"),  # Custom legend labels
    values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("grey", 0.8))
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text = element_text(size = 20),    # Axis labels
    strip.text = element_text(size = 20),   # Facet strip text
    legend.text = element_text(size = 20)   # Legend text
  )

# Print the plot
print(spearman_plot)

CorMat_plot <- # Create the correlation matrix using the magma palette with lighter extremes
  ggcorrmat(
    data = Median_data %>% filter(parasitism=="Parasite"),
    type = "robust",
    label = TRUE,
    grouping.var = parasitism,
    p.adjust.method = "fdr",
    digits = 3
  ) +
  # Apply the viridis magma palette with adjusted scaling
  scale_fill_viridis_c(
    option = "magma", 
    name = "Winsorized pearson's \nCorrelation",
    limits = c(-1, 1),              
    begin = 0.25,                    
    end = 1                
  ) +
  theme(legend.position = "none")+
  # Create the correlation matrix using the magma palette with lighter extremes
  ggcorrmat(
    data = Median_data %>% filter(parasitism=="Host"),
    type = "robust",
    label = TRUE,
    grouping.var = parasitism,
    size = 0,
    p.adjust.method = "fdr",
    digits = 3
  ) +
  # Apply the viridis magma palette with adjusted scaling
  scale_fill_viridis_c(
    option = "magma", 
    name = "Winsorized pearson's \nCorrelation",
    limits = c(-1, 1),              
    begin = 0.25,                    
    end = 1                
  ) +
  theme(
    legend.position = "right",       # Adjust legend position
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
CorMat_plot
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
  model2 <- rlm(as.formula(paste(trait1, "~", trait2)),
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
  # Fit the robust regression model for the whole dataset
  
  plot <- Median_data %>%
    ggplot(aes(x = .data[[trait2]], y = .data[[trait1]])) +  # Map column names dynamically
    geom_point(aes(color = parasitism), size = 3) +
    
    # Group-specific robust regression lines
    stat_smooth(
      method = function(formula, data, weights = weight) 
        rlm(formula, data, weights = weight, method = "MM",maxit = 100),
      aes(color = parasitism), size = 3
    ) +
    
    # Add dashed line for the whole dataset regression
    geom_abline(
      intercept = coef(model2)[1], 
      slope = coef(model2)[2], 
      color = "blue", 
      linetype = "dashed", 
      size = 1
    ) +
    
    # Customize colors for groups
    scale_colour_manual(values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("black", 0.7))) +
    scale_fill_manual(values = c("Parasite" = alpha("red", 0.7), "Host" = alpha("black", 0.7))) +
    
    # Parasite slope annotation
    annotate(
      "text", 
      x = max(Median_data[[trait2]]) * 0.5,  # Position annotation
      y = max(Median_data[[trait1]]) * 0.25, 
      label = paste0("Parasite slope: ", signif(parasite_slope,2)),
      hjust = 0, vjust = 1, color = "red", size = 6
    ) +
    
    # Host slope annotation
    annotate(
      "text", 
      x = max(Median_data[[trait2]]) * 0.5,  # Position annotation
      y = max(Median_data[[trait1]]) * 0.15, 
      label = paste0("Host slope: ", signif(host_slope, 2)),
      hjust = 0, vjust = 1, color = "black", size = 6
    ) +
    
    # p-value annotation with dynamic display
    annotate(
      "text", 
      x = max(Median_data[[trait2]]) * 0.5,  # Position annotation
      y = max(Median_data[[trait1]]) * 0.05, 
      label = ifelse(p_value < 0.001, "p-value < 0.001", paste0("p-value: ", round(p_value, 3))),
      hjust = 0, vjust = 1, color = "black", size = 6
    ) +
    
    # Styling and theme adjustments
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 15),  # Axis titles
      axis.text = element_text(size = 14)
    ) +
    
    xlim(c(0, max(Median_data[[trait2]]) * 1.2)) +
    ylim(c(0, max(Median_data[[trait1]]) * 1.35))  # Adjust margins to fit annotations outside
  
  # Print the plot
  print(plot)
  
  # Print the plot
  print(plot)
  regression_plots[[pair]] <- plot
}
spur <- c(11,13,15,19)
regression_plots <- regression_plots[-spur]
# Combine all plots and save
regression_board <- lapply(regression_plots, function(p) p + theme(aspect.ratio = 1))
# Combine the plots into rows using valid list subsetting
regression_board <- patchwork::wrap_plots(regression_board, ncol = 3)


ggsave(filename = here("outputs", "figs", "trait_regression_plot.png"),regression_board,
       dpi = 600, units = "in", height = 28, width = 24)
ggsave(filename = here("outputs", "figs", "spearman_plot.png"),spearman_plot,
       dpi = 600, units = "in", height = 28, width = 28)
ggsave(filename = here("outputs","figs","CorrMat_plot.png"),CorMat_plot,dpi=600,units = "in",height = 9,width = 12)
