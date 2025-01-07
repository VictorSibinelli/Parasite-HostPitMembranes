# Load necessary libraries
library(smatr)
library(ggplot2)

# Prepare data
numeric_vars <- Median_data[, sapply(Median_data, is.numeric)]
pairs <- combn(names(numeric_vars), 2, simplify = FALSE)
results <- list()

# SMA Analysis Loop
for (pair in pairs) {
  var_x <- pair[1]
  var_y <- pair[2]
  
  # Prepare data for SMA
  combined_data <- data.frame(
    x = numeric_vars[[var_x]],
    y = numeric_vars[[var_y]],
    group = Median_data$parasitism
  )
  combined_data <- na.omit(combined_data)  # Remove missing values
  
  if (length(unique(combined_data$group)) < 2) {
    cat("Skipping pair:", var_x, "vs", var_y, "- Not enough groups\n")
    next
  }
  
  # Perform SMA and slope comparison
  sma_result <- sma(y ~ x * group, data = combined_data)
  slope_test <- slope.com(combined_data$y, combined_data$x, combined_data$group)
  
  # Store results
  results[[paste(var_x, "vs", var_y)]] <- list(
    slope_test = slope_test,
    sma_result = sma_result
  )
}

# Plot Results Loop
for (comparison in names(results)) {
  pair_result <- results[[comparison]]
  sma_result <- pair_result$sma_result
  slope_test <- pair_result$slope_test
  
  group_slopes <- slope_test$bs
  slopes <- group_slopes["slope", ]
  lower_CI <- group_slopes["lower.CI.lim", ]
  upper_CI <- group_slopes["upper.CI.lim", ]
  
  r2_values <- c(
    "Host" = as.numeric(sma_result$r2[1]),
    "Parasite" = as.numeric(sma_result$r2[2])
  )
  p_values <- signif(slope_test$p, 3)
  combined_data <- sma_result$data
  
  # Generate ggplot
  plot <- ggplot(combined_data, aes(x = x, y = y)) +
    geom_point(data = combined_data[combined_data$group == "Host", ], aes(color = "Host"), size = 3) +
    geom_point(data = combined_data[combined_data$group == "Parasite", ], aes(color = "Parasite"), size = 3) +
    geom_smooth(data = combined_data[combined_data$group == "Host", ], aes(color = "Host"),
                method = "lm", se = TRUE, formula = y ~ x, linetype = "solid", size = 1) +
    geom_smooth(data = combined_data[combined_data$group == "Parasite", ], aes(color = "Parasite"),
                method = "lm", se = TRUE, formula = y ~ x, linetype = "solid", size = 1) +
    scale_color_manual(values = c("Host" = "black", "Parasite" = "red"), name = "Group") +
    labs(
      title = paste("SMA:", comparison),
      subtitle = paste(
        "Host Slope: ", round(slopes["Host"], 2), 
        " (95% CI: [", round(lower_CI["Host"], 2), ", ", round(upper_CI["Host"], 2), "]), R² = ", round(r2_values["Host"], 2),
        "\nParasite Slope: ", round(slopes["Parasite"], 2), 
        " (95% CI: [", round(lower_CI["Parasite"], 2), ", ", round(upper_CI["Parasite"], 2), "]), R² = ", round(r2_values["Parasite"], 2),
        "\np-value (SMA test): ", p_values
      ),
      x = strsplit(comparison, " vs ")[[1]][1],
      y = strsplit(comparison, " vs ")[[1]][2]
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      legend.position = "bottom"
    )
  
  # Display plot
  print(plot)
}

