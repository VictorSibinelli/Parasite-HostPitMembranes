
pca_median <- Median_data[3:11]
pc <- prcomp(pca_median, center = TRUE, scale. = TRUE)

# Print PCA results
print(pc)
summary(pc)

# Extract and sort the loadings for PC1
loadings<- pc$rotation[order(abs(pc$rotation[, "PC1"]), decreasing = TRUE), ]

# Print the sorted loadings
print(loadings)


S# Add variable names

# Create the plot
p <- ggbiplot(pc, 
              obs.scale = 1, 
              var.scale = 1, 
              groups = Median_data$parasitism, 
              ellipse = TRUE, 
              circle = TRUE, 
              ellipse.prob = 0.68) +
  scale_colour_manual(
    values = c("Parasite" = alpha("red", 0.8), "Host" = alpha("grey", 0.8))  # Apply alpha here
  ) +
  scale_fill_manual(
    values = c("Parasite" = alpha("red", 0.5), "Host" = alpha("grey", 0.5))  # Apply alpha here
  ) +
  # Adjust species name size (smaller text)
  geom_text_repel(aes(label = Median_data$ssp), size = 2, max.overlaps = 50) + # Bigger text for arrows
  xlab(paste0("PC1 (", round(summary(pc)$importance[2, 1] * 100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(summary(pc)$importance[2, 2] * 100, 1), "% variance)")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "PCA Biplot: Parasites vs Hosts",
    color = "Parasitism"
  )

print(p)

