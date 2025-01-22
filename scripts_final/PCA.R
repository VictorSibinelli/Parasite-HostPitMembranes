# Load necessary libraries
library(here)
library(FactoMineR)
library(ggbiplot)
library(tidyverse)
library(GGally)
library(ggrepel)
library(factoextra)
library(patchwork)  # For arranging multiple plots
library(viridis)  # For viridis color scale

# Load data
Median_data <- read.csv(here("data", "processed", "Median_data.csv"))

# --- Step 1: Define short names for species ---
short_names <- c(
  "Psittacanthus robustus" = "P. robustus",
  "Vochysia thyrsoidea" = "V. thyrsoidea",
  "Phoradendron perrotettii" = "P. perrotettii",
  "Tapirira guianensis" = "T. guianensis",
  "Struthanthus rhynchophyllus" = "S. rhynchophyllus",
  "Tipuana tipu" = "T. tipu",
  "Viscum album" = "V. album",
  "Populus nigra" = "P. nigra"
)

# Replace full names with short names in the data
Median_data$ssp_short <- short_names[Median_data$ssp]

# Assign different point shapes for each species
shapes <- c(
  "P. robustus" = 0,  # Circle
  "V. thyrsoidea" = 15,  # Triangle
  "P. perrotettii" = 1,  # Diamond
  "T. guianensis" = 16,  # Circle (larger)
  "S. rhynchophyllus" = 2,  # Square
  "T. tipu" = 17,  # Plus
  "V. album" = 5,  # Cross
  "P. nigra" = 18   # Star
)

# --- Step 2: Perform PCA ---
pca_result <- PCA(Median_data[, 3:11], scale.unit = TRUE, graph = FALSE)

# --- Step 3: PCA Summary ---
# Extract and print eigenvalues
summary_pc <- pca_result$eig  # Eigenvalues
print(summary_pc)

# Extract and sort the loadings for PC1
loadings <- as.data.frame(pca_result$var$coord)
loadings <- loadings[order(abs(loadings[, "Dim.1"]), decreasing = TRUE), ]
print("Sorted Loadings for PC1:")
print(loadings)

# --- Step 4: Scree Plot ---
scree_plot <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 100)) +
  labs(title = "Scree Plot", x = "Principal Components", y = "% of Variance Explained") +
  theme_minimal()
print(scree_plot)

# --- Step 5: Variable Contributions Plot ---
variable_plot <- fviz_pca_var(
  pca_result,
  col.var = "contrib",                         # Color by contribution
  gradient.cols = viridis::magma(10)[3:9],          # Use viridis color scale
  repel = TRUE                                 # Avoid overlapping text
) +
  labs(title = "PCA: Variable Contributions") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )
variable_plot

variable_plot2 <- fviz_pca_var(
  pca_result,axes = c(2,3),
  col.var = "contrib",                         # Color by contribution
  gradient.cols = viridis::magma(10)[3:9],          # Use viridis color scale
  repel = TRUE                                 # Avoid overlapping text
) +
  labs(title = "PCA: Variable Contributions") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )
print(variable_plot2)


filtered_data <- Median_data %>%
  filter(parasitism == "Parasite") 

p1 <- fviz_pca_var(
  PCA(filtered_data[3:12], scale.unit = TRUE, graph = FALSE),axes = c(1,2),
  col.var = "contrib",                         # Color by contribution
  gradient.cols = viridis::magma(10)[3:9],          # Use viridis color scale
  repel = TRUE                                 # Avoid overlapping text
) +
  labs(title = "PCA: Variable Contributions") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

filtered_data <- Median_data %>%
  filter(parasitism == "Host") 
p2 <- fviz_pca_var(
  PCA(filtered_data[3:12], scale.unit = TRUE, graph = FALSE),axes = c(1,2),
  col.var = "contrib",                         # Color by contribution
  gradient.cols = viridis::magma(10)[3:9],          # Use viridis color scale
  repel = TRUE                                 # Avoid overlapping text
) +
  labs(title = "PCA: Variable Contributions") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )
p1+p2


# --- Step 6: PCA Biplot ---
# Perform PCA with prcomp for ggbiplot compatibility
pc <- prcomp(Median_data[, 3:11], center = TRUE, scale. = TRUE)

# Flip signs for PC1 and PC2 as needed
pc$rotation[, c("PC1", "PC2")] <- -pc$rotation[, c("PC1", "PC2")]
pc$x[, c("PC1", "PC2")] <- -pc$x[, c("PC1", "PC2")]

# Create the ggbiplot
biplot <- ggbiplot(pc, 
                   obs.scale = 1, 
                   var.scale = 1, 
                   groups = Median_data$parasitism, 
                   ellipse = TRUE, 
                   circle = TRUE, 
                   ellipse.prob = 0.68,
                   point.size = 0) +
  scale_colour_manual(
    values = c("Parasite" = alpha("red", 0.8), "Host" = alpha("grey", 0.8))  # Apply alpha blending
  ) +
  scale_fill_manual(
    values = c("Parasite" = alpha("red", 0.5), "Host" = alpha("grey", 0.5))  # Apply alpha blending
  ) +
  geom_point(aes(shape = Median_data$ssp_short,colour = Median_data$parasitism), size = 3) +  # Different point types for each species
  scale_shape_manual(values = shapes) +  
  xlab(paste0("Dim1 (", round(summary(pc)$importance[2, 1] * 100, 1), "% variance)")) +
  ylab(paste0("Dim2 (", round(summary(pc)$importance[2, 2] * 100, 1), "% variance)")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "PCA Biplot: Parasites vs Hosts",
    color = "Parasitism",
    shape = "Species"
  )

biplot2 <- ggbiplot(pc, 
         choices = c(2, 3),
         obs.scale = 1, 
         var.scale = 1, 
         groups = Median_data$parasitism, 
         ellipse = TRUE, 
         circle = TRUE, 
         ellipse.prob = 0.68,
         point.size = 0) +
  scale_colour_manual(
    values = c("Parasite" = alpha("red", 0.8), "Host" = alpha("grey", 0.8))
  ) +
  scale_fill_manual(
    values = c("Parasite" = alpha("red", 0.5), "Host" = alpha("grey", 0.5))
  ) +
  geom_point(aes(shape = Median_data$ssp_short, colour = Median_data$parasitism), size = 3) +  # Adjust size as needed
  scale_shape_manual(values = shapes) +  # Define 'shapes' appropriately
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  xlab(paste0("Dim2 (", round(summary(pc)$importance[2, 2] * 100, 1), "% variance)")) +
  ylab(paste0("Dim3 (", round(summary(pc)$importance[2, 3] * 100, 1), "% variance)")) +
  labs(
    title = "PCA Biplot: Parasites vs Hosts",
    color = "Parasitism",
    shape = "Species"
  )

print(biplot)
print(biplot2)
# --- Combine Plots into a Board ---
board <- biplot + variable_plot
board2 <- biplot2 + variable_plot2

ggsave(filename = here("outputs","figs","PCA_plot.png"),board,width = 20,height = 10,units = "in",dpi = 600)
ggsave(filename = here("outputs","figs","PCA_plot2.png"),board2,width = 20,height = 10,units = "in",dpi = 600)
ggsave(filename = here("outputs","figs","PCA_scree_plot.png"),scree_plot,dpi = 600)
write.csv(loadings,file = here("outputs","tables","PCA_loadings.csv"))
