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
  "P. robustus" = 16,  # Circle
  "V. thyrsoidea" = 17,  # Triangle
  "P. perrotettii" = 18,  # Diamond
  "T. guianensis" = 19,  # Circle (larger)
  "S. rhynchophyllus" = 15,  # Square
  "T. tipu" = 3,  # Plus
  "V. album" = 4,  # Cross
  "P. nigra" = 8   # Star
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
                   ellipse.prob = 0.68) +
  scale_colour_manual(
    values = c("Parasite" = alpha("red", 0.8), "Host" = alpha("grey", 0.8))  # Apply alpha blending
  ) +
  scale_fill_manual(
    values = c("Parasite" = alpha("red", 0.5), "Host" = alpha("grey", 0.5))  # Apply alpha blending
  ) +
  geom_point(aes(shape = Median_data$ssp_short,colour = Median_data$parasitism), size = 2) +  # Different point types for each species
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

# --- Combine Plots into a Board ---
board <- biplot + variable_plot

ggsave(filename = here("outputs","figs","PCA_plot.png"),board,width = 20,height = 10,units = "in",dpi = 600)
ggsave(filename = here("outputs","figs","PCA_scree_plot.png"),scree_plot,dpi = 600)
