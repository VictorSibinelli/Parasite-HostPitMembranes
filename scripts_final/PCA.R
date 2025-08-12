# Load necessary libraries
library(here)
library(FactoMineR)
library(ggbiplot)
library(tidyverse)
library(ggrepel)
library(factoextra)
library(patchwork)  # For arranging multiple plots
library(viridis)  # For viridis color scale
source(here("scripts","Functions.R"))
# Load data
Median_data <- read.csv(here("data", "processed", "Median_data.csv"))
Mean_data <- read.csv(here("data", "processed", "Mean_data.csv"))
PitMembrane_data<- read.csv(here("data", "processed", "PitMembrane_data.csv"))
head(PitMembrane_data)


# Perform bootstrapping and extract medians for Tpm and pcd
proxy <- PitMembrane_data %>%
  group_by(ssp) %>% 
  summarise(
    Tpm =median(Tpm,na.rm=T),  # Bootstrap medians for Tpm
    Pcd = median(pcd,na.rm=T)
  ) %>%
  unnest(c(Tpm, Pcd))
# View results
print(proxy)

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
trait_names <- c(
  "VesselDiameter" = "D",
  "TopVesselDiameter" = "Dtop",
  "HydraulicDiameter" = "Dh",
  "VesselDensity" = "VD",
  "VesselFraction" = "Fv",
  "Kmax" = "Kmax",
  "Wthickness" = "Tvw",
  "PitOpening" = "Dpa",
  "PitDiameter" = "Dpit",
  "PitFraction" = "Fp",
  "pcd" = "Hpit",
  "Tpm" = "Tpm"
)



Median_data <- merge(Median_data,proxy, by="ssp")
# Replace full names with short names in the data
Median_data$ssp_short <- short_names[Median_data$ssp]
colnames(Median_data) <- ifelse(
  colnames(Median_data) %in% names(trait_names),
  trait_names[colnames(Median_data)],
  colnames(Median_data)  # keep unchanged if not in trait_names
)


# --- Step 2: Perform PCA ---
pca_result <- PCA(
  Median_data %>% dplyr::select(where(is.numeric)),
  scale.unit = TRUE,
  graph = FALSE
)


relevel_factors(ls())

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
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
  )
variable_plot

variable_plot2 <- fviz_pca_var(
  pca_result,axes = c(2,3),
  col.var = "contrib",                         # Color by contribution
  gradient.cols = viridis::magma(10)[3:9],          # Use viridis color scale
  repel = TRUE                                 # Avoid overlapping text
)  +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
  )
print(variable_plot2)



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
  labs(shape = "Species")

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
    shape = "Species"
  )

print(biplot)
print(biplot2)
# --- Combine Plots into a Board ---
board <- (biplot + variable_plot)/(p1+p2)
print(board)

# Ensure consistent aspect ratio for all plots
biplot <- biplot + theme(aspect.ratio = 1)
variable_plot <- variable_plot + theme(aspect.ratio = 1)

# Arrange plots with specific layout
bi_board <- (biplot+variable_plot)+
  plot_layout(guides = "collect")

# Display the aligned board
print(bi_board)


# Ensure consistent aspect ratio for all plots
biplot2 <- biplot2 + theme(aspect.ratio = 1)
variable_plot2 <- variable_plot2 + theme(aspect.ratio = 1)


# Arrange plots with specific layout
board2 <- (biplot2 + variable_plot2) + plot_layout(guides = "collect")


print(board2)
ggsave(filename = here("outputs","figs","PCA_biplot.png"),bi_board,width = 12,height = 10,units = "in",dpi = 600)
ggsave(filename = here("outputs","figs","PCA_plot2.png"),board2,width =12,height = 10,units = "in",dpi = 600)
ggsave(filename = here("outputs","figs","PCA_scree_plot.png"),scree_plot,dpi = 600)
write.csv(loadings,file = here("outputs","tables","PCA_loadings.csv"))
