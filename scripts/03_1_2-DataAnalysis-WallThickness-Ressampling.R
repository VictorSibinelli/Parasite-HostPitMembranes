######################################################################
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.1.2 - Data Analysis Resampling - Vessel walls
######################################################################

# Clear environment
rm(list = ls())  # Remove all objects from memory

# Load required libraries, scripts, and data
library(here)
source(here("scripts", "00-library.R"))
source(here("scripts", "00-Functions.R"))
wdata <- read.csv(here("data", "processed", "wdata.csv"))
wdata_clean <- read.csv(here("data", "processed", "wdata_clean.csv"))

# --------------------------------------------------------
# Calculate mean wall thickness per label and species
# --------------------------------------------------------
WT_avg <- wdata %>%
  group_by(label, ssp) %>%
  summarise(
    wthickness = mean(wthickness, na.rm = TRUE),
    .groups = "drop"
  )

# Add parasitism column to both datasets
parasitic_species <- c("Psittacanthus robustus", "Phoradendron perrotettii", 
                       "Struthanthus rhynchophyllus", "Viscum album")

WT_avg <- WT_avg %>%
  mutate(parasitism = ifelse(ssp %in% parasitic_species, "Parasite", "Host"))
wdata <- wdata %>%
  mutate(parasitism = ifelse(ssp %in% parasitic_species, "Parasite", "Host"))

# --------------------------------------------------------
# Reorder factor levels for all data frames containing 'ssp'
# --------------------------------------------------------
species_levels <- c("Psittacanthus robustus", "Vochysia thyrsoidea",
                    "Phoradendron perrotettii", "Tapirira guianensis",
                    "Struthanthus rhynchophyllus", "Tipuana tipu",
                    "Viscum album", "Populus nigra")

for (df_name in ls()) {
  df <- get(df_name)
  if ("ssp" %in% colnames(df)) {
    df$ssp <- factor(df$ssp, levels = species_levels)
    assign(df_name, df)
  }
}

# Define species pairs for pairwise comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

# --------------------------------------------------------
# Permutation test setup
# --------------------------------------------------------
set.seed(42)
iterations <- 10000  # Number of iterations for resampling

# Observed mean wall thickness difference by parasitism
WT_obs <- tapply(WT_avg$wthickness, WT_avg$parasitism, mean)
WT_obsdiff <- WT_obs["Parasite"] - WT_obs["Host"]

# Resample and calculate differences across iterations
WT_resample <- t(replicate(iterations, 
                           shuffle_means(x = WT_avg, cols = "wthickness", cat = "parasitism")))
WT_resample[1, ] <- WT_obs  # Set first row to observed values

WT_diff <- WT_resample[, "Parasite"] - WT_resample[, "Host"]

# Calculate p-value from resampling distribution
WT_pvalue <- mean(abs(WT_diff) >= abs(WT_obsdiff))

# Plot density of differences and mark observed difference and 95% CI
plot(density(WT_diff), main = "Density of Differences - Permutation Test")
abline(v = WT_obsdiff, col = "red", lwd = 2)
abline(v = quantile(WT_diff, c(0.025, 0.975)), lwd = 2)
text(x = 0.5, y = 0.5, paste("p-value =", WT_pvalue))

# --------------------------------------------------------
# Bootstrap test for resampling with replacement
# --------------------------------------------------------
bootstrap_result <- t(replicate(iterations, 
                                shuffle_means(x = WT_avg, cols = "wthickness", cat = "parasitism", 
                                              rcol = TRUE, rcat = TRUE)))

bootstrap_result[1, ] <- WT_obs
boot_diff <- bootstrap_result[, "Parasite"] - bootstrap_result[, "Host"]
CI_95 <- quantile(boot_diff, c(0.025, 0.975))

# Plot bootstrap density with observed difference and CI
plot(density(boot_diff), main = "Density of Differences - Bootstrap Test")
abline(v = WT_obsdiff, col = "red", lwd = 2)
abline(v = CI_95, col = "black", lwd = 2)

# --------------------------------------------------------
# Pairwise species comparisons
# --------------------------------------------------------
WTssp_obs <- tapply(WT_avg$wthickness, WT_avg$ssp, mean)
WTssp_diff_obs <- c(WTssp_obs["Psittacanthus robustus"] - WTssp_obs["Vochysia thyrsoidea"],
                    WTssp_obs["Phoradendron perrotettii"] - WTssp_obs["Tapirira guianensis"],
                    WTssp_obs["Struthanthus rhynchophyllus"] - WTssp_obs["Tipuana tipu"],
                    WTssp_obs["Viscum album"] - WTssp_obs["Populus nigra"])

# Initialize storage for pairwise bootstrap results
pair_results <- list()

# Loop through each species pair for bootstrapping
for (pair in species_pairs) {
  subset_data <- subset(WT_avg, ssp %in% pair)
  boot_resample <- t(replicate(iterations, 
                               shuffle_means(x = subset_data, cols = "wthickness", cat = "ssp", rcol = TRUE)))
  pair_results[[paste(pair, collapse = "_vs_")]] <- boot_resample
}

# Compute difference matrix for each pair
diff_matrix <- do.call(cbind, lapply(pair_results, function(x) x[,1] - x[,2]))
diff_matrix[1, ] <- WTssp_diff_obs

# Calculate p-values for species pairs
ssp_pvalues <- apply(diff_matrix, MARGIN = 2, function(x) {
  mean(abs(x) >= abs(x[1]))
})

# Plot densities for each pairwise comparison
for (i in 1:ncol(diff_matrix)) {
  plot(density(diff_matrix[, i]), main = paste("Density Plot - Pair", i))
  abline(v = quantile(diff_matrix[, i], c(0.025, 0.975)))
  abline(v = diff_matrix[1, i], col = "red")
  text(x = 0.5, y = 0.5, paste("p-value =", ssp_pvalues[i]))
}

# --------------------------------------------------------
# Generate density plots with 95% CI for species pairs
# --------------------------------------------------------
for (pair in species_pairs) {
  subset_data <- boot_sspWT_long %>% filter(species %in% pair)
  
  p <- ggplot(subset_data, aes(x = wthickness, fill = parasitism)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Parasite" = "red", "Host" = "grey")) +
    geom_vline(data = subset(subset_data, parasitism == "Parasite"), aes(xintercept = lower), linetype = "dashed", color = "red") +
    geom_vline(data = subset(subset_data, parasitism == "Parasite"), aes(xintercept = upper), linetype = "dashed", color = "red") +
    geom_vline(data = subset(subset_data, parasitism == "Host"), aes(xintercept = lower), linetype = "dashed", color = "black") +
    geom_vline(data = subset(subset_data, parasitism == "Host"), aes(xintercept = upper), linetype = "dashed", color = "black") +
    labs(title = paste("Density Plot for", pair[1], "and", pair[2]),
         x = "Vessel Wall Thickness", y = "Density") +
    scale_x_continuous(limits = c(min(subset_data$wthickness) - 0.1, max(subset_data$wthickness) + 0.1)) +
    theme_classic()
  
  print(p)
}

