######################################################################
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.1.2 - Data Analysis Resampling - Vessel walls
######################################################################

# Clear environment
rm(list = ls())  # Remove all objects from memory

library(here)
source(here("scripts", "00-library.R"))  # Load packages
source(here("scripts", "00-Functions.R"))  # Load custom functions

# Load data
wdata <- read.csv(here("data", "processed", "wdata.csv"))


# --------------------------------------------------------
# Calculate mean wall thickness per label and species
# --------------------------------------------------------
WT_avg <- wdata %>%
  group_by(label, ssp) %>%
  summarise(
    wthickness = mean(wthickness, na.rm = TRUE),
    .groups = "drop"  # Drop grouping after summarising
  )

# Add parasitism information based on species name
WT_avg <- WT_avg %>%
  mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  ))

# Apply the same parasitism classification to wdata
wdata <- wdata %>%
  mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  ))

# --------------------------------------------------------
# Factor level reordering across all data frames
# --------------------------------------------------------
dataframes <- ls()  # List all objects in the environment

for (df_name in dataframes) {
  df <- get(df_name)  # Get the object by name
  
  if ("ssp" %in% colnames(df)) {  # If the object contains the 'ssp' column
    df$ssp <- factor(df$ssp, levels = c(
      "Psittacanthus robustus", "Vochysia thyrsoidea",
      "Phoradendron perrotettii", "Tapirira guianensis",
      "Struthanthus rhynchophyllus", "Tipuana tipu",
      "Viscum album", "Populus nigra"
    ))
    assign(df_name, df)  # Save the modified data frame
  }
  rm(df, df_name)  # Remove temporary variables
}

# Define species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

# --------------------------------------------------------
# Permutation test setup
# --------------------------------------------------------
set.seed(42)  # Set seed for reproducibility
iterations <- 100000 # Number of iterations for resampling

# Calculate observed difference in mean wall thickness by parasitism
WT_obs <- tapply(WT_avg$wthickness, WT_avg$parasitism, mean)
WT_obsdiff <- WT_obs[1] - WT_obs[2]  # Observed difference

# #####Perform resampling to create null distribution
WT_resample <- read.csv(file = here("data","processed","ressampled","WallThickness_ressampled.csv"))
# WT_resample <- t(replicate(iterations,
#                            shuffle_means(x = WT_avg, cols = "wthickness", cat = "parasitism")))


WT_resample[1, ] <- WT_obs  # First row: observed values

WT_diff <- WT_resample[, 1] - WT_resample[, 2]  # Differences between groups

fwrite(WT_resample,file = here("data","processed","ressampled","WallThickness_ressampled.csv")) 
# --------------------------------------------------------
# Calculate p-value from resampling distribution
# --------------------------------------------------------
WT_pvalue <- sum(abs(WT_diff) >=  abs(WT_diff[1]))/ length(WT_diff)

# Plot resampling density and mark the observed difference

plot(density(WT_diff), main = "Resampling Density Plot")
abline(v = WT_diff[1], col = "red")
abline(v = quantile(WT_diff, c(0.025, 0.975)), col = "black", lty = 2)  # 95% CI
text(x = 0.5, y = 0.5, paste("p-value =", WT_pvalue))


# --------------------------------------------------------
# Bootstrap test for resampling with replacement
# --------------------------------------------------------

# bootstrap_result <- t(replicate(iterations,
#                                  shuffle_means(x = WT_avg, cols = "wthickness", 
#                                                cat = "parasitism", rcol = TRUE, rcat = TRUE)))
# head(bootstrap_result)


bootstrap_result[1, ] <- WT_obs  # Observed values in the first row
fwrite(bootstrap_result,file=here("data","processed","ressampled","WallThickness_bootstrap.csv"))


bootstrap_result <- read.csv(file = here("data","processed","ressampled","WallThickness_bootstrap.csv"))
boot_diff <- bootstrap_result[, 1] - bootstrap_result[, 2]  # Difference in bootstrap results
CI_95 <- quantile(boot_diff, c(0.025, 0.975))  # 95% confidence interval

# Plot bootstrap density and mark observed difference and CI

plot(density(boot_diff), main = "Bootstrap Density Plot")
abline(v = WT_obsdiff, col = "red", lwd = 2)
abline(v = CI_95, col = "black", lwd = 2)
bootp <- sum(boot_diff >= abs(boot_diff[1])) / length(boot_diff)
text(x = 0.5, y = 0.5, paste("p-value =", bootp))


# --------------------------------------------------------
# Pair-wise species comparisons
# --------------------------------------------------------
WTssp_obs <- tapply(WT_avg$wthickness, WT_avg$ssp, mean)  # Observed species means
WTssp_diff_obs <- c(
  WTssp_obs[1] - WTssp_obs[2],  # P. robustus - V. thyrsoidea
  WTssp_obs[3] - WTssp_obs[4],  # P. perrotettii - T. guianensis
  WTssp_obs[5] - WTssp_obs[6],  # S. rhynchophyllus - T. tipu
  WTssp_obs[7] - WTssp_obs[8]   # V. album - P. nigra
)

# Initialize a list to store results for each pair
pair_results <- list()
iterations <- 50000
# Loop through each species pair and perform resampling
for (pair in species_pairs) {
  subset_data <- subset(WT_avg, ssp %in% pair)
  boot_resample <- t(replicate(iterations, shuffle_means(x = subset_data, cols = "wthickness", cat = "ssp", rcol = TRUE)))
  pair_results[[paste(pair, collapse = "_vs_")]] <- boot_resample
}

# Calculate the differences for each pair
diff_matrix <- do.call(cbind, lapply(pair_results, function(x) {
  x[, 1] - x[, 2]
}))

diff_matrix[1, ] <- WTssp_diff_obs  # Set observed differences in the first row

fwrite(diff_matrix,file=here("data","processed","ressampled","WT_ssp_Diffbootstrap.csv"))

ssp_pvalues <- apply(diff_matrix, MARGIN = 2, function(x) {
  sum(abs(x) >= abs(x[1])) / length(x)
})

# --------------------------------------------------------
# Density plots for each species pair
# --------------------------------------------------------

for (i in 1:length(species_pairs)) {
  plot(density(diff_matrix[, i]), main = paste("Density Plot for", paste(species_pairs[[i]], collapse = " vs ")))
  abline(v = quantile(diff_matrix[, i], c(0.025, 0.975)), col = "black", lty = 2)  # 95% CI
  abline(v = diff_matrix[1, i], col = "red")
  text(x = 0.5, y = 0.5, paste("p-value =", ssp_pvalues[i]))
}


# --------------------------------------------------------
# Bootstrap CIs for each species
# --------------------------------------------------------
boot_sspWT <- tapply(WT_avg$wthickness, WT_avg$ssp, function(x) {
  replicate(100, mean(sample(x, replace = TRUE), na.rm = TRUE))  # 100 bootstrap means
})

boot_sspWT_df <- as.data.frame(do.call(cbind, boot_sspWT))
colnames(boot_sspWT_df) <- names(boot_sspWT)

ssp_95ci <- as.data.frame(t(apply(boot_sspWT_df, MARGIN = 2, function(x) {
  quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
})))
colnames(ssp_95ci) <- c("lower", "upper")

# Reshape boot_sspWT from wide to long format
boot_sspWT_long <- boot_sspWT_df %>%
  pivot_longer(cols = everything(), 
               names_to = "species", 
               values_to = "wthickness") %>%
  left_join(WT_avg %>% select(ssp, parasitism) %>% distinct(), 
            by = c("species" = "ssp"))  # Joining on species and ssp



WT_MonteCarlo <- data.frame(
  Groups = c("Parasite x Hosts", "P. robustus x V. Thyrsoidea",
             "P. perrotettii x T. guianensis", "S. rhynchophyllus x T. tipu",
             "V. album x P. nigra"),
  Observed_diff = numeric(5),
  Pvalue = numeric(5),
  stringsAsFactors = FALSE
)
WT_MonteCarlo$Observed_diff <- c(WT_obsdiff,WTssp_diff_obs)
WT_MonteCarlo$Pvalue <- c(bootp,ssp_pvalues)
fwrite(WT_MonteCarlo,file=here("outputs","tables","WT_MonteCarlo_results.csv"))
