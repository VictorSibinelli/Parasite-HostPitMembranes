######################################################################
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.1.2 - Data Analysis Resampling - Vessel walls
######################################################################

# Clear the environment
rm(list = ls())  # Remove all objects from memory

# Load required libraries and custom functions
library(here)
source(here("scripts", "00-library.R"))  # Load packages
source(here("scripts", "Functions.R"))    # Load custom functions

# Load data
wdata <- read.csv(here("data", "processed", "wdata.csv"))


# --------------------------------------------------------
# Calculate mean wall thickness per label and species
# --------------------------------------------------------
WT_EV <- wdata %>%
  group_by(label, ssp,parasitism) %>%
  summarise(
    wthickness = median(wthickness, na.rm = TRUE),
    .groups = "drop"  # Drop grouping after summarising
  )


# --------------------------------------------------------
# Factor level reordering across all data frames
relevel_factors(ls())

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
iterations <- 1000000  # Number of iterations for resampling

# Calculate observed difference in mean wall thickness by parasitism
WT_obs <- tapply(WT_EV$wthickness, WT_EV$parasitism, median)
WT_obsdiff <- WT_obs[1] - WT_obs[2]  # Observed difference

# Perform resampling to create null distribution --> If it is the first time runing this script,
#uncoment lines 54-59
# WT_resample <- t(replicate(iterations,
#                             shuffle_means(x = WT_EV, cols = "wthickness", cat = "parasitism")))
# 
# 
# WT_resample[1, ] <- WT_obs  # First row: observed values
# fwrite(WT_resample, file = here("data", "processed", "ressampled", "WallThickness_ressampled.csv")) 

WT_resample <- as.matrix(read.csv(file = here(
  "data", "processed", "ressampled", "WallThickness_ressampled.csv")))

# Calculate differences between groups
WT_diff <- WT_resample[, 1] - WT_resample[, 2]  

# --------------------------------------------------------
# Calculate p-value from resampling distribution
# --------------------------------------------------------
WT_pvalue <- sum(abs(WT_diff) >= abs(WT_diff[1])) / length(WT_diff)

# Plot resampling density and mark the observed difference
plot(density(WT_diff), main = "Resampling Density Plot")
abline(v = WT_diff[1], col = "red")  # Observed difference line
abline(v = quantile(WT_diff, c(0.025, 0.975)), col = "black", lty = 2)  # 95% CI
text(x = 0.45, y = 0.5,cex=0.9, paste("p-value =", round(WT_pvalue,digits = 3)))


# --------------------------------------------------------
# Bootstrap test for resampling with replacement
# --------------------------------------------------------
bootstrap_result <- t(replicate(iterations,
                                shuffle_means(x = WT_EV, cols = "wthickness", 
                                              cat = "parasitism", rcol = TRUE, rcat = TRUE)))

# Store observed values in the first row
bootstrap_result[1, ] <- WT_obs  
fwrite(bootstrap_result, file = here("data", "processed", "ressampled", "WallThickness_bootstrap.csv"))

# Load bootstrap results and calculate differences
bootstrap_result <- read.csv(file = here("data", "processed", "ressampled", "WallThickness_bootstrap.csv"))
boot_diff <- bootstrap_result[, 1] - bootstrap_result[, 2]  # Difference in bootstrap results
CI_95 <- quantile(boot_diff, c(0.025, 0.975))  # 95% confidence interval

# Plot bootstrap density and mark observed difference and CI
plot(density(boot_diff), main = "Bootstrap Density Plot")
abline(v = WT_obsdiff, col = "red", lwd = 2)  # Observed difference line
abline(v = CI_95, col = "black", lwd = 2)  # 95% CI lines
bootp <- sum(boot_diff >= abs(boot_diff[1])) / length(boot_diff)  # Calculate p-value
text(x = 0.5, y = 0.5, paste("p-value =", bootp))
# --------------------------------------------------------
# Pair-wise species comparisons
# --------------------------------------------------------
WTssp_obs <- tapply(WT_EV$wthickness, WT_EV$ssp, mean)  # Observed species means
WTssp_diff_obs <- c(
  WTssp_obs[1] - WTssp_obs[2],  # P. robustus - V. thyrsoidea
  WTssp_obs[3] - WTssp_obs[4],  # P. perrotettii - T. guianensis
  WTssp_obs[5] - WTssp_obs[6],  # S. rhynchophyllus - T. tipu
  WTssp_obs[7] - WTssp_obs[8]   # V. album - P. nigra
)

# Initialize a list to store results for each pair
pair_results <- list()
iterations <- 50  # Number of iterations for resampling

# Loop through each species pair and perform resampling
for (pair in species_pairs) {
  subset_data <- subset(WT_EV, ssp %in% pair)  # Subset data for the current pair
  boot_resample <- t(replicate(iterations, 
                               shuffle_means(x = subset_data, 
                                             cols = "wthickness", 
                                             cat = "ssp", 
                                             rcol = TRUE)))
  pair_results[[paste(pair, collapse = "_vs_")]] <- boot_resample  # Store results
}

# Calculate the differences for each pair
diff_matrix <- do.call(cbind, lapply(pair_results, function(x) {
  x[, 1] - x[, 2]  # Calculate differences
}))

diff_matrix[1, ] <- WTssp_diff_obs  # Set observed differences in the first row

# Write the differences to CSV
fwrite(diff_matrix, file = here("data", "processed", "ressampled", "WT_ssp_Diffbootstrap.csv"))

# Calculate p-values for each species pair
ssp_pvalues <- apply(diff_matrix, MARGIN = 2, function(x) {
  sum(abs(x) >= abs(x[1])) / length(x)  # Calculate p-value
})

# --------------------------------------------------------
# Density plots for each species pair
# --------------------------------------------------------
for (i in 1:length(species_pairs)) {
  tryCatch({
    # Create the density plot
    plot(density(diff_matrix[, i]), 
         main = paste("Density Plot for", paste(species_pairs[[i]], collapse = " vs ")))
    
    # Add quantile lines for 95% CI
    abline(v = quantile(diff_matrix[, i], c(0.025, 0.975)), col = "black", lty = 2)  
    
    # Add a vertical line for the observed value
    abline(v = diff_matrix[1, i], col = "red")
    
    # Add text for p-value
    text(x = 0.5, y = 0.5, paste("p-value =", ssp_pvalues[i]))
    
  }, error = function(e) {
    message(paste("An error occurred for species pair", 
                  paste(species_pairs[[i]], collapse = " vs "), 
                  ":", e$message))
  })
}

# --------------------------------------------------------
# Bootstrap CIs for each species
# --------------------------------------------------------
boot_sspWT <- tapply(WT_EV$wthickness, WT_EV$ssp, function(x) {
  replicate(100, mean(sample(x, replace = TRUE), na.rm = TRUE))  # Generate 100 bootstrap means
})

# Convert to data frame
boot_sspWT_df <- as.data.frame(do.call(cbind, boot_sspWT))
colnames(boot_sspWT_df) <- names(boot_sspWT)

# Calculate 95% confidence intervals
ssp_95ci <- as.data.frame(t(apply(boot_sspWT_df, MARGIN = 2, function(x) {
  quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)  # Calculate CI
})))
colnames(ssp_95ci) <- c("lower", "upper")

# Reshape boot_sspWT from wide to long format
boot_sspWT_long <- boot_sspWT_df %>%
  pivot_longer(cols = everything(), 
               names_to = "species", 
               values_to = "wthickness") %>%
  left_join(WT_EV %>% select(ssp, parasitism) %>% distinct(), 
            by = c("species" = "ssp"))  # Join with species and parasitism info

# Prepare results for Monte Carlo
WT_MonteCarlo <- data.frame(
  Groups = c("Parasite x Hosts", "P. robustus x V. Thyrsoidea",
             "P. perrotettii x T. guianensis", "S. rhynchophyllus x T. tipu",
             "V. album x P. nigra"),
  Observed_diff = numeric(5),
  Pvalue = numeric(5),
  stringsAsFactors = FALSE
)

WT_MonteCarlo$Observed_diff <- c(WT_obsdiff, WTssp_diff_obs)  # Add observed differences
WT_MonteCarlo$Pvalue <- c(bootp, ssp_pvalues)  # Add p-values
fwrite(WT_MonteCarlo, file = here("outputs", "tables", "WT_MonteCarlo_results.csv"))
