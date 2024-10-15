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
WT_obs <- tapply(WT_EV$wthickness, WT_EV$parasitism, mean)
WT_obsdiff <- WT_obs[1] - WT_obs[2]  # Observed difference



# --------------------------------------------------------
# Bootstrap test for resampling with replacement
# --------------------------------------------------------

 #bootstrap_result <- replicate(iterations,
                               shuffle_means(x = WT_EV, cols = "wthickness", 
                                             cat = "parasitism", rcol = TRUE, rcat = TRUE),
                               simplify = FALSE)
 
 bootstrap_result <- do.call(rbind, lapply(bootstrap_result, function(x) as.numeric(x)))
 
 colnames(bootstrap_result) <- levels(WT_EV$parasitism)
 # Store observed values in the first row
 bootstrap_result[1, ] <- WT_obs  
fwrite(bootstrap_result, file = here("data", "processed", "ressampled", "WallThickness_bootstrap.csv"))

# Load bootstrap results and calculate differences
bootstrap_result <- read.csv(file = here("data", "processed", "ressampled", "WallThickness_bootstrap.csv"))
boot_diff <- bootstrap_result[, 1] - bootstrap_result[, 2]  # Difference in bootstrap results

# Plot bootstrap density and mark observed difference and CI
plot(density(boot_diff), main = "Bootstrap Density Plot")
abline(v = WT_obsdiff, col = "red", lwd = 2)  # Observed difference line
abline(v = quantile(boot_diff, c(0.025, 0.975)), col = "black", lwd = 2)  # 95% CI lines
bootp <- sum(boot_diff >= abs(WT_obsdiff)) / length(boot_diff)  # Calculate p-value
text(x = 0.5, y = 0.5, paste("p-value =", bootp))


para_boot <- replicate(n = 5, {
  WT_EV %>%
    filter(parasitism == "Parasite") %>% 
    pull(wthickness) %>%  # Extract the wthickness column
    sample(replace = TRUE) %>% mean()  # Sample with replacement and calculate mean
})
host_boot <- replicate(n = 5, {
  WT_EV %>%
    filter(parasitism == "Host") %>% 
    pull(wthickness) %>%  # Extract the wthickness column
    sample(replace = TRUE) %>% mean()  # Sample with replacement and calculate mean
})
host_CI <- host_boot %>% quantile(probs=c(0.025,0.975))
para_CI <- para_boot %>% quantile(probs=c(0.025,0.975))



# --------------------------------------------------------
# Pair-wise species comparisons
# # --------------------------------------------------------
 WTssp_obs <- tapply(WT_EV$wthickness, WT_EV$ssp, mean)  # Observed species means
 WTssp_diff_obs <- c(
   WTssp_obs[1] - WTssp_obs[2],  # P. robustus - V. thyrsoidea
   WTssp_obs[3] - WTssp_obs[4],  # P. perrotettii - T. guianensis
   WTssp_obs[5] - WTssp_obs[6],  # S. rhynchophyllus - T. tipu
   WTssp_obs[7] - WTssp_obs[8]   # V. album - P. nigra
 )
#
# # Initialize a list to store results for each pair
 pair_results <- list()
 iterations <- 200  # Number of iterations for resampling

 # Loop through each species pair and perform resampling
# for (pair in species_pairs) {
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
 head(diff_matrix)
 diff_matrix[1, ] <- WTssp_diff_obs  # Set observed differences in the first row

# Write the differences to CSV
fwrite(diff_matrix, file = here("data", "processed", "ressampled", "WT_ssp_Diffbootstrap.csv"))

diff_matrix <- read.csv(file = here("data", "processed", "ressampled", "WT_ssp_Diffbootstrap.csv"))
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
  replicate(iterations, mean(sample(x, replace = TRUE), na.rm = TRUE))  
})

# Convert to data frame
boot_sspWT <- as.data.frame(do.call(cbind, boot_sspWT))
colnames(boot_sspWT) <- names(boot_sspWT)

fwrite(boot_sspWT, file = here("data", "processed", "ressampled", "WT_ssp_bootstrap.csv"))

# Calculate 95% confidence intervals
ssp_95ci <- as.data.frame(t(apply(boot_sspWT, MARGIN = 2, function(x) {
  quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)  # Calculate CI
})))

CI95 <- rbind(para_CI,host_CI,ssp_95ci)
row.names(CI95)[1:2] <- c("Parasite","Host")
colname(CI95) <- c("Lower","Upper")
CI95$Estimate <- na.omit(c(WT_obs, WTssp_obs))

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
fwrite(CI95, file = here("outputs", "tables", "WT_MonteCarlo_CI95.csv"))

