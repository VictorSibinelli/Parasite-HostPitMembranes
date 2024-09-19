######################################################################
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.1.2 - Data Analysis Resampling - Vessel walls
######################################################################


# Clear environment
rm(list = ls()) # Remove all objects from memory

# Load pdata and packages
library(here) 
source(here("scripts", "00-library.R")) 
source(here("scripts","00-Functions.R"))
wdata <- read.csv(here("data", "processed", "wdata.csv"))
wdata_clean <- read.csv(here("data", "processed", "wdata_clean.csv"))

# --------------------------------------------------------
# Calculate mean wall thickness per label and species
# --------------------------------------------------------
WT_avg <- wdata %>%
  group_by(label, ssp) %>%
  summarise(
    wthickness = mean(wthickness, na.rm = TRUE),
    .groups = "drop" # Drop grouping after summarising
  )

# Add parasitism column to the dataset
WT_avg <- WT_avg %>%
  mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  ))
wdata <- wdata %>%
  mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  ))
# --------------------------------------------------------
# Factor level reordering across all data frames
# --------------------------------------------------------
dataframes <- ls() # Get the list of data frame names
for (df_name in dataframes) {
  df <- get(df_name) # Retrieve the data frame by its name
  
  if ("ssp" %in% colnames(df)) { # If 'ssp' column exists
    df$ssp <- factor(df$ssp, levels = c(
      "Psittacanthus robustus", "Vochysia thyrsoidea",
      "Phoradendron perrotettii", "Tapirira guianensis",
      "Struthanthus rhynchophyllus", "Tipuana tipu",
      "Viscum album", "Populus nigra"
    ))
    assign(df_name, df) # Save the modified data frame
  }
  rm(df, df_name) # Remove temporary variables
}

# --------------------------------------------------------
# Permutation test setup
# --------------------------------------------------------
iterations <- 10000 # Number of iterations for the resampling procedure

# Observed difference in mean wall thickness by parasitism
WT_obs <- tapply(WT_avg$wthickness, WT_avg$parasitism, mean)
WT_obsdiff <- WT_obs[1] - WT_obs[2] # Observed difference

# Resample and calculate differences across iterations
WT_resample <- t(replicate(iterations,
                           shuffle_means(x = WT_avg, cols = "wthickness", cat = "parasitism")))
WT_resample[1, ] <- WT_obs # Set first row to observed values

WT_diff <- WT_resample[, 1] - WT_resample[, 2] # Difference between groups

# --------------------------------------------------------
# Calculate p-value from resampling distribution
# --------------------------------------------------------
WT_pvalue <- sum(WT_diff >= WT_diff[1] | WT_diff <= (WT_diff[1] * -1)) / length(WT_diff)

# Plot density of differences and mark observed difference
plot(density(WT_diff))
abline(v = WT_diff[1], col = "red")
abline(v = WT_diff[1] * -1, col = "red")
abline(v=quantile(WT_diff,c(0.025,0.975)))
text(x = 0.5, y = 0.5, paste("p-value =", WT_pvalue))

# --------------------------------------------------------
# Bootstrap test for resampling with replacement
# --------------------------------------------------------
bootstrap_result <- t(replicate(iterations,
                                shuffle_means(x = WT_avg, cols = "wthickness", 
                                              cat = "parasitism", rcol = TRUE, rcat = TRUE)))
head(bootstrap_result)
bootstrap_result[1, ] <- WT_obs

boot_diff <- bootstrap_result[, 1] - bootstrap_result[, 2] # Difference in bootstrap
CI_95 <- quantile(boot_diff, c(0.025, 0.975)) # 95% confidence interval

# Plot bootstrap density and mark observed difference and CI
plot(density(boot_diff))
abline(v = WT_obsdiff, col = "red", lwd = 2)
abline(v = CI_95, col = "blue", lwd = 2)


# --------------------------------------------------------
# Pair-wise species comparisons
# --------------------------------------------------------

# Resample based on species
Wssp_resample <- t(replicate(iterations,
                             shuffle_means(x = WT_avg, cols = "wthickness", cat = "ssp")))
WTssp_obs <- tapply(WT_avg$wthickness, WT_avg$ssp, mean) # Observed species means
Wssp_resample[1, ] <- WTssp_obs # Set first row to observed values

# Bootstrap resampling for species pairs
boot_ssp <- t(replicate(iterations,
                        shuffle_means(x = WT_avg, cols = "wthickness", cat = "ssp", 
                                      rcol = TRUE, rcat = TRUE)))
boot_ssp[1, ] <- WTssp_obs

# --------------------------------------------------------
# Psittacanthus robustus vs Vochysia thyrsoidea
# --------------------------------------------------------
PrvsVt <- Wssp_resample[, "Vochysia thyrsoidea"]-Wssp_resample[, "Psittacanthus robustus"]

# Calculate p-value for the difference
PrVtpvalue <- sum(PrvsVt >= PrvsVt[1] | PrvsVt <= (PrvsVt[1] * -1)) / length(PrvsVt)

# Plot result
plot(density(PrvsVt))
abline(v = PrvsVt[1], col = "red")
abline(v = PrvsVt[1] * -1, col = "red")
abline(v=quantile(PrvsVt,c(0.025,0.975)))
text(x = 0.5, y = 0.5, paste("p-value =", PrVtpvalue))

# Bootstrap and CI for this pair
boot_PrvsVt <-  boot_ssp[, "Vochysia thyrsoidea"]-boot_ssp[, "Psittacanthus robustus"]
pr_CI_95 <- quantile(boot_PrvsVt, c(0.025, 0.975), na.rm = TRUE)

# Plot bootstrap results
plot(density(boot_PrvsVt, na.rm = TRUE))
abline(v = boot_PrvsVt[1], col = "red", lwd = 2)
abline(v = pr_CI_95, col = "blue", lwd = 2)

# --------------------------------------------------------
# Phoradendron perrotettii vs Tapirira guianensis
# --------------------------------------------------------
PpvsTg <- Wssp_resample[, "Tapirira guianensis"]- Wssp_resample[, "Phoradendron perrotettii"]

# Calculate p-value for the difference
PpvTgpvalue <- sum(PpvsTg >= PpvsTg[1] | PpvsTg <= (PpvsTg[1] * -1)) / length(PpvsTg)

# Plot result
plot(density(PpvsTg))
abline(v = PpvsTg[1], col = "red")
abline(v = PpvsTg[1] * -1, col = "red")
abline(v=quantile(PpvsTg,c(0.025,0.975)))
text(x = 0.5, y = 0.5, paste("p-value =", PpvTgpvalue))

# Bootstrap and CI for this pair
boot_PpvsTg <-  boot_ssp[, "Tapirira guianensis"]-boot_ssp[, "Phoradendron perrotettii"]
pp_CI_95 <- quantile(boot_PpvsTg, c(0.025, 0.975), na.rm = TRUE)

# Plot bootstrap results
plot(density(boot_PpvsTg, na.rm = TRUE))
abline(v = boot_PpvsTg[1], col = "red", lwd = 2)
abline(v = pp_CI_95, col = "blue", lwd = 2)

# --------------------------------------------------------
# Struthanthus rhynchophyllus vs Tipuana tipu
# --------------------------------------------------------
SrvsTt <- Wssp_resample[, "Tipuana tipu"]-Wssp_resample[, "Struthanthus rhynchophyllus"]

# Calculate p-value for the difference
SrTtpvalue <- sum(SrvsTt >= SrvsTt[1] | SrvsTt <= (SrvsTt[1] * -1)) / length(SrvsTt)

# Plot result
plot(density(SrvsTt))
abline(v = SrvsTt[1], col = "red")
abline(v = SrvsTt[1] * -1, col = "red")
abline(v=quantile(SrvsTt,c(0.025,0.975)))
text(x = 0.5, y = 0.5, paste("p-value =", SrTtpvalue))

# Bootstrap and CI for this pair
boot_SrvsTt <- boot_ssp[, "Tipuana tipu"]- boot_ssp[, "Struthanthus rhynchophyllus"]
sr_CI_95 <- quantile(boot_SrvsTt, c(0.025, 0.975), na.rm = TRUE)

# Plot bootstrap results
plot(density(boot_SrvsTt, na.rm = TRUE))
abline(v = boot_SrvsTt[1], col = "red", lwd = 2)
abline(v = sr_CI_95, col = "blue", lwd = 2)

# --------------------------------------------------------
# Viscum album vs Populus nigra
# --------------------------------------------------------
VavsPn <- Wssp_resample[, "Populus nigra"]-Wssp_resample[, "Viscum album"] 

# Calculate p-value for the difference
VaPnpvalue <- sum(VavsPn >= VavsPn[1] | VavsPn <= (VavsPn[1] * -1)) / length(VavsPn)

# Plot result
plot(density(VavsPn))
abline(v = VavsPn[1], col = "red")
abline(v = VavsPn[1] * -1, col = "red")
abline(v=quantile(VavsPn,c(0.025,0.975)))
text(x = 0.5, y = 0.5, paste("p-value =", VvPnpvalue))

# Bootstrap and CI for this pair
boot_VavsPn <-  boot_ssp[, "Populus nigra"]-boot_ssp[, "Viscum album"]
va_CI_95 <- quantile(boot_VavsPn, c(0.025, 0.975), na.rm = TRUE)

# Plot bootstrap results
plot(density(boot_VavsPn, na.rm = TRUE))
abline(v = boot_VavsPn[1], col = "red", lwd = 2)
abline(v = va_CI_95, col = "blue", lwd = 2)

WT_MonteCarlo <- data.frame(
  Groups = c("Parasite x Hosts", "P. robustus x V. Thyrsoidea",
             "P. perrotettii x T. guianensis", "S. rhynchophyllus x T. tipu",
             "V. album x P. nigra"),
  Observed_diff = numeric(5),
  Pvalue = numeric(5),
  BootstrapCI=numeric(5),
  stringsAsFactors = FALSE
)

ci95s <- list(round(CI_95,pr_CI_95,pp_CI_95,sr_CI_95,va_CI_95),digits=3)
ci95_combined <- sapply(
  lapply(ci95s, function(ci) round(ci, digits = 3)), function(ci) paste(ci[1], ci[2], sep = "-"))

WT_MonteCarlo$Observed_diff <- c(WT_obsdiff,PrvsVt[1],PpvsTg[1],SrvsTt[1],VavsPn[1])
WT_MonteCarlo$Pvalue <- c(WT_pvalue,PrVtpvalue,PpvTgpvalue,SrTtpvalue,VaPnpvalue)
WT_MonteCarlo$BootstrapCI <- ci95_combined
# End of script


hostw <- subset(wdata,parasitism=="Host")
host_boot <- replicate(iterations,
                                shuffle_means(x = hostw, cols = "wthickness", 
                                              cat = "parasitism", rcol = TRUE))

plot(density(host_boot,na.rm = T))
abline(v=quantile(host_boot, c(0.025, 0.975), na.rm = TRUE),col="red")
abline(v=mean(hostw$wthickness))


paraw <- subset(wdata,parasitism=="Parasite")
para_boot <- replicate(iterations,
                       shuffle_means(x = paraw, cols = "wthickness", 
                                     cat = "parasitism", rcol = TRUE))

plot(density(para_boot,na.rm = T))
abline(v=quantile(para_boot, c(0.025, 0.975), na.rm = TRUE),col="red")
abline(v=mean(paraw$wthickness))

# Convert the bootstrapped results to a long format
host_df <- data.frame(wthickness = as.vector(host_boot), group = "Host")
para_df <- data.frame(wthickness = as.vector(para_boot), group = "Parasite")

# Combine the two data frames
parahost_boot <- rbind(host_df, para_df)

# Create density plot
ggplot(data = combined_df, aes(x = wthickness, fill = group)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = quantile(para_boot, c(0.025, 0.975), na.rm = TRUE),
             color = "black", linetype = "dashed",size=0.8) +
  geom_vline(xintercept = quantile(host_boot, c(0.025, 0.975), na.rm = TRUE),
             color = "black", linetype = "dashed",size=0.8) +
  scale_fill_manual(values = c("Host" = "grey", "Parasite" = "red")) + # Set colors
  labs(title = "Density of Wall Thickness for Hosts and Parasites",
       x = "Wall Thickness",
       y = "Density") +
  theme_minimal()
