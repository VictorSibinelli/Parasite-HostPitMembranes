######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.2 - Data Analysis - Vessels
#################################################################
library(here)
source(here("scripts", "02_2-TestAssumptions-Vessels.R"))
rm(list = ls())

#load data
HydraulicData <- read.csv(here("data", "processed", "HydraulicData.csv"))
vdata <- read.csv(here("data", "processed", "vdata.csv"))
vadata <- read.csv(here("data", "processed", "vadata.csv"))



# List of species pairs
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu")  #,
  #c("Viscum album", "Populus nigra")
)


Hydraulic_means <- HydraulicData %>%
  group_by(indiv) %>%
  summarise(
    ssp = first(ssp),  # Or use any other summary function if needed
    HydraulicDiameter = mean(HydraulicDiameter, na.rm = TRUE),
    vdensity = mean(vdensity, na.rm = TRUE),
    Kmax = mean(Kmax, na.rm = TRUE),
    parasitism = first(parasitism),  # This assumes the same for each individual
    .groups = "drop"
  )



shuffle_calculation <- function(x, cols, cat) {
  # Copy the data frame
  shuffled <- x
  
  # Shuffle numeric columns (cols) using apply
  shuffled[, cols] <- apply(shuffled[, cols], 2, sample)
  
  # Shuffle the categorical column
  shuffled[, cat] <- sample(shuffled[[cat]])
  
  # Initialize a numeric vector to store the results
  results <- numeric(length(cols))
  
  # Loop through each numeric column to calculate differences in means
  for (i in seq_along(cols)) {
    # Calculate the difference in means for the current numeric column
    results[i] <- diff(tapply(shuffled[[cols[i]]], shuffled[[cat]], mean))
  }
  
  # Assign names to the results vector based on the column names
  names(results) <- cols
  return(results)
}


# 
# # Function to perform one iteration of shuffling and calculating differences
# shuffle_calculation <- function() {
#   # Shuffle numeric columns and the parasitism column
#   shuffled <- Hydraulic_means
#   shuffled[, 3:5] <- apply(shuffled[, 3:5], 2, sample)  # Shuffle numeric columns
#   shuffled$parasitism <- sample(shuffled$parasitism)   # Shuffle parasitism column
#   
#   # Calculate differences in medians
#   VD <- diff(tapply(shuffled$HydraulicDiameter, shuffled$parasitism, median))
#   Vdens <- diff(tapply(shuffled$vdensity, shuffled$parasitism, median))
#   K <- diff(tapply(shuffled$Kmax, shuffled$parasitism, median))
#   
#   return(c(VD, Vdens, K))
# }


# Function to perform one bootstrap resample and calculate differences
bootstrap_calculation <- function() {
  # Resample with replacement
  resampled_data <- Hydraulic_means[sample(nrow(Hydraulic_means), replace = TRUE), ]
  
  # Calculate differences in medians for the resampled data
  VD_bootstrap <- diff(tapply(resampled_data$HydraulicDiameter, resampled_data$parasitism, median))
  Vdens_bootstrap <- diff(tapply(resampled_data$vdensity, resampled_data$parasitism, median))
  K_bootstrap <- diff(tapply(resampled_data$Kmax, resampled_data$parasitism, median))
  
  return(c(VD_bootstrap, Vdens_bootstrap, K_bootstrap))
}

# Number of iterations for the permutation test
iterations <- 100


# Use replicate to vectorize the permutations and calculations
Hydraulic_Results <- t(replicate(iterations, shuffle_calculation()))
# Assign column names to the result matrix
colnames(Hydraulic_Results) <- colnames(Hydraulic_means)[3:5]

#Calculate observed differences
VD_obs <- diff(
  tapply(Hydraulic_means$HydraulicDiameter, Hydraulic_means$parasitism, median))

Vdens_obs <- diff(
  tapply(Hydraulic_means$vdensity, Hydraulic_means$parasitism, median))

K_obs <- diff(
  tapply(Hydraulic_means$Kmax, Hydraulic_means$parasitism, median))

#add observed diferences to the dataframe
Hydraulic_Results[1,] <- c(VD_obs,Vdens_obs,K_obs) 


#vessel diameter
plot(density(Hydraulic_Results[,1]))
abline(v = Hydraulic_Results[1,1], col = "red")
abline(v = Hydraulic_Results[1,1]*-1, col = "red")

#vessel density
plot(density(Hydraulic_Results[,2]))
abline(v = Hydraulic_Results[1,2], col = "red")
abline(v = Hydraulic_Results[1,2]*-1, col = "red")


#Kmax
plot(density(Hydraulic_Results[,3]))
abline(v = Hydraulic_Results[1,3], col = "red")
abline(v = Hydraulic_Results[1,3]*-1, col = "red")


#calculate p-values
HD_pvalue <- sum(Hydraulic_Results[,1] <= Hydraulic_Results[1,1]|
                  Hydraulic_Results[,1] >= (Hydraulic_Results[1,1]*-1))/
  length(Hydraulic_Results[,1])

Vdens_pvalue <- sum(Hydraulic_Results[,2] >= Hydraulic_Results[1,2]|
                  Hydraulic_Results[,2] <= (Hydraulic_Results[1,2]*-1))/
  length(Hydraulic_Results[,2])

Kmax_pvalue <- sum(Hydraulic_Results[,3] <= Hydraulic_Results[1,3]|
                  Hydraulic_Results[,3] >= (Hydraulic_Results[1,3]*-1))/
  length(Hydraulic_Results[,3])



# Use replicate to perform the bootstrap resampling
Bootstrap_Results <- t(replicate(iterations, bootstrap_calculation()))


# Add observed differences to the results
Bootstrap_Results <- rbind(c(VD_obs, Vdens_obs, K_obs), Bootstrap_Results)


# Calculate bootstrap confidence intervals
CI_95 <- apply(Bootstrap_Results[-1, ], 2, quantile, probs = c(0.025, 0.975)) 
colnames(CI_95) <- colnames(Hydraulic_means)[3:5]

# Calculate bootstrap p-values
HD_pvalue_bootstrap <- mean(abs(Bootstrap_Results[-1, 1]) >= abs(Bootstrap_Results[1, 1]))
Vdens_pvalue_bootstrap <- mean(abs(Bootstrap_Results[-1, 2]) >= abs(Bootstrap_Results[1, 2]))
Kmax_pvalue_bootstrap <- mean(abs(Bootstrap_Results[-1, 3]) >= abs(Bootstrap_Results[1, 3]))

# Display bootstrap p-values and confidence intervals
list(
  HD_pvalue_bootstrap = HD_pvalue_bootstrap,
  Vdens_pvalue_bootstrap = Vdens_pvalue_bootstrap,
  Kmax_pvalue_bootstrap = Kmax_pvalue_bootstrap,
  CI_95 = CI_95
)
# Assign column names to the result matrix
colnames(Bootstrap_Results) <- colnames(Hydraulic_means)[3:5]
(HD_CI <- quantile(Bootstrap_Results[,1], prob = c(0.025, 0.975)))
(VDens_CI <- quantile(Bootstrap_Results[,2], prob = c(0.025, 0.975)))
(Kmax_CI <- quantile(Bootstrap_Results[,3], prob = c(0.025, 0.975)))

# Plot Bootstrap Distributions

plot(density(Bootstrap_Results[, 1]), main = "Bootstrap - Vessel Diameter", xlab = "Difference in Medians")
abline(v = VD_obs, col = "red", lwd = 2)
abline(v=CI_95[,1],col="blue",lwd=2)

plot(density(Bootstrap_Results[, 2]), main = "Bootstrap - Vessel Density", xlab = "Difference in Medians")
abline(v = Vdens_obs, col = "red", lwd = 2)
abline(v=CI_95[,2],col="blue",lwd=2)

plot(density(Bootstrap_Results[, 3]), main = "Bootstrap - Kmax", xlab = "Difference in Medians")
abline(v = K_obs, col = "red", lwd = 2)
abline(v=CI_95[,3],col="blue",lwd=2)

fwrite(Hydraulic_Results,file=here("outputs","tables","Hydraulic_MonteCarlo_results.csv"))

