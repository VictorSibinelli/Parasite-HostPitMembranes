######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.2 - Data Analysis - Vessels
#################################################################
library(here)
source(here("scripts", "02_2-TestAssumptions-Vessels.R"))
rm(list = ls())
source(here("scripts", "00-Functions.R"))
#load data
HydraulicData <- read.csv(here("data", "processed", "HydraulicData.csv"))
vdata <- read.csv(here("data", "processed", "vdata.csv"))
vadata <- read.csv(here("data", "processed", "vadata.csv"))



# List of species pairs
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu")  ,
  c("Viscum album", "Populus nigra")
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


# Number of iterations for the permutation test
iterations <- 10
set.seed(42)
# Specify the columns for which you want to calculate bootstrap values
vars <- colnames(Hydraulic_means)[3:5]

# Initialize a list to store results for each variable
bootstrap_results <- list()

# Loop through each variable to calculate bootstrap values
for (v in vars) {
  name <- paste(v, "bootstrap values", sep = " ")
  
  # Store the results in the list using the variable name as the key
  bootstrap_results[[name]] <- t(replicate(iterations,
                                           shuffle_means(Hydraulic_means, cols = v, cat = "parasitism",rcol = T)))
}

Obs_values <-  Hydraulic_means %>% group_by(parasitism) %>% 
  summarize(HydraulicDiameter = mean(HydraulicDiameter, na.rm = TRUE),
            vdensity = mean(vdensity, na.rm = TRUE),
            Kmax = mean(Kmax, na.rm = TRUE),
            parasitism = first(parasitism),  # This assumes the same for each individual
            .groups = "drop")


HD_bootstrap <- bootstrap_results[[1]]
VesselDensity_bootstrap <-bootstrap_results[[2]] 
Kmax_bootstrap <- bootstrap_results[[3]]


#adding observed value
HD_bootstrap[1,] <- t(Obs_values[,2])
VesselDensity_bootstrap[1,] <- t(Obs_values[,3])
Kmax_boostrap[1,] <- t(Obs_values[,4])

# Create a list of the bootstrap matrices
boot_frames <- list(HD_bootstrap = HD_bootstrap, 
                    VesselDensity_bootstrap = VesselDensity_bootstrap, 
                    Kmax_bootstrap = Kmax_bootstrap)

# Calculate differences using lapply and bind the results into a data frame
boot_diffs <- lapply(boot_frames, function(matrix) {
  matrix[, 1] - matrix[, 2]
})

# Combine the differences into a single data frame
boot_diffs <- do.call(cbind, boot_diffs)
colnames(boot_diffs) <- names(boot_frames) 

## Create a vector of column names
column_names <- colnames(boot_diffs_df)

# Loop through each column in the boot_diffs data frame using lapply
lapply(seq_along(column_names), function(i) {
  x <- boot_diffs_df[, i]  # Get the current column data
  col_name <- column_names[i]  # Get the current column name
  
  # Create a density plot
  plot(density(x, na.rm = TRUE), main = "", xlab = "Differences", ylab = "Density")
  
  # Add quantile lines
  abline(v = quantile(x, c(0.025, 0.975), na.rm = TRUE), lwd = 2, col = "black")
  
  # Add a vertical line for the observed value
  abline(v = x[1], col = "red", lwd = 2)
  
  # Calculate p-value and display it on the plot
  p_value <- sum(x >= abs(x[1]), na.rm = TRUE) / length(x)
  text(x = mean(x, na.rm = TRUE), y = max(density(x)$y, na.rm = TRUE) * 0.9, paste("p-value =", round(p_value, 4)))
  
  # Title for the plot using the current column name
  title(main = col_name)
})




# Initialize an empty list to store the long-format data
long_data_list <- list()

# Loop through the matrices, making each long and adding to the list
for (name in names(boot_frames)) {
  
  # Access each matrix by name
  matrix_data <- boot_frames[[name]]
  
  # Convert the matrix to a data frame and make it long format
  long_data <- as.data.frame(matrix_data) %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(cols = c("Host", "Parasite"), names_to = "parasitism", values_to = name) # Use the matrix name as the value column
  
  # Store the long data frame in the list
  long_data_list[[name]] <- long_data
}

# Combine all the long-format data frames into a single data frame
long_df <- reduce(long_data_list, full_join, by = c("iteration", "parasitism"))



# Ensure long_df contains the necessary columns
colnames(long_df)

# Create a dataframe for confidence intervals
ci_data <- data.frame(
  parasitism = rep(c("Host", "Parasite"), times = 3),
  variable = rep(c("HD_bootstrap", "VesselDensity_bootstrap", "Kmax_bootstrap"), each = 2),
  lower = c(HD_CI95[, 1], VD_CI95[, 1], Kmax_CI95[, 1]),
  upper = c(HD_CI95[, 2], VD_CI95[, 2], Kmax_CI95[, 2])
)

# List of variable names to loop over
vars <- c("HD_bootstrap", "VesselDensity_bootstrap", "Kmax_bootstrap")

# Loop through variables and create density plots
for (v in vars) {
  # Extract the data for parasitism and the current variable
  data <- long_df[, c("parasitism", v)]
  
  # Dynamically assign the variable name for x
  var_name <- v
  
  # Create the plot using aes()
  p <- ggplot(data, aes(x = .data[[var_name]], fill = parasitism)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Parasite" = "red", "Host" = "grey")) +
    # Add dashed lines for confidence intervals from ci_data
    geom_vline(data = ci_data[ci_data$variable == var_name & ci_data$parasitism == "Parasite", ],
               aes(xintercept = lower), linetype = "dashed", color = "red") +
    geom_vline(data = ci_data[ci_data$variable == var_name & ci_data$parasitism == "Parasite", ],
               aes(xintercept = upper), linetype = "dashed", color = "red") +
    geom_vline(data = ci_data[ci_data$variable == var_name & ci_data$parasitism == "Host", ],
               aes(xintercept = lower), linetype = "dashed", color = "black") +
    geom_vline(data = ci_data[ci_data$variable == var_name & ci_data$parasitism == "Host", ],
               aes(xintercept = upper), linetype = "dashed", color = "black") +
    scale_x_continuous(limits = c(min(data[[var_name]]) * 0.2, max(data[[var_name]]) * 1.5))+
    theme_classic() +
    labs(title = paste("Density Plot for", var_name))
  
  # Print the plot
  print(p)
}


# The results are now stored in the bootstrap_results list, accessible by name


hydrHydraulic_Results <- t(replicate(iterations,
                                     shuffle_means(Hydraulic_means,cols = v,cat = "parasitism")))












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

