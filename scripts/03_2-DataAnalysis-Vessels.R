######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.2 - Data Analysis - Vessels
#################################################################
library(here)
rm(list = ls())

source(here("scripts", "Functions.R"))
#load data
HydraulicData <- read.csv(here("data", "processed", "HydraulicData.csv"))
vdata <- read.csv(here("data", "processed", "vdata.csv"))
vadata <- read.csv(here("data", "processed", "vadata.csv"))



# List of species pairs
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu") ,
  c("Viscum album", "Populus nigra")
)


Hydraulic_EV <- HydraulicData %>%
  group_by(indiv) %>%
  summarise(
    ssp = first(ssp),  # Or use any other summary function if needed
    HydraulicDiameter = median(HydraulicDiameter, na.rm = TRUE),
    vdensity = median(vdensity, na.rm = TRUE),
    Kmax = median(Kmax, na.rm = TRUE),
    parasitism = first(parasitism),  # This assumes the same for each individual
    .groups = "drop"
  )

relevel_factors(ls())
# Number of iterations for the permutation test
iterations <- 500000
set.seed(42)
# Specify the columns for which you want to calculate bootstrap values
vars <- colnames(Hydraulic_EV)[3:5]

# Initialize a list to store results for each variable
bootstrap_results <- list()

# # Loop through each variable to calculate bootstrap values
for (v in vars) {

  # Generate descriptive name for storing bootstrap values
  bootstrap_results[[v]] <- t(replicate(iterations,
                                        shuffle_means(
                                          Hydraulic_EV,
                                          cols = v,
                                          cat = "parasitism",
                                          rcol = TRUE)))

  # Assign column names based on parasitism levels
  colnames(bootstrap_results[[v]]) <- levels(Hydraulic_EV$parasitism)
}

lapply(bootstrap_results,head)


Obs_values <-  Hydraulic_EV %>% group_by(parasitism) %>%
  summarize(HydraulicDiameter = mean(HydraulicDiameter, na.rm = TRUE),
            vdensity = mean(vdensity, na.rm = TRUE),
            Kmax = mean(Kmax, na.rm = TRUE),
            parasitism = first(parasitism),  # This assumes the same for each individual
            .groups = "drop")

bootstrap_results[[1]][1,] <- t(Obs_values[,2])
bootstrap_results[[2]][1,] <- t(Obs_values[,3])
bootstrap_results[[3]][1,] <- t(Obs_values[,4])


for (t in seq_along(bootstrap_results)){
  # Create the file name for each table, appending "_bootstrap" before the file extension
  table_name <- paste(names(bootstrap_results)[t], "_boot", ".csv", sep = "")

  # Use fwrite to save each table to the specified directory
  fwrite(bootstrap_results[[t]],
         file = here("data", "processed", "ressampled", table_name))

  # Print a success message to the console
  cat(paste(table_name, "successfully saved\n"))
}


# Read the CSV files into a list with specified names
bootstrap_results <- list(
  HydraulicDiameter = read.csv(here("data", "processed", "ressampled", "HydraulicDiameter_boot.csv")),
  VesselDensity = read.csv(here("data", "processed", "ressampled", "vdensity_boot.csv")),
  Kmax = read.csv(here("data", "processed", "ressampled", "Kmax_boot.csv"))
)


# Calculate differences using lapply and bind the results into a data frame

boot_diffs <- bootstrap_results %>% lapply(FUN = function(x){
  x[,1]-x[,2]
}) %>% do.call(what=cbind)


## Create a vector of column names
column_names <- colnames(boot_diffs)
# Loop through each column in the boot_diffs data frame using lapply
sapply(seq_along(column_names), function(i) {
  x <- boot_diffs[, i]  # Get the current column data
  col_name <- column_names[i]  # Get the current column name
  
  # Create a density plot
  plot(density(x, na.rm = TRUE), main = "", xlab = "Differences", ylab = "Density")
  
  # Add quantile lines
  abline(v = quantile(x, c(0.025, 0.975), na.rm = TRUE), lwd = 2, col = "black")
  
  # Add a vertical line for the observed value
  abline(v = x[1], col = "red", lwd = 2)
  
  # Calculate p-value and display it on the plot
  p_value <- sum(abs(x) >= abs(x[1]), na.rm = TRUE) / length(x)
  text(x = mean(x, na.rm = TRUE), y = max(density(x)$y, na.rm = TRUE) * 0.9, paste("p-value =",p_value))
  
  # Title for the plot using the current column name
  title(main = col_name)
})

# Calculate p-values
p_values <- apply(boot_diffs, 2, function(x) {
  sum(abs(x) >= abs(x[1]), na.rm = TRUE) / length(x)
}) %>% as.data.frame() %>% t()
row.names(p_values) <- "Parasite vs Host"
print(p_values)


# Function to shuffle and calculate means for each variable (Host)
host_boot <- replicate(n = iterations, {
  sapply(vars, function(var) {
    Hydraulic_EV %>%
      subset(parasitism == "Host") %>%
      shuffle_means(cols = var, cat = "parasitism", rcol = TRUE)
  })
}, simplify = TRUE) %>% t()


# Function to shuffle and calculate means for each variable (Parasite)
para_boot <- replicate(n = iterations, {
  sapply(vars, function(var) {
    Hydraulic_EV %>%
      subset(parasitism == "Parasite") %>%
      shuffle_means(cols = var, cat = "parasitism", rcol = TRUE)
  })
}, simplify = TRUE) %>% t()


fwrite(host_boot,file = here("data","processed","ressampled","Host_Vessels_boot.csv"))
fwrite(para_boot,file = here("data","processed","ressampled","Parasite_Vessels_boot.csv"))

para_boot <- fread(here("data","processed","ressampled","Parasite_Vessels_boot.csv"))
colnames(para_boot)[2] <- "VesselDensity"
host_boot <- fread(here("data","processed","ressampled","Host_Vessels_boot.csv"))
colnames(host_boot)[2] <- "VesselDensity"
# Calculate CI95 for each variable (using 2.5th and 97.5th percentiles)
host_CI95 <- t(apply(host_boot, 2, function(x) {
  quantile(x, c(0.025, 0.975), na.rm = TRUE)
})) %>% as.data.frame()

para_CI95 <- t(apply(para_boot, 2, function(x) {
  quantile(x, c(0.025, 0.975), na.rm = TRUE)
})) %>% as.data.frame()

# Add the observed values to each as a new column
para_CI95 <- cbind(para_CI95, Obs = as.numeric(Obs_values[1, 2:4]))
host_CI95 <- cbind(host_CI95, Obs = as.numeric(Obs_values[2, 2:4]))
colnames(para_CI95) <- c("Lower","Upper","Mean")
colnames(host_CI95) <- c("Lower","Upper","Mean")
host_CI95$ssp <- "Host"
para_CI95$ssp <- "Parasite"


#########pair wise comparisson

ssp_obs <- apply(Hydraulic_EV[,vars],2,function(x){
  tapply(x,Hydraulic_EV$ssp,mean)
})


iterations <- 250000
pair_boot=list()
# Iterate through each species pair


for (pair in species_pairs) {
  # Subset the data for the current species pair
  subset_data <- subset(Hydraulic_EV, ssp %in% pair)

  # Initialize an empty result data frame
  result <- data.frame()

  # Loop through each variable in vars
  for (v in vars) {
    # Generate bootstrap samples for the current variable
    resample_data <- as.data.frame(t(replicate(iterations, shuffle_means(subset_data, cols = v, cat = "ssp", rcol = TRUE))))

    # Add a bootstrap identifier
    resample_data$bootstrap_id <- 1:nrow(resample_data)

    # Transform to long format
    resample_data_long <- resample_data %>%
      pivot_longer(cols = -bootstrap_id,  # Exclude the identifier column
                   names_to = "ssp",       # New column for species names
                   values_to = "value")    # New column for values

    # Rename the value column to the variable name
    resample_data_long <- resample_data_long %>%
      rename(!!v := value)  # Using the variable name as the column name for values

    # Combine results
    if (nrow(result) == 0) {
      result <- resample_data_long  # If result is empty, initialize it with the first variable's data
    } else {
      result <- left_join(result, resample_data_long, by = c("bootstrap_id", "ssp"))  # Join by bootstrap_id and ssp
    }
  }

  pair_boot[[paste(pair,collapse = "vs")]] <- result
}

names(pair_boot)
lapply(pair_boot,head)
# Loop through each table in bootstrap_results and save it
for (t in seq_along(pair_boot)) {
  # Create the file name for each table, appending "_bootstrap" before the file extension
  table_name <- paste(names(pair_boot)[t], "_Vessels_bootstrap", ".csv", sep = "")

  # Use fwrite to save each table to the specified directory
  fwrite(pair_boot[[t]],
         file = here("data", "processed", "ressampled", table_name))

  # Print a success message to the console
  cat(paste(table_name, "successfully saved\n"))
}

# Create a list of species pairs with their respective file names
pair_boot2 <- list(
  "Psittacanthus robustusvsVochysia thyrsoidea"= fread(here("data", "processed", "ressampled","Psittacanthus robustusvsVochysia thyrsoidea_Vessels_bootstrap.csv")),
  "Phoradendron perrotettiivsTapirira guianensis" =fread(here("data", "processed", "ressampled","Phoradendron perrotettiivsTapirira guianensis_Vessels_bootstrap.csv")),
  "Struthanthus rhynchophyllusvsTipuana tipu" = fread(here("data", "processed", "ressampled","Struthanthus rhynchophyllusvsTipuana tipu_Vessels_bootstrap.csv")),
  "Viscum albumvsPopulus nigra"= fread(here("data", "processed", "ressampled","Viscum albumvsPopulus nigra_Vessels_bootstrap.csv"))
)


# Initialize an empty list to store the results
ssp_diffs_list <- vector("list", length = length(pair_boot2))

ssp_diffs_list <- sapply(pair_boot2, simplify = FALSE, function(x) {
  # Skip if the data frame is empty
  if (nrow(x) == 0) {
    warning("Skipping empty data frame.")
    return(NULL)  # Return NULL for empty data frames
  }
  
  # Split the data frame by 'ssp'
  split_data <- split(x, x$ssp)
  
  # Check if there are exactly 2 species
  if (length(split_data) == 2) {
    # Calculate the difference between the first and second species
    diff_result <- split_data[[1]][, ..vars, with = FALSE] - split_data[[2]][, ..vars, with = FALSE]
    return(as.data.frame(diff_result))  # Convert to data frame
  } else {
    # Print a warning message for pairs that do not have exactly two species
    warning(paste("Expected 2 species, but found:", length(split_data)))
    return(NULL)  # Return NULL if there are not exactly 2 species
  }
})
# Print the resulting list to verify
print(lapply(ssp_diffs_list, head))





## Check pair_obs_diff structure and make sure they are data frames
pair_obs_diff <- setNames(list(
  as.data.frame(t(as.matrix(ssp_obs[1, ] - ssp_obs[2, ]))),  # Convert to data frame
  as.data.frame(t(as.matrix(ssp_obs[3, ] - ssp_obs[4, ]))),  # Convert to data frame
  as.data.frame(t(as.matrix(ssp_obs[5, ] - ssp_obs[6, ]))),  # Convert to data frame
  as.data.frame(t(as.matrix(ssp_obs[7, ] - ssp_obs[8, ])))   # Convert to data frame
), species_pairs)

# Assign the named vectors to the first row of each element in ssp_diffs_list
ssp_diffs_list$`Psittacanthus robustusvsVochysia thyrsoidea` [1,]<- pair_obs_diff[[1]]
ssp_diffs_list$`Phoradendron perrotettiivsTapirira guianensis`[1, ] <- pair_obs_diff[[2]]
ssp_diffs_list$`Struthanthus rhynchophyllusvsTipuana tipu`[1, ] <- pair_obs_diff[[3]]
#ssp_diffs_list$`Viscum album vs Populus nigra`[1, ] <- pair_obs_diff[[4]]

# Remove any NULL elements from ssp_diffs_list before calculating p-values
valid_diffs_list <- Filter(Negate(is.null), ssp_diffs_list)

# Calculate p-values for each element in the valid ssp_diffs_list
ssp_pvalues <- sapply(valid_diffs_list, simplify = TRUE, function(x) {
  apply(x, 2, function(y) {
    sum(abs(y) >= abs(y[1]), na.rm = TRUE) / length(y)
  })
}) %>% t() %>% as.data.frame()
colnames(ssp_pvalues)[2] <- "VesselDensity"
# Print the results to verify

P_values <- rbind(p_values,ssp_pvalues)
rownames(P_values)[1] <- "Parasite vs Host"
print(P_values)

# Initialize a list to store bootstrap results
ssp_CI95 <- vector("list", iterations)

# Bootstrap sampling and calculation of CI95
# Using replicate to perform bootstrap sampling
iterations <- 100000
ssp_CI95 <- replicate(
  iterations,
  Hydraulic_EV %>%
    group_by(ssp) %>%
    summarise(
      HydraulicDiameter = mean(sample(HydraulicDiameter, replace = TRUE)),
      VesselDensity = mean(sample(vdensity, replace = TRUE)),
      Kmax = mean(sample(Kmax, replace = TRUE)),
      .groups = 'drop'
    ),
  simplify = FALSE  # Ensure the result is a list
)


# Combine all bootstrap samples into a data frame
ssp_CI95 <- bind_rows(ssp_CI95, .id = "iteration")

# Create a list to hold CI95 for each variable
CI95 <- list()
variables <- c("HydraulicDiameter", "VesselDensity", "Kmax")

# Calculate CI95 for each variable and store in the list
for (var in variables) {
  CI95[[var]] <- ssp_CI95 %>%
    group_by(ssp) %>%
    summarise(
      Lower = quantile(get(var), 0.025),
      Mean = mean(get(var)),
      Upper = quantile(get(var), 0.975),
      .groups = 'drop'
    )
}

# Loop through each variable and update CI95
for (var in variables) {
  CI95[[var]] <- bind_rows(
    para_CI95 %>% filter(rownames(para_CI95) == var) %>% mutate(ssp = "Parasite", .before = 1),
    host_CI95 %>% filter(rownames(host_CI95) == var) %>% mutate(ssp = "Host", .before = 1),
    CI95[[var]] %>% select(ssp, Lower, Mean, Upper)  # Reorder to ensure ssp is first
  ) %>% select(ssp, Lower, Mean, Upper)
  
  rownames(CI95[[var]]) <- NULL  # Remove rownames
}

# Convert `ssp_obs` to a tibble and add the `parasitism` (or `Group`) column
ssp_obs <- as_tibble(ssp_obs, rownames = "parasitism")

# Bind rows while keeping numeric columns intact
C_obs <- rbind(Obs_values, ssp_obs) %>% 
  rename(Group = parasitism)

fwrite(C_obs,here("outputs","tables","Vessel_obs.csv"))
# Display the updated CI95 list
print(CI95)
print(P_values)

fwrite(P_values,here("outputs","tables","Vessels_MonteCarlo_Pvalues.csv"))
# Save each data frame in CI95 list to separate CSV files
lapply(names(CI95), function(name) {
  fwrite(CI95[[name]], file = here("outputs", "tables", paste0("Vessels_MonteCarlo_CI95_", name, ".csv")))
})

print("Data analysis - vessels completed. Head to Graphics vessels")

