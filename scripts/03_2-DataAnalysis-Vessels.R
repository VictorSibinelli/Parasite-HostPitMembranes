######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.2 - Data Analysis - Vessels
#################################################################
library(here)
#source(here("scripts", "02_2-TestAssumptions-Vessels.R"))
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
    if ("parasitism" %in% colnames(df)) {  # If the object contains the 'ssp' column
      df$parasitism <- factor(df$parasitism, levels = c("Parasite","Host"))
      assign(df_name, df)  # Save the modified data frame
    }
    rm(df, df_name)  # Remove temporary variables
  }}
# Number of iterations for the permutation test
iterations <- 100
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
                                           shuffle_means(
                                             Hydraulic_means, cols = v, cat = "parasitism",rcol = T)))
}

lapply(bootstrap_results,head)

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
Kmax_bootstrap[1,] <- t(Obs_values[,4])

# Create a list of the bootstrap matrices
boot_frames <- list(HD_bootstrap = HD_bootstrap, 
                    VesselDensity_bootstrap = VesselDensity_bootstrap, 
                    Kmax_bootstrap = Kmax_bootstrap)

# Calculate differences using lapply and bind the results into a data frame
boot_diffs <- lapply(boot_frames, function(matrix) {
  matrix[, 1] - matrix[, 2]
}) %>% do.call(what=cbind,)

colnames(boot_diffs) <- names(boot_frames) 




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

p_values <- apply(boot_diffs,2,function(x){
  sum(abs(x) >= abs(x[1]), na.rm = TRUE) / length(x)
})


# Function to shuffle and calculate means for each variable (Host)
host_boot <- replicate(n = iterations, {
  sapply(vars, function(var) {
    Hydraulic_means %>%
      subset(parasitism == "Host") %>%
      shuffle_means(cols = var, cat = "parasitism", rcol = TRUE)
  })
}, simplify = TRUE)

# Convert results to a data frame for each variable
host_boot <- as.data.frame(t(host_boot))
colnames(host_boot) <- vars

# Convert to a matrix (if needed)
host_boot <- as.matrix(host_boot)

# Function to shuffle and calculate means for each variable (Parasite)
para_boot <- replicate(n = iterations, {
  sapply(vars, function(var) {
    Hydraulic_means %>%
      subset(parasitism == "Parasite") %>%
      shuffle_means(cols = var, cat = "parasitism", rcol = TRUE)
  })
}, simplify = TRUE)

# Convert results to a data frame for each variable
para_boot <- as.data.frame(t(para_boot))
colnames(para_boot) <- vars

# Convert to a matrix (if needed)
para_boot <- as.matrix(para_boot)

# Calculate CI95 for each variable (using 2.5th and 97.5th percentiles)
host_CI95 <- t(apply(host_boot, 2, function(x) {
  quantile(x, c(0.025, 0.975), na.rm = TRUE)
}))

para_CI95 <- t(apply(para_boot, 2, function(x) {
  quantile(x, c(0.025, 0.975), na.rm = TRUE)
}))

# Add the observed values to each as a new column
para_CI95 <- cbind(para_CI95, Obs = as.numeric(Obs_values[1, 2:4]))
host_CI95 <- cbind(host_CI95, Obs = as.numeric(Obs_values[2, 2:4]))



#########pair wise comparisson
# Initialize an empty list to store results for each species pair
ssp_obs <- apply(Hydraulic_means[,vars],2,function(x){
  tapply(x,Hydraulic_means$ssp,mean)
})



pair_boot=list()
# Iterate through each species pair
for (pair in species_pairs) {
  # Subset the data for the current species pair
  subset_data <- subset(Hydraulic_means, ssp %in% pair)
  
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

# Initialize an empty list to store the results
ssp_diffs_list <- vector("list", length = length(pair_boot))

# Calculate differences between species in pair_boot
ssp_diffs_list <- sapply(pair_boot, simplify = FALSE, function(x) {
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
    diff_result <- split_data[[1]][, vars] - split_data[[2]][, vars]
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
ssp_diffs_list$`Psittacanthus robustusvsVochysia thyrsoidea`[1, ] <- pair_obs_diff[[1]]
ssp_diffs_list$`Phoradendron perrotettiivsTapirira guianensis`[1, ] <- pair_obs_diff[[2]]
ssp_diffs_list$`Struthanthus rhynchophyllusvsTipuana tipu`[1, ] <- pair_obs_diff[[3]]

# Remove any NULL elements from ssp_diffs_list before calculating p-values
valid_diffs_list <- Filter(Negate(is.null), ssp_diffs_list)

# Calculate p-values for each element in the valid ssp_diffs_list
ssp_pvalues <- sapply(valid_diffs_list, simplify = TRUE, function(x) {
  apply(x, 2, function(y) {
    sum(abs(y) >= abs(y[1]), na.rm = TRUE) / length(y)
  })
}) %>% t()

# Print the results to verify
print(ssp_pvalues)


# Initialize an empty list to store results
ssp_CI95 <- list()

# Calculate confidence intervals for each species
for (y in pair_boot) {
  # Use tryCatch to handle potential errors
  tryCatch({
    if (ncol(y) >= 5 && "ssp" %in% colnames(y)) {
      # Split the data by 'ssp' and calculate confidence intervals
      split_data <- split(y[, 3:5], y$ssp)
      for (ssp_name in names(split_data)) {
        ci_df <- as.data.frame(t(apply(split_data[[ssp_name]], 2, function(col) {
          quantile(col, c(0.025, 0.975), na.rm = TRUE)
        })))
        # Add variable names and reorder columns
        ci_df <- cbind(variable = rownames(ci_df), ci_df)
        colnames(ci_df) <- c("variable", "lower", "upper")
        rownames(ci_df) <- NULL
        
        # Combine results for the current species
        if (is.null(ssp_CI95[[ssp_name]])) {
          ssp_CI95[[ssp_name]] <- ci_df  # Initialize if NULL
        } else {
          ssp_CI95[[ssp_name]] <- rbind(ssp_CI95[[ssp_name]], ci_df)  # Combine results
        }
      }
    } else {
      warning("Data frame does not have expected structure or 'ssp' column is missing.")
    }
  }, error = function(e) {
    # Print a message and continue if there's an error
    message("Error processing one of the data frames: ", conditionMessage(e))
  })
}

# Print the resulting list of data frames
print(ssp_CI95)


print("Data analysis - vessels completed. Head to Graphics vessels")


