######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.1.2 - Data Analysis Ressampling- Vessel walls
#################################################################
library(here)
source(here("scripts", "02_1-TestAssumptions-WallThickness.R"))
rm(list=ls())

wdata <- read.csv(here("data", "processed", "wdata.csv"))
wdata_clean <- read.csv(here("data", "processed", "wdata_clean.csv"))



# List of data frame names
dataframes <- ls()
# Relevel the factors for each data frame
for (df_name in dataframes) {
  df <- get(df_name) # Get the data frame by name
  
  if ("ssp" %in% colnames(df)) { # Check if 'ssp' column exists
    df$ssp <- factor(df$ssp, levels = c(
      "Psittacanthus robustus", "Vochysia thyrsoidea",
      "Phoradendron perrotettii", "Tapirira guianensis",
      "Struthanthus rhynchophyllus", "Tipuana tipu",
      "Viscum album", "Populus nigra"
    ))
    assign(df_name, df) # Assign the modified data frame back to its name
  }
  rm(df, df_name) # remove duplicated dataframe
}

# Data subsets by species
species_data <- list(
  "Psittacanthus robustus" = wdata_clean[wdata_clean$ssp == "Psittacanthus robustus", ],
  "Vochysia thyrsoidea" = wdata_clean[wdata_clean$ssp == "Vochysia thyrsoidea", ],
  "Phoradendron perrotettii" = wdata_clean[wdata_clean$ssp == "Phoradendron perrotettii", ],
  "Tapirira guianensis" = wdata_clean[wdata_clean$ssp == "Tapirira guianensis", ],
  "Struthanthus rhynchophyllus" = wdata_clean[wdata_clean$ssp == "Struthanthus rhynchophyllus", ],
  "Tipuana tipu" = wdata_clean[wdata_clean$ssp == "Tipuana tipu", ],
  "Viscum album" = wdata_clean[wdata_clean$ssp == "Viscum album", ],
  "Populus nigra" = wdata_clean[wdata_clean$ssp == "Populus nigra", ]
)

# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

# Calculate means for each label and retain all columns
WT_avg <- wdata %>%
  group_by(label, ssp) %>%
  summarise(
    wthickness = mean(wthickness, na.rm = TRUE),
    .groups = "drop"  # Drop the grouping after summarizing
  )
WT_avg <- WT_avg %>%
  mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", 
               "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", 
               "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  ))

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


bootstrap_calculation <- function(x, cols, cat) {
  # Resample the data frame with replacement
  resampled_data <- x[sample(nrow(x), replace = TRUE), ]
  
  # Initialize a numeric vector to store the results
  results <- numeric(length(cols))
  
  # Loop through each numeric column to calculate differences in means
  for (i in seq_along(cols)) {
    column_name <- cols[i]
    group_means <- tapply(resampled_data[[column_name]], resampled_data[[cat]], mean)
    
    
    # Check if we have more than one group
    if (length(group_means) > 1) {
      results[i] <- diff(group_means)
    } else {
      results[i] <- NA  # Assign NA if not enough groups for diff
    }
  }
  
  # Return the results as a named vector
  names(results) <- cols
  return(results)
}


# Number of iterations for the permutation test
iterations <- 100000

WT_obs <- diff(tapply(WT_avg$wthickness,WT_avg$parasitism, FUN = mean))

WT_resample <- t(replicate(iterations,
                           shuffle_calculation(
                             x=WT_avg,cols ="wthickness",cat = "parasitism"))) %>%
  matrix(ncol = 1)
colnames(WT_resample) <- "wthickness"
WT_resample[1] <- WT_obs
#calculate p-values
WT_pvalue <- sum( WT_resample <=  WT_resample[1]|
                    WT_resample >= ( WT_resample[1]*-1))/
  length(WT_resample)
#vessel density
plot(density(WT_resample))
abline(v = WT_resample[1], col = "red")
abline(v = WT_resample[1]*-1, col = "red")
text(x=0.5,y=0.5,paste("p-value=",WT_pvalue))



# Perform the bootstrap calculation multiple times and store results
Bootstrap_Results <- replicate(
  iterations, bootstrap_calculation(x = WT_avg, cols = "wthickness", cat = "parasitism"))

# Convert to a matrix and assign column names
Bootstrap_Results <- matrix(Bootstrap_Results, ncol = 1) %>% na.rm()
colnames(Bootstrap_Results) <- "wthickness"
Bootstrap_Results[1] <- WT_obs
# View the result
View(Bootstrap_Results)


# Calculate bootstrap confidence intervals
CI_95 <- apply(Bootstrap_Results, 2, quantile, probs = c(0.025, 0.975)) 

plot(density(Bootstrap_Results), main = "Bootstrap - WThickness", xlab = "Difference in means")
abline(v = VD_obs, col = "red", lwd = 2)
abline(v=CI_95[,1],col="blue",lwd=2)

fwrite(Hydraulic_Results,file=here("outputs","tables","Hydraulic_MonteCarlo_results.csv"))

