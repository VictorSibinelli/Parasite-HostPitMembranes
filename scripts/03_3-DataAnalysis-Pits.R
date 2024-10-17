######################################################################
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.3 - Data Analysis - Pit membranes - Bootstrap
#################################################################

# Load script for testing assumptions of pit data and remove any existing objects from the environment
source(here("scripts", "02_3-TestAssumptions-Pit.R"))
rm(list = ls())  # Clear environment
source(here("scripts", "Functions.R"))  # Load custom functions

# Load pre-processed pit membrane data
pitdata <- fread(here("data", "processed", "pitdata.csv")) %>% as.tibble()  # Load main pit data
pitOdata <- fread(here("data", "processed", "pitOdata.csv")) %>% as.tibble()# Load pit opening data


library(dplyr)

# Load your data (assuming these lines are already included)
pitdata <- fread(here("data", "processed", "pitdata.csv")) %>% as.tibble()  # Load main pit data
pitOdata <- fread(here("data", "processed", "pitOdata.csv")) %>% as.tibble()  # Load pit opening data

# Keep only the first occurrence of each 'ssp' (assuming you want the first based on the order in the data)


# Summarize data for each species and parasite/host
Pit_EV <- pitOdata %>%
  group_by(label) %>%
  summarise(
    ssp = first(ssp),                # Extract first species (ssp) name for each label
    parasitism = first(parasitism),   # Extract parasitism status (Parasite/Host)
    PitDiameter = median(PitDiameter, na.rm = TRUE),   # Calculate median PitDiameter
    PitOpening = median(PitOpening, na.rm = TRUE),     # Calculate median PitOpening
    .groups = "drop"
  )

# List of species pairs for comparisons
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),   # Pair 1
  c("Phoradendron perrotettii", "Tapirira guianensis"), # Pair 2
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),     # Pair 3
  c("Viscum album", "Populus nigra")                    # Pair 4
)
# Define the variables of interest for bootstrapping
vars <- colnames(Pit_EV[4:5])  # For Pit_EV: PitDiameter and PitOpening
vars2 <- colnames(pitdata[, c(7, 11)])  # For pitdata: columns 7, 10, and 11

# Relevel factors for analysis
relevel_factors(ls())  # Reorder factors for categorical variables

# Set the number of bootstrap iterations and random seed for reproducibility
it <- 2  # Number of bootstrap replicates
set.seed(42)  # Set seed for consistent random sampling

# Bootstrap sampling function
bootstrap_sampling <- function(data, index, vars) {
  replicate(n = it, {
    sapply(vars, function(var) {
      data %>%
        subset(parasitism == index) %>%  # Subset data for the specified parasitism
        shuffle_means(cols = var, cat = "parasitism", rcol = TRUE)  # Shuffle and calculate means
    })
  }, simplify = TRUE) %>% t()  # Transpose the result
}

# Bootstrap for Host and Parasite species in Pit_EV
host_boot <- bootstrap_sampling(Pit_EV, "Host", vars)
para_boot <- bootstrap_sampling(Pit_EV, "Parasite", vars)

# Bootstrap for Host and Parasite species in pitdata
host_boot2 <- bootstrap_sampling(pitdata, "Host", vars2)
para_boot2 <- bootstrap_sampling(pitdata, "Parasite", vars2)

# Bootstrap for each species (ssp) in Pit_EV
ssp_boot <- setNames(lapply(levels(Pit_EV$ssp), function(current_ssp) {
  t(replicate(it, sapply(vars, function(var) {
    shuffle_means(Pit_EV[Pit_EV$ssp == current_ssp, ], cols = var, cat = "ssp", rcol = TRUE)
  })))
}), levels(Pit_EV$ssp))  # Set species names as list element names

# Bootstrap for each species (ssp) in pitdata
ssp_boot2 <- setNames(lapply(levels(pitdata$ssp), function(current_ssp) {
  t(replicate(it, sapply(vars2, function(var) {
    shuffle_means(pitdata[pitdata$ssp == current_ssp, ], cols = var, cat = "ssp", rcol = TRUE)
  })))
}), levels(pitdata$ssp))  # Set species names as list element names

# Combine bootstrap results into a single list
CI_boot <- list(
  Parasite = cbind(para_boot, para_boot2),  # Combine parasite results from both data sources
  Host = cbind(host_boot, host_boot2)       # Combine host results from both data sources
)

# Combine species-specific results into the combined list
CI_boot <- c(combined_boot, lapply(levels(Pit_EV$ssp), function(current_ssp) {
  cbind(
    ssp_boot[[current_ssp]],
    ssp_boot2[[current_ssp]]
  )
}))

# Name each species element appropriately
names(CI_boot)[3:length(CI_boot)] <- levels(Pit_EV$ssp)

boot_results <- list()



