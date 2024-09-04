######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.3 - Data Analysis - Pit membranes
#################################################################

source(here("scripts", "02_3-TestAssumptions-Pit.R"))
rm(list = ls())

# Load data
pitdata_clean <- fread(here("data", "processed", "pitdata_clean.csv"))
pitOdata_clean <- fread(here("data", "processed", "pitO_clean.csv"))

# Relevel the factors
if ("ssp" %in% colnames(pitdata_clean)) {
  pitdata_clean$ssp <- factor(pitdata_clean$ssp, levels = c(
    "Psittacanthus robustus", "Vochysia thyrsoidea",
    "Phoradendron perrotettii", "Tapirira guianensis",
    "Struthanthus rhynchophyllus", "Tipuana tipu",
    "Viscum album", "Populus nigra"
  ))
}

if ("ssp" %in% colnames(pitOdata_clean)) {
  pitOdata_clean$ssp <- factor(pitOdata_clean$ssp, levels = c(
    "Psittacanthus robustus", "Vochysia thyrsoidea",
    "Phoradendron perrotettii", "Tapirira guianensis",
    "Struthanthus rhynchophyllus", "Tipuana tipu",
    "Viscum album", "Populus nigra"
  ))
}

# Initialize result data frames
PitMembrane_results <- data.frame(Parasite = character(),
                                  Host = character(),
                                  ParasiteMean = numeric(),
                                  HostMean = numeric(),
                                  MeanDifference = numeric(),
                                  pvalue = numeric(),
                                  stringsAsFactors = FALSE)

pcd_results <- data.frame(Parasite = character(),
                          Host = character(),
                          ParasiteMean = numeric(),
                          HostMean = numeric(),
                          MeanDifference = numeric(),
                          pvalue = numeric(),
                          stringsAsFactors = FALSE)

pitDiameter_results <- data.frame(Parasite = character(),
                                  Host = character(),
                                  ParasiteMean = numeric(),
                                  HostMean = numeric(),
                                  MeanDifference = numeric(),
                                  pvalue = numeric(),
                                  stringsAsFactors = FALSE)

pitOpening_results <- data.frame(Parasite = character(),
                                 Host = character(),
                                 ParasiteMean = numeric(),
                                 HostMean = numeric(),
                                 MeanDifference = numeric(),
                                 pvalue = numeric(),
                                 stringsAsFactors = FALSE)


# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

# Perform t-tests for each measure
for (pair in species_pairs) {
  
  parasite_data <- pitdata_clean[pitdata_clean$ssp == pair[1], ]
  host_data <- pitdata_clean[pitdata_clean$ssp == pair[2], ]
  
  # Handle errors using tryCatch
  tryCatch({
    # Perform t-test for pit membrane thickness
    ttest <- t.test(parasite_data$pitavg, host_data$pitavg, var.equal = FALSE)
    PitMembrane_results <- rbind(PitMembrane_results, data.frame(
      Parasite = pair[1],
      Host = pair[2],
      ParasiteMean = ttest$estimate[1] * 1000,
      HostMean = ttest$estimate[2] * 1000,
      MeanDifference = diff(ttest$estimate) * 1000,
      pvalue = ttest$p.value,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    message("Error with pair ", pair[1], " vs ", pair[2], ": ", e$message)
  })
}

# Perform t-tests for Pit chamber depth
for (pair in species_pairs) {
  
  parasite_data <- pitdata_clean[pitdata_clean$ssp == pair[1], ]
  host_data <- pitdata_clean[pitdata_clean$ssp == pair[2], ]
  
  tryCatch({
    # Perform t-test for pit chamber depth
    ttest <- t.test(parasite_data$pcd, host_data$pcd, var.equal = FALSE)
    pcd_results <- rbind(pcd_results, data.frame(
      Parasite = pair[1],
      Host = pair[2],
      ParasiteMean = ttest$estimate[1] * 1000,
      HostMean = ttest$estimate[2] * 1000,
      MeanDifference = diff(ttest$estimate) * 1000,
      pvalue = ttest$p.value,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    message("Error with pair ", pair[1], " vs ", pair[2], ": ", e$message)
  })
}

# Perform t-tests for Pit Diameter
for (pair in species_pairs) {
  
  parasite_data <- pitOdata_clean[pitOdata_clean$ssp == pair[1], ]
  host_data <- pitOdata_clean[pitOdata_clean$ssp == pair[2], ]
  
  tryCatch({
    # Perform t-test for pit diameter
    ttest <- t.test(parasite_data$PitDiameter, host_data$PitDiameter, var.equal = FALSE)
    pitDiameter_results <- rbind(pitDiameter_results, data.frame(
      Parasite = pair[1],
      Host = pair[2],
      ParasiteMean = ttest$estimate[1],
      HostMean = ttest$estimate[2],
      MeanDifference = diff(ttest$estimate),
      pvalue = ttest$p.value,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    message("Error with pair ", pair[1], " vs ", pair[2], ": ", e$message)
  })
}

# Perform t-tests for Pit Opening Diameter
for (pair in species_pairs) {
  
  parasite_data <- pitOdata_clean[pitOdata_clean$ssp == pair[1], ]
  host_data <- pitOdata_clean[pitOdata_clean$ssp == pair[2], ]
  
  tryCatch({
    # Perform t-test for pit opening diameter
    ttest <- t.test(parasite_data$PitOpening, host_data$PitOpening, var.equal = FALSE)
    pitOpening_results <- rbind(pitOpening_results, data.frame(
      Parasite = pair[1],
      Host = pair[2],
      ParasiteMean = ttest$estimate[1],
      HostMean = ttest$estimate[2],
      MeanDifference = diff(ttest$estimate),
      pvalue = ttest$p.value,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    message("Error with pair ", pair[1], " vs ", pair[2], ": ", e$message)
  })
}

PitMembrane_results
pcd_results
pitDiameter_results
pitOpening_results


# Save results to CSV
fwrite(PitMembrane_results, file = here("outputs", "tables", "pit_membrane_diff.csv"))
fwrite(pcd_results, file = here("outputs", "tables", "pcd_results.csv"))
fwrite(pitDiameter_results, file = here("outputs", "tables", "PitDiameter_results.csv"))
fwrite(pitOpening_results, file = here("outputs", "tables", "PitOpening_results.csv"))

