######################################################################
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 01 - Data Wrangling
######################################################################

# Load required libraries and custom functions
library(here)
source(here("scripts", "00-library.R"))
library(tidyverse)
library(data.table)

# Load and clean Wall data
Wall_data <- read.csv(here("data", "raw", "VesselWall.csv"))
colnames(Wall_data)[3] <- "WallThickness"
Wall_data <- within(Wall_data, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)  # Correct misspelling
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
  ssp <- gsub("Populus nigra ", "Populus nigra", ssp)  # Remove extra space
})
Wall_data$WallThickness <- ifelse(Wall_data$ssp != "Populus nigra", 
                                  Wall_data$WallThickness * 1.43, 
                                  Wall_data$WallThickness)  # Adjust thickness for non-"Populus nigra"
Wall_data$WallThickness <- Wall_data$WallThickness / 2  # Convert thickness unit
Wall_data$WallThickness[Wall_data$WallThickness < 0.2] <- NA  # Replace small values with NA

# Load and calculate vessel density
VDensity_data <- read.csv(here("data", "raw", "VesselCount.csv"))[-9]
VDensity_data$VesselDensity <- (VDensity_data$Count / (VDensity_data$Total.Area / VDensity_data$X.Area)) * 10000
VDensity_data <- within(VDensity_data, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
  Label <- gsub("Phoradendeon", "Phoradendron", Label)
})

# Load and clean Vessel Diameter data
VesselDiameter_data <- read.csv(here("data", "raw", "VesselDiameter.csv"))
VesselDiameter_data <- within(VesselDiameter_data, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  ssp <- gsub("Tapirira guianensis ", "Tapirira guianensis", ssp)  # Correct spacing
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
  Label <- gsub("Phoradendeon", "Phoradendron", Label)
})
VesselDiameter_data$VesselDiameter <- sqrt(VesselDiameter_data$Area * VesselDiameter_data$Circ.) * 2  # Calculate diameter

# Load and clean Pit Fraction data
PitFraction_data <- read.csv(here("data", "raw", "PitFraction.csv"))[-c(3:5)]
PitFraction_data <- within(PitFraction_data, {
  indiv <- gsub("Tapirira guianensis 1 ", "Tapirira guianensis 1", indiv)  # Correct spacing
  ssp <- gsub("perrottetti", "perrotettii", ssp)  # Correct species name
  indiv <- gsub("perrottetti", "perrotettii", indiv)
})

# Load and clean Pit Opening data
PitOp <- read.csv(here("data", "raw", "PitOpening.csv"))[-3]
PitOp <- within(PitOp, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
})
PitOp <- PitOp %>%  group_by(ssp, indiv) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

# Load and clean Pit Diameter data
PitDi <- read.csv(here("data","raw","PitDiameter.csv"))[-3]
PitDi <- within(PitDi, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
})
PitDi <- PitDi %>% group_by(ssp, indiv) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

# Merge Pit Diameter and Opening data
PitDiOp_data <- merge(
  PitDi,
  PitOp,
  by = c("ssp", "indiv", "row_id"), all = TRUE)[-3]
colnames(PitDiOp_data)[3:4] <- c("PitDiameter", "PitOpening")
PitDiOp_data$PitOpening <- ifelse(PitDiOp_data$ssp != "Populus nigra", 
                                  PitDiOp_data$PitOpening * 1.43, 
                                  PitDiOp_data$PitOpening)  # Adjust opening for non-"Populus nigra"

rm(PitDi, PitOp)  # Remove intermediate data

# Load and process Pit Membrane data
PitMembrane_data <- read.csv(here("data", "raw", "PitMembrane.csv"))[-8]
PitMembrane_data$Tpm <- rowMeans(PitMembrane_data[, 2:6], na.rm = TRUE) * 1000
PitMembrane_data$pcd <- PitMembrane_data$pcd*1000
PitMembrane_data$ssp <- gsub("^(.)", "\\U\\1", PitMembrane_data$ssp, perl = TRUE)
PitMembrane_data <- within(PitMembrane_data, {
  ssp <- gsub("Vochisya thirsoidea", "Vochysia thyrsoidea", ssp)  # Correct misspelling
  ssp <- gsub("perrotetti", "perrotettii", ssp)
})
PitMembrane_data <- PitMembrane_data[,c(1,7,8)]
# Calculate mean thickness

# Initialize Hydraulic Data
HydraulicData <- data.frame(
  ssp = NA,
  indiv = NA,
  Label = unique(VesselDiameter_data$Label),
  HydraulicDiameter = NA,
  stringsAsFactors = FALSE
)

# Compute Hydraulic Diameter for each unique label
for (i in seq_along(unique(VesselDiameter_data$Label))) {
  label_value <- unique(VesselDiameter_data$Label)[i]
  filter_data <- VesselDiameter_data[VesselDiameter_data$Label == label_value, ]
  
  HydraulicData[i, "HydraulicDiameter"] <- (sum(filter_data$VesselDiameter^4) / 
                                              length(filter_data$VesselDiameter))^(1/4)  # Calculate hydraulic diameter
  
  HydraulicData[i, "Label"] <- as.character(filter_data$Label[1])  # Assign label
  HydraulicData[i, "ssp"] <- as.character(filter_data$ssp[1])      # Assign species
  HydraulicData[i, "indiv"] <- as.character(filter_data$indiv[1])  # Assign individual
  rm(filter_data)
}

# Merge Hydraulic Data with Vessel Density data
HydraulicData <- merge(HydraulicData, VDensity_data[, c("Label", "VesselDensity", "X.Area")], by = "Label", all.x = TRUE)
colnames(HydraulicData)[6] <- "VesselFraction"

# Constants for hydraulic conductivity calculation
n <- 1.002e-9          # Viscosity of water in MPa·s
pw <- 998.2            # Density of water in kg/m³

# Calculate Kmax
HydraulicData$Kmax <- ((pi * pw) / (n * 128)) *      # Constants
  (HydraulicData$VesselDensity * 1e6) *             # Vessel density in vessels/m²
  ((HydraulicData$HydraulicDiameter * 1e-6)^4)      # Vessel diameter in meters (to the 4th power)


# Create and merge summarized data frames for median values
Median_data <- list(
  PitDiOp_data %>%
    group_by(indiv, ssp) %>%
    summarise(
      PitOpening = median(PitOpening, na.rm = TRUE),
      PitDiameter = median(PitDiameter, na.rm = TRUE),
      .groups = "drop"
    ),
  
  Wall_data %>%
    group_by(indiv, ssp) %>%
    summarise(
      Wthickness = median(WallThickness, na.rm = TRUE),
      .groups = "drop"
    ),
  VesselDiameter_data %>%
    group_by(indiv, ssp)  %>%
    summarise(
VesselDiameter = median(VesselDiameter, na.rm = TRUE),
      .groups = "drop"
    ),
  
  VesselDiameter_data %>%
    group_by(indiv, ssp) %>%
    mutate(threshold = quantile(VesselDiameter, 0.9, na.rm = TRUE)) %>%
    filter(VesselDiameter >= threshold) %>%
    summarise(
      TopVesselDiameter = median(VesselDiameter, na.rm = TRUE),
      .groups = "drop"
    ),
  
  PitFraction_data %>%
    group_by(indiv, ssp) %>%
    summarise(
      PitFraction = median(PitFraction, na.rm = TRUE),
      .groups = "drop"
    ),
  
  HydraulicData %>%
    group_by(indiv, ssp) %>%
    summarise(
      HydraulicDiameter = median(HydraulicDiameter, na.rm = TRUE),
      VesselDensity = median(VesselDensity, na.rm = TRUE),
      VesselFraction = median(VesselFraction, na.rm = TRUE),
      Kmax = median(Kmax, na.rm = TRUE),
      .groups = "drop"
    )
) %>% 
  reduce(full_join, by = c("ssp", "indiv"))




# Create and merge summarized data frames for mean values
Mean_data <- list(
  PitDiOp_data %>%
    group_by(indiv, ssp) %>%
    summarise(
      PitOpening = mean(PitOpening, na.rm = TRUE),
      PitDiameter = mean(PitDiameter, na.rm = TRUE),
      .groups = "drop"
    ),
  Wall_data %>%
    group_by(indiv, ssp) %>%
    summarise(Wthickness = mean(WallThickness, na.rm = TRUE), .groups = "drop"),
  VesselDiameter_data %>%
    group_by(indiv, ssp)  %>%
    summarise(
      VesselDiameter = mean(VesselDiameter, na.rm = TRUE),
      .groups = "drop"
    ),
  
  VesselDiameter_data %>%
    group_by(indiv, ssp) %>%
    mutate(threshold = quantile(VesselDiameter, 0.9, na.rm = TRUE)) %>%
    filter(VesselDiameter >= threshold) %>%
    summarise(
      TopVesselDiameter = mean(VesselDiameter, na.rm = TRUE),
      .groups = "drop"
    ),
  PitFraction_data %>%
    group_by(indiv, ssp) %>%
    summarise(PitFraction = mean(PitFraction, na.rm = TRUE), .groups = "drop"),
  HydraulicData %>%
    group_by(indiv, ssp) %>%
    summarise(
      HydraulicDiameter = mean(HydraulicDiameter, na.rm = TRUE),
      VesselDensity = mean(VesselDensity, na.rm = TRUE),
      VesselFraction = mean(VesselFraction, na.rm = TRUE),
      Kmax = mean(Kmax, na.rm = TRUE),
      .groups = "drop"
    )
) %>% 
  reduce(full_join, by = c("ssp", "indiv"))




# Add parasitism column to all data frames
dataframes <- Filter(function(x) is.data.frame(get(x)), ls())
for (df_name in dataframes) {
  df <- get(df_name)
  df <- df %>% mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", 
               "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", 
               "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  ))
  assign(df_name, df)
  rm(df)
}


# Save each data frame as a CSV
for (df_name in dataframes) {
  df <- get(df_name)

  fwrite(df, file = here("data", "processed", paste0(df_name, ".csv")))
  cat("Data frame", df_name, "saved successfully.\n")
}

# Clean up environment
rm(list = ls())

#######################################################
print("Data Wrangling complete")
