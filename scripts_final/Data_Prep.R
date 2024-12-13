######################################################################
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 01 - Data Wrangling
######################################################################

# Load required libraries and custom functions
library(here)
source(here("scripts", "00-library.R"))
rm(list = ls())  # Clear environment

Wall_data <- read.csv(here("data", "raw", "VesselWall.csv"))
colnames(Wall_data)[3] <- "WallThickness"
Wall_data$WallThickness <- ifelse(Wall_data$ssp != "Populus nigra", 
                                  Wall_data$WallThickness * 1.43, 
                                  Wall_data$WallThickness)
Wall_data$WallThickness <- Wall_data$WallThickness/2
Wall_data <- within(Wall_data, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
})


VDensity_data <- read.csv(here("data", "raw", "VesselCount.csv"))[-9]
VDensity_data$VesselDendity <- (VDensity_data$Count / (VDensity_data$Total.Area / VDensity_data$X.Area)) * 10000
VDensity_data <- within(VDensity_data, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
  Label <- gsub("Phoradendeon", "Phoradendron", Label)
})




VesselDiameter_data <- read.csv(here("data", "raw", "VesselDiameter.csv"))

VesselDiameter_data <- within(VesselDiameter_data, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
  Label <- gsub("Phoradendeon", "Phoradendron", Label)
})
VesselDiameter_data$VesselDiameter <- sqrt(VesselDiameter_data$Area*VesselDiameter_data$Circ.)*2


PitFraction_data <- read.csv(here("data", "raw", "PitFraction.csv"))[-c(3:5)]
PitFraction_data <- within(PitFraction_data, {
  ssp <- gsub("perrottetti", "perrotettii", ssp)
  indiv <- gsub("perrottetti", "perrotettii", indiv)
})



PitOp <- read.csv(here("data", "raw", "PitOpening.csv"))[-3]
PitOp <- within(PitOp, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
})
PitOp <- PitOp %>%  group_by(ssp, indiv) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

PitDi <- read.csv(here("data","raw","PitDiameter.csv"))[-3]
PitDi <- within(PitDi, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
})

PitDi <- PitDi %>% group_by(ssp, indiv) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

PitDiOp_data <- merge(
  PitDi,
  PitOp,
  by = c("ssp", "indiv", "row_id"),all=T)[-3]

colnames(PitDiOp_data)[3:4] <- c("PitDiameter","PitOpening")

PitDiOp_data$PitOpening <- ifelse(PitDiOp_data$ssp != "Populus nigra", 
                                  PitDiOp_data$PitOpening * 1.43, 
                                  PitDiOp_data$PitOpening)

rm(PitDi, PitOp)




PitMembrane_data <- read.csv(here("data", "raw", "PitMembrane.csv"))[-8]
PitMembrane_data$Tpm <- rowMeans(PitMembrane_data[,2:6],na.rm = T)*1000



dataframes <- Filter(function(x) is.data.frame(get(x)), ls())



for (df_name in dataframes) {
  df <- get(df_name)  # Get the data frame by its name
  
  # Modify the data frame with parasitism column
  df <- df %>% mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", 
               "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", 
               "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  ))
  
  assign(df_name, df)  # Save the modified data frame back to the environment
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
