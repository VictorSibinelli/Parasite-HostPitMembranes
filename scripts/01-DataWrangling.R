######################################################################
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 01 - Data Wrangling
######################################################################

# Load required libraries and custom functions
library(here)
source(here("scripts", "00-library.R"))
rm(list = ls())  # Clear environment

# Read in raw data files
vdata <- read.csv(here("data", "raw", "VesselsDiameter.csv"))
vadata <- read.csv(here("data", "raw", "VesselDensity.csv"))
wdata <- read.csv(here("data", "raw", "2xWallThickness.csv"), sep = ";")
pitdata <- read.csv(here("data", "raw", "Pits.csv"), sep = ";")[, 1:8]
pitOdata <- read.csv(here("data", "raw", "PitOpeningData.csv"), sep = ";")

#### Modifying Data Frames

## vdata - Vessel Diameter Data
# Correct species names and reformatting
vdata$indiv <- paste(vdata$ssp, vdata$indiv, sep = " ")
colnames(vdata)[3] <- "pic"
vdata[739:747, "pic"] <- "Phoradendron perrotettii 2 - GC4841 x10 005-1.tif"

vdata$pic <- gsub("(^[a-z])", "\\U\\1", vdata$pic, perl = TRUE)  # Capitalize first letter
vdata <- within(vdata, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
  pic <- gsub("Phoradendeon", "Phoradendron", pic)
})
vdata$VesselDiameter <- 2 * sqrt(vdata$Area / pi)  # Calculate vessel diameter

## vadata - Vessel Density Data
# Correct species names and calculate vessel density
vadata$vdensity <- vadata$Count / (vadata$Total.Area / vadata$X.Area) * 10000
colnames(vadata)[2] <- "indiv"
vadata <- within(vadata, {
  ssp <- gsub("Phoradendeon", "Phoradendron", ssp)
  indiv <- gsub("Phoradendeon", "Phoradendron", indiv)
  pic <- gsub("Phoradendeon", "Phoradendron", pic)
})


## wdata - Wall Thickness Data
# Calculate single wall thickness
wdata$wthickness <- wdata$Length / 2

## pitdata - Pit Data
# Calculate average pit membrane thickness at edges and center
pitdata$peavg <- rowMeans(pitdata[, 5:6], na.rm = TRUE)
pitdata$pcavg <- rowMeans(pitdata[, 2:4], na.rm = TRUE)
pitdata$pitavg <- rowMeans(pitdata[, 2:6], na.rm = TRUE)

## pitOdata - Pit Opening Data
# Adjust pit opening diameter
pitOdata$PitDiameter <- read.csv(here("data", "raw", "PitDiameter.csv"), sep = ";")$PitDiameter / 10
colnames(pitOdata)[4] <- "PitOpening"

#### Save Processed Data Frames
# Filter only data frames from the environment
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
