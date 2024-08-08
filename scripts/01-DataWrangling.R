######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 01 - Data wrangling
#################################################################


# loading packages used run the analysis (you can also run the 00-library script manually)
library(here)
source(here("scripts", "00-library.R"))
rm(list = ls())

# reading data
vdata <- read.csv(here("data", "raw", "VesselsDiameter.csv"), sep = ";")
vadata <- read.csv(here("data", "raw", "VesselArea2.csv"), sep = ",")
wdata <- read.csv(here("data", "raw", "2xWallThickness.csv"), sep = ";")
pitdata <- read.csv(here("data", "raw", "Pits.csv"), sep = ";")
pitOdata <- read.csv(here("data", "raw", "PitOpeningData.csv"), sep = ",")

#### modifying data frames


## vdata

# calculated equivalent area circle diameter based on vessel area
vdata$VD <- 2 * (sqrt(vdata$Area / pi))
# correcting indiv
vdata$label <- paste(vdata$ssp, vdata$label, sep = " ")

## vadata

# calculated vessel density based on vessel count,total area and vessel area fraction in mm2
vadata$vdensity <- vadata$Count / (vadata$Total.Area / vadata$X.Area) * 10000
# Inserting missing name
vadata$ssp[33:63] <- "Psittacanthus robustus"

## wdata
# calculate single wall thickness
wdata$wthickness <- wdata$Length / 2

# pitdata

pitdata$peavg <- rowMeans(pitdata[, 5:6], na.rm = T) # calculating pit membrane thickness average at the edges
pitdata$pcavg <- rowMeans(pitdata[, 2:4], na.rm = T) # calculating pit membrane thickness average at the center
pitdata$pitavg <- rowMeans(pitdata[, 2:6], na.rm = T) # calculating pit membrane thickness average


#### saving dataframes

# List of data frame names
dataframes <- ls()

# Loop over each data frame name
for (df_name in dataframes) {
  df <- get(df_name) # Get the data frame by name

  # Construct file path
  dfile <- paste0(df_name, ".csv")
  file_path <- here("data", "processed", dfile)

  # Save data frame as CSV using fwrite
  fwrite(df, file = file_path)

  cat("Data frame", df_name, "saved successfully to", file_path, "\n")
}
rm(list = ls())

#######################################################
