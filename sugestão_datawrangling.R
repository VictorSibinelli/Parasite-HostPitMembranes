# 01-DataWrangling.R SUGESTÂO
# Linha 78 em diante alterada, para tornar o loop mais eficiente

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
vdata <- read.csv(here("data", "raw", "VesselsDiameter.csv"), sep = ",")
vadata <- read.csv(here("data", "raw", "VesselCount.csv"), sep = ",")
wdata <- read.csv(here("data", "raw", "2xWallThickness.csv"), sep = ";")
pitdata <- read.csv(here("data", "raw", "Pits.csv"), sep = ";")[,1:8]
pitOdata <- read.csv(here("data", "raw", "PitOpeningData.csv"), sep = ";")


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


#### modifying data frames

## vdata

# calculated equivalent area circle diameter based on vessel area
vdata$ssp[is.na(vdata$ssp)] <- "Phoradendron perrotettii"
vdata$VesselDiameter <- 2 * (sqrt(vdata$Area/pi))

# correcting indiv
vdata$indiv <- paste(vdata$ssp, vdata$indiv, sep = " ")

## vadata

# calculated vessel density based on vessel count,total area and vessel area fraction in mm2
vadata$vdensity <- vadata$Count / (vadata$Total.Area / vadata$X.Area) * 10000

# Inserting missing name
vadata$ssp[1:32] <- "Phoradendron perrotettii"
vadata$ssp[33:63] <- "Psittacanthus robustus"

## wdata
# calculate single wall thickness
wdata$wthickness <- wdata$Length / 2

# pitdata

pitdata$peavg <- rowMeans(pitdata[, 5:6], na.rm = T) # calculating pit membrane thickness average at the edges
pitdata$pcavg <- rowMeans(pitdata[, 2:4], na.rm = T) # calculating pit membrane thickness average at the center
pitdata$pitavg <- rowMeans(pitdata[, 2:6], na.rm = T) # calculating pit membrane thickness average

#pitdata
pitOdata$PitDiameter <- read.csv(here("data", "raw", "PitDiameter.csv"), sep = ";")$PitDiameter/10
colnames(pitOdata)[4] <- "PitOpening"


#### saving dataframes
# List of all objects in the environment
dataframes <- Filter(function(x) is.data.frame(get(x)), ls())

# Loop over each object
for (df_name in dataframes) {
  df <- get(df_name) # Get the object by name
  
  dfile <- paste0(df_name, ".csv")
  file_path <- here("data", "processed", dfile)
  
  # Save data frame as CSV using fwrite
  fwrite(df, file = file_path)
  
  cat("Data frame", df_name, "saved successfully to", file_path, "\n")
  
}

# Remove all objects from the environment
rm(list = ls())



#######################################################
