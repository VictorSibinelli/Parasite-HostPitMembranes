######################################################################
#
#Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
#13/07/2024
# Scrpt 01 - Data wrangling
#################################################################


#loading packages used uin the analysis (you can also run the 00-library script manualy)
source("scripts/00-library.R")
rm(list=ls())

#reading data
vdata <- read.csv(here("data","raw","VesselsDiameter.csv"))
vadata <- read.csv(here("data","raw","VesselArea.csv"), sep=";")
wdata <- read.csv(here("data","raw","2xWallThickness.csv"), sep=";")
pitdata <- read.csv(here("data","raw","Pits.csv"),sep=";")


# List of data frame names
dataframes <- ls()

# Relevel the factors for each data frame
for (df_name in dataframes) {
  df <- get(df_name)  # Get the data frame by name
  
  if ("ssp" %in% colnames(df)) {  # Check if 'ssp' column exists
    df$ssp <- factor(df$ssp, levels = c("Psittacanthus robustus", "Vochysia thyrsoidea",
                                        "Phoradendeon perrotettii", "Tapirira guianensis",
                                        "Struthanthus rhynchophyllus", "Tipuana tipu",
                                        "Viscum album", "Populus nigra"))
    assign(df_name, df)  # Assign the modified data frame back to its name
  }
}
 rm(df) #remove duplicated dataframe
 
 
#### modifiying dataframes
 
 
##vdata
 
# calculated equivalent area circle diameter based on vessel area
vdata$VD <- 2*(sqrt(vdata$Area/pi)) 
#correcting indiv
vdata$indiv <- paste(vdata$ssp,vdata$indiv, sep=" ")

##vadata

#calculated vessel density based on vessel count,total area and vessel area fraction in mm2
vadata$vdensity <- vadata$Count/(vadata$Total.Area/vadata$X.Area)*10000
#corrcting missing name
vadata$ssp[33:63] <- "Psittacanthus robustus"


#calculate single wall thickness
wdata$wthickness <- wdata$Length/2


####saving dataframes 


# Loop over each data frame name
for (df_name in dataframes) {
  df <- get(df_name)  # Get the data frame by name
  
  # Construct file path
  dfile <- paste0(df_name, ".csv")
  file_path <- here("data", "processed", dfile)
  
  # Save data frame as CSV using fwrite
  fwrite(df, file = file_path)
  
  cat("Data frame", df_name, "saved successfully to", file_path, "\n")
}
rm(df)
