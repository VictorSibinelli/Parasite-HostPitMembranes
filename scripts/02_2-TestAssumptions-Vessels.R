######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Scrpt 02.3 - Test assumptions test - Vessels
#################################################################
# load packages
library(here)
source(here("scripts", "01-DataWrangling.R"))


# load data
vdata <- read.csv(here("data", "processed", "vdata.csv"))
vadata <- read.csv(here("data", "processed", "vadata.csv"))


for (i in unique(vdata$ssp)) {
  plot(density(vdata$VesselDiameter[vdata$ssp == i]),
       main = paste("Density of VesselDiameter for", i),
       xlab = "Vessel Diameter", ylab = "Density")
}

# Create a data frame instead of a matrix
HydraulicData <- data.frame(
  pic = unique(vdata$Label),
  ssp= NA,
  Indiv = NA,
  HydraulicDiameter = NA,
  stringsAsFactors = FALSE
)

# Calculate Hydraulic Diameter
for (i in seq_along(unique(vdata$Label))) {
  label_value <- unique(vdata$Label)[i]
  filter_data <- vdata %>%
    filter(Label == label_value)
  
  HydraulicData[i, "HydraulicDiameter"] <- (sum(filter_data$VesselDiameter^4) / 
                                                   length(filter_data$VesselDiameter))^(1/4)
  HydraulicData[i, c("pic", "ssp", "Indiv")] <- filter_data[1, c("Label", "ssp", "indiv")]
}

HydraulicData <- merge(HydraulicData, vadata[, c("pic", "vdensity")], by = "pic", all.x = TRUE)

# η is the viscosity index of water (1.002 × 10−9 MPa s at 20°C)
n <- 1.002 * 10^-9

#ρw is the density of water at 20°C (998.2*kg*m−3 at 20°C)
pw <- 998.2

#Max hydraulic conductivity
HydraulicData$Kmax <- ((pi*pw)/(n*128)) * #constants
  (HydraulicData$vdensity*1e+6) * #Vessel density in vessels/ meter2 (converted from vessel/mm2)
  ((HydraulicData$HydraulicDiameter*10e-6)^4) #Vessel diameter in meters to the power of 4
