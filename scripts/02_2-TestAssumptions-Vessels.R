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
Hydraulic_Diameter <- data.frame(
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
  
  Hydraulic_Diameter[i, "HydraulicDiameter"] <- (sum(filter_data$VesselDiameter^4) / 
                                                   length(filter_data$VesselDiameter))^(1/4)
  Hydraulic_Diameter[i, c("pic", "ssp", "Indiv")] <- filter_data[1, c("Label", "ssp", "indiv")]
}

Hydraulic_Diameter


df_merged <- merge(Hydraulic_Diameter, vadata[, c("pic", "vdensity")], by = "pic", all.x = TRUE)

