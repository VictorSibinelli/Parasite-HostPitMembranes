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


# Create a data frame for Hydraulic data
HydraulicData <- data.frame(
  ssp = NA,
  indiv = NA,
  pic = unique(vdata$pic),
  HydraulicDiameter = NA,
  stringsAsFactors = F
)

# Ensure the 'ssp' column is a factor and use its value correctly
for (i in seq_along(unique(vdata$pic))) {
  label_value <- unique(vdata$pic)[i]
  
  # Filter data for the current 'pic'
  filter_data <- vdata[vdata$pic == label_value, ]
  
  # Compute Hydraulic Diameter
  HydraulicData[i, "HydraulicDiameter"] <- (sum(filter_data$VesselDiameter^4) / 
                                              length(filter_data$VesselDiameter))^(1/4)
  
  # Explicitly assign 'pic', 'ssp', 'indiv' values
  HydraulicData[i, "pic"] <- as.character(filter_data$pic[1])       # Convert to character if needed
  HydraulicData[i, "ssp"] <- as.character(filter_data$ssp[1])       # Convert to character to avoid factor level issue
  HydraulicData[i, "indiv"] <- as.character(filter_data$indiv[1])   # Convert to character if needed
  rm(filter_data)
}
# Update  to include the `parasitism` variable
HydraulicData <- HydraulicData %>%
  mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", 
               "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", 
               "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  )) %>% drop_na()

HydraulicData <- merge(HydraulicData, vadata[, c("pic", "vdensity")], by = "pic", all.x = TRUE)

# η is the viscosity index of water (1.002 × 10−9 MPa s at 20°C)
n <- 1.002 * 10^-9

#ρw is the density of water at 20°C (998.2*kg*m−3 at 20°C)
pw <- 998.2

#Max hydraulic conductivity
HydraulicData$Kmax <- ((pi*pw)/(n*128)) * #constants
  (HydraulicData$vdensity*1e+6) * #Vessel density in vessels/ meter2 (converted from vessel/mm2)
  ((HydraulicData$HydraulicDiameter*10e-6)^4) #Vessel diameter in meters to the power of 4




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



for (i in unique(vdata$ssp)) {
  plot(density(vdata$VesselDiameter[vdata$ssp == i]),
       main = paste("Density of VesselDiameter for", i),
       xlab = "Vessel Diameter", ylab = "Density")
}

for (i in unique(vadata$ssp)) {
  plot(density(vadata$vdensity[vadata$ssp == i]),
       main = paste("Density of Vessel Density for", i),
       xlab = "Vessel Diameter", ylab = "Density")
}

for (i in unique(HydraulicData$ssp)) {
  plot(density(HydraulicData$Kmax[HydraulicData$ssp == i]),
       main = paste("Kmax", i),
       xlab = "Vessel Diameter", ylab = "Density")
}

HydraulicData %>% ggplot(aes(x=ssp,y=HydraulicDiameter))+
  geom_boxplot()


HydraulicData %>% ggplot(aes(x=ssp,y=vdensity))+
  geom_boxplot()

HydraulicData %>% ggplot(aes(x=ssp,y=Kmax))+
  geom_boxplot()
#Non normality and no homocedasticity for all traits
