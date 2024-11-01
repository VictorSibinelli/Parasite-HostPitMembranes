######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 02.3 - Test assumptions test - Vessels
#####################################################################

# Load packages
library(here)
source(here("scripts", "functions.R"))
# Load data
vdata <- read.csv(here("data", "processed", "vdata.csv"))
vadata <- read.csv(here("data", "processed", "vadata.csv"))

output_dir_figs <- here("outputs", "figs", "assumptions")

# Create a data frame for Hydraulic data
HydraulicData <- data.frame(
  ssp = NA,
  indiv = NA,
  pic = unique(vdata$pic),
  HydraulicDiameter = NA,
  stringsAsFactors = FALSE
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

# Merge with vessel density data
HydraulicData <- merge(HydraulicData, vadata[, c("pic", "vdensity","X.Area")], by = "pic", all.x = TRUE)

# Update to include the `parasitism` variable
HydraulicData <- HydraulicData %>%
  mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", 
               "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", 
               "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  )) %>%
  drop_na()
colnames(HydraulicData)[6] <- "VesselFraction"
# Constants for hydraulic conductivity calculation
n <- 1.002 * 10^-9  # Viscosity index of water (MPa·s at 20°C)
pw <- 998.2         # Density of water at 20°C (kg/m³)

# Max hydraulic conductivity calculation
HydraulicData$Kmax <- ((pi * pw) / (n * 128)) *  # Constants
  (HydraulicData$vdensity * 1e+6) *               # Vessel density in vessels/m² (converted from vessels/mm²)
  ((HydraulicData$HydraulicDiameter * 1e-6)^4)    # Vessel diameter in meters to the power of 4



####relevel factors
relevel_factors(ls())
head(HydraulicData)

# Set up plotting for Vessel Diameter Density
png(filename = file.path(output_dir_figs, "VesselDiameterDensity.png"), 
    width = 6400, height = 4800, res = 600)  # Increase image size for bigger graphs
par(mfrow = c(2, 4), oma = c(4, 4, 4, 2), mar = c(2, 2, 2, 2))

for (i in unique(vdata$ssp)) {
 p <-  plot(density(vdata$VesselDiameter[vdata$ssp == i]),
       main = paste(i),
       xlab = "Vessel Diameter", 
       ylab = "Density")
  mtext("Vessel Diameter", side = 3, line = 0, outer = TRUE)
  p
  print(p)
}
dev.off()

# Set up plotting for Vessel Density Density
png(filename = file.path(output_dir_figs, "VesselDensityDensity.png"), 
    width = 6400, height = 4800, res = 600)  # Increase image size for bigger graphs
par(mfrow = c(2, 4), oma = c(4, 4, 4, 2), mar = c(2, 2, 2, 2))

for (i in unique(vadata$ssp)) {
 p <-  plot(density(vadata$vdensity[vadata$ssp == i]),
       main = paste("Density of Vessel Density for", i),
       xlab = "Vessel Density", 
       ylab = "Density")
  mtext("Vessel Density", side = 3, line = 0, outer = TRUE)
  print(p)
}
dev.off()

# Set up plotting for Kmax Density
png(filename = file.path(output_dir_figs, "KmaxDensity.png"), 
    width = 6400, height = 4800, res = 600)  # Increase image size for bigger graphs
par(mfrow = c(2, 4), oma = c(4, 4, 4, 2), mar = c(2, 2, 2, 2))

for (i in unique(HydraulicData$ssp)) {
  p <- plot(density(HydraulicData$Kmax[HydraulicData$ssp == i]),
       main = paste("Kmax", i),
       xlab = "Vessel Diameter", 
       ylab = "Density")
  mtext("Kmax", side = 3, line = 0, outer = TRUE)
  print(p)
}
dev.off()

# Boxplots for HydraulicDiameter, vdensity, and Kmax
h <- HydraulicData %>% ggplot(aes(x = ssp, y = HydraulicDiameter)) +
  geom_boxplot()
ggsave(here(output_dir_figs, "VesselDiameterVariance.png"), plot = h, dpi = 600, width = 10, height = 7)
print(h)

h <- HydraulicData %>% ggplot(aes(x = ssp, y = vdensity)) +
  geom_boxplot()
ggsave(here(output_dir_figs, "VesselDensityVariance.png"), plot = h, dpi = 600, width = 10, height = 7)
print(h)

h <- HydraulicData %>% ggplot(aes(x = ssp, y = Kmax)) +
  geom_boxplot()
ggsave(here(output_dir_figs, "KmaxVariance.png"), plot = h, dpi = 600, width = 10, height = 7)
print(h)

h <- HydraulicData %>% ggplot(aes(x = ssp, y = X.Area)) +
  geom_boxplot()
ggsave(here(output_dir_figs, "VesselFractionArea.png"), plot = h, dpi = 600, width = 10, height = 7)
print(h)

# Non-normality and no homoscedasticity for all traits
write.csv(HydraulicData, file = file.path(here("data", "processed", "HydraulicData.csv")), row.names = FALSE)


cat("Vessel Traits Test assumptions Complete 
    \nSummary:
    \n1. Great deviation from normality in Vessel diameter
    \n2. Small deviation from normality in Vessel Density and Kmax
    \n3.Major heteroscedasticity observed in all vessel traits")
