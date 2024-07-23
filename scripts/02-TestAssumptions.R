######################################################################
#
#Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
#13/07/2024
# Scrpt 02 - Test assumptions test
#################################################################
#load packages
source(here("scripts","01-DataWrangling.R"))



#load data
vdata <- read.csv(here("data","processed","vdata.csv"))
vadata <- read.csv(here("data","processed","vadata.csv"))
wdata <- read.csv(here("data","processed","wdata.csv"))
pitdata <- read.csv(here("data","processed","pitdata.csv"))

################################################################################


##############################################################################
#Pits analysis
################################################################################

# Create an empty dataframe to store results
pitresults <- data.frame(
  species = character(),
  pitavg_p_value = numeric(),
  pcavg_p_value = numeric(),
  peavg_p_value = numeric(),
  pcdvg_p_value = numeric(),
  stringsAsFactors = FALSE
)

# Assuming `species` is a vector containing unique species names from pitdata$ssp
species <- unique(pitdata$ssp)

#normality test
#######################################################
# Loop through each species and perform the Shapiro-Wilk test for each column
for (ssp in species) {
  # Subset data for the current species
  subset_data <- pitdata[pitdata$ssp == ssp, ]
  
  # Perform Shapiro-Wilk test for each column
  shapiro_results <- sapply(subset_data[, c("pitavg", "pcavg", "peavg","pcd")], function(x) {
    shapiro_test <- shapiro.test(x)
    shapiro_test$p.value
  })
  
  # Append the results to the dataframe
  pitresults <- rbind(pitresults, data.frame(
    species = ssp,
    pitavg_p_value = shapiro_results["pitavg"],
    pcavg_p_value = shapiro_results["pcavg"],
    peavg_p_value = shapiro_results["peavg"],
    pcd_p_value = shapiro_results["pcd"]
  ))
}

# Print the final results dataframe
print(pitresults)


####removing outliers
#################################################################


####pitavg
outlier<- boxplot.stats(pitdata$pitavg[pitdata$ssp=="Viscum album"])$out #Find outliers
x1 <- outlier[which.max(abs(outlier - median(pitdata$pitavg[pitdata$ssp=="Viscum album"],na.rm = T)))] #find most extreme outlier
pitdata$pitavg[which(pitdata$pitavg == x1)] <- NA
shapiro.test(pitdata$pitavg[pitdata$ssp == "Viscum album"])###re test
plot(density(pitdata$pitavg[pitdata$ssp == "Viscum album"], na.rm = TRUE))

###remove another 
outlier<- boxplot.stats(pitdata$pitavg[pitdata$ssp=="Viscum album"])$out #Find outliers
x1 <- outlier[which.max(abs(outlier - median(pitdata$pitavg[pitdata$ssp=="Viscum album"],na.rm = T)))] #find most extreme outlier
pitdata$pitavg[which(pitdata$pitavg == x1)] <- NA
shapiro.test(pitdata$pitavg[pitdata$ssp == "Viscum album"])###re test
plot(density(pitdata$pitavg[pitdata$ssp == "Viscum album"], na.rm = TRUE))

###remove another 
outlier<- boxplot.stats(pitdata$pitavg[pitdata$ssp=="Viscum album"])$out #Find outliers
x1 <- outlier[which.max(abs(outlier - median(pitdata$pitavg[pitdata$ssp=="Viscum album"],na.rm = T)))] #find most extreme outlier
pitdata$pitavg[which(pitdata$pitavg == x1)] <- NA
shapiro.test(pitdata$pitavg[pitdata$ssp == "Viscum album"])###re test
plot(density(pitdata$pitavg[pitdata$ssp == "Viscum album"], na.rm = TRUE))

###normality achieved with 3 points removed

##############################################################
#pcavg
plot(density(pitdata$pcavg[pitdata$ssp == "Viscum album"], na.rm = TRUE))

outlier<- boxplot.stats(pitdata$pcavg[pitdata$ssp=="Viscum album"])$out #Find outliers
x1 <- outlier[which.max(abs(outlier - median(pitdata$pcavg[pitdata$ssp=="Viscum album"],na.rm = T)))] #find most extreme outlier
pitdata$pcavg[which(pitdata$pcavg == x1)] <- NA
shapiro.test(pitdata$pcavg[pitdata$ssp == "Viscum album"])###re test

outlier<- boxplot.stats(pitdata$pcavg[pitdata$ssp=="Viscum album"])$out #Find outliers
x1 <- outlier[which.max(abs(outlier - median(pitdata$pcavg[pitdata$ssp=="Viscum album"],na.rm = T)))] #find most extreme outlier
pitdata$pcavg[which(pitdata$pcavg == x1)] <- NA
shapiro.test(pitdata$pcavg[pitdata$ssp == "Viscum album"])###re test

outlier<- boxplot.stats(pitdata$pcavg[pitdata$ssp=="Viscum album"])$out #Find outliers
x1 <- outlier[which.max(abs(outlier - median(pitdata$pcavg[pitdata$ssp=="Viscum album"],na.rm = T)))] #find most extreme outlier
pitdata$pcavg[which(pitdata$pcavg == x1)] <- NA
shapiro.test(pitdata$pcavg[pitdata$ssp == "Viscum album"])###re test

###normality achieved with 3 points removed

####################################################################################
#peavg

plot(density(pitdata$peavg[pitdata$ssp == "Vochysia thyrsoidea"], na.rm = TRUE))

outlier<- boxplot.stats(pitdata$peavg[pitdata$ssp=="Vochysia thyrsoidea"])$out #Find outliers
x1 <- outlier[which.max(abs(outlier - median(pitdata$peavg[pitdata$ssp=="Vochysia thyrsoidea"],na.rm = T)))] #find most extreme outlier
pitdata$peavg[which(pitdata$peavg == x1)] <- NA
shapiro.test(pitdata$peavg[pitdata$ssp == "Vochysia thyrsoidea"])###re test

######normality achievide with 1 point removed



############################################
#pcd

plot(density(pitdata$pcd[pitdata$ssp == "Populus nigra"], na.rm = TRUE))

outlier<- boxplot.stats(pitdata$pcd[pitdata$ssp=="Populus nigra"])$out #Find outliers
x1 <- outlier[which.max(abs(outlier - median(pitdata$pcd[pitdata$ssp=="Populus nigra"],na.rm = T)))] #find most extreme outlier
pitdata$pcd[which(pitdata$pcd == x1)] <- NA
shapiro.test(pitdata$pcd[pitdata$ssp == "Populus nigra"])###re test

#homogeneity of variance test
###############################################

ggplot(pitdata, aes(x = ssp, y = pitavg, fill = ssp)) +
  geom_boxplot() +
  labs(title = "Boxplot of pitavg Across Species",
       x = "Species", y = "pitavg") +
  theme_minimal()


###Graph shows heterocedasticity between groups

ggplot(pitdata, aes(x = ssp, y = peavg, fill = ssp)) +
  geom_boxplot() +
  labs(title = "Boxplot of peavg Across Species",
       x = "Species", y = "peavg") +
  theme_minimal()

ggplot(pitdata, aes(x = ssp, y = pcavg, fill = ssp)) +
  geom_boxplot() +
  labs(title = "Boxplot of pcavg Across Species",
       x = "Species", y = "pcavg") +
  theme_minimal()


ggplot(pitdata, aes(x = ssp, y = pcd, fill = ssp)) +
  geom_boxplot() +
  labs(title = "Boxplot of pcd Across Species",
       x = "Species", y = "pcd") +
  theme_minimal()

###heterocedasticity present in all pit traits



fwrite(pitdata,file=here("data","processed","pitdata_clean.csv"))





