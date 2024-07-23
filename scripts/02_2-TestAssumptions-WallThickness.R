######################################################################
#
#Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
#13/07/2024
# Scrpt 02.2 - Test assumptions test
#################################################################
#load packages
source(here("scripts","01-DataWrangling.R"))



#load data
vdata <- read.csv(here("data","processed","vdata.csv"))
vadata <- read.csv(here("data","processed","vadata.csv"))
wdata <- read.csv(here("data","processed","wdata.csv"))
pitdata <- read.csv(here("data","processed","pitdata.csv"))

#wdata
###########

# Assuming wdata is your data frame and it contains columns 'ssp' and 'wthickness'
for (ssp in unique(wdata$ssp)) {
  # Subset data for the current species
  subset_data <- wdata[wdata$ssp == ssp, ]
  
  # Perform Shapiro-Wilk test for the current species
  shapiro_test <- shapiro.test(subset_data$wthickness)
  print(paste("Shapiro-Wilk test for species", ssp, ": p-value =", shapiro_test$p.value))
  
  # Plot the density of 'wthickness' for the current species
  p <- ggplot(subset_data, aes(x = wthickness)) +
    geom_density() +
    ggtitle(paste("Density plot for species", ssp)) +
    xlab("Wthickness") +
    ylab("Density")
  
  print(p)
}
#####################

# Function to check normality by removing one data point at a time
check_normality <- function(data) {
  # Perform initial Shapiro-Wilk test
  shapiro_test <- shapiro.test(data$wthickness)
  p_value <- shapiro_test$p.value
  
  # Continue removing one data point at a time until normality is reached or only one data point is left
  while (p_value < 0.05 && nrow(data) > 1) {
    # Find the data point whose removal improves normality the most
    best_p_value <- 0
    best_data <- data
    
    for (i in 1:nrow(data)) {
      temp_data <- data[-i, ]
      temp_shapiro_test <- shapiro.test(temp_data$wthickness)
      
      if (temp_shapiro_test$p.value > best_p_value) {
        best_p_value <- temp_shapiro_test$p.value
        best_data <- temp_data
      }
    }
    
    # Update data and p_value
    data <- best_data
    p_value <- best_p_value
  }
  
  return(data)
}

# Assuming wdata is your data frame and it contains columns 'ssp' and 'wthickness'
wdata_clean <- list()

for (ssp in unique(wdata$ssp)) {
  # Subset data for the current species
  subset_data <- wdata[wdata$ssp == ssp, ]
  
  # Check normality and clean data
  clean_data <- check_normality(subset_data)
  
  # Append cleaned data to the list
  wdata_clean[[ssp]] <- clean_data
  
  # Plot the density of 'wthickness' for the current species
  p <- ggplot(clean_data, aes(x = wthickness)) +
    geom_density() +
    ggtitle(paste("Density plot for species", ssp)) +
    xlab("Wthickness") +
    ylab("Density")
  
  print(p)
}

# Combine all cleaned data into a single dataframe
wdata_clean <- bind_rows(wdata_clean)

# View the cleaned dataframe
head(wdata_clean)

fwrite(wdata_clean,file=here("data","processed","wdata_clean"))
