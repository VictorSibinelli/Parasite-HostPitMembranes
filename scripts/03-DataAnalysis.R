######################################################################
#
#Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
#13/07/2024
# Scrpt 03 - Data Analysis
#################################################################

#load data
vdata <- read.csv(here("data","processed","vdata.csv"))
vadata <- read.csv(here("data","processed","vadata.csv"))
wdata <- read.csv(here("data","processed","wdata.csv"))
pitdata <- read.csv(here("data","processed","pitdata.csv"))


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
  rm(df,df_name) #remove duplicated dataframe
}



##sub seting species
psittacanthus_data <-pitdata[pitdata$ssp == "Psittacanthus robustus", ]
vochysia_data <-pitdata[pitdata$ssp == "Vochysia thyrsoidea", ]
phoradendron_data <-pitdata[pitdata$ssp == "Phoradendron perrotettii", ]
tapirira_data <-pitdata[pitdata$ssp == "Tapirira guianensis", ]
struthanthus_data <-pitdata[pitdata$ssp == "Struthanthus rhynchophyllus", ]
tipuana_data <-pitdata[pitdata$ssp == "Tipuana tipu", ]
viscum_data <-pitdata[pitdata$ssp == "Viscum album", ]
populus_data <-pitdata[pitdata$ssp == "Populus nigra", ]