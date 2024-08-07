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
wdata <- read.csv(here("data", "processed", "wdata.csv"))
pitdata <- read.csv(here("data", "processed", "pitdata.csv"))



