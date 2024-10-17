######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.2 - Data Analysis - Vessels LME
#################################################################
library(here)
source(here("scripts", "02_2-TestAssumptions-Vessels.R"))
rm(list = ls())

source(here("scripts", "Functions.R"))
#load data
HydraulicData <- read.csv(here("data", "processed", "HydraulicData.csv"))
vdata <- read.csv(here("data", "processed", "vdata.csv"))
vadata <- read.csv(here("data", "processed", "vadata.csv"))

 install.packages("lqmm")
library(lqmm)
 