#Overview

This repository contains all data, scripts and outputs (including Supplementary data) used in the paper "XXX". *Still in progress*


All data were collected from mature wood samples, either with light or electronic miscroscopy techniques (for more details, check Materials and Methods sections of the article at "xxxxx".

##Repository Structure

*data: contains all data collect needed to run the analysis in csv. and a .txt files with metadata. All primary measurements were acquired with ImageJ software
	* raw: contains original data for the analysis 
		*2xWallThickness: Intervascular double wall thickness
		*VesselsDiameter: Individual vessel data
		*VesselArea: Vessel frequency and fraction data
		*Pits: Intervascular wall pits linear data
		*PitsArea: Intervascular wall pit área data

	* processed: contains raw .csv after data wrangling 

*scripts: contains scripts used to run the analysis
	*00-library: checks if required packages are installed (and installs if they are not) and loads them. Additionally, menages packages version controls (if desired) using the package groundhog.
	*01-data_wrangling: data manipulation
	*02-analysis: performs tthe statistical analysis used in the paper
	*03-graphics: builds the figures presente int he article

*outputs: Contains the outputs from the analysis and graphics scripts
	*figs: contains the figures
	*tables: contains the tables


