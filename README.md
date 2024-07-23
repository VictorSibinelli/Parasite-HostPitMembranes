# Overview

This repository contains all data, scripts, and outputs (including Supplementary data) used in the paper "XXX". *Still in progress*.

All data were collected from mature wood samples, either with light or electronic microscopy techniques (for more details, check the Materials and Methods sections of the article at "xxxxx").

## Repository Structure

### data
Contains all data needed to run the analysis in `.csv` and `.txt` files with metadata. All primary measurements were acquired with ImageJ software.

#### raw
Contains original data for the analysis:
- `2xWallThickness`: Intervascular double wall thickness
- `VesselsDiameter`: Individual vessel data
- `VesselArea`: Vessel frequency and fraction data
- `Pits`: Intervascular wall pits linear data
- `PitsArea`: Intervascular wall pit area data

#### processed
Contains raw `.csv` files after data wrangling.

### scripts
Contains scripts used to run the analysis:
- `00-library`: Checks if required packages are installed (and installs if they are not) and loads them. Additionally, manages package version controls (if desired) using the package groundhog.
- `01-data_wrangling`: Data manipulation.
- `02-testAssumption`: Performs the statistical analysis used in the paper.
- `03-Analysis`: Builds the figures presented in the article.
- `04-Graphics`: Builds the figures presented in the article.

### outputs
Contains the outputs from the analysis and graphics scripts:
- `figs`: Contains the figures.
- `tables`: Contains the tables.


