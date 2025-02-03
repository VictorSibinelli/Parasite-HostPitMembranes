# Load required libraries
library(ggplot2)
library(dplyr)
library(here)
library(factoextra)

# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data <- read.csv(here("data", "processed", "VesselDiameter_data.csv")) 
Hydraulic_data <- read.csv(here("data", "processed", "HydraulicData.csv"))
PitFraction_data <- read.csv(here("data", "processed", "PitFraction_data.csv"))
PitDiOp_data <- read.csv(here("data", "processed", "PitDiOp_data.csv"))

source(here("scripts", "Functions.R"))

relevel_factors(ls())

# List of individual pairs for testing
indiv_test <- list(
  c("Psittacanthus robustus 1", "Vochysia thyrsoidea 1"),
  c("Psittacanthus robustus 2", "Vochysia thyrsoidea 2"),
  c("Psittacanthus robustus 3", "Vochysia thyrsoidea 3"),
  c("Phoradendron perrotettii 1", "Tapirira guianensis 3"),
  c("Phoradendron perrotettii 2", "Tapirira guianensis 2"),
  c("Phoradendron perrotettii 3", "Tapirira guianensis 3"),
  c("Struthanthus rhynchophyllus 1", "Tipuana tipu 1"),
  c("Struthanthus rhynchophyllus 2", "Tipuana tipu 2"),
  c("Struthanthus rhynchophyllus 3", "Tipuana tipu 3"),
  c("Viscum album 1", "Populus nigra 1"),
  c("Viscum album 2", "Populus nigra 2"),
  c("Viscum album 3", "Populus nigra 3")
)


# Running tests for all traits
Wall_ind_result <- perform_t_tests(Wall_data, "WallThickness", indiv_test)
Wall_ind_result$trait <- "Tw"

VesselD_ind_result <- perform_t_tests(VesselDiameter_data, "VesselDiameter", indiv_test)
VesselD_ind_result$trait <- "D"

VesselDTop_ind <- perform_t_tests(
  VesselDiameter_data %>%
    group_by(ssp, indiv) %>%
    filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
    ungroup(),
  "VesselDiameter",
  indiv_test
)
VesselDTop_ind$trait <- "Dtop"

Hd_ind_results <- perform_t_tests(Hydraulic_data, "HydraulicDiameter", indiv_test)
Hd_ind_results$trait <- "Dh"

VD_ind_results <- perform_t_tests(Hydraulic_data, "VesselDensity", indiv_test)
VD_ind_results$trait <- "VD"

Fv_ind_results <- perform_t_tests(Hydraulic_data, "VesselFraction", indiv_test)
Fv_ind_results$trait <- "Fv"

Kmax_ind_results <- perform_t_tests(Hydraulic_data, "Kmax", indiv_test)
Kmax_ind_results$trait <- "Kmax"

Fpit_ind_results <- perform_t_tests(PitFraction_data, "PitFraction", indiv_test)
Fpit_ind_results$trait <- "Fpit"

Dpo_ind_results <- perform_t_tests(PitDiOp_data, "PitOpening", indiv_test)
Dpo_ind_results$trait <- "Dpo"

Dpit_ind_results <- perform_t_tests(PitDiOp_data, "PitDiameter", indiv_test)
Dpit_ind_results$trait <- "Dpit"

# Combine all results
ind_results <- rbind(
  Wall_ind_result,
  Hd_ind_results,
  VD_ind_results,
  VesselD_ind_result,
  VesselDTop_ind,  # Fixed: Previously missing in rbind()
  Fv_ind_results,
  Kmax_ind_results,
  Fpit_ind_results,
  Dpo_ind_results,
  Dpit_ind_results
)
ind_results$pvalue <- ifelse(ind_results$pvalue<0.0001,"<0.0001",round(ind_results$pvalue,digits = 5))
pv_wide <- ind_results %>% select(indiv_tested,pvalue,trait) %>% pivot_wider(names_from = trait,values_from = pvalue, names_prefix = "Pvalue_")
diff_wide <- ind_results %>% select(indiv_tested,difference_observed,trait) %>% pivot_wider(names_from = trait,values_from = difference_observed, names_prefix = "Diff_")
Rdiff_wide <- ind_results %>% select(indiv_tested,Rdiff,trait) %>% pivot_wider(names_from = trait,values_from = Rdiff, names_prefix = "Effect_")
ind_wide <- merge(diff_wide,pv_wide,by="indiv_tested") %>% merge(Rdiff_wide,by="indiv_tested")
colnames(ind_wide)

plot(density(abs(ind_results$Rdiff)))

# Extract column names
diff_cols <- grep("^Diff_", colnames(ind_wide), value = TRUE) %>% sort()
pval_cols <- grep("^Pvalue_", colnames(ind_wide), value = TRUE) %>% sort()
effect_cols <- grep("^Effect_", colnames(ind_wide), value = TRUE) %>% sort()

# Interleave Diff_ and Pvalue_ columns
reordered_cols <- as.vector(rbind(diff_cols, pval_cols,effect_cols))

# Final column order: indiv_tested first, then interleaved Diff_ and Pvalue_ columns
final_col_order <- c("indiv_tested", reordered_cols)

# Reorder the dataframe
ind_wide <- ind_wide[, final_col_order]

ind_wide
cite("DHARMa")
library(glm)
citation("glmmTMB")
