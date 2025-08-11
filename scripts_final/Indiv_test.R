# Load necessary libraries
library(tidyverse)
library(here)
source(here("scripts", "Functions.R"))

# Load the datasets
PitMembrane_data <- read.csv(here("data", "processed", "PitMembrane_data.csv"))
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data <- read.csv(here("data", "processed", "VesselDiameter_data.csv")) 
Hydraulic_data <- read.csv(here("data", "processed", "HydraulicData.csv"))
PitFraction_data <- read.csv(here("data", "processed", "PitFraction_data.csv"))
PitDiOp_data <- read.csv(here("data", "processed", "PitDiOp_data.csv"))

# Add unique identifier for each individual in PitMembrane_data
PitMembrane_data$indiv <- paste(PitMembrane_data$ssp, 1, sep = " ")

# Relevel factors in the dataset
relevel_factors(ls())

# Define list of individual pairs for testing
indiv_test <- list(
  c("Psittacanthus robustus 1", "Vochysia thyrsoidea 1"),
  c("Psittacanthus robustus 2", "Vochysia thyrsoidea 2"),
  c("Psittacanthus robustus 3", "Vochysia thyrsoidea 3"),
  c("Phoradendron perrotettii 1", "Tapirira guianensis 1"),
  c("Phoradendron perrotettii 2", "Tapirira guianensis 2"),
  c("Phoradendron perrotettii 3", "Tapirira guianensis 3"),
  c("Struthanthus rhynchophyllus 1", "Tipuana tipu 1"),
  c("Struthanthus rhynchophyllus 2", "Tipuana tipu 2"),
  c("Struthanthus rhynchophyllus 3", "Tipuana tipu 3"),
  c("Viscum album 1", "Populus nigra 1"),
  c("Viscum album 2", "Populus nigra 2"),
  c("Viscum album 3", "Populus nigra 3")
)

perform_t_tests <- function(data, trait_col, indiv_test) {
  # Ensure 'indiv' column exists
  if (!"indiv" %in% colnames(data)) {
    stop(paste("Error: 'indiv' column not found in dataset for trait:", trait_col))
  }
  
  # Initialize a list to store the results
  results_list <- list()
  
  for (ind in indiv_test) {
    # Filter the data for the current individual pair
    filtered_data <- data %>% filter(indiv %in% ind)
    
    # Check if there are exactly two distinct individuals in the pair
    if (length(unique(filtered_data$indiv)) == 2) {
      group1 <- filtered_data[[trait_col]][filtered_data$indiv == ind[1]]
      group2 <- filtered_data[[trait_col]][filtered_data$indiv == ind[2]]
      
      # Proceed if both groups have more than one data point
      if (length(group1) > 1 & length(group2) > 1) {
        result <- t.test(group1, group2, var.equal = FALSE)
        
        indiv_tested <- paste(ind[1], "vs", ind[2])
        difference_observed <- result$estimate[1] - result$estimate[2]
        CI95 <- paste0("(", round(result$conf.int[1], 2), ", ", round(result$conf.int[2], 2), ")")
        pvalue <- result$p.value
        mean_group2 <- mean(group2, na.rm = TRUE)  # Compute mean of group 2
        
        # Compute Rdiff safely (avoid division by zero)
        Rdiff <- ifelse(mean_group2 != 0, (difference_observed / mean_group2) * 100, NA)
        
        # Add the results to the list, include 'trait'
        results_list[[length(results_list) + 1]] <- data.frame(
          indiv_tested = indiv_tested,
          difference_observed = round(difference_observed, 3),
          CI95 = CI95,
          pvalue = round(pvalue, 5),
          Rdiff = round(Rdiff, 3),
          trait = trait_col,  # Ensure the 'trait' column is included
          stringsAsFactors = FALSE
        )
        
      } else {
        cat("\nSkipping pair:", ind, "- Not enough data for valid t-test.\n")
      }
    } else {
      cat("\nSkipping pair:", ind, "- Insufficient or invalid data.\n")
    }
  }
  
  # Combine the list of results into a single data frame
  results_df <- bind_rows(results_list)
  
  return(results_df)
}


# Run t-tests for each trait
test_results <- list()

test_results[["Wall"]] <- perform_t_tests(Wall_data, "WallThickness", indiv_test)
test_results[["VesselDiameter"]] <- perform_t_tests(VesselDiameter_data, "VesselDiameter", indiv_test)
test_results[["VesselDTop"]] <- perform_t_tests(
  VesselDiameter_data %>% group_by(ssp, indiv) %>% filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>% ungroup(),
  "VesselDiameter",
  indiv_test
)
test_results[["VesselDTop"]]$trait <- "VesselDTop"
test_results[["HydraulicDiameter"]] <- perform_t_tests(Hydraulic_data, "HydraulicDiameter", indiv_test)
test_results[["VesselDensity"]] <- perform_t_tests(Hydraulic_data, "VesselDensity", indiv_test)
test_results[["VesselFraction"]] <- perform_t_tests(Hydraulic_data, "VesselFraction", indiv_test)
test_results[["Kmax"]] <- perform_t_tests(Hydraulic_data, "Kmax", indiv_test)
test_results[["PitFraction"]] <- perform_t_tests(PitFraction_data, "PitFraction", indiv_test)
test_results[["PitOpening"]] <- perform_t_tests(PitDiOp_data, "PitOpening", indiv_test)
test_results[["PitDiameter"]] <- perform_t_tests(PitDiOp_data, "PitDiameter", indiv_test)
test_results[["PitMembrane_pcd"]] <- perform_t_tests(PitMembrane_data, "pcd", indiv_test)
test_results[["PitMembrane_Tpm"]] <- perform_t_tests(PitMembrane_data, "Tpm", indiv_test)



# Combine all test results into one dataframe
ind_results <- bind_rows(test_results)


# Format p-values
ind_results$pvalue <- ifelse(ind_results$pvalue < 0.0001, "<0.0001", round(ind_results$pvalue, 5))

# Reshape the data into wide format for easier interpretation
pv_wide <- ind_results %>%
  select(indiv_tested, pvalue, trait) %>%
  pivot_wider(names_from = trait, values_from = pvalue, names_prefix = "Pvalue_")

diff_wide <- ind_results %>%
  select(indiv_tested, difference_observed, trait) %>%
  pivot_wider(names_from = trait, values_from = difference_observed, names_prefix = "Diff_")

Rdiff_wide <- ind_results %>%
  select(indiv_tested, Rdiff, trait) %>%
  pivot_wider(names_from = trait, values_from = Rdiff, names_prefix = "Effect_")

# Merge wide-format data
ind_wide <- merge(diff_wide, pv_wide, by = "indiv_tested") %>%
  merge(Rdiff_wide, by = "indiv_tested")

# Reorder the columns: Diff_ columns first, then Pvalue_ columns, then Effect_ columns
diff_cols <- grep("^Diff_", colnames(ind_wide), value = TRUE)
pval_cols <- grep("^Pvalue_", colnames(ind_wide), value = TRUE)
effect_cols <- grep("^Effect_", colnames(ind_wide), value = TRUE)

# Interleave Diff_, Pvalue_, and Effect_ columns
reordered_cols <- c(diff_cols, pval_cols, effect_cols)

# Final column order: indiv_tested first, followed by interleaved columns
final_col_order <- c("indiv_tested", reordered_cols)

# Reorder the dataframe according to the final column order
ind_wide <- ind_wide[, final_col_order]

# Output the final table
colnames(ind_wide)

write.csv(ind_results,file = here("outputs","tables","Indiv_Welches.csv"))
          
          