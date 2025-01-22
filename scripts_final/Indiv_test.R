###individual test
# Load data
Wall_data <- read.csv(here("data", "processed", "Wall_data.csv"))
VesselDiameter_data<- read.csv(here("data", "processed", "VesselDiameter_data.csv")) 
Hydraulic_data<- read.csv(here("data", "processed", "HydraulicData.csv"))
PitFraction_data<- read.csv(here("data", "processed", "PitFraction_data.csv"))
PitDiOp_data<- read.csv(here("data", "processed", "PitDiOp_data.csv"))

source(here("scripts", "Functions.R"))

relevel_factors(ls())

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

perform_t_tests <- function(data, trait_col, indiv_test) {
  # Initialize a results data frame
  results_df <- data.frame(
    indiv_tested = character(),
    difference_observed = numeric(),
    CI95 = character(),
    pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each pair of individuals
  for (ind in indiv_test) {
    # Filter the data for the current pair
    filtered_data <- data %>% filter(indiv %in% ind)
    
    # Ensure there are exactly two groups
    if (length(unique(filtered_data$indiv)) == 2) {
      # Extract the values of the specified trait for each group
      group1 <- filtered_data[[trait_col]][filtered_data$indiv == ind[1]]
      group2 <- filtered_data[[trait_col]][filtered_data$indiv == ind[2]]
      
      # Perform Welch's t-test
      result <- t.test(group1, group2, var.equal = FALSE)
      
      # Extract relevant information
      indiv_tested <- paste(ind[1], "vs", ind[2])
      difference_observed <- result$estimate[1] - result$estimate[2]
      CI95 <- paste0("(", round(result$conf.int[1], 2), ", ", round(result$conf.int[2], 2), ")")
      pvalue <- ifelse(result$p.value < 0.001, "<0.001", round(result$p.value, 5))
      
      # Append the results to the results data frame
      results_df <- rbind(
        results_df,
        data.frame(
          indiv_tested = indiv_tested,
          difference_observed = round(difference_observed, 3),
          CI95 = CI95,
          pvalue = round(pvalue, 5),
          stringsAsFactors = FALSE
        )
      )
    } else {
      cat("\nSkipping pair:", ind, "- Insufficient or invalid data.\n")
    }
  }
  
  return(results_df)
}

Wall_ind_result <- perform_t_tests(Wall_data,"WallThickness",indiv_test)
Wall_ind_result$trait <- "Tw"
VesselD_ind_result <- perform_t_tests(VesselDiameter_data,"VesselDiameter",indiv_test)
VesselD_ind_result$trait <- "D"
VesselDTop_ind <- perform_t_tests(VesselDiameter_data %>% group_by(ssp,indiv) %>%
                                    filter(VesselDiameter >= quantile(VesselDiameter, 0.9)) %>%
                                    ungroup(),"VesselDiameter",indiv_test)
VesselDTop_ind$trit <- "Dtop"
Hd_ind_results <- perform_t_tests(Hydraulic_data,"HydraulicDiameter",indiv_test)
Hd_ind_results$trait <- "Dh"
VD_ind_results <- perform_t_tests(Hydraulic_data,"VesselDensity",indiv_test)
VD_ind_results$trait <- "VD"
Fv_ind_results <- perform_t_tests(Hydraulic_data,"VesselFraction",indiv_test)
Fv_ind_results$trait <- "Fv"
Kmax_ind_results <- perform_t_tests(Hydraulic_data,"Kmax",indiv_test)
Kmax_ind_results$trait <- "Kmax"
Fpit_ind_results <- perform_t_tests(PitFraction_data,"PitFraction",indiv_test)
Fpit_ind_results$trait <- "Fpit"
Dpo_ind_results <- perform_t_tests(PitDiOp_data,"PitOpening",indiv_test)
Dpo_ind_results$trait <- "Dpo"
Dpit_ind_results <- perform_t_tests(PitDiOp_data,"PitDiameter",indiv_test)
Dpit_ind_results$trait <- "Dpit"

ind_results <- rbind(Wall_ind_result,
                     Hd_ind_results,
                     VD_ind_results,
                     VesselD_ind_result,
                     VesselD_ind_result,
                     Fv_ind_results,
                     Kmax_ind_results,
                     Fpit_ind_results,
                     Dpo_ind_results,
                     Dpit_ind_results
                     )
