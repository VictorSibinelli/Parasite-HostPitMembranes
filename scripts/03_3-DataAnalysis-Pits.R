######################################################################
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.3 - Data Analysis - Pit membranes - Bootstrap
#################################################################
library(here)
# Load script for testing assumptions of pit data and remove any existing objects from the environment
#source(here("scripts", "02_3-TestAssumptions-Pit.R"))
rm(list = ls())  # Clear environment
source(here("scripts","00-library.R"))
source(here("scripts", "Functions.R"))  # Load custom functions

# Load your data (assuming these lines are already included)
pitdata <- fread(here("data", "processed", "pitdata.csv"))[,c(1,7,8,11)] %>% as.tibble()  # Load main pit data
pitOdata <- fread(here("data", "processed", "pitOdata.csv")) %>% as.tibble()  # Load pit opening data


# Summarize data for each species and parasite/host
Pit_EV <- pitOdata %>%
  group_by(label) %>%
  summarise(
    ssp = first(ssp),                # Extract first species (ssp) name for each label
    parasitism = first(parasitism),   # Extract parasitism status (Parasite/Host)
    PitDiameter = median(PitDiameter, na.rm = TRUE),   # Calculate median PitDiameter
    PitOpening = median(PitOpening, na.rm = TRUE),     # Calculate median PitOpening
    .groups = "drop"
  )

# List of species pairs for comparisons
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),   # Pair 1
  c("Phoradendron perrotettii", "Tapirira guianensis"), # Pair 2
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),     # Pair 3
  c("Viscum album", "Populus nigra")                    # Pair 4
)

# Summarize Pit_EV by parasitism
Obs1 <- Pit_EV %>%
  group_by(parasitism) %>%
  summarize(
    Grouping = first(parasitism),  # Rename parasitism to Grouping
    PitDiameter = mean(PitDiameter, na.rm = TRUE),
    PitOpening = mean(PitOpening, na.rm = TRUE),
    .groups = "drop"
  ) %>%select(Grouping, PitDiameter, PitOpening) %>%
  rbind(Pit_EV %>%
          group_by(ssp) %>%
          summarize(
            Grouping = first(ssp),  # Rename ssp to Grouping
            PitDiameter = mean(PitDiameter, na.rm = TRUE),
            PitOpening = mean(PitOpening, na.rm = TRUE),
            .groups = "drop")%>%select(Grouping, PitDiameter, PitOpening) 
  )

# Summarize pitdata by parasitism
Obs2 <- pitdata %>%
  group_by(parasitism) %>%
  summarize(
    Grouping = first(parasitism),  # Rename parasitism to Grouping
    Pcd = mean(pcd, na.rm = TRUE),
    Tpm = mean(pitavg, na.rm = TRUE),
    .groups = "drop"
  ) %>%select(Grouping, Pcd, Tpm) %>%
  rbind(pitdata %>%
          group_by(ssp) %>%
          summarize(
            Grouping = first(ssp),  # Rename ssp to Grouping
            Pcd = mean(pcd, na.rm = TRUE),
            Tpm = mean(pitavg, na.rm = TRUE),
            .groups = "drop") %>% select(Grouping, Pcd, Tpm)
  )

# Merge the two summarized data frames by the Grouping column
Obs_values <- merge(Obs1, Obs2, by = "Grouping", all = TRUE) %>% as.tibble() %>%
  arrange(Grouping)  # Change all=TRUE to all.x=TRUE for left join

rm(list = c("Obs1", "Obs2", "pitOdata"))


# Define the variables of interest for bootstrapping
vars <- colnames(Pit_EV[4:5])  # For Pit_EV: PitDiameter and PitOpening
vars2 <- colnames(pitdata[, c(2, 4)])  # 

# Relevel factors for analysis
relevel_factors(ls())  # Reorder factors for categorical variables

# Set the number of bootstrap iterations and random seed for reproducibility
it <- 1000# Number of bootstrap replicates
set.seed(42)  # Set seed for consistent random sampling

# Bootstrap sampling function
bootstrap_sampling <- function(data, index, vars) {
  replicate(n = it, {
    sapply(vars, function(var) {
      data %>%
        subset(parasitism == index) %>%  # Subset data for the specified parasitism
        shuffle_means(cols = var, cat = "parasitism", rcol = TRUE)  # Shuffle and calculate means
    })
  }, simplify = TRUE) %>% t()  # Transpose the result
}

# Bootstrap for Host and Parasite species in Pit_EV
host_boot <- bootstrap_sampling(Pit_EV, "Host", vars)
para_boot <- bootstrap_sampling(Pit_EV, "Parasite", vars)

# Bootstrap for Host and Parasite species in pitdata
host_boot2 <- bootstrap_sampling(pitdata, "Host", vars2)
para_boot2 <- bootstrap_sampling(pitdata, "Parasite", vars2)

# Bootstrap for each species (ssp) in Pit_EV
ssp_boot <- setNames(lapply(levels(Pit_EV$ssp), function(current_ssp) {
  t(replicate(it, sapply(vars, function(var) {
    shuffle_means(Pit_EV[Pit_EV$ssp == current_ssp, ], cols = var, cat = "ssp", rcol = TRUE)
  })))
}), levels(Pit_EV$ssp))  # Set species names as list element names

# Bootstrap for each species (ssp) in pitdata
ssp_boot2 <- setNames(lapply(levels(pitdata$ssp), function(current_ssp) {
  t(replicate(it, sapply(vars2, function(var) {
    shuffle_means(pitdata[pitdata$ssp == current_ssp, ], cols = var, cat = "ssp", rcol = TRUE)
  })))
}), levels(pitdata$ssp))  # Set species names as list element names

# Combine bootstrap results into a single list
CI_boot <- list(
  Parasite = cbind(para_boot, para_boot2),  # Combine parasite results from both data sources
  Host = cbind(host_boot, host_boot2)       # Combine host results from both data sources
)

# Combine species-specific results into the combined list
CI_boot <- c(CI_boot, lapply(levels(Pit_EV$ssp), function(current_ssp) {
  cbind(
    ssp_boot[[current_ssp]],
    ssp_boot2[[current_ssp]]
  )
}))
# Name each species element appropriately
names(CI_boot)[3:length(CI_boot)] <- levels(Pit_EV$ssp)

rm(list = c("host_boot","host_boot2","para_boot","para_boot2","ssp_boot","ssp_boot2"))

# Create a list of datasets and their corresponding variable sets
data_list <- list(
  list(data = Pit_EV, vars = vars),
  list(data = pitdata, vars = vars2)
)

boot_results <- list()  # Initialize an empty list to store results

# Loop over each dataset and its corresponding variables
for (dataset in data_list) {
  for (v in dataset$vars) {
    boot_results[[v]] <- t(replicate(it, 
                                     shuffle_means(
                                       dataset$data, 
                                       cols = v, 
                                       cat = "parasitism", 
                                       rcol = TRUE)))
    
    # Assign column names based on parasitism levels
    colnames(boot_results[[v]]) <- levels(dataset$data$parasitism)
  }
}

lapply(boot_results,head)

# Initialize list to store results
Ssp_bootstrap_result <- list()  

for (data_item in data_list) {
  dataset <- data_item$data  # Extract the dataset
  vars <- data_item$vars      # Extract the corresponding variables
  
  for (pair in species_pairs) {
    subset_data <- subset(dataset, ssp %in% pair)  # Subset for the current species pair
    
    for (v in vars) {
      # Generate bootstrap samples for the current variable
      pair_boot <- t(replicate(it, 
                               shuffle_means(subset_data, cols = v, cat = "ssp", rcol = TRUE)))
      
      if (ncol(pair_boot) == length(pair)) {
        colnames(pair_boot) <- pair
        
        # Combine results
        Ssp_bootstrap_result[[v]] <- if (is.null(Ssp_bootstrap_result[[v]])) {
          pair_boot
        } else {
          cbind(Ssp_bootstrap_result[[v]], pair_boot)
        }
      }
    }
  }
}


lapply(Ssp_bootstrap_result,head)


boot_diff <- lapply(boot_results,function(x){
  apply(x,1,diff)
}) %>% do.call(what=cbind)

# Calculate differences for each variable and species pair
ssp_boot_diff <- lapply(Ssp_bootstrap_result, function(x) {
  # Calculate differences for each species pair, handling missing pairs
  pair_diff_df <- as.data.frame(do.call(cbind, lapply(species_pairs, function(pair) {
    if (all(pair %in% colnames(x))) {
      x[, pair[1]] - x[, pair[2]]  # Calculate differences
    } else {
      NA_real_  # Return NA if the pair is not found
    }
  })))
  
  # Name columns with species pair names
  colnames(pair_diff_df) <- sapply(species_pairs, function(pair) paste(pair, collapse = " vs "))
  
  return(pair_diff_df)
})

# Name the list elements with variable names
names(ssp_boot_diff) <- names(Ssp_bootstrap_result)

# Resulting list of data frames
lapply(ssp_boot_diff,head)

###########################################################
####testing

CI95 <- CI_boot %>%
  lapply(as.data.frame) %>%
  lapply(function(df) {
    # Replace non-numeric values with NA
    df[] <- lapply(df, function(col) {
      if (!is.numeric(col)) {
        return(rep(NA, length(col)))  # Replace non-numeric columns with NA
      } else {
        return(col)  # Keep numeric columns as they are
      }
    })
    # Apply quantile function to each column with NA handling
    apply(df, 2, function(q) {
      quantile(q, probs = c(0.025, 0.975), na.rm = TRUE)
    })
  })

head(boot_diff)
boot_diff[1,] <- as.matrix(Obs_values[1,2:5]-Obs_values[2,2:5])
PH_pvalue <- apply(boot_diff,2,function(x){
  sum(abs(x) >= abs(x[1]), na.rm = TRUE) / length(x)
})


# Calculate differences and combine as rows
ssp_obs_diffs <- rbind(
  Obs_values[3, 2:5] - Obs_values[4, 2:5],  # Difference between rows 3 and 4
  Obs_values[5, 2:5] - Obs_values[6, 2:5],  # Difference between rows 5 and 6
  Obs_values[7, 2:5] - Obs_values[8, 2:5],  # Difference between rows 7 and 8
  Obs_values[9, 2:5] - Obs_values[10, 2:5]  # Difference between rows 9 and 10
)
ssp_obs_diffs <- as.data.frame(ssp_obs_diffs)

colnames(ssp_obs_diffs) <- colnames(Obs_values)[2:5]
row.names(ssp_obs_diffs) <-colnames(ssp_boot_diff[[1]])

for (i in seq_along(ssp_boot_diff)) {
  x <- t(ssp_obs_diffs)
  ssp_boot_diff[[i]][1, ] <- x[i,]  # Assign first column of ssp_obs_diffs as the first row
}

Ssp_pvalue <- lapply(ssp_boot_diff, function(x) {
  apply(x, 2, function(y) {
    sum(abs(y) >= abs(y[1]), na.rm = TRUE) / length(y)
  })
}) %>% do.call(what=cbind) %>% rbind(PH_pvalue)

