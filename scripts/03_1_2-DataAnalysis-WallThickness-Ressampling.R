######################################################################
#
# Victor Sibinelli (victor.sibinelli@usp.br / sibinelli95@gmail.com)
# 13/07/2024
# Script 03.1.2 - Data Analysis Ressampling- Vessel walls
#################################################################
library(here)
source(here("scripts", "02_1-TestAssumptions-WallThickness.R"))
rm(list=ls())
wdata <- read.csv(here("data", "processed", "wdata.csv"))
wdata_clean <- read.csv(here("data", "processed", "wdata_clean.csv"))



# Data subsets by species
species_data <- list(
  "Psittacanthus robustus" = wdata_clean[wdata_clean$ssp == "Psittacanthus robustus", ],
  "Vochysia thyrsoidea" = wdata_clean[wdata_clean$ssp == "Vochysia thyrsoidea", ],
  "Phoradendron perrotettii" = wdata_clean[wdata_clean$ssp == "Phoradendron perrotettii", ],
  "Tapirira guianensis" = wdata_clean[wdata_clean$ssp == "Tapirira guianensis", ],
  "Struthanthus rhynchophyllus" = wdata_clean[wdata_clean$ssp == "Struthanthus rhynchophyllus", ],
  "Tipuana tipu" = wdata_clean[wdata_clean$ssp == "Tipuana tipu", ],
  "Viscum album" = wdata_clean[wdata_clean$ssp == "Viscum album", ],
  "Populus nigra" = wdata_clean[wdata_clean$ssp == "Populus nigra", ]
)


# List of species pairs for comparison
species_pairs <- list(
  c("Psittacanthus robustus", "Vochysia thyrsoidea"),
  c("Phoradendron perrotettii", "Tapirira guianensis"),
  c("Struthanthus rhynchophyllus", "Tipuana tipu"),
  c("Viscum album", "Populus nigra")
)

# Calculate means for each label and retain all columns
WT_avg <- wdata %>%
  group_by(label, ssp) %>%
  summarise(
    wthickness = mean(wthickness, na.rm = TRUE),
    .groups = "drop"  # Drop the grouping after summarizing
  )
WT_avg <- WT_avg %>%
  mutate(parasitism = case_when(
    ssp %in% c("Psittacanthus robustus", 
               "Phoradendron perrotettii", 
               "Struthanthus rhynchophyllus", 
               "Viscum album") ~ "Parasite",
    TRUE ~ "Host"
  ))

# List of data frame names
dataframes <- ls()
# Relevel the factors for each data frame
for (df_name in dataframes) {
  df <- get(df_name) # Get the data frame by name
  
  if ("ssp" %in% colnames(df)) { # Check if 'ssp' column exists
    df$ssp <- factor(df$ssp, levels = c(
      "Psittacanthus robustus", "Vochysia thyrsoidea",
      "Phoradendron perrotettii", "Tapirira guianensis",
      "Struthanthus rhynchophyllus", "Tipuana tipu",
      "Viscum album", "Populus nigra"
    ))
    assign(df_name, df) # Assign the modified data frame back to its name
  }
  rm(df, df_name) # remove duplicated dataframe
}

# Number of iterations for the permutation test
iterations <- 1000




WT_obs <- tapply(WT_avg$wthickness,WT_avg$parasitism, FUN = mean)
WT_resample <- t(replicate(iterations,
                           shuffle_means(
                             x=WT_avg,cols ="wthickness",cat = "parasitism")))
WT_obs
head(WT_resample)
WT_resample[1,] <- WT_obs
WT_diff <- WT_resample[,1]-WT_resample[,2]

#calculate p-values
WT_pvalue <- sum( WT_diff <=  WT_diff[1]|
                    WT_diff >= ( WT_diff[1]*-1))/
  length(WT_resample)
#vessel density
plot(density(WT_diff))
abline(v = WT_diff[1], col = "red")
abline(v = WT_diff[1]*-1, col = "red")
text(x=0.5,y=0.5,paste("p-value=",WT_pvalue))

bootstrap_result <- t(replicate(iterations,
                                shuffle_means(
                                  x=WT_avg,cols ="wthickness",cat = "parasitism",rcol=T)))
boot_diff <- bootstrap_result[,1]-bootstrap_result[,2]
CI_95 <- quantile(boot_diff,c(0.025,0.975))






###pair-wise comparison
Wssp_resample <- t(replicate(iterations,
                           shuffle_means(
                             x=WT_avg,cols ="wthickness",cat = "ssp")))
WTssp_obs <- tapply(WT_avg$wthickness,WT_avg$ssp,mean)
WTssp_obs
head(Wssp_resample)
Wssp_resample[1,] <- WTssp_obs

PrvsVt <- Wssp_resample[,"Psittacanthus robustus"]-Wssp_resample[,"Vochysia thyrsoidea"]
#calculate p-values
PrVtpvalue <- sum( PrvsVt <= PrvsVt[1]|
                      PrvsVt >= ( PrvsVt[1]*-1))/
  length(PrvsVt)
plot(density(PrvsVt))
abline(v = PrvsVt[1], col = "red")
abline(v = PrvsVt[1]*-1, col = "red")
text(x=0.5,y=0.5,paste("p-value=",PrVtpvalue))

PpvsTg <- Wssp_resample[,"Phoradendron perrotettii"]-Wssp_resample[,"Tapirira guianensis"]
#calculate p-values
PvTgpvalue <- sum( PpvsTg <= PpvsTg[1]|
                     PpvsTg >= ( PpvsTg[1]*-1))/
  length(PpvsTg)
plot(density(PpvsTg))
abline(v = PpvsTg[1], col = "red")
abline(v = PpvsTg[1]*-1, col = "red")
text(x=0.5,y=0.5,paste("p-value=",PvTgpvalue))

StvsTt <- Wssp_resample[,"Struthanthus rhynchophyllus"]-Wssp_resample[,"Tipuana tipu"]
#calculate p-values

StTtpvalue <- sum( StvsTt <= StvsTt[1]|
                     StvsTt >= ( StvsTt[1]*-1))/
  length(StvsTt)
plot(density(StvsTt))
abline(v = StvsTt[1], col = "red")
abline(v = StvsTt[1]*-1, col = "red")
text(x=0.5,y=0.5,paste("p-value=",StTtpvalue))

VavsPn <- Wssp_resample[,"Viscum album"]-Wssp_resample[,"Vopulus nigra"]
#calculate p-values
VaPnpvalue <- sum( VavsPn <= VavsPn[1]|
                     VavsPn >= ( VavsPn[1]*-1))/
  length(VavsPn)
plot(density(VavsPn))
abline(v = VavsPn[1], col = "red")
abline(v = VavsPn[1]*-1, col = "red")
text(x=0.5,y=0.5,paste("p-value=",VaPnpvalue))







WTssp_pvalue
#vessel density
plot(density(WTssp_resample[,1]))
abline(v = WTssp_resample[1,1], col = "red")
abline(v = WTssp_resample[1,1]*-1, col = "red")
text(x=0.5,y=0.5,paste("p-value=",WTssp_pvalue[1]))

plot(density(WTssp_resample[,2]))
abline(v = WTssp_resample[1,2], col = "red")
abline(v = WTssp_resample[1,2]*-1, col = "red")
text(x=0.5,y=0.5,paste("p-value=",WTssp_pvalue[2]))

plot(density(WTssp_resample[,3]))
abline(v = WTssp_resample[1,3], col = "red")
abline(v = WTssp_resample[1,3]*-1, col = "red")
text(x=0.5,y=0.5,paste("p-value=",WTssp_pvalue[3]))


