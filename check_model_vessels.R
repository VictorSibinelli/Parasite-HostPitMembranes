models_Vdiam <- list()
for (pair in species_pairs) {
  subset_data <- subset(VesselDiameter_data, ssp %in% pair)
  
  m1 <- lme(
    VesselDiameter ~ ssp,
    random = ~ 1 | ssp / indiv,
    data = subset_data,
    control = list(maxIter = 150, msMaxIter = 150),
    weights = varIdent(form = ~ 1 | ssp),
    method = "ML"
  )
  models_Vdiam[[paste0(pair,collapse=" x ")]] <- m1
}  
for (i in length(models_Vdiam)) {
  # Extract the model for the current pair
  model <- models_Vdiam[[i]]
  # Residual plot
  residplot(model, newwd = FALSE)
  title(sub = names(models_Vdiam)[i])
  
  # Model check plot
  print(check_model(model, show_dots = FALSE))
  
  # Cook's Distance plot
  CookD(model, idn = 20, newwd = FALSE)
  abline(h = 4 / nrow(subset_data), col = "red")
  title(sub = names(models_Vdiam)[i])
}
nrow(m1)
