####covariance
install.packages("ggbiplot")
library(GGally)
library(ggstatsplot)
library(devtools)
library(ggbiplot)

ggpairs(Median_data,columns = 3:11,aes(colour = parasitism),
        lower = list(continuous="smooth"),
        upper=list(continuous=wrap("cor",method="pearson")))
ggpairs(Mean_data,columns = 3:11,aes(colour = parasitism),
        lower = list(continuous="smooth"),
        upper=list(continuous=wrap("cor",method="pearson")))
library(ggstatsplot)

# Create grouped correlation matrix without the `title` argument
grouped_corr_plot <- grouped_ggcorrmat(
  data = Median_data,
  type = "p",                 # Non-parametric correlation
  grouping.var = parasitism,   # Grouping variable
  method = "pearson",         # Spearman rank correlation
  label = TRUE                 # Display correlation coefficients
)

# Add a title after the plot is created
grouped_corr_plot + labs(title = "Grouped Correlation Matrix by Parasitism")



# Perform PCAhttp://127.0.0.1:25125/graphics/3a69a903-5638-4b79-b36b-54fb268756f3.png
pca_median <- Median_data[3:11]
pc <- prcomp(pca_median, center = TRUE, scale. = TRUE)

# Print PCA results
print(pc)
summary(pc)

# Create the PCA biplot
library(ggbiplot)

p <- ggbiplot(pc, 
              obs.scale = 1, 
              var.scale = 1, 
              groups = Median_data$parasitism, 
              ellipse = TRUE, 
              circle = TRUE, 
              ellipse.prob = 0.68)

# Add point labels (names) from Median_data$ssp
p <- p + 
  geom_text(aes(label = Median_data$ssp), 
            vjust = -0.5, 
            hjust = 0.5, 
            size = 3, 
            color = "black")

# Customize and display the plot
p <- p + 
  theme_minimal() + 
  ggtitle("PCA Biplot with Point Labels") +
  theme(legend.position = "right")

print(p)

