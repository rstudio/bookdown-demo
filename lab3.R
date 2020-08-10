# Load the R package that we need for this analysis
library(eRm)

library(readr) # Import the data from your computer
transreas <- read_csv("transreas.csv")


# Trim the data
Di_Rash_data <- transreas[,c(-1,-2)]
head(Di_Rash_data) # Take a look

# Run the Dichotomous Data
Di_Model <- RM(Di_Rash_data)
summary(Di_Model)


