#OverSampling (SMOTE)

# Setting File Location
# From Menu: Session --> Set Working Directory --> To Source File Location demelisiniz
#Clear environment
rm(list=ls())

# Reading data as a dataframe
data <- read.table("cancer.txt", header = FALSE, sep = "\t")

# Changing column names of dataframe
colnames(data) <- c("sample_code_number", 
                    "clump_thickness", 
                    "uniformity_of_cell_size", 
                    "uniformity_of_cell_shape", 
                    "marginal_adhesion", 
                    "single_epithelial_cell_size", 
                    "bare_nuclei", 
                    "bland_chromatin", 
                    "normal_nucleoli", 
                    "mitosis", 
                    "classes")

# Removing first columns since it is not necessary for calculations
data<-data[,-1]

# Convert Target variable type to factor
data$classes <- as.factor(data$classes)

# Rename the factor variable in target column
# Type of Tumor: 2-benign, 4- malignant
levels(data$classes)
levels(data$classes) <- c("benign", "malignant")
levels(data$classes)

# Target variable distribution
table(data$classes)

# Malignant category is defined as reference value (This means we changed level order)
data$classes <- relevel(data$classes, ref = "malignant")

# Target variable distribution
table(data$classes)


# Summary of dataset
summary(data)

# Installing ggplot2 package
install.packages("ggplot2")
library(ggplot2)
ggplot(data, aes(x = classes, fill = classes)) +
  geom_bar()


# Analyzes of imbalanced data set
# Data Preparation
install.packages("caret")
install.packages("e1071")

library(caret)
set.seed(42)
trainingset <- createDataPartition(data$classes, p = 0.7, list = FALSE)

train_data <- data[trainingset, ]
test_data  <- data[-trainingset, ]

# Train the model
set.seed(42)
model_knn <- train(classes ~ ., data = train_data, method = "knn")

#Make predictions
predictions <- predict(model_knn, newdata = test_data, type = "raw")
predictionPerformance <- confusionMatrix(predictions, test_data$classes)

# OVER-SAMPLING
# To make oversampling trcontrol parameter is added
set.seed(42)
model_knn_over <- train(classes ~ ., data = train_data, method = "knn", trControl=trainControl(sampling = "up"))

#Predict test dataset
predictionOver <- predict(model_knn_over, newdata = test_data, type = "raw")
predictionOverPerformance <- confusionMatrix(predictionOver, test_data$classes)

# SMOTE Oversampling
install.packages("DMwR")
library(DMwR)

# Only trcontrol parameter is changed
set.seed(42)
model_knn_smote <- train(classes ~ ., data = train_data, method = "knn", trControl=trainControl(sampling = "smote"))

#Predictions:
smotePredictions <- predict(model_knn_smote, newdata = test_data, type = "raw")
smotePerformance <- confusionMatrix(smotePredictions, test_data$classes)

