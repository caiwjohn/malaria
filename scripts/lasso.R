library(tidyverse)
library(useful)
library(ggplot2)
library("glmnet")
library(lmtest)
library(textclean)
library(ROCR)
library(randomForest)

#####
# PREP DATA FOR PREDICTION
#####
# Load
expr<- readRDS("../clean_data/modelData.RDS")

# Build training and validation sets
model_dat <- expr %>%
  as_tibble() %>%
  group_by(class)

# Shuffle the order of samples
model_dat <- model_dat[sample(nrow(model_dat)),]

# OPTIONAL: shuffle class labels
#model_dat<- transform(model_dat, numeric_class= sample(numeric_class)) %>%
#  as_tibble() %>%
#  group_by(numeric_class)
  
# See how many of each class we have
length(which(model_dat$numeric_class==0))
length(which(model_dat$numeric_class==1))

# Select training and validation
train<- sample_n(model_dat, 200)
test<- model_dat[which(!(model_dat$GenotypeID %in% train$GenotypeID)),]

# Isolate predictors and responses
train_resp<- as.factor(as.matrix(train[,c("numeric_class")]))
train_pred<- as.matrix(train[,8:ncol(train)])
test_resp<- as.factor(as.matrix(test[,c("numeric_class")]))
test_pred<- as.matrix(test[,8:ncol(test)])

#####
# LASSO PRED
#####
# Determine optimal lambda value
fit<- cv.glmnet(train_pred, train_resp, family = "binomial", 
                type.measure= "class", nfolds = 10)

# Plot lambda optimization
plot(fit)

# Find corresponding minimum misclassification
mse.min <- fit$cvm[fit$lambda == fit$lambda.min]

# View coefficients
myCoefs<- coef(fit, s = "lambda.min")

# Assemble into df
features <- data.frame(
  features = myCoefs@Dimnames[[1]][which(myCoefs != 0 )],
  coefs    = myCoefs[ which(myCoefs != 0 ) ]
)

# OPTIONAL: Save list of selected features
features<- features[-c(1),]
write.table(features$features, file ="../viz_results/lasso_sel_transcripts.txt" , 
            sep="\t", row.names=FALSE, quote = FALSE)

# Predict classes for held-out data
model_pred<- as.numeric(predict(fit, newx = test_pred, s = "lambda.min", type="class"))

# Compute accuracy
acc<- length(which(model_pred==test_resp))/length(model_pred)






















