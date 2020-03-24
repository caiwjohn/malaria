library(tidyverse)
library(useful)
library(ggfortify)
library(mgsub)
theme_set(theme_gray(base_size = 18))

#####
# PCA
#####
# Load data
expr<- readRDS("../clean_data/expression.RDS")
samples<- readRDS("../clean_data/model_data.RDS")

# Merge meta into expression data
## IMPORTANT: this drops extra samples in expr
model<- merge(samples, expr)

#####
# Auto decomposition
#####
# Isolate expression data and compute
num<- model[,8:ncol(model)]
pca_results<- prcomp(num, scale=F, center= T, tol= 0.01)

# Visualize PCA
autoplot(pca_results, data= model, label= F, colour= "class") +
  labs(colour="IC50", title= "Sample PCA") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle= element_text(hjust=0.5), 
        legend.text = element_text(size=12), legend.title = element_text(size = 12)) +
  scale_colour_discrete(name="IC50 Status")

#####
# After LASSO Complete
#####
# Reduce transcripts to those selected from prediction
red_trans<- expr[,which(colnames(expr) %in% features$features)]

# Repeat with reduced set
pca_results <- red_trans %>%
  prcomp(scale=F, center= T, tol= 0.01)

# Visualize PCA
autoplot(pca_results, data= model, label= F, colour= "class") +
  labs(colour="IC50", title= "Sample PCA", subtitle= "LASSO Selected Transcripts") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle= element_text(hjust=0.5), 
        legend.text = element_text(size=12), legend.title = element_text(size = 12)) +
  scale_colour_discrete(name="IC50 Status")
