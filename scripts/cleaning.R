library(tidyverse)
library(mgsub)
library(naniar)
library(VIM)

# Read in data
meta<- read.table("../raw_data/metadata.txt", sep="\t")
exprRaw<- read.table("../raw_data/gse_matrix.txt", sep="\t", header=T, nrow=11005)

# Tidy expression data
rownames(exprRaw)<- exprRaw$ID_REF
exprRaw<- as.data.frame(t(exprRaw))
exprRaw$GenotypeID<- rownames(exprRaw)
expr <- exprRaw %>%
  select(GenotypeID, everything()) %>%
  slice(-1)

# Tidy meta data
colnames(meta)<- c("GenotypeID", "SampleID", "IC50", "Timepoint", "Origin")
meta$IC50<- as.numeric(meta$IC50)
meta$Timepoint<- mgsub(meta$Timepoint, 
                       c("timepoint: prior to artemisinin combination therapy (ACT)",
                         "timepoint: 0h post-sampling from patient",
                         "timepoint: 16h post-sampling from patient",
                         "timepoint: 24h post-sampling from patient",
                         "timepoint: 32h post-sampling from patient",
                         "timepoint: 40h post-sampling from patient",
                         "timepoint: 8h post-sampling from patient"),
                       c("prior", "0h", "16h", "24h", "32h", "40h", "8h"))
meta$Origin<- gsub("geographic origin: ", "", meta$Origin)

# Inspect distribution of IC50 values
ggplot(data= meta, aes(IC50)) +
  geom_bar()

# Create high (above 3rd quartile) and low (below 1st quartile) groups
summary(meta$IC50)
low<- meta[which(meta$IC50 < 265),]
high<- meta[which(meta$IC50> 753),]
high$class<- "High"
low$class<- "Low"
high$numeric_class<- 1
low$numeric_class<- 0

# Merge into one df for investigating
model_dat<- rbind(high, low)

# Fix class types
model_dat$class<- as.factor(model_dat$class)
temp<- as.data.frame(apply(expr[,2:ncol(expr)], 2, as.numeric))
expr_clean<- cbind(expr$GenotypeID, temp)
colnames(expr_clean)[1]<- c("GenotypeID")

#####
# Data Imputation
#####
# Check meta has no missing values
gg_miss_var(model_dat)

# Calculate missing values
res<-summary(aggr(expr_clean[,-c(1)], sortVar=TRUE))

# Extract all column names
# String cut at first "_"
# Reduce to unique list
# Loop over prefixes
  # Extract column index in full df for those prefixes
  ## "starts_with()" dplyr command
  # Extract columns from full df
  # Test for no missing values
    # None missng, contine to next prefix
  # Test for missing value across all columns
    # If true, remove all columns from full df
  # Test for missing value in some columns
    # If true fill with median of other columns
    # Replace columns in full df with new imputed columns

cols<- colnames(expr_clean[,c(-1)])
prefix<- sapply(strsplit(cols, "_", fixed= T),`[`, 1)
prefix<- unique(prefix)
for(curr in prefix){
  # Extract columns and indices
  probes<- select(expr_clean, starts_with(curr))
  probe_index<- which(colnames(expr_clean) %in% colnames(probes))
  
  # Test for all samples complete
  if(all(complete.cases(probes))){
    next
  }
  # Test if any samples missing across whole probe, if True, drop
  else if(any(rowSums(is.na(probes)) == ncol(probes))){
    expr_clean<- expr_clean[,-probe_index]
  }
  # If full probe set not missing, impute
  else{
    # Which samples incomplete
    missing<- which(complete.cases(probes)==F)
    medians<- apply(probes[missing,], 1, median, na.rm=T)
    
    # Replace missing value with median
    for(counter in 1:length(missing)){
      colMiss<- which(is.na(probes[missing[counter],]))
      probes[missing[counter], colMiss]<- medians[counter]
    }
    
    # Replace original columns with probes
    expr_clean[,probe_index]<- probes
  }
}

# Save object for modelling
saveRDS(model_dat, "../clean_data/model_data.RDS")
saveRDS(expr_clean, "../clean_data/expression.RDS")









