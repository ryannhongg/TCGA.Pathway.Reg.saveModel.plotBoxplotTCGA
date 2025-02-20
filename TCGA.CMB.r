library(stringr)
library(ggpubr)
library(rstatix)
library(umap)
library(ggrepel)
library(randomForest)
library(caret)
library(gtools)
library(tidyr)
library(dplyr)
library(openxlsx)

# Clear previous variables
rm(list = ls())

###############################################################################################################
# Define Random Forest Cross Validation Function
###############################################################################################################
randomForest_CrossVal <- function(data, formula, fold=5) {
  num_sample = dim(data)[1]
  foldsize = floor(num_sample / fold)
  expidx = matrix(fold, nrow=num_sample, ncol=1)
  tmp = rep(1:fold, length.out = num_sample)
  expidx[c(1:length(tmp))] = tmp
  expidx = as.data.frame(expidx)
  colnames(expidx) = c('foldID')
  
  results_holdout = c()
  for (fid in c(1:fold)) {
    data_train = data[expidx$foldID != fid, ]
    data_test = data[expidx$foldID == fid, ]
    
    rf <- randomForest(formula, data=data_train, proximity=FALSE)
    saveRDS(rf, file=paste(result_root_sub, '/', 'rfmodel-fold', fid, '.rds', sep=''))
    
    print(rf)
    data_test$lbl_rf_prediction = predict(rf, data_test)
    results_holdout = rbind(results_holdout, data_test)
    
    importance_scores = importance(rf)
    feature_importance = data.frame(Feature=rownames(importance_scores), Importance=importance_scores[, 1])
    
    importance_plot = ggplot(feature_importance, aes(x = reorder(Feature, Importance), y = Importance)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Features", y = "Importance", title = "Feature Importance in Random Forest Model")
    
    pdf(paste(result_root_sub, '/', 'importance-fold', fid, '.pdf', sep=''), width=45, height=10)
    print(importance_plot)
    dev.off()
  }
  return(results_holdout)
}

###############################################################################################################
# Load Data
###############################################################################################################
result_root = './Results'
dir.create(result_root)
versioned_dir = paste(result_root, '/TCGA_KIRC_CMB_Analysis', sep='')  # New folder for this analysis
dir.create(versioned_dir, showWarnings = FALSE)

# Load the data
cmb_data = read.xlsx('./Data/TCGA.KIRC.73CMB.Manuscript.xlsx', rowNames = T)
gene_data = read.csv('./Data/TCGA_KIRC_gene_data_primary_average.csv', row.names = 1)
# In this case, assuming there are only TCGA and CMB columns for features
cmb_names = colnames(cmb_data)
cmb_gene_data = merge(cmb_data, gene_data, by="row.names")
rownames(cmb_gene_data) = cmb_gene_data$Row.names

###############################################################################################################
# Settings
###############################################################################################################
features = colnames(gene_data)
nfold = 10

#features = c("MTOR")
summary_data = matrix(NA, length(features), 2)
rownames(summary_data) = features
colnames(summary_data) = c("r", "p")
summary_data = as.data.frame(summary_data)
start_gene = "TMEM133"
start_index = which(features == start_gene)

# Loop over features starting from the "TMEM133" gene
for (ii in c(start_index:length(features))) {
  
  # Define the feature name and the directory for results
  fname = features[ii]
  print(fname)
  
  # Define the formula for the random forest model
  formula_str = paste(fname, '~', paste(cmb_names, collapse = '+'), sep='')
  
  # Preprocess data: only rows where feature and CMB columns are complete
  cmb_tme_data_complete = cmb_gene_data[complete.cases(cmb_gene_data[, c(fname, cmb_names)]), c(fname, cmb_names)]
  if (dim(cmb_tme_data_complete)[1] < nfold) {
    next
  }
  result_root_sub = paste(versioned_dir, '/', fname, '_nfold', nfold, sep='')
  dir.create(result_root_sub, showWarnings = FALSE)
  
  # Run the random forest cross-validation function
  raman_data_pigment_complete_inuse = randomForest_CrossVal(cmb_tme_data_complete, as.formula(formula_str), nfold)
  
  # Scatter plot between predictions and actual values
  data_plot = data.frame(lbl_prediction = raman_data_pigment_complete_inuse$lbl_rf_prediction,
                         lbl_groundtruth = raman_data_pigment_complete_inuse[, fname])
  
  scatter_plot = ggplot(data_plot, aes(lbl_prediction, lbl_groundtruth)) +
    geom_point() +
    stat_smooth(method = 'lm', formula = y ~ x) +
    ylab(fname) +
    xlab('Random Forest Prediction') +
    stat_cor(method = 'spearman')
  
  pdf(paste(result_root_sub, '/scatterplot_', fname, '_lblPrediction.pdf', sep=''))
  print(scatter_plot)
  dev.off()
}

