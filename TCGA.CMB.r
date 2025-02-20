#!/usr/bin/env Rscript
# library(glmnet)
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

rm(list = ls())

###############################################################################################################
# define RF CV function
###############################################################################################################
randomForest_CrossVal <- function(data, formula, fold=5){
  # config the experiments
  num_sample = dim(data)[1]
  foldsize = floor(num_sample/fold)
  expidx = matrix(fold, nrow=num_sample,ncol=1)
  tmp = rep(c(1:fold),foldsize)
  expidx[c(1:length(tmp))]=tmp
  # expidx = permute(expidx) # may want to permute the idx
  expidx = as.data.frame(expidx)
  colnames(expidx)=c('foldID')
  results_holdout = c()
  for(fid in c(1:fold)){
    data_train = data[expidx$foldID!=fid,]
    data_test = data[expidx$foldID==fid,]
    rf <- randomForest(formula, data=data_train, proximity=FALSE)
    saveRDS(rf,file=paste(result_root_sub,'/','rfmodel-fold',fid,'.rds',sep = ''))
    print(rf)
    data_test$lbl_rf_prediction <- predict(rf, data_test)
    results_holdout = rbind(results_holdout,data_test)
    importance_scores = importance(rf)
    feature_importance = data.frame(Feature=rownames(importance_scores), Importance = importance_scores[,1])
    importance_plot = ggplot(feature_importance, aes(x = reorder(Feature, Importance), y = Importance)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Features", y = "Importance", title = "Feature Importance in Random Forest Model")
    pdf(paste(result_root_sub,'/','importance-fold',fid,'.pdf',sep = ''), width=45, height=10)
    print(importance_plot)
    dev.off()
  }
  return(results_holdout)
}

###############################################################################################################
# load data
###############################################################################################################
result_root = './Results'
dir.create(result_root)
cmb_data = read.xlsx('./Data/TCGA.Prostate.CMB.xlsx',rowNames = T)
path_data = read.xlsx('./Data/sig_res_sig_heatmap.xlsx',rowNames = T)
cmb_tme_data = merge(cmb_data,path_data,by='row.names')
row.names(cmb_tme_data)=cmb_tme_data$Row.names
cmb_tme_data$Row.names=NULL

cmb_names = colnames(cmb_data)[1:13]
tme_names = colnames(path_data)
cmbs = cmb_names
result_root = paste(result_root,'/CMB_Path_Regression_ModelSaved_TCGABoxplots','/',sep='')
dir.create(result_root)

###############################################################################################################
# settings
###############################################################################################################

features=tme_names

##############################################################################################################################
#fit lasso on continouse variables (dose, timepoint)
##############################################################################################################################
nfold = 5
# write.csv(raman_data_pigment_complete,paste(result_root,'/',performer,'_phase',phase,'_',modality,'-complete','.csv',sep = ''))

for(ii in c(1:length(features))){
  
  fname = features[ii]
  result_root_sub = paste(result_root,'/',fname,'_nfold',nfold,sep = '')
  dir.create(result_root_sub)
  
  ## random forest
  formula_str = paste(fname,'~',paste(cmb_names,collapse = '+'),sep = '')
  
  cmb_tme_data_complete = cmb_tme_data[complete.cases(cmb_tme_data[,c(fname,cmb_names)]), c(fname,cmb_names,'TreatmentResistancePredictionCategory')]
  raman_data_pigment_complete_inuse = randomForest_CrossVal(cmb_tme_data_complete,as.formula(formula_str),nfold)
  
  ## plot scatter plot
  
  data_plot = c()
  data_plot$lbl_prediction = raman_data_pigment_complete_inuse$lbl_rf_prediction
  data_plot$lbl_groundtruth = raman_data_pigment_complete_inuse[,fname]
  data_plot = as.data.frame(data_plot)
  
  scatter_plot = ggplot(data_plot,aes(lbl_prediction,lbl_groundtruth))+
    geom_point()+
    stat_smooth(method = 'lm',formula = y~x)+
    ylab(fname)+
    xlab('Randomforest_prediction')+
    stat_cor(method='spearman')
  pdf(paste(result_root_sub,'/scatterplot_',fname,'_lblPrediction.pdf',sep = ''))
  print(scatter_plot)
  dev.off()
  tiff(paste(result_root_sub,'/scatterplot_',fname,'_lblPrediction.tif',sep = ''))
  print(scatter_plot)
  dev.off()
  
  ## plot boxplot in TCGA between INR-like and ER-like groups
  bplot = ggboxplot(raman_data_pigment_complete_inuse, x = "TreatmentResistancePredictionCategory", y = fname,
                    outlier.colour="black",
                    fill = "TreatmentResistancePredictionCategory")+
    scale_color_manual(values=c("#EFC000FF", "#0073C2FF")) +
    scale_fill_manual(values=c("#EFC000FF", "#0073C2FF")) +
    stat_compare_means()
  pdf(paste(result_root_sub,'/','TCGA_PRAD_',fname,'_ERvsINR.pdf',sep = ''))
  print(bplot)
  dev.off()
}

