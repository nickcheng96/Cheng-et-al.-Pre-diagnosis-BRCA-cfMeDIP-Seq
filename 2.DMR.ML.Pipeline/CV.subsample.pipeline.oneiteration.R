#!/usr/bin/env Rscript
#libraries to load (you will need to install DESeq2 and sva)
library(plyr)
library(ROCR)
library(edgeR)
library(rstatix)
library(netDx)
library(DESeq2)
library(sva)
#####functions to use later######
autosome_filt = function(medips.count_df) {
  tmp = gsub(':.*','',rownames(medips.count_df))
  return_df = medips.count_df[-which(tmp %in% c('chrX','chrY')),]
  return(return_df)
}
train_test_partition_cv = function(year2_diagnosis_samples, splits = 10, seed = 10) { 
  set.seed(seed)
  year2_diagnosis_samples$diagnosis_time_group =year2_diagnosis_samples$Diagnosis_Time
  tmp_list = split(year2_diagnosis_samples, year2_diagnosis_samples$Diagnosis_Time) #split into list according to diagnosis time group
  tmp_list_train = lapply(tmp_list, function(x) {
    yourData<-x[sample(nrow(x)),]
    #Create 10 equally size folds
    folds <- cut(seq(1,nrow(yourData)),breaks=splits,labels=FALSE)
    yourData$fold = folds
    return(yourData)
  })
  return_df = do.call('rbind',tmp_list_train)
  return(return_df)
}

matrix_setup = function(all_count_matrix, train_samples, features, sample_information_df) {
  train_df = all_count_matrix[rownames(all_count_matrix) %in% features,]
  train_df = train_df[,colnames(train_df) %in% train_samples]
  train_df_t = data.frame(t(train_df), check.names = F)
  group_dat = sample_information_df[sample_information_df$GRP_Id %in% train_samples,]
  group_dat = group_dat[order(match(group_dat$GRP_Id,rownames(train_df_t))),]
  train_df_t$GRP_Id = rownames(train_df_t)
  return_df = merge(group_dat[,c('GRP_Id','group')], train_df_t, by = 'GRP_Id', all.y = T)
  rownames(return_df) = return_df$GRP_Id
  return(return_df[,-1])
}


auc_calc = function(prediction_table) {
  tmp = prediction_table
  tmp = tmp[order(-tmp$methylation_score),]
  labels = c(tmp$reported)
  pred = prediction(predictions = c(tmp$methylation_score) ,labels =  labels, label.ordering = c(2,1))
  perf_AUC=performance(pred,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  return(AUC)
}
rfe_calc = function(train_set, size_estimate = 100) {
  control <- rfeControl(functions=rfFuncs, method="cv", number=10)
  size = min(size_estimate, ncol(train_set[,-1]))
  rfe_model <- rfe(train_set[,-1], train_set$group, sizes=size, rfeControl=control, metric = "Accuracy")
  rfe_features = rfe_model$variables
  targ_size = min(size, ncol(train_set[,-1]))
  rfe_features = rfe_features[rfe_features$Variables == size,]
  rfe_features = rfe_features[order(-rfe_features$Overall),]
  features = unique(rfe_features$var)[1:targ_size]
  return(features)
}

auc_calc = function(prediction_table) {
  tmp = prediction_table
  tmp = tmp[order(-tmp$methylation_score),]
  labels = c(tmp$reported)
  pred = prediction(predictions = c(tmp$methylation_score) ,labels =  labels, label.ordering = c(2,1))
  perf_AUC=performance(pred,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  return(AUC)
}

####else#####

cpg_count = readRDS('cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
cpg_count = cpg_count[cpg_count$count > 5,] #selecting windows with at least 6 or mor CpG sites

combined_tissue_summary = readRDS('blood_wgbs_summary.RDS')
blood_wgbs_windows = as.character(combined_tissue_summary[combined_tissue_summary$Blood_Cells >= 0.4,'window'])

fantom_regulatory_information = readRDS('fantom_regulatory_information_300.RDS')
regulatory = fantom_regulatory_information[fantom_regulatory_information$CpG_Region != 'null' | fantom_regulatory_information$Regulatory != 'None','window']


#setting file directories
medips.count_dir='/Path/to/counts/' #directory containing sample/region counts
dmrtable_dir='/Path/to/save/DMR/tables/'
dir.create(dmrtable_dir, recursive = T)

setwd(dmrtable_dir) #setting work directory
medips.count_df = readRDS(paste0(medips.count_dir, 'brpa.extended_v2_medips.count.RDS')) #reading in sample/region count matrix 
all_train_test_samples = readRDS('brca.sample.information.RDS')
####
medips.count_df.filt = medips.count_df[,colnames(medips.count_df) %in% c(all_train_test_samples$GRP_Id)] #removing samples from count matrix to only keep samples of interest
all_train_test_samples = all_train_test_samples[all_train_test_samples$GRP_Id %in% colnames(medips.count_df.filt),]
all_train_test_samples  = all_train_test_samples[order(match(all_train_test_samples$GRP_Id, colnames(medips.count_df.filt))),] #match order of samples between count matrix and data information df


#####DMR calling######
set.seed(seedno)
year2_diagnosis_samples_splits= train_test_partition_cv(all_train_test_samples, splits = 10, seedno)  #selecting case samples only to split into training and test sets

dmr_list =list()

for (fold in 1:10) {
  group_data_filt_matched_samples = year2_diagnosis_samples_splits[year2_diagnosis_samples_splits$fold != fold,]
  medips.count_df.filt = medips.count_df[,colnames(medips.count_df) %in% c(group_data_filt_matched_samples$GRP_Id)] #filtering medips matrix for train samples
  
  group_data_filt_matched_samples  = group_data_filt_matched_samples[order(match(group_data_filt_matched_samples$GRP_Id, colnames(medips.count_df.filt))),] #matching order for sample information df with count matrix
  
  #group_data_filt_matched_samples = matchcontrols(year2_diagnosis_train_samples, group_col = 'group', cases = 'breast_cancer', controls = 'control', quant = 0.1,id.var = 'GRP_Id') #selecting matched controls for training cases
  medips.count_df.filt = medips.count_df[,colnames(medips.count_df) %in% c(group_data_filt_matched_samples$GRP_Id)] #filtering medips matrix for 
  
  #converting covariates to factors
  group_data_filt_matched_samples$Batch = factor(group_data_filt_matched_samples$Batch, levels = unique(group_data_filt_matched_samples$Batch))
  group_data_filt_matched_samples$Sex = factor(group_data_filt_matched_samples$Sex, levels = c('Female','Male'))
  
  #
  
  #filtering regions
  region_sum_counts = rowSums(medips.count_df.filt)
  medips.count_df.filt = medips.count_df.filt[which(region_sum_counts >0),]
  medips.count_df.filt = medips.count_df.filt[rownames(medips.count_df.filt) %in% cpg_count$window,]
  medips.count_df.filt = medips.count_df.filt[!rownames(medips.count_df.filt) %in% blood_wgbs_windows,]
  medips.count_df.filt = medips.count_df.filt[rownames(medips.count_df.filt) %in% regulatory,]
  #
  dds <- DESeqDataSetFromMatrix(countData = medips.count_df.filt,
                                colData = group_data_filt_matched_samples[,-c(1:2)],
                                design= ~ Batch + age_group + group) #can add gender here but we're only looking at female samples currently
  dds$condition = group_data_filt_matched_samples$group
  
  #surrogate variable analysis
  #setwd(wkdir)
  dds <- estimateSizeFactors(dds) #estimating size factors for tmm normalization
  dat <- counts(dds, normalized=TRUE) #getting tmm-normalised counts
  dat = dat[rowSums(dat) > 0,] #removing rows with no signals
  
  coldata = data.frame(colData(dds)[,c('group','age_group','Batch')]) #again can add gender but we're only looking at female samples
  mod <- model.matrix(~ Batch, coldata) 
  mod0 <- model.matrix(~ 1, coldata)
  sv_count = num.sv(dat, mod) #calculating number of significant surrogate variables
  svseq <- svaseq(dat, mod, mod0, n.sv=sv_count) #modelling batch effect with sva
  # output:
  ddssva <- dds
  
  #updating model with surrogate variables
  if (sv_count == 1) {
    ddssva$SV1 <- svseq$sv[,1]
    design(ddssva) <- ~ SV1  + age_group + group
  } else if (sv_count == 2) {
    ddssva$SV1 <- svseq$sv[,1]
    ddssva$SV2 <- svseq$sv[,2]
    design(ddssva) <- ~ SV1 + SV2 + age_group + group
  } else if (sv_count == 3) {
    ddssva$SV1 <- svseq$sv[,1]
    ddssva$SV2 <- svseq$sv[,2]
    ddssva$SV3 <- svseq$sv[,3]
    design(ddssva) <- ~ SV1 + SV2 + SV3 + age_group + group
  }
  
 ddssva <- DESeq(ddssva) #differential methylation analysis 
  ressva <- results(ddssva) #generating results table 
  ressva = ressva[complete.cases(ressva$pvalue),] #removing NAs 
  ressva$logFC = ressva$log2FoldChange
  dmr_list[[fold]] = ressva
}

saveRDS(year2_diagnosis_samples_splits,paste0(dmrtable_dir, 'sample_split_read.readno','_year.all_seed.',seedno,'.RDS'))
saveRDS(dmr_list, paste0(dmrtable_dir, 'dmr_table_list_read.readno','_year.all_seed.',seedno,'.RDS'))

matrix = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/medips_files/300bp/medips.extended/breast_cpg.reg.blood.filt_read10_mans.gerdfilt_deseq2.norm.matrix.RDS') #normalised matrix


samples_df_folds = readRDS(paste0(dmrtable_dir, 'sample_split_read.',10,'_year.all_seed.',seedno,'.RDS'))
dmr_fold_list = readRDS(paste0(dmrtable_dir, 'dmr_table_list_read.',10,'_year.all_seed.',seedno,'.RDS'))

#rf
for (fold in 1:length(dmr_fold_list)) {
  train_samples = samples_df_folds[samples_df_folds$fold != fold,]
  test_samples =  samples_df_folds[samples_df_folds$fold == fold,]
  dmr_table = dmr_fold_list[[fold]]
  dmr_table$window = rownames(dmr_table)
  
  #####random forest approach#####
  savedir='/Path/to/saving/files/'
  dir.create(savedir, recursive = T)
  
  
  train_set = train_samples
  train_set$group = factor(train_set$group, levels = c('breast_cancer','control'))
  test_set = test_samples
  test_set$group = factor(test_set$group, levels = c('breast_cancer','control'))
  all_count_matrix = matrix
  train_test_matrix = all_count_matrix[,colnames(all_count_matrix) %in% c(train_set$GRP_Id, test_set$GRP_Id)]
  matrix.mean = mean(rowMeans(train_test_matrix))
  
  sample_information_df = rbind(train_set,test_set)
  feature_filt = data.frame(dmr_table,stringsAsFactors = F)

  control_train = train_set[train_set$group == 'control',]
  control_train_matrix = all_count_matrix[,colnames(all_count_matrix) %in% control_train$GRP_Id]
  
  tmp = rowMeans(control_train_matrix)
  tmp1 = tmp[tmp < matrix.mean]
  cancer_train = train_set[train_set$group != 'control',]
  cancer_train_matrix = all_count_matrix[,colnames(all_count_matrix) %in% cancer_train$GRP_Id]
  feature_filt = feature_filt[rownames(feature_filt) %in% names(tmp1),] #
  
  hyper_sig_features = feature_filt[feature_filt$logFC > 0 & feature_filt$pvalue < 0.05,]
  hyper_sig_features = hyper_sig_features[order(hyper_sig_features$pvalue),]
  hypo_sig_features = feature_filt[feature_filt$logFC < 0 & feature_filt$pvalue < 0.05,]
  hypo_sig_features = hypo_sig_features[order(hypo_sig_features$pvalue),]
  combined_sig_features = rbind(hyper_sig_features, hypo_sig_features)
  combined_sig_features = combined_sig_features[order(combined_sig_features$pvalue),]
  
  
  
  #all_ids = all_ids[!all_ids %in% ]
  feature_set = seq(100,200,25)
  
  #tmp =do.call('rbind',lapply(rf_list, function(x)x[1,] ))
  #machine learning training with filtered features
  rf_list = list()
  
  #####450k array filter - beta#####
  
  rf_list_hyper = list()
  rf_list_balanced = list()
  glm_list_balanced = list()
  glm_list_hyper = list()
  
  library(caret)
  library(ROCR)
  
  #hyper only test
  for (feature_number in 1:length(feature_set)) {
    feature_length = feature_set[feature_number]
    dirsave = paste0(savedir,'/rf/',feature_length,'.hyper/performance/')
    dir.create(dirsave, recursive = T)
    features = hyper_sig_features$window
    train_samples = train_set$GRP_Id
    test_samples = test_set$GRP_Id
    train_set_matrix = matrix_setup(all_count_matrix, train_samples, features, train_set) 
    test_set_matrix = matrix_setup(all_count_matrix, test_samples, features, test_set)

    
    targ_features = hyper_sig_features$window[1:feature_length]
    
    train_set_matrix = train_set_matrix[,colnames(train_set_matrix) %in% c('group',targ_features)]
    test_set_matrix = test_set_matrix[,colnames(test_set_matrix) %in% c('group',targ_features)]
    control = trainControl(method = 'repeatedcv', number = 10, repeats = 10, search = 'random', classProbs = TRUE, summaryFunction = twoClassSummary)
    mtry = 100
    tunegrid = expand.grid(.mtry = 100)
    metric = 'Accuracy'
    
    ntrees = c(50,250,500,1000)
    ntree_list = list()
    acc = c()
    for (tree in 1:length(ntrees)) {
      rf_model = train(group ~ ., data = train_set_matrix, method = 'rf', metric = metric, tuneGrid = tunegrid, trControl = control, ntrees = ntrees[tree])
      predictions = predict(rf_model, test_set_matrix[,-1])
      predictions_prob = predict(rf_model, test_set_matrix[,-1], type = 'prob')
      prediction_table = data.frame(GRP_Id= rownames(test_set_matrix), predictions = predictions, reported = test_set_matrix$group, methylation_score = predictions_prob[,colnames(predictions_prob) != 'control'], seed = seedno)
      prediction_table = prediction_table[order(-prediction_table$methylation_score),]
      prediction_table$features = length(targ_features)
      prediction_table$auroc = auc_calc(prediction_table)
      print(confusionMatrix(predictions, test_set_matrix$group) )
      prediction_table$trees = ntrees[tree]
      prediction_table$array_filt = 'balanced'
      prediction_table$model = 'RF'
      train_performance = getTrainPerf(rf_model)
      prediction_table$TrainROC = train_performance$TrainROC
      prediction_table$TrainSens = train_performance$TrainSens
      prediction_table$TrainSpec = train_performance$TrainSpec
      pred_table = (confusionMatrix(predictions, test_set_matrix$group))
      prediction_table$Sensitivity_overall = pred_table$byClass[1]
      prediction_table$Specificity_overall = pred_table$byClass[2]
      prediction_table = merge(prediction_table, sample_information_df[,c('GRP_Id','Diagnosis_Time')], all.x = T)
      
      prediction_table = time_breakdown_performance(prediction_table, sample_information_df, control_report = F, cancer_group = 'breast_cancer')
      ntree_list[[tree]] = prediction_table
      acc = c(acc,prediction_table$auroc[1])
    }
    return_df = ntree_list[[which(acc == max(acc))[1]]]
    saveRDS(return_df, paste0(dirsave, 'rf_top.hyper_read.10','_seed.',seedno,'_fold.',fold,'.RDS'))
    saveRDS(features, paste0(dirsave, 'rf_top.hyper_read.10','_seed.',seedno,'_fold.',fold,'.features.RDS'))
  }
  
  
  
}
