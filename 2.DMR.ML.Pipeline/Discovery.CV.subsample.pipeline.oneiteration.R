#!/usr/bin/env Rscript
library(plyr)
library(ROCR)
library(edgeR)
library(rstatix)
library(DESeq2)
library(sva)
library(matrixTests)
library(caret)

#####functions to use later######
autosome_filt = function(medips.count_df) {
  tmp =gsub(':.*','',medips.count_df)
  return_df = medips.count_df[!tmp %in% c('chrY','chrX','chrM')]
  return(return_df)
}

tpr.fpr.calc = function(x){
  tmp1 = x
  tmp1$f = as.integer(ifelse(tmp1$reported == 'control', 0, 1))
  tmp1$f.str = tmp1$reported
  
  tmp1 = tmp1[order(-tmp1$methylation_score),]
  case_no = nrow(tmp1[tmp1$reported!='control',])
  control_no = nrow(tmp1[tmp1$reported=='control',])
  auc.df =data.frame(matrix(nrow = 0, ncol=3))
  for(l in 1:nrow(tmp1)) {
    x = tmp1[1:l,]
    case_cum = nrow(x[x$reported!='control',])
    control_cum = nrow(x[x$reported=='control',])
    tpr = case_cum/case_no
    fpr = control_cum/control_no
    return.tmp = data.frame(l,tpr,fpr)
    auc.df = rbind(auc.df, return.tmp)
  }
  base = data.frame(l = 0, tpr = 0, fpr = 0)
  return(auc.df)
  #auc.df = auc.df[auc.df$fpr > 0,]
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

auc_calc = function(prediction_table,labels = c('control','cancer')) {
  tmp = prediction_table
  tmp = tmp[order(-tmp$methylation_score),]
  #labels = c(tmp$reported)
  pred = prediction(predictions = c(tmp$methylation_score) ,labels =  tmp$reported, labels)
  perf_AUC=performance(pred,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  return(AUC)
}


####else#####
cpg_count = readRDS('hg38_cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
cpg_count = cpg_count[cpg_count$count > 5,] #selecting windows with at least 6 or mor CpG sites

regulatory = readRDS('hg38.window300.promoter.enhancer.RDS')
cgi=readRDS('hg38.cig.info.300.RDS')
blood.wgbs = readRDS('hg38.wbc.ihec.300.mean.RDS')
blood.remove = unique(unlist(lapply(blood.wgbs[1:3], function(x) x[x$pct > 0.1 & x$cell_mean > 0.25,'window'] )))
blacklist = readRDS('hg38-blacklist.v2_window300.RDS')
repeat.element = readRDS('hg38.repeat.window.300.filter.elements.cut.RDS')

#dmr
dmrtable_dir='/dmr/output/dir/'

dir.create(dmrtable_dir, recursive = T)
setwd(dmrtable_dir)
#
medips.count_df = readRDS('brca.count_matrix.RDS') #reading in sample/region count matrix 
all_train_test_samples = readRDS('sample.information.RDS')
all_train_test_samples = all_train_test_samples[all_train_test_samples$Set == 'Discovery',]
all_train_test_samples = all_train_test_samples[all_train_test_samples$GRP_Id %in% colnames(medips.count_df.filt),]
medips.count_df = medips.count_df[,all_train_test_samples$GRP_Id]
#####DMR calling######
seedno=10
set.seed(seedno)
diagnosis_samples_splits= train_test_partition_cv(all_train_test_samples, splits = 10, seedno)  #selecting case samples only to split into training and test sets

deseq2.list = list()
for (fold in 1:10) {
  group_data_filt_matched_samples = diagnosis_samples_splits[diagnosis_samples_splits$fold != fold,]
  medips.count_df.filt = medips.count_df[,colnames(medips.count_df) %in% c(group_data_filt_matched_samples$GRP_Id)] #filtering medips matrix for train samples
  
  group_data_filt_matched_samples  = group_data_filt_matched_samples[order(match(group_data_filt_matched_samples$GRP_Id, colnames(medips.count_df.filt))),] #matching order for sample information df with count matrix
  
  medips.count_df.filt = medips.count_df[,colnames(medips.count_df) %in% c(group_data_filt_matched_samples$GRP_Id)] #filtering medips matrix for 
  
  #converting covariates to factors
  group_data_filt_matched_samples$Batch = factor(group_data_filt_matched_samples$Batch, levels = unique(group_data_filt_matched_samples$Batch))
  group_data_filt_matched_samples$Sex = factor(group_data_filt_matched_samples$Sex, levels = c('Female','Male'))
  group_data_filt_matched_samples$group = factor(group_data_filt_matched_samples$group, levels = c('control','cancer'))
  #filtering regions
  region_sum_counts = rowSums(medips.count_df.filt)
  medips.count_df.filt = medips.count_df.filt[which(region_sum_counts >0),]
  medips.count_df.filt = medips.count_df.filt[rownames(medips.count_df.filt) %in% cpg_count$window,]
  medips.count_df.filt = medips.count_df.filt[!rownames(medips.count_df.filt) %in% blacklist$window,]
  medips.count_df.filt = medips.count_df.filt[rownames(medips.count_df.filt) %in% c(regulatory$window, cgi$window,repeat.element$window),]
  medips.count_df.filt = autosome_filt(medips.count_df.filt)
  #
  dds <- DESeqDataSetFromMatrix(countData = medips.count_df.filt,
                                colData = group_data_filt_matched_samples[,-c(1:2)],
                                design= ~ Batch + age_group + group) 
  dds$condition = group_data_filt_matched_samples$group
  
  #surrogate variable analysis
  dds <- estimateSizeFactors(dds) #estimating size factors for tmm normalization
  dat <- counts(dds, normalized=TRUE) #getting tmm-normalised counts
  dat = dat[rowSums(dat) > 0,] #removing rows with no signals
  mod <- model.matrix(~ Batch, coldata) 
  mod0 <- model.matrix(~ 1, coldata)
  sv_count = num.sv(dat, mod) #calculating number of significant surrogate variables
  sv_count=2
  svseq <- svaseq(dat, mod, mod0, n.sv=sv_count) #modelling batch effect with sva
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
  deseq2.list[[length(deseq2.list) + 1]] = data.frame(ressva)
}

saveRDS(deseq2.list,paste0('brca.predx.deseq2.seed.',seedno,'dmr.RDS'))
saveRDS(year2_diagnosis_samples_splits,paste0('brca.predx.sample.split.seed.',seedno,'.RDS'))

######ML testing####
library(gridExtra)
library(caret)
library(ROCR)
library(pROC)
matrix = readRDS('deseq2.normalized.matrix.RDS')
deseq2.list = readRDS(paste0('brca.predx.deseq2.seed.',seedno,'dmr.RDS'))
samples_df_folds = readRDS(paste0('brca.predx.sample.split.seed.',seedno,'.RDS'))



for (fold in 1:10) {
  
  train.samples = samples_df_folds[samples_df_folds$fold != fold,]
  test.samples = samples_df_folds[samples_df_folds$fold == fold,]
  train.matrix = matrix[,train.samples$GRP_Id]
  test.matrix = matrix[,test.samples$GRP_Id]
  
  regions = rownames(train.matrix)
  ressva = deseq2.list[[fold]]
  
  
  deseq2.sig = ressva[ressva$pvalue < 0.05 & ressva$log2FoldChange > 0.5 ,]
  deseq2.sig = deseq2.sig[order(deseq2.sig$pvalue),]
  targ.features = rownames(deseq2.sig)
  feature_size = seq(50,500,50)
  plot.list = list()
  
  model.performance.rf.rank = function(feature_size, type = 'deseq2',targ.features,savedir) {
    performance.summary = list()
    performance.summary.logreg = list()
    
    for (f in feature_size) {
      feature_length = f
      dirsave = paste0(savedir,type,'/rf/',feature_length,'.hyper/performance/')
      dir.create(dirsave, recursive = T)
      
      #
      features =  autosome_filt(targ.features)[1:min(f,length(targ.features))]#
      train_set_matrix = matrix_setup(matrix, train.samples$GRP_Id, features, train.samples) 
      test_set_matrix = matrix_setup(matrix, test.samples$GRP_Id, features, test.samples)
      
      
      targ_features =features
      train_set_matrix = train_set_matrix[,colnames(train_set_matrix) %in% c('group',targ_features)]
      test_set_matrix = test_set_matrix[,colnames(test_set_matrix) %in% c('group',targ_features)]
      
      #rf
      control = trainControl(method = 'repeatedcv', number = 10, repeats = 10, search = 'random', classProbs = TRUE, summaryFunction = twoClassSummary)
      tunegrid = expand.grid(.mtry = seq(5:50))
      metric = 'Accuracy'
      
      ntrees = c(250,500,750,1000)
      ntree_list = list()
      acc = c()
      train_set_matrix$group = factor(train_set_matrix$group, levels = c('control','cancer'))
      rf_model = train(group ~ ., data = train_set_matrix, method = 'rf', metric = metric, tuneGrid = tunegrid, trControl = control, ntrees = ntrees)
      predictions = predict(rf_model, test_set_matrix[,-1])
      predictions_prob = predict(rf_model, test_set_matrix[,-1], type = 'prob')
      prediction_table = data.frame(GRP_Id= rownames(test_set_matrix), predictions = predictions, reported = test_set_matrix$group, methylation_score = predictions_prob[,colnames(predictions_prob) != 'control'], model = 'rf')
      prediction_table = prediction_table[order(-prediction_table$methylation_score),]
      prediction_table$features = length(targ_features)
      prediction_table$auroc = auc_calc(prediction_table,labels = c('control','cancer'))
      prediction_table$model = 'RF'
      prediction_table$type = type
      saveRDS(prediction_table, paste0(dirsave, 'rf','_seed.',seedno,'_fold.',fold,'.RDS'))
      
      performance.summary[[length(performance.summary) + 1]] = prediction_table
      
    }
    return(performance.summary)
  }
  
  targ.filt = function(x) {
    a = x[x %in% cpg_count$window]
    a = a[!a %in% blacklist$window]
    a = a[a %in% regulatory$window | a %in% cgi$window | a %in% repeat.element$window]
    a = a[!a %in% blood.remove]
    
    return(a)
  }
  
  savedir='/output/for/ml/performance/'
  
  dir.create(savedir,recursive = T)
  setwd(savedir)
  targ.features = rownames(deseq2.sig)
  targ.features = targ.filt(targ.features)
  feature_size = seq(25,500,25)
  model.performance.rf.rank(feature_size, type = 'deseq2.nb',targ.features,savedir)
  
  
  #
}




