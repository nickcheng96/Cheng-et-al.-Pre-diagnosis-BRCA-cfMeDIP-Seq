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


feature_selection = function(dmr_table, tissue_filter_windows, remove_prob, keep_prob, option = optionno, cpg_count, train_samples, all_count_matrix, sample_information_df, size= 100) {
  matrix_setup = function(all_count_matrix, train_samples, features, sample_information_df) {
    train_df = all_count_matrix[rownames(all_count_matrix) %in% features,]
    train_df = train_df[,colnames(train_df) %in% train_samples]
    train_df_t = data.frame(t(train_df), check.names = F)
    group_dat = sample_information_df[sample_information_df$GRP_Id %in% train_samples,]
    group_dat = group_dat[order(match(rownames(train_df_t), group_dat$GRP_Id)),]
    return_df = data.frame(cbind(group = group_dat$group, train_df_t), check.names = F)
    return(return_df)
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
  
  
  filtered_dmr_table = autosome_filt(dmr_table)
  cpg_filt = cpg_count[cpg_count$count > 5,]
  filtered_dmr_table = filtered_dmr_table[rownames(filtered_dmr_table) %in% cpg_filt$window,]
  #filtering for probes to keep
  if(keep_prob == 'None') {
    filtered_dmr_table = filtered_dmr_table
  } else {
    filtered_dmr_table = filtered_dmr_table[rownames(filtered_dmr_table) %in% tissue_filter_windows[[1]],]
  }
  
  if (remove_prob == 'None') {
    filtered_dmr_table = filtered_dmr_table
  } else {
    filtered_dmr_table = filtered_dmr_table[!rownames(filtered_dmr_table) %in% tissue_filter_windows[[2]],]
  }
  
  if (option == 1) { #all significant DMRs
    sig_dmr = filtered_dmr_table[filtered_dmr_table$pvalue < 0.05 & abs(filtered_dmr_table$log2FoldChange) > 0.7,]
    sig_dmr = rownames(sig_dmr[order(sig_dmr$pvalue),])
    if (length(sig_dmr) > 1) {
      train_set = matrix_setup(all_count_matrix, train_samples, sig_dmr, sample_information_df)
      features = rfe_calc(train_set, size = size)
    } else {
      features = sig_dmr
    }
  } else if (option == 2) { #hyper DMRs only 
    sig_dmr = filtered_dmr_table[filtered_dmr_table$pvalue < 0.05 & (filtered_dmr_table$log2FoldChange) < -0.7,]
    sig_dmr = rownames(sig_dmr[order(sig_dmr$pvalue),])
    if (length(sig_dmr) > 1) {
      train_set = matrix_setup(all_count_matrix, train_samples, sig_dmr, sample_information_df)
      features = rfe_calc(train_set, size = size)
    } else {
      features = sig_dmr
    }
  } else if (option == 3) { #top 150 hypo and hyper
    hyper_sig_dmr = filtered_dmr_table[filtered_dmr_table$pvalue < 0.05 & (filtered_dmr_table$log2FoldChange) > 0.7,]
    hypo_sig_dmr = filtered_dmr_table[filtered_dmr_table$pvalue < 0.05 & (filtered_dmr_table$log2FoldChange) < -0.7,]
    hyper_sig_dmr = rownames(hyper_sig_dmr[order( hyper_sig_dmr$pvalue),])
    hypo_sig_dmr = rownames(hypo_sig_dmr[order(hypo_sig_dmr$pvalue),])
    
    if (length(hyper_sig_dmr) > 1 ) {
      hyper_train_set = matrix_setup(all_count_matrix, train_samples, hyper_sig_dmr, sample_information_df)
      hyper_features = rfe_calc(hyper_train_set, size = size/2)
    } else {
      hyper_features = hyper_sig_dmr
    }
    if (length(hypo_sig_dmr) > 1 ) {
      hypo_train_set = matrix_setup(all_count_matrix, train_samples, hypo_sig_dmr, sample_information_df)
      hypo_features = rfe_calc(hypo_train_set, size = size/2)
    } else {
      hypo_features = hypo_sig_dmr
    }
    features = c(hyper_features, hypo_features)
  }
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

cpg_count = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
cpg_count = cpg_count[cpg_count$count > 5,] #selecting windows with at least 6 or mor CpG sites

tissuedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/tissue_reference/'
combined_tissue_summary = readRDS(paste0(tissuedir,'blood_wgbs_summary.RDS'))
blood_wgbs_windows = as.character(combined_tissue_summary[combined_tissue_summary$Blood_Cells >= 0.4,'window'])

fantom_regulatory_information = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/medips_regulatory_window/fantom_regulatory_information_300.RDS')
regulatory = fantom_regulatory_information[fantom_regulatory_information$CpG_Region != 'null' | fantom_regulatory_information$Regulatory != 'None','window']


#setting file directories
medips.count_dir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/medips_files/300bp/medips.extended/' #directory containing sample/region counts
dmrtable_dir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/cancer_analysis/dmr_analysis/batch1_8.extended/all_samples_dmr/breast.all.bloodfilt.deseq2/read10.cv/dmr_tables.ensemble.90.10.gerd.filt.v3/'
dir.create(dmrtable_dir, recursive = T)
#/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/cancer_analysis/auc_ensembl/close_diagnosis/batch1_8/trimq_control/ensemble/breast/
setwd(dmrtable_dir) #setting wkdir

medips.count_dir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/medips_files/300bp/medips.extended/' #directory containing sample/region counts
medips.count_df = readRDS(paste0(medips.count_dir, 'brpa.extended_v2_medips.count.RDS')) #reading in sample/region count matrix 
all_train_test_samples = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/participant_data/hq.brca.sample.information.RDS')
####
medips.count_df.filt = medips.count_df[,colnames(medips.count_df) %in% c(all_train_test_samples$GRP_Id)] #removing samples from count matrix to only keep samples of interest
all_train_test_samples = all_train_test_samples[all_train_test_samples$GRP_Id %in% colnames(medips.count_df.filt),]
all_train_test_samples  = all_train_test_samples[order(match(all_train_test_samples$GRP_Id, colnames(medips.count_df.filt))),] #match order of samples between count matrix and data information df
#gerd = c('AIX_0035','AIX_0264','AIX_0252','AIX_0152','AIX_0058','AIX_0157')
#group_data_filt = group_data_filt[!group_data_filt$GRP_Id %in% gerd,]

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
  #saveRDS(dat, '/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/medips_files/300bp/medips.extended/breast_cpg_filt_read5_batch1_8_deseq2.bloodfilt.dgnormmatrix.RDS')
  dat = dat[rowSums(dat) > 0,] #removing rows with no signals
  
  coldata = data.frame(colData(dds)[,c('group','age_group','Batch')]) #again can add gender but we're only looking at female samples
  mod <- model.matrix(~ Batch, coldata) 
  mod0 <- model.matrix(~ 1, coldata)
  #sv_count = num.sv(dat, mod) #calculating number of significant surrogate variables
  sv_count=2
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
  
  #vsd <- vst(ddssva, blind=T, nsub = 10000)
  #saveRDS(assay(vsd), '/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/medips_files/300bp/medips.extended/breast_cpg_filt_read10_batch1_11_vst.all.RDS')
  ddssva <- DESeq(ddssva) #differential methylation analysis 
  ressva <- results(ddssva) #generating results table 
  ressva = ressva[complete.cases(ressva$pvalue),] #removing NAs 
  ressva$logFC = ressva$log2FoldChange
  dmr_list[[fold]] = ressva
}

saveRDS(year2_diagnosis_samples_splits,paste0(dmrtable_dir, 'sample_split_read.readno','_year.all_seed.',seedno,'.RDS'))
saveRDS(dmr_list, paste0(dmrtable_dir, 'dmr_table_list_read.readno','_year.all_seed.',seedno,'.RDS'))
