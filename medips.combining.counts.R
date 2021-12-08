####regular counts#####
#library
library(DESeq2)
library(plyr)
library(data.table)
#####functions#####
combined_counts_dt = function(medips.count_list, method = 'medips') {
  if (method == 'medips') {
    return_list = lapply(medips.count_list, function(x) {
      tmp = x
      tmp$window = paste0(tmp[,1],':',tmp[,2],'-',tmp[,3])
      tmp = tmp[tmp[,5] != 0 ,]
      return(tmp[,c(6,5)])
    })
    
    sample = names(medips.count_list)
    sample = gsub('_CM.*','',sample)
    datatable_list = list()
    for (i in 1:length(return_list)) {
      colnames(return_list[[i]])[2] = sample[i]
      datatable_list[[i]] = data.table(return_list[[i]], key = 'window')
    }
    
    for (i in 1:length(datatable_list)) {
      if (i == 1) {
        return_df = datatable_list[[i]]
      } else {
        return_df = merge(return_df, datatable_list[[i]], by = 'window', all = T)
      }
    }
    #rownames(return_df) = return_df[,1]
    return_df[is.na(return_df)] <- 0 
    return_df = data.frame(return_df, check.names = F, stringsAsFactors = F)
    rownames(return_df) = return_df$window
    return(return_df)
  }
}

#####new setup#####
#reference files
cpg_count = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
cpg_count5 = cpg_count[cpg_count$count > 5,] #selecting windows with at least 6 or mor CpG sites

tissuedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/tissue_reference/'
combined_tissue_summary = readRDS(paste0(tissuedir,'blood_wgbs_summary.RDS'))
blood_wgbs_windows = as.character(combined_tissue_summary[combined_tissue_summary$Blood_Cells >= 0.4,'window'])

fantom_regulatory_information = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/medips_regulatory_window/fantom_regulatory_information_300.RDS')
regulatory = fantom_regulatory_information[fantom_regulatory_information$CpG_Region != 'null' | fantom_regulatory_information$Regulatory != 'None','window']

#generating raw counts
medips.count_dir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/batchall/novaseq/novaseq_umitools_output_trimq/6_medips_qc/window_300/'
savedir = '/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/medips_files/300bp/medips.extended/'
dir.create(savedir,recursive = T)
medips.count_files = list.files(path = medips.count_dir, pattern = '*_medips_count.bed')
medips.count_list = lapply(medips.count_files, function(x) {
  return(read.table(paste0(medips.count_dir, x),
                    header = F,
                    stringsAsFactors = F,
                    quote = "",
                    sep = '\t'))
} )

names(medips.count_list) = gsub('_medips_count.bed','', medips.count_files)
#combining to single matrix
medips.count_df = combined_counts_dt(medips.count_list, method = 'medips')
saveRDS(medips.count_df, paste0(savedir, 'brca.raw.medips.count_matrix.RDS'))

medips.count_df.filtered = autosome_filt(medips.count_df)
medips.count_df.filtered = medips.count_df.filtered[rownames(medips.count_df.filtered) %in% cpg_count5$window,] #retaining windows with > 5cpg sites
medips.count_df.filtered = medips.count_df.filtered[!rownames(medips.count_df.filtered) %in% blood_wgbs_windows,] #removing regions methylated in PBLs
medips.count_df.filtered = medips.count_df.filtered[rownames(medips.count_df.filtered) %in% regulatory,] #retaining windows located in CpG dense regions, or regulatory/enhancer sites

#library size normalization with deseq2
qc_dir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/medips_files/500bp/qc_metrics/'
qc_data_df = readRDS(paste0(qc_dir, 'batch1_11_qc_summary.extended.RDS'))
qc_data_df = qc_data_df[qc_data_df$GRP_Id %in% colnames(medips.count_df),]
dds <- DESeqDataSetFromMatrix(countData = medips.count_df[,-1],colData = qc_data_df, design = ~Batch  ) #can add gender here but we're only looking at female
dds <- estimateSizeFactors(dds) #estimating size factors for tmm normalization
dat <- counts(dds, normalized=TRUE) #getting tmm-normalised counts
dat = dat[rowSums(dat) > 0,] 
medips.count_df.deseq = dat
#deseq2 normalization
saveRDS(medips.count_df.deseq, paste0(savedir, 'brca.deseq2.medips.count_matrix.RDS'))

#variance stabilizing transformation
vsd <- vst(dds, blind=T, nsub = 10000)
saveRDS(assay(vsd), paste0(savedir,'brca.vst.medips.count_matrix.RDS'))



