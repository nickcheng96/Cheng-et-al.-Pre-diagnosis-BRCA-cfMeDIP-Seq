######functions#####
library(ROCR)
library(caret)
library(plyr)
library(pROC)
library(plotROC)
library(survival)
library(survminer)
library(survAUC)
library(survivalROC)
library(ncvreg)
library(parallel)
library(roperators)
library(foreach)
library(doParallel)
library(survMarkerTwoPhase)
library(RCurl)
library(cvAUC)
library(DESeq2)
library(plyr)

#weighted incidence
all_cause_mortality_CAN<-structure(list(x = c(0, 2.5, 7, 12, 17, 22, 27, 32, 37, 42, 47, 
                                              52, 57, 62, 67, 72, 77, 82, 87, 90), qx = c(0.0047, 2e-04, 1e-04, 
                                                                                          1e-04, 5e-04, 8e-04, 0.001, 0.0011, 0.0013, 0.0015, 0.0023, 0.0036, 
                                                                                          0.0056, 0.0088, 0.0134, 0.0216, 0.0346, 0.061, 0.1089, 0.2163
                                              ), lx = c(1e+05, 97650, 97552.35, 97503.573825, 97454.8220380875, 
                                                        97211.1849829923, 96822.3402430603, 96338.228541845, 95808.3682848649, 
                                                        95185.6138910133, 94471.7217868307, 93385.2969862821, 91704.361640529, 
                                                        89136.6395145942, 85214.6273759521, 79505.2473417633, 70918.6806288528, 
                                                        58649.7488800613, 40761.5754716426, 18566.8976273332), dx = c(470, 
                                                                                                                      19.53, 9.755235, 9.7503573825, 48.7274110190437, 77.7689479863938, 
                                                                                                                      96.8223402430603, 105.97205139603, 124.550878770324, 142.77842083652, 
                                                                                                                      217.284960109711, 336.187069150616, 513.544425186962, 784.402427728429, 
                                                                                                                      1141.87600683776, 1717.31334258209, 2453.78634975831, 3577.63468168374, 
                                                                                                                      4438.93556886188, 4016.01995679217), qx.1 = c(0.004, 1e-04, 1e-04, 
                                                                                                                                                                    1e-04, 2e-04, 3e-04, 4e-04, 5e-04, 6e-04, 9e-04, 0.0015, 0.0023, 
                                                                                                                                                                    0.0037, 0.0058, 0.0087, 0.0142, 0.0232, 0.0416, 0.0768, 0.179
                                                                                                                      ), lx.1 = c(1e+05, 98000, 97951, 97902.0245, 97853.07348775, 
                                                                                                                                  97755.2204142622, 97608.5875836408, 97413.3704084736, 97169.8369824524, 
                                                                                                                                  96878.327471505, 96442.3749978832, 95719.0571853991, 94618.288027767, 
                                                                                                                                  92867.8496992533, 90174.682057975, 86252.0833884531, 80128.1854678729, 
                                                                                                                                  70833.3159535996, 56099.9862352509, 34557.5915209146), dx.1 = c(400, 
                                                                                                                                                                                                  9.8, 9.7951, 9.79020245, 19.57061469755, 29.3265661242787, 39.0434350334563, 
                                                                                                                                                                                                  48.7066852042368, 58.3019021894714, 87.1904947243545, 144.663562496825, 
                                                                                                                                                                                                  220.153831526418, 350.087665702738, 538.633528255669, 784.519733904382, 
                                                                                                                                                                                                  1224.77958411603, 1858.97390285465, 2946.66594366975, 4308.47894286727, 
                                                                                                                                                                                                  6185.80888224371)), class = "data.frame", row.names = c(50436L, 
                                                                                                                                                                                                                                                          50442L, 50448L, 50454L, 50460L, 50466L, 50472L, 50478L, 50484L, 
                                                                                                                                                                                                                                                          50490L, 50496L, 50502L, 50508L, 50514L, 50520L, 50526L, 50532L, 
                                                                                                                                                                                                                                                          50538L, 50544L, 50550L))

age_incidence<-structure(list(Age.group = c("0 to 04", "05 to 09", "10 to 14", 
                                            "15 to 19", "20 to 24", "25 to 29", "30 to 34", "35 to 39", "40 to 44", 
                                            "45 to 49", "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74", 
                                            "75 to 79", "80 to 84", "85 to 89", "90+"), PanB = c(0, 0, 0, 
                                                                                                 0, 0, 0, 0.8, 1.3, 1.6, 5, 9.7, 15.1, 27.1, 40.9, 56.5, 71.4, 
                                                                                                 75, 77, 56.2), PanM = c(0, 0, 0, 0, 0, 0.5, 0.5, 1.1, 2.2, 5.9, 
                                                                                                                         9.8, 19.2, 31.1, 46.7, 66.7, 79.2, 82, 85.1, 64.4), PanF = c(0, 
                                                                                                                                                                                      0, 0, 0.6, 0, 0, 1, 1, 1.1, 4.2, 10.2, 11.1, 23.3, 35.4, 47, 
                                                                                                                                                                                      65.8, 71, 71.6, 55.7), LungB = c(0, 0, 0.3, 0, 0.5, 0.5, 1.3, 
                                                                                                                                                                                                                       1.9, 4.7, 10, 32.1, 71.6, 132.3, 208.3, 302.4, 387.2, 373.1, 
                                                                                                                                                                                                                       302.4, 188), LungM = c(0, 0, 0, 0, 0.5, 0.5, 1, 2.2, 3.9, 9.6, 
                                                                                                                                                                                                                                              28.8, 65.8, 134.1, 217, 314.4, 423.3, 439.9, 384.8, 300.7), LungF = c(0, 
                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0.6, 0.5, 1.5, 2.1, 5.4, 10.5, 35.4, 77.3, 130, 200.8, 
                                                                                                                                                                                                                                                                                                                    292.3, 357.4, 322, 248.3, 139.3), BreastB = c(0, 0, 0, 0, 0.8, 
                                                                                                                                                                                                                                                                                                                                                                  4.8, 13.6, 28, 55, 84, 101, 115, 144, 175, 208, 188, 196.7, 212, 
                                                                                                                                                                                                                                                                                                                                                                  175), BreastM = c(0, 0, 0, 0, 0, 0, 0, 0.5, 0.6, 1.1, 0.5, 1, 
                                                                                                                                                                                                                                                                                                                                                                                    2.8, 4.1, 5.5, 6.8, 6, 3.4, 0), BreastF = c(0, 0, 0, 0.6, 1.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                9.3, 26.8, 55, 109, 166, 200, 227, 280, 335, 393, 346, 345.6, 
                                                                                                                                                                                                                                                                                                                                                                                                                                349, 250.7), ProsB = c(0, 0, 0, 0, 0, 0, 0, 0, 2.7, 9.8, 40.6, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                       93, 161.6, 265.1, 296.3, 278.6, 228.5, 167.4, 110.2), ProsM = c(0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       0, 0, 0, 0, 0, 0, 0, 5.6, 19.2, 81.6, 187.6, 330.4, 547.9, 618.7, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       602.1, 517.9, 422.2, 365.1)), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         -19L))



#functions
auc_calc.seed=   function(plot_list) {
  tmp_list = lapply(plot_list, function(x) x[order(-(x$methylation_score)),] )
  auc_all = lapply(tmp_list, function(dat)  {
    prauc <- function(dat) {
      x <- dat@x.values[[1]]
      y <- dat@y.values[[1]]
      idx <- which(is.nan(y))
      if (any(idx)) {
        x <- x[-idx]
        y <- y[-idx]
      }
      return(pracma::trapz(x, y))
    }
    dat$STATUS = dat$reported
    dat$PRED_CLASS = dat$predictions
    if(nrow(dat[dat$STATUS == st[1],]) > 1 & nrow(dat[dat$STATUS == st[2],]) > 1) {
      pred <- ROCR::prediction(dat$methylation_score, dat$STATUS == 
                                 st[1])
      c1 <- st[1]
      tp <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS == 
                  c1)
      tn <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS != 
                  c1)
      fp <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS != 
                  c1)
      fn <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS == 
                  c1)
      curRoc <- ROCR::performance(pred, "tpr", "fpr")
      curPr <- ROCR::performance(pred, "prec", "rec")
      auroc <- performance(pred, "auc")@y.values[[1]]
      aupr <- prauc(curPr)
      return_df = data.frame(auroc = auroc, aupr = aupr, seed = dat$seed[1])
      return(return_df)
    } else {
      return(NULL)
    }
    
    
  })
  auc_all_df = do.call('rbind',auc_all)
  return(auc_all_df)
}
auc_calc = function(plot_list) {
  tmp_list = lapply(plot_list, function(x) x[order(-(x$methylation_score)),] )
  auc_all = lapply(tmp_list, function(dat)  {
    prauc <- function(dat) {
      x <- dat@x.values[[1]]
      y <- dat@y.values[[1]]
      idx <- which(is.nan(y))
      if (any(idx)) {
        x <- x[-idx]
        y <- y[-idx]
      }
      return(pracma::trapz(x, y))
    }
    dat$STATUS = dat$reported
    dat$PRED_CLASS = dat$predictions
    if(nrow(dat[dat$STATUS == st[1],]) > 1 & nrow(dat[dat$STATUS == st[2],]) > 1) {
      pred <- ROCR::prediction(dat$methylation_score, dat$STATUS == 
                                 st[1])
      c1 <- st[1]
      tp <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS == 
                  c1)
      tn <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS != 
                  c1)
      fp <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS != 
                  c1)
      fn <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS == 
                  c1)
      curRoc <- ROCR::performance(pred, "tpr", "fpr")
      curPr <- ROCR::performance(pred, "prec", "rec")
      auroc <- performance(pred, "auc")@y.values[[1]]
      aupr <- prauc(curPr)
      return_df = data.frame(auroc = auroc, aupr = aupr)
      return(return_df)
    } else {
      return(NULL)
    }
    
    
  })
  auc_all_df = do.call('rbind',auc_all)
  return(auc_all_df)
}
confidence_interval <- function(vector) {
  ordered = vector[order(vector)]
  lower_bound = round(length(ordered)*0.025)
  upper_bound = round(length(ordered)*0.975)
  # Confidence interval as a vector
  result <- c("lower" = ordered[lower_bound], "upper" = ordered[upper_bound])
  return(result)
}
roc_plot_function.only = function(plot_list, file_name = 'performance_auc',dir = figdir, st,auc_summary.df,sample_info, core.number,foldno = 2) {
  incr.count = function(x) {

    tmp.df = x
    tmp.df = tmp.df[order(tmp.df$fpr, tmp.df$tpr),]
    tmp.df$fpr = round(tmp.df$fpr, digits =3)
    tmp.df$tpr = round(tmp.df$tpr, digits =3)

    return_df = NULL
    for (i in 2:nrow(tmp.df)) {

      start_fpr = tmp.df$fpr[i-1]
      end_fpr = tmp.df$fpr[i]
      
      start_tpr = tmp.df$tpr[i-1]
      end_tpr = tmp.df$tpr[i]
      
      fpr.seq = seq(start_fpr,end_fpr, 0.001)
      if(length(fpr.seq) > 1) {
        tpr.seq = seq(start_tpr,end_tpr,by = (end_tpr -start_tpr)/(length(fpr.seq)-1))
      } else{
        tpr.seq = end_tpr
      }
      
      tmp.return.df = data.frame(tpr_mean = tpr.seq, fpr_mean = fpr.seq, seed = x$seed[1])
      return_df = rbind(return_df, tmp.return.df)
    }
    return(return_df)
  }
  fu.times = 365*c(0)
  for (fu.ind in fu.times) {
    sample_info.filt = sample_info
    targ.ids = unique(do.call('rbind',plot_list)$GRP_Id)
    targ.ids= targ.ids[targ.ids %in% sample_info.filt$GRP_Id]
    
    tmp_list = lapply(plot_list, function(x){
      tmp = x[x$GRP_Id %in% targ.ids,]
      tmp = tmp[-tmp$methylation_score,]
      
      tmp$STATUS  = tmp$reported
      tmp$PRED_CLASS = tmp$predictions
      return(tmp)
    }  )
    auc_all = lapply(tmp_list, function(dat)  {
      prauc <- function(dat) {
        x <- dat@x.values[[1]]
        y <- dat@y.values[[1]]
        idx <- which(is.nan(y))
        if (any(idx)) {
          x <- x[-idx]
          y <- y[-idx]
        }
        return(pracma::trapz(x, y))
      }
      if(nrow(dat[dat$STATUS == st[1],]) > 1 & nrow(dat[dat$STATUS == st[2],]) > 1) {
        pred <- ROCR::prediction(dat$methylation_score, dat$STATUS == 
                                   st[1])
        c1 <- st[1]
        tp <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS == 
                    c1)
        tn <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS != 
                    c1)
        fp <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS != 
                    c1)
        fn <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS == 
                    c1)
        curRoc <- ROCR::performance(pred, "tpr", "fpr")
        curPr <- ROCR::performance(pred, "prec", "rec")
        auroc <- performance(pred, "auc")@y.values[[1]]
        aupr <- prauc(curPr)
        return_df = data.frame(auroc = auroc, aupr = aupr)
        return(return_df)
      } else {
        return(NULL)
      }
      
      
    })
    
    auc_all_df = do.call('rbind',auc_all)
    
    average.auc = NULL
    

    tmp = lapply(plot_list, function(x) {
      tmp1 = merge(x[,c('GRP_Id','methylation_score','reported','predictions')], sample_info[,c('GRP_Id','Age','Diagnosis_Time')], by = 'GRP_Id')
      return(tmp1)
    } )
    seeds = unlist(lapply(seq(1:c(length(plot_list)/foldno)),function(x) rep(x, 10) ))
    for (i in 1:length(tmp)) { 
      if (nrow(tmp[[i]]) > 0) { 
        tmp[[i]]$seed = seeds[i]}
    }
    tmp.auc = auc_calc.seed(tmp)
    tmp.seed.average = ddply(tmp.auc, 'seed',numcolwise(mean))
    
    perf.df = NULL
    j = 1
    perf.list= lapply(seq(1:200), function(j) {
      set.seed(j)
      for (i in 1:nrow(tmp.seed.average)) {
        selected.samples = tmp.seed.average[sample(x = seq(1:nrow(tmp.seed.average)), i),]
        tmp.df = data.frame(n = i, auroc_mean = mean(selected.samples$auroc),seed = j)
        perf.df = rbind(perf.df, tmp.df)
      }
      return(perf.df)
    }  )
    
    perf.iteration.df = do.call('rbind',perf.list)
    perf.iteration.df.split = split(perf.iteration.df, perf.iteration.df$n)
    perf.iteration.df.split = lapply(perf.iteration.df.split, function(x) data.frame(n = x$n[1], mean = mean(x$auroc_mean), ci.l = confidence_interval(x$auroc_mean)[1], ci.u = confidence_interval(x$auroc_mean)[2]))
    combined.iteration.perf = do.call('rbind',perf.iteration.df.split)
    tissue = st[1]
    pdf(paste0(tissue,"_perfiteration_AUC_CV.pdf"),height = 4, width = 4)
    plot = ggplot(combined.iteration.perf, aes(ymin = ci.l,ymax=ci.u,x = n),fill = '#5A7684') + 
      geom_ribbon(fill = '#619B8A',alpha = 0.5) + 
      geom_line(aes(y = mean),size = 1, col = '#273F38') +
      theme_bw()+
      scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = seq(0,1,0.1)) +
      scale_x_continuous(limits = c(min(combined.iteration.perf$n),max(combined.iteration.perf$n )), expand = c(0, 0), breaks = seq(min(combined.iteration.perf$n),max(combined.iteration.perf$n),50)) +
      theme(text = element_text(size=6),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = 'bottom',
            legend.title = element_blank(),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold")) +
      xlab('Number of repeats') + 
      ylab('Mean AUC') + ggtitle(round(combined.iteration.perf$mean[nrow(combined.iteration.perf)],digits = 3))
    print(plot)
 
    dev.off()
    
    
    #ROC curve binary
    if (nrow(auc_all_df) > 1) {
      ci_auc = round(confidence_interval(auc_summary.df$average_auroc), digits =3)
      auc_summary = auc_summary.df$average_auroc
      mean.auc.all = round(mean(auc_summary),digits =3)
      median_auc = round(auc_summary[order(auc_summary)][round(length(auc_summary)/2)],digits = 3)
      median_auc_ind = auc_summary.df[round(auc_summary.df$average_auroc, digits = 3) == median_auc,'seed'][1]
      
      print(median_auc)
      auc_list = list()
      auc_list.raw = list()
      auc_list.combined = list()
      auc.mean.scores = list()
      for (i in 1:(length(tmp_list)/foldno)) { #length(tmp_list )
        print(i)
        colors = rev(hcl.colors(50, palette='Purple-Yellow'))
        start_ind = 1+ ((i-1)*foldno)
        end_ind = i*foldno
        cv.average = auc_all[start_ind:end_ind]
        cv.average = do.call('rbind',cv.average)
        cv.average.df = data.frame(auc = mean(cv.average$auroc), seed = i)
        auc.mean.scores[[i]] = cv.average.df
        tmp = tmp_list[start_ind:end_ind]
        tmp.combined= do.call('rbind',tmp)
        tmp.combined = tmp.combined[order(-tmp.combined$methylation_score),]
        
        tpr.fpr.calc = function(x){
          tmp1 = x
          tmp1$f = as.integer(ifelse(tmp1$STATUS == st[2], 0, 1))
          tmp1$f.str = tmp1$reported
          
          tmp1 = tmp1[order(-tmp1$methylation_score),]
          case_no = nrow(tmp1[tmp1$STATUS!=st[2],])
          control_no = nrow(tmp1[tmp1$STATUS==st[2],])
          auc.df =data.frame(matrix(nrow = 0, ncol=3))
          for(l in 1:nrow(tmp1)) {
            x = tmp1[1:l,]
            case_cum = nrow(x[x$STATUS!=st[2],])
            control_cum = nrow(x[x$STATUS==st[2],])
            tpr = case_cum/case_no
            fpr = control_cum/control_no
            return.tmp = data.frame(l,tpr,fpr)
            auc.df = rbind(auc.df, return.tmp)
          }
          return(auc.df)
          #auc.df = auc.df[auc.df$fpr > 0,]
        }
        auc.df_list.combined  = tpr.fpr.calc(tmp.combined)
        auc.df_list.combined$seed = i
        auc.df_list.combined$median = F
        auc.df_list.combined = incr.count(auc.df_list.combined)
        auc = function(tmp.combined) {
          tmp = tmp.combined
          tmp = tmp[order(-tmp$methylation_score),]
          labels = as.character(tmp$reported)
          pred = prediction(predictions = c(tmp$methylation_score) ,labels =  labels, label.ordering = rev(st))
          perf_AUC=performance(pred,"auc") #Calculate the AUC value
          AUC=perf_AUC@y.values[[1]]
          return(AUC)
        }
        
        auc.df_list.combined$auc= auc(tmp.combined)
        auc.df_list = lapply(tmp, function(x){
          tmp1 = x
          tmp1$f = as.integer(ifelse(tmp1$STATUS == st[2], 0, 1))
          tmp1$f.str = tmp1$reported
          
          tmp1 = tmp1[order(-tmp1$methylation_score),]
          case_no = nrow(tmp1[tmp1$STATUS!=st[2],])
          control_no = nrow(tmp1[tmp1$STATUS==st[2],])
          auc.df =data.frame(matrix(nrow = 0, ncol=3))
          if (case_no > 0) {
            for (l in 1:nrow(tmp1)) {
              x = tmp1[1:l,]
              case_cum = nrow(x[x$STATUS!=st[2],])
              control_cum = nrow(x[x$STATUS==st[2],])
              tpr = case_cum/case_no
              fpr = control_cum/control_no
              return.tmp = data.frame(l,tpr,fpr)
              auc.df = rbind(auc.df, return.tmp)
              
            }
          } else {
            auc.df= NULL
          }
          auc.df$seed = i
          tmp.df= data.frame(0,0,0,i)
          colnames(tmp.df) = colnames(auc.df)
          auc.df = rbind(tmp.df, auc.df)
          return(auc.df)
          #auc.df = auc.df[auc.df$fpr > 0,]
        } )
        auc.df_list = lapply(auc.df_list, function(x) incr.count(x))
        auc.df_list = lapply(auc.df_list, function(x) ddply(x, c('seed','fpr_mean'), numcolwise(max)))
        auc.df_list = lapply(auc.df_list, function(x) {
          t = x
          t$fpr_mean = as.character(t$fpr_mean)
          return(t)
        })
        null.ind = sapply(auc.df_list, is.null)
        auc.df_list = auc.df_list[!null.ind]
        auc.df = auc.df_list[[1]][,c('fpr_mean','tpr_mean')]
        for (m in 2:length(auc.df_list)) {
          auc.df = merge(auc.df, auc.df_list[[m]][,c('fpr_mean','tpr_mean')], by = 'fpr_mean', all = T)
        }
        auc.df$fpr_mean = as.numeric(auc.df$fpr_mean)
        auc.df$tpr_mean = rowMeans(auc.df[,-1],na.rm = T)
        auc.df = auc.df[,c('fpr_mean','tpr_mean')]
        auc.df = rbind( data.frame(tpr_mean = 0 , fpr_mean = 0),auc.df)
        
        auc.df = auc.df[order(auc.df$fpr_mean),]
        auc.df.new = auc.df
        auc.df.raw = auc.df
        auc.df.raw$seed = i
        auc.df.raw$auc = cv.average.df$auc
        auc.df$seed = i
        auc.df$median = F
        auc.df$auc = cv.average.df$auc
        
        auc_list[[i]] = auc.df
        auc_list.raw[[i]] = auc.df.raw
        auc_list.combined[[i]] = auc.df_list.combined
      }
      auc_list.combined.incr = mclapply(auc_list.combined, function(x) {
        tmp1 = ddply(x,c('fpr_mean','seed','auc'),numcolwise(max))
        tmp1 = incr.count(tmp1)
        return(tmp1)
      } ,mc.cores = core.number)
      auc.cv.combined.df= do.call('rbind',auc_list.combined.incr)
      auc.combined.df = do.call('rbind',auc_list.combined) #
      mean_auc.combined = ddply(auc.combined.df[,c('tpr_mean','fpr_mean')], c('fpr_mean'), numcolwise(mean))
      auc.df = do.call('rbind',auc_list) #indv fold performnace average per cv iteration
      
      #
      auc.scores = do.call('rbind',auc.mean.scores)
      ci_auc = round(confidence_interval(auc.scores$auc), digits =3)
      auc_summary = auc.scores$auc
      mean.auc.all = round(mean(auc_summary),digits =3)
      median_auc = round(auc_summary[order(auc_summary)][round(length(auc_summary)/2)],digits = 3)
      median_auc_ind = auc_summary.df[round(auc_summary.df$average_auroc, digits = 3) == median_auc,'seed'][1]
      #
      mean_auc = ddply(auc.df[,c('tpr_mean','fpr_mean')], c('fpr_mean'), numcolwise(mean))
      mean_auc$seed = 0
      mean_auc$auc = mean(auc.df$auc)
      mean_auc = mean_auc[order(mean_auc[,2],mean_auc[,3]),]
      mean_auc$median = T
      auc.df.raw = do.call('rbind',auc_list.raw)
      mean_auc.raw = ddply(auc.df.raw[,c('tpr_mean','fpr_mean')], c('fpr_mean'), numcolwise(mean))
      auc.df.raw$median = F
      mean_auc.raw$seed = 0
      mean_auc.raw$auc = mean(auc.df$auc)
      mean_auc.raw = mean_auc.raw[order(mean_auc.raw[,2],mean_auc.raw[,3]),]
      mean_auc.raw$median = T
      auc.df = rbind(auc.df, mean_auc)
      auc.df.raw = rbind(auc.df.raw, mean_auc.raw)
      uni_seed = unique(auc.df[,c('seed','median','auc')])
      col = sapply(uni_seed$auc, function(x){
        auc_round = round(x*100, digits = 0)
        col = colors[max(round(auc_round)-50, 1)]
        return(col)
      } )
      col[median_auc_ind] = 'black'
      auc.df$seed = as.character(auc.df$seed)
      png(paste0(dir,'/',file_name,'_average.smooth.ROC.v2.png'), height = 500, width = 500, type = 'cairo')
      plot =  ggplot(auc.df, aes(x=fpr_mean, y = tpr_mean, col = seed, group = seed))+ theme_bw() + 
        geom_vline(xintercept = 0.05, linetype = "longdash") + 
        geom_line(alpha = 0.3)+
        theme_bw() + 
        ylab('True Positive Rate') + 
        xlab('False Positive Rate') + 
        ggtitle(paste0('AUROC (95% CI): ', median_auc, ' (', ci_auc[1], '-',ci_auc[2],')')) +
        theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),legend.position = "none")+
        scale_color_manual(values = col)#
      plot = plot + geom_line(data = auc.df[auc.df$median == T,], color = 'black', aes(x = fpr_mean, y = tpr_mean)) + geom_vline(xintercept = 0.05, linetype = "longdash")
      print(plot)
      dev.off()
      
      
      

      pdf(paste0(dir,'/',file_name,'.futime.',fu.ind,'_average.smooth.raw.ROC.v2.pdf'), height = 5, width = 5)
      
      auc.df.raw$seed = as.character(auc.df.raw$seed)
      plot =  ggplot(auc.df.raw, aes(x=fpr_mean, y = tpr_mean, col = seed, group = seed))+ theme_bw() + 
        geom_vline(xintercept = 0.05, linetype = "longdash") + 
        geom_line(alpha = 0.3)+
        theme_bw() + 
        ylab('True Positive Rate') + 
        xlab('False Positive Rate') + 
        ggtitle(paste0('AUROC (95% CI): ', median_auc, ' (', ci_auc[1], '-',ci_auc[2],')')) +
        theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),legend.position = "none")+
        scale_color_manual(values = col)#
      plot = plot + geom_line(data = auc.df.raw[auc.df.raw$median == T,], color = 'black', aes(x = fpr_mean, y = tpr_mean)) #+ geom_vline(xintercept = 0.05, linetype = "longdash")
      
      print(plot)
      dev.off()
      
      test = ddply(auc.df.raw[,c('tpr_mean','fpr_mean','seed')],c('fpr_mean','seed'), numcolwise(mean))
      test$tpr_mean = round(test$tpr_mean, digits = 3 )
      test$fpr_mean = round(test$fpr_mean, digits = 3 )
      test.count = ddply(test[,c('tpr_mean','fpr_mean')],c('fpr_mean'), numcolwise(sum))
      test.avg = ddply(test[,c('tpr_mean','fpr_mean')], c('fpr_mean'), numcolwise(mean,na.rm = T))
      test.ciu = ddply(test[,c('tpr_mean','fpr_mean')], c('fpr_mean'), function(x) (quantile(x$tpr_mean,.975,na.rm=T) ))
      test.cil = ddply(test[,c('tpr_mean','fpr_mean')], c('fpr_mean'), function(x) (quantile(x$tpr_mean,.025,na.rm=T) ))
      
      combined.ci = merge(test.ciu, test.cil, by = 'fpr_mean')
      colnames(test.ciu)[2] = 'tpr_mean'
      colnames(test.cil)[2] = 'tpr_mean'
      colnames(combined.ci) = c('fpr_mean','ci.upper','ci.lower')
      # mean score curve#
      mean_scores = do.call('rbind',plot_list)
      mean_scores = ddply(mean_scores[,c('GRP_Id','methylation_score','predictions')], 'GRP_Id', numcolwise(mean))
      mean_scores = mean_scores[order(-mean_scores$methylation_score),]
      mean_scores = merge(mean_scores,combined_list[[1]][,c('GRP_Id','reported')],by = 'GRP_Id')
      mean_scores$reported = factor(as.character(mean_scores$reported), levels = rev(st))
      tmp.list = list(mean_scores)
      mean_scores = mean_scores[order(-mean_scores$methylation_score),]
      auc = auc_calc(tmp.list)
      auc = data.frame(average_auroc = auc[1], seed = 0)
      auc$subtype ='all'
      
      colnames(auc) = c('average_auroc','seed')
      mean_scores.list = list(mean_scores)
      tmp_list.mean = lapply(mean_scores.list, function(x){
        tmp = x[order(-x$methylation_score),]
        tmp$STATUS  = tmp$reported
        tmp$PRED_CLASS = tmp$predictions
        return(tmp)
      }  )
      auc_all.mean = lapply(tmp_list.mean, function(dat)  {
        prauc <- function(dat) {
          x <- dat@x.values[[1]]
          y <- dat@y.values[[1]]
          idx <- which(is.nan(y))
          if (any(idx)) {
            x <- x[-idx]
            y <- y[-idx]
          }
          return(pracma::trapz(x, y))
        }
        if(nrow(dat[dat$STATUS == st[1],]) > 1) {
          pred <- ROCR::prediction(dat$methylation_score, dat$STATUS == 
                                     st[1])
          c1 <- st[1]
          tp <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS == 
                      c1)
          tn <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS != 
                      c1)
          fp <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS != 
                      c1)
          fn <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS == 
                      c1)
          curRoc <- ROCR::performance(pred, "tpr", "fpr")
          curPr <- ROCR::performance(pred, "prec", "rec")
          auroc <- performance(pred, "auc")@y.values[[1]]
          aupr <- prauc(curPr)
          return_df = data.frame(auroc = auroc, aupr = aupr)
          return(return_df)
        } else {
          return(NULL)
        }
        
        
      })
      
      auc_all_df = do.call('rbind',auc_all.mean)
      
      auc_summary = auc$average_auroc
      median_auc = round(auc$average_auroc[1], digits = 3)
      median_auc_ind = auc[round(auc$average_auroc, digits = 3) == median_auc,'seed'][1]
      #median test score aucs
      tmp1= tmp_list.mean[[1]]
      tmp1$f = as.integer(ifelse(tmp1$STATUS == st[2], 0, 1))
      tmp1$f.str = tmp1$reported
      
      tmp1 = tmp1[order(-tmp1$methylation_score),]
      mean.test.auc = auc_calc(list(tmp1))
      case_no = nrow(tmp1[tmp1$STATUS!=st[2],])
      control_no = nrow(tmp1[tmp1$STATUS==st[2],])
      auc.df =data.frame(matrix(nrow = 0, ncol=3))
      for(l in 1:nrow(tmp1)) {
        x = tmp1[1:l,]
        case_cum = nrow(x[x$STATUS!=st[2],])
        control_cum = nrow(x[x$STATUS==st[2],])
        tpr = case_cum/case_no
        fpr = control_cum/control_no
        return.tmp = data.frame(l,tpr,fpr)
        auc.df = rbind(auc.df, return.tmp)
      }
      
      auc.df = rbind(data.frame(l=0,tpr=0,fpr=0),auc.df)
      
      
      test.avg = test.avg[order(test.avg$fpr_mean,test.avg$tpr_mean),]
      plot =  ggplot(combined.ci, aes(x=fpr_mean, ymax = ci.upper,ymin = ci.lower))+ theme_bw() + 

        scale_x_continuous(limits = c(0,1.01), expand= c(0,0),breaks = seq(0,1,0.1))+
        scale_y_continuous(limits = c(0,1), expand= c(0,0),breaks = seq(0,1,0.1))+
        geom_ribbon(alpha = 0.5,fill = '#619B8A')+
        theme_bw() + 
        ylab('True Positive Rate') + 
        xlab('False Positive Rate') + 
        theme(text = element_text(size=6),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),legend.position = "none")+
        scale_color_manual(values = col)#
      test.avg$ci.upper = test.avg$tpr_mean
      test.avg$ci.lower = test.avg$tpr_mean
      sens.90 = round(test.avg[test.avg$fpr_mean == 0.1, 'tpr_mean'],digits = 3)
      sens.90.ciu = round(combined.ci[combined.ci$fpr_mean == 0.1, 'ci.upper'],digits = 3)
      sens.90.cil = round(combined.ci[combined.ci$fpr_mean == 0.1, 'ci.lower'],digits = 3)
      sens.95 = round(test.avg[test.avg$fpr_mean == 0.05, 'tpr_mean'],digits = 3)
      sens.95.ciu = round(combined.ci[combined.ci$fpr_mean == 0.05, 'ci.upper'],digits = 3)
      sens.95.cil = round(combined.ci[combined.ci$fpr_mean == 0.05, 'ci.lower'],digits = 3)
      
      
      plot = plot + geom_line(data = test.avg, aes(x=fpr_mean, y = tpr_mean), size = 1, col = '#273F38') +   ggtitle(paste0(mean.auc.all, ' (', ci_auc[1], '-',ci_auc[2],')',' sens.90 = ', sens.90,'(',sens.90.cil,'-',sens.90.ciu,') \n sens.95 = ',sens.95,'(',sens.95.cil,'-',sens.95.ciu,')')) 
    
      pdf(paste0(dir,'/',file_name,'.futime.',fu.ind,'_average.smooth.raw.ROC2.nomedtest.v2.pdf'), height = 4, width = 4)
      
      print(plot)
      dev.off()
      
      auc.df$fpr_mean = auc.df$fpr
      auc.df$tpr_mean = auc.df$tpr
      auc.df$ci.upper = auc.df$tpr_mean
      auc.df$ci.lower = auc.df$ci.upper
      plot = plot +  geom_line(data =auc.df, aes(x=fpr_mean, y = tpr_mean), size = 1, col = '#612940') # +   ggtitle(paste0('AUROC (95% CI): ', median_auc, ' (', ci_auc[1], '-',ci_auc[2],')',' median.test auc = ', round(mean.test.auc[1],digits = 3))) 
      
      pdf(paste0(dir,'/',file_name,'.futime.',fu.ind,'_average.smooth.raw.ROC2.v2.pdf'), height = 4, width = 4)
      
      print(plot)
      dev.off()
      
      
      ###combined cv#####
      auc_list.combined.incr.df = do.call('rbind',auc_list.combined.incr)
      test = auc_list.combined.incr.df
    
      test.avg = ddply(test[,c('tpr_mean','fpr_mean')], c('fpr_mean'), numcolwise(mean,na.rm = T))
      test.ciu = ddply(test[,c('tpr_mean','fpr_mean')], c('fpr_mean'), function(x) (quantile(x$tpr_mean,.975,na.rm=T) ))
      test.cil = ddply(test[,c('tpr_mean','fpr_mean')], c('fpr_mean'), function(x) (quantile(x$tpr_mean,.025,na.rm=T) ))
      
      combined.ci = merge(test.ciu, test.cil, by = 'fpr_mean')
      colnames(test.ciu)[2] = 'tpr_mean'
      colnames(test.cil)[2] = 'tpr_mean'
      colnames(combined.ci) = c('fpr_mean','ci.upper','ci.lower')
      # mean score curve#
    
      mean_scores = ddply(mean_scores[,c('GRP_Id','methylation_score','reported')], c('GRP_Id','reported'), numcolwise(mean))
      mean_scores = mean_scores[order(-mean_scores$methylation_score),]

      mean_scores$reported = factor(as.character(mean_scores$reported), levels = rev(st))
      tmp.list = list(mean_scores)
      mean_scores = mean_scores[order(-mean_scores$methylation_score),]
      auc = auc_calc(tmp.list)
      auc = data.frame(average_auroc = auc[1], seed = 0)
      auc$subtype ='all'
      
      colnames(auc) = c('average_auroc','seed')
      mean_scores = list(mean_scores)
      tmp_list.mean = lapply(mean_scores, function(x){
        tmp = x[order(-x$methylation_score),]
        tmp$STATUS  = tmp$reported
        tmp$PRED_CLASS = tmp$predictions
        return(tmp)
      }  )
      auc_all.mean = lapply(tmp_list.mean, function(dat)  {
        prauc <- function(dat) {
          x <- dat@x.values[[1]]
          y <- dat@y.values[[1]]
          idx <- which(is.nan(y))
          if (any(idx)) {
            x <- x[-idx]
            y <- y[-idx]
          }
          return(pracma::trapz(x, y))
        }
        if(nrow(dat[dat$STATUS == st[1],]) > 1) {
          pred <- ROCR::prediction(dat$methylation_score, dat$STATUS == 
                                     st[1])
          c1 <- st[1]
          tp <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS == 
                      c1)
          tn <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS != 
                      c1)
          fp <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS != 
                      c1)
          fn <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS == 
                      c1)
          curRoc <- ROCR::performance(pred, "tpr", "fpr")
          curPr <- ROCR::performance(pred, "prec", "rec")
          auroc <- performance(pred, "auc")@y.values[[1]]
          aupr <- prauc(curPr)
          return_df = data.frame(auroc = auroc, aupr = aupr)
          return(return_df)
        } else {
          return(NULL)
        }
        
        
      })
      
      auc_all_df = do.call('rbind',auc_all.mean)
      
      auc_summary = auc$average_auroc
      median_auc = round(auc$average_auroc[1], digits = 3)
      median_auc_ind = auc[round(auc$average_auroc, digits = 3) == median_auc,'seed'][1]
      #median test score aucs
      tmp1= tmp_list.mean[[1]]
      tmp1$f = as.integer(ifelse(tmp1$STATUS == st[2], 0, 1))
      tmp1$f.str = tmp1$reported
      
      tmp1 = tmp1[order(-tmp1$methylation_score),]
      mean.test.auc = auc_calc(list(tmp1))
      case_no = nrow(tmp1[tmp1$STATUS!=st[2],])
      control_no = nrow(tmp1[tmp1$STATUS==st[2],])
      auc.df =data.frame(matrix(nrow = 0, ncol=3))
      for(l in 1:nrow(tmp1)) {
        x = tmp1[1:l,]
        case_cum = nrow(x[x$STATUS!=st[2],])
        control_cum = nrow(x[x$STATUS==st[2],])
        tpr = case_cum/case_no
        fpr = control_cum/control_no
        return.tmp = data.frame(l,tpr,fpr)
        auc.df = rbind(auc.df, return.tmp)
      }
      
      auc.df = rbind(data.frame(l=0,tpr=0,fpr=0),auc.df)
      
      
      test.avg = test.avg[order(test.avg$fpr_mean,test.avg$tpr_mean),]
      plot =  ggplot(combined.ci, aes(x=fpr_mean, ymax = ci.upper,ymin = ci.lower))+ theme_bw() +
        scale_x_continuous(limits = c(0,1.01), expand= c(0,0),breaks = seq(0,1,0.1))+
        scale_y_continuous(limits = c(0,1), expand= c(0,0),breaks = seq(0,1,0.1))+
        geom_ribbon(alpha = 0.5,fill = '#619B8A')+
        theme_bw() + 
        ylab('True Positive Rate') + 
        xlab('False Positive Rate') + 
        theme(text = element_text(size=6),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),legend.position = "none")+
        scale_color_manual(values = col)#
      test.avg$ci.upper = test.avg$tpr_mean
      test.avg$ci.lower = test.avg$tpr_mean
      sens.90 = round(test.avg[test.avg$fpr_mean == 0.1, 'tpr_mean'],digits = 3)
      sens.90.ciu = round(combined.ci[combined.ci$fpr_mean == 0.1, 'ci.upper'],digits = 3)
      sens.90.cil = round(combined.ci[combined.ci$fpr_mean == 0.1, 'ci.lower'],digits = 3)
      sens.95 = round(test.avg[test.avg$fpr_mean == 0.05, 'tpr_mean'],digits = 3)
      sens.95.ciu = round(combined.ci[combined.ci$fpr_mean == 0.05, 'ci.upper'],digits = 3)
      sens.95.cil = round(combined.ci[combined.ci$fpr_mean == 0.05, 'ci.lower'],digits = 3)
      
      
      plot = plot + geom_line(data = test.avg, aes(x=fpr_mean, y = tpr_mean), size = 1, col = '#273F38') +   ggtitle(paste0(mean.auc.all, ' (', ci_auc[1], '-',ci_auc[2],')',' sens.90 = ', sens.90,'(',sens.90.cil,'-',sens.90.ciu,') \n sens.95 = ',sens.95,'(',sens.95.cil,'-',sens.95.ciu,')')) 
      pdf(paste0(dir,'/',file_name,'.futime.',fu.ind,'_average.smooth.raw.CV.AUC.pdf'), height = 4, width = 4)
      
      print(plot)
      dev.off()
      
      
    } 
  }
}
sens.calc= function(plot_list, file_name = 'performance_auc',dir = figdir, st,auc_summary.df,sample_info, core.number,foldno = 2) {
  incr.count = function(x) {

    tmp.df = x
    tmp.df = tmp.df[order(tmp.df$fpr, tmp.df$tpr),]
    tmp.df$fpr = round(tmp.df$fpr, digits =3)
    tmp.df$tpr = round(tmp.df$tpr, digits =3)

    return_df = NULL
    for (i in 2:nrow(tmp.df)) {

      start_fpr = tmp.df$fpr[i-1]
      end_fpr = tmp.df$fpr[i]
      
      start_tpr = tmp.df$tpr[i-1]
      end_tpr = tmp.df$tpr[i]
      
      fpr.seq = seq(start_fpr,end_fpr, 0.001)
      if(length(fpr.seq) > 1) {
        tpr.seq = seq(start_tpr,end_tpr,by = (end_tpr -start_tpr)/(length(fpr.seq)-1))
      } else{
        tpr.seq = end_tpr
      }
      
      tmp.return.df = data.frame(tpr_mean = tpr.seq, fpr_mean = fpr.seq, seed = x$seed[1])
      return_df = rbind(return_df, tmp.return.df)
    }
    return(return_df)
  }
  fu.times = 365*c(0)
  tissue<-st[1]
  for (fu.ind in fu.times) {
    
    
    #performance across iteratison
    tmp = lapply(plot_list, function(x) {
      tmp1 = merge(x[,c('GRP_Id','methylation_score','reported','predictions')], sample_info[,c('GRP_Id','Age','Diagnosis_Time')], by = 'GRP_Id')

      return(tmp1)
    } )
    seeds = unlist(lapply(seq(1:c(length(plot_list)/foldno)),function(x) rep(x, 10) ))
    for (i in 1:length(tmp)) { 
      if (nrow(tmp[[i]]) > 0) { 
        tmp[[i]]$seed = seeds[i]}
    }
    tmp.auc = auc_calc.seed(tmp)
    tmp.seed.average = ddply(tmp.auc, 'seed',numcolwise(mean))
    
    perf.df = NULL
    j = 1
    perf.list= lapply(seq(1:200), function(j) {
      set.seed(j)
      for (i in 1:nrow(tmp.seed.average)) {
        selected.samples = tmp.seed.average[sample(x = seq(1:nrow(tmp.seed.average)), i),]
        tmp.df = data.frame(n = i, auroc_mean = mean(selected.samples$auroc),seed = j)
        perf.df = rbind(perf.df, tmp.df)
      }
      return(perf.df)
    }  )
    
    perf.iteration.df = do.call('rbind',perf.list)
    perf.iteration.df.split = split(perf.iteration.df, perf.iteration.df$n)
    perf.iteration.df.split = lapply(perf.iteration.df.split, function(x) data.frame(n = x$n[1], mean = mean(x$auroc_mean), ci.l = confidence_interval(x$auroc_mean)[1], ci.u = confidence_interval(x$auroc_mean)[2]))
    combined.iteration.perf = do.call('rbind',perf.iteration.df.split)
    pdf(paste0(tissue,"_perfiteration_AUC_CV.pdf"),height = 4, width = 4)
    plot = ggplot(combined.iteration.perf, aes(ymin = ci.l,ymax=ci.u,x = n),fill = '#5A7684') + 
      geom_ribbon(fill = '#619B8A',alpha = 0.5) + 
      geom_line(aes(y = mean),size = 1, col = '#273F38') +
      theme_bw()+
      scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = seq(0,1,0.1)) +
      scale_x_continuous(limits = c(min(combined.iteration.perf$n),max(combined.iteration.perf$n )), expand = c(0, 0), breaks = seq(min(combined.iteration.perf$n),max(combined.iteration.perf$n),50)) +
      theme(text = element_text(size=6),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = 'bottom',
            legend.title = element_blank(),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold")) +
      xlab('Number of repeats') + 
      ylab('Mean AUC') + ggtitle(round(combined.iteration.perf$mean[nrow(combined.iteration.perf)],digits = 3))# + geom_line(data = concordance.df.mediantest, aes(x = year, y = cil), size = 1, col = '#612940')
    print(plot)
    #+ scale_fill_manual(values = cpg_region_colors)+ guides(colour = guide_legend(override.aes = list(size=4))) + ggtitle(gtitle)
    dev.off()
    
    
    
    #all auc
    sample_info.filt = sample_info
    targ.ids = unique(do.call('rbind',plot_list)$GRP_Id)
    targ.ids= targ.ids[targ.ids %in% sample_info.filt$GRP_Id]
    
    tmp_list = lapply(plot_list, function(x){
      tmp = x[x$GRP_Id %in% targ.ids,]
      tmp = tmp[-tmp$methylation_score,]
      
      tmp$STATUS  = tmp$reported
      tmp$PRED_CLASS = tmp$predictions
      return(tmp)
    }  )
    auc_all = lapply(tmp_list, function(dat)  {
      prauc <- function(dat) {
        x <- dat@x.values[[1]]
        y <- dat@y.values[[1]]
        idx <- which(is.nan(y))
        if (any(idx)) {
          x <- x[-idx]
          y <- y[-idx]
        }
        return(pracma::trapz(x, y))
      }
      if(nrow(dat[dat$STATUS == st[1],]) > 1 & nrow(dat[dat$STATUS == st[2],]) > 1) {
        pred <- ROCR::prediction(dat$methylation_score, dat$STATUS == 
                                   st[1])
        c1 <- st[1]
        tp <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS == 
                    c1)
        tn <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS != 
                    c1)
        fp <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS != 
                    c1)
        fn <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS == 
                    c1)
        curRoc <- ROCR::performance(pred, "tpr", "fpr")
        curPr <- ROCR::performance(pred, "prec", "rec")
        auroc <- performance(pred, "auc")@y.values[[1]]
        aupr <- prauc(curPr)
        return_df = data.frame(auroc = auroc, aupr = aupr)
        return(return_df)
      } else {
        return(NULL)
      }
      
      
    })
    
    auc_all_df = do.call('rbind',auc_all)
    
    average.auc = NULL
    
    #ROC curve binary
    if (nrow(auc_all_df) > 1) {
      ci_auc = round(confidence_interval(auc_summary.df$average_auroc), digits =3)
      auc_summary = auc_summary.df$average_auroc
      mean.auc.all = round(mean(auc_summary),digits =3)
      median_auc = round(auc_summary[order(auc_summary)][round(length(auc_summary)/2)],digits = 3)
      median_auc_ind = auc_summary.df[round(auc_summary.df$average_auroc, digits = 3) == median_auc,'seed'][1]
      
      print(median_auc)
      auc_list = list()
      auc_list.raw = list()
      auc_list.combined = list()
      auc.mean.scores = list()
      for (i in 1:(length(tmp_list)/foldno)) {
        print(i)
        colors = rev(hcl.colors(50, palette='Purple-Yellow'))
        start_ind = 1+ ((i-1)*foldno)
        end_ind = i*foldno
        cv.average = auc_all[start_ind:end_ind]
        cv.average = do.call('rbind',cv.average)
        cv.average.df = data.frame(auc = mean(cv.average$auroc), seed = i)
        auc.mean.scores[[i]] = cv.average.df
        tmp = tmp_list[start_ind:end_ind]
        tmp.combined= do.call('rbind',tmp)
        tmp.combined = tmp.combined[order(-tmp.combined$methylation_score),]
        
        tpr.fpr.calc = function(x){
          tmp1 = x
          tmp1$f = as.integer(ifelse(tmp1$STATUS == st[2], 0, 1))
          tmp1$f.str = tmp1$reported
          
          tmp1 = tmp1[order(-tmp1$methylation_score),]
          case_no = nrow(tmp1[tmp1$STATUS!=st[2],])
          control_no = nrow(tmp1[tmp1$STATUS==st[2],])
          auc.df =data.frame(matrix(nrow = 0, ncol=3))
          for(l in 1:nrow(tmp1)) {
            x = tmp1[1:l,]
            case_cum = nrow(x[x$STATUS!=st[2],])
            control_cum = nrow(x[x$STATUS==st[2],])
            tpr = case_cum/case_no
            fpr = control_cum/control_no
            return.tmp = data.frame(l,tpr,fpr)
            auc.df = rbind(auc.df, return.tmp)
          }
          return(auc.df)
          #auc.df = auc.df[auc.df$fpr > 0,]
        }
        auc.df_list.combined  = tpr.fpr.calc(tmp.combined)
        auc.df_list.combined$seed = i
        auc.df_list.combined$median = F
        auc.df_list.combined = incr.count(auc.df_list.combined)
        auc = function(tmp.combined) {
          tmp = tmp.combined
          tmp = tmp[order(-tmp$methylation_score),]
          labels = as.character(tmp$reported)
          pred = prediction(predictions = c(tmp$methylation_score) ,labels =  labels, label.ordering = rev(st))
          perf_AUC=performance(pred,"auc") #Calculate the AUC value
          AUC=perf_AUC@y.values[[1]]
          return(AUC)
        }
        
        auc.df_list.combined$auc= auc(tmp.combined)
        auc.df_list = lapply(tmp, function(x){
          tmp1 = x
          tmp1$f = as.integer(ifelse(tmp1$STATUS == st[2], 0, 1))
          tmp1$f.str = tmp1$reported
          
          tmp1 = tmp1[order(-tmp1$methylation_score),]
          case_no = nrow(tmp1[tmp1$STATUS!=st[2],])
          control_no = nrow(tmp1[tmp1$STATUS==st[2],])
          auc.df =data.frame(matrix(nrow = 0, ncol=3))
          if (case_no > 0 & control_no > 0 & sum(case_no,control_no) > 2 ) {
            for (l in 1:nrow(tmp1)) {
              x = tmp1[1:l,]
              case_cum = nrow(x[x$STATUS!=st[2],])
              control_cum = nrow(x[x$STATUS==st[2],])
              tpr = case_cum/case_no
              fpr = control_cum/control_no
              return.tmp = data.frame(l,tpr,fpr)
              auc.df = rbind(auc.df, return.tmp)
            }
            auc.df$seed = i
            tmp.df= data.frame(0,0,0,i)
            colnames(tmp.df) = colnames(auc.df)
            auc.df = rbind(tmp.df, auc.df)
          } else {
            tmp.df= NULL
          }
          
          
          return(auc.df)

        } )
        auc.df_list = auc.df_list[sapply(auc.df_list,nrow) > 0]
        auc.df_list = lapply(auc.df_list, function(x) incr.count(x))
        auc.df_list = lapply(auc.df_list, function(x) ddply(x, c('seed','fpr_mean'), numcolwise(max)))
        auc.df_list = lapply(auc.df_list, function(x) {
          t = x
          t$fpr_mean = as.character(t$fpr_mean)
          return(t)
        })
        null.ind = sapply(auc.df_list, is.null)
        auc.df_list = auc.df_list[!null.ind]
        auc.df = auc.df_list[[1]][,c('fpr_mean','tpr_mean')]
        if (length(auc.df_list) > 1) {
          for (m in 2:length(auc.df_list)) {
            auc.df = merge(auc.df, auc.df_list[[m]][,c('fpr_mean','tpr_mean')], by = 'fpr_mean', all = T)
          }
          auc.df$tpr_mean = rowMeans(auc.df[,-1],na.rm = T)
          
        }
        
        auc.df$fpr_mean = as.numeric(auc.df$fpr_mean)
        auc.df = auc.df[,c('fpr_mean','tpr_mean')]
        auc.df = rbind( data.frame(tpr_mean = 0 , fpr_mean = 0),auc.df)
        
        auc.df = auc.df[order(auc.df$fpr_mean),]
        auc.df.new = auc.df
        auc.df.raw = auc.df
        auc.df.raw$seed = i
        auc.df.raw$auc = cv.average.df$auc
        auc.df$seed = i
        auc.df$median = F
        auc.df$auc = cv.average.df$auc
        
        auc_list[[i]] = auc.df
        auc_list.raw[[i]] = auc.df.raw
        auc_list.combined[[i]] = auc.df_list.combined
      }
      auc_list.combined.incr = mclapply(auc_list.combined, function(x) {
        tmp1 = ddply(x,c('fpr_mean','seed','auc'),numcolwise(max))
        tmp1 = incr.count(tmp1)
        return(tmp1)
      } ,mc.cores = core.number)
      auc.cv.combined.df= do.call('rbind',auc_list.combined.incr)
      auc.combined.df = do.call('rbind',auc_list.combined) #
      mean_auc.combined = ddply(auc.combined.df[,c('tpr_mean','fpr_mean')], c('fpr_mean'), numcolwise(mean))
      auc.df = do.call('rbind',auc_list) #indv fold performnace average per cv iteration
      
      #
      auc.scores = do.call('rbind',auc.mean.scores)
      ci_auc = round(confidence_interval(auc.scores$auc), digits =3)
      auc_summary = auc.scores$auc
      mean.auc.all = round(mean(auc_summary),digits =3)
      median_auc = round(auc_summary[order(auc_summary)][round(length(auc_summary)/2)],digits = 3)
      median_auc_ind = auc_summary.df[round(auc_summary.df$average_auroc, digits = 3) == median_auc,'seed'][1]
      #
      mean_auc = ddply(auc.df[,c('tpr_mean','fpr_mean')], c('fpr_mean'), numcolwise(mean))
      mean_auc$seed = 0
      mean_auc$auc = mean(auc.df$auc)
      mean_auc = mean_auc[order(mean_auc[,2],mean_auc[,3]),]
      mean_auc$median = T
      auc.df.raw = do.call('rbind',auc_list.raw)
      mean_auc.raw = ddply(auc.df.raw[,c('tpr_mean','fpr_mean')], c('fpr_mean'), numcolwise(mean))
      auc.df.raw$median = F
      mean_auc.raw$seed = 0
      mean_auc.raw$auc = mean(auc.df$auc)
      mean_auc.raw = mean_auc.raw[order(mean_auc.raw[,2],mean_auc.raw[,3]),]
      mean_auc.raw$median = T
      auc.df = rbind(auc.df, mean_auc)
      auc.df.raw = rbind(auc.df.raw, mean_auc.raw)
      uni_seed = unique(auc.df[,c('seed','median','auc')])
      col = sapply(uni_seed$auc, function(x){
        auc_round = round(x*100, digits = 0)
        col = colors[max(round(auc_round)-50, 1)]
        return(col)
      } )
      col[median_auc_ind] = 'black'
      auc.df$seed = as.character(auc.df$seed)
      #

      test = ddply(auc.df.raw[,c('tpr_mean','fpr_mean','seed')],c('fpr_mean','seed'), numcolwise(mean))
      test$tpr_mean = round(test$tpr_mean, digits = 3 )
      test$fpr_mean = round(test$fpr_mean, digits = 3 )
      test.count = ddply(test[,c('tpr_mean','fpr_mean')],c('fpr_mean'), numcolwise(sum))
      test.avg = ddply(test[,c('tpr_mean','fpr_mean')], c('fpr_mean'), numcolwise(mean,na.rm = T))
      test.ciu = ddply(test[,c('tpr_mean','fpr_mean')], c('fpr_mean'), function(x) (quantile(x$tpr_mean,.975,na.rm=T) ))
      test.cil = ddply(test[,c('tpr_mean','fpr_mean')], c('fpr_mean'), function(x) (quantile(x$tpr_mean,.025,na.rm=T) ))
      
      combined.ci = merge(test.ciu, test.cil, by = 'fpr_mean')
      colnames(test.ciu)[2] = 'tpr_mean'
      colnames(test.cil)[2] = 'tpr_mean'
      colnames(combined.ci) = c('fpr_mean','ci.upper','ci.lower')
      # mean score curve#
      mean_scores = do.call('rbind',plot_list)
      mean_scores = ddply(mean_scores[,c('GRP_Id','methylation_score')], 'GRP_Id', numcolwise(mean))
      mean_scores = mean_scores[order(-mean_scores$methylation_score),]
      mean_scores = merge(mean_scores,combined_list[[1]][,c('GRP_Id','reported')],by = 'GRP_Id')
      mean_scores$reported = factor(as.character(mean_scores$reported), levels = rev(st))
      tmp.list = list(mean_scores)
      mean_scores = mean_scores[order(-mean_scores$methylation_score),]
      auc = auc_calc(tmp.list)
      auc = data.frame(average_auroc = auc[1], seed = 0)
      auc$subtype ='all'
      
      colnames(auc) = c('average_auroc','seed')
      mean_scores.list = list(mean_scores)
      tmp_list.mean = lapply(mean_scores.list, function(x){
        tmp = x[order(-x$methylation_score),]
        tmp$STATUS  = tmp$reported
        tmp$PRED_CLASS = tmp$predictions
        return(tmp)
      }  )
      auc_all.mean = lapply(tmp_list.mean, function(dat)  {
        prauc <- function(dat) {
          x <- dat@x.values[[1]]
          y <- dat@y.values[[1]]
          idx <- which(is.nan(y))
          if (any(idx)) {
            x <- x[-idx]
            y <- y[-idx]
          }
          return(pracma::trapz(x, y))
        }
        if(nrow(dat[dat$STATUS == st[1],]) > 1) {
          pred <- ROCR::prediction(dat$methylation_score, dat$STATUS == 
                                     st[1])
          c1 <- st[1]
          tp <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS == 
                      c1)
          tn <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS != 
                      c1)
          fp <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS != 
                      c1)
          fn <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS == 
                      c1)
          curRoc <- ROCR::performance(pred, "tpr", "fpr")
          curPr <- ROCR::performance(pred, "prec", "rec")
          auroc <- performance(pred, "auc")@y.values[[1]]
          aupr <- prauc(curPr)
          return_df = data.frame(auroc = auroc, aupr = aupr)
          return(return_df)
        } else {
          return(NULL)
        }
        
        
      })
      
      auc_all_df = do.call('rbind',auc_all.mean)
      
      auc_summary = auc$average_auroc
      median_auc = round(auc$average_auroc[1], digits = 3)
      median_auc_ind = auc[round(auc$average_auroc, digits = 3) == median_auc,'seed'][1]
      #median test score aucs
      tmp1= tmp_list.mean[[1]]
      tmp1$f = as.integer(ifelse(tmp1$STATUS == st[2], 0, 1))
      tmp1$f.str = tmp1$reported
      
      tmp1 = tmp1[order(-tmp1$methylation_score),]
      mean.test.auc = auc_calc(list(tmp1))
      case_no = nrow(tmp1[tmp1$STATUS!=st[2],])
      control_no = nrow(tmp1[tmp1$STATUS==st[2],])
      auc.df =data.frame(matrix(nrow = 0, ncol=3))
      for(l in 1:nrow(tmp1)) {
        x = tmp1[1:l,]
        case_cum = nrow(x[x$STATUS!=st[2],])
        control_cum = nrow(x[x$STATUS==st[2],])
        tpr = case_cum/case_no
        fpr = control_cum/control_no
        return.tmp = data.frame(l,tpr,fpr)
        auc.df = rbind(auc.df, return.tmp)
      }
      
      test$tpr_mean = round(test$tpr_mean, digits = 3 )
      test$fpr_mean = round(test$fpr_mean, digits = 3 )
      test.sens.95.spec = test[test$fpr_mean == '0.05',c('seed','tpr_mean')]
      colnames(test.sens.95.spec) = c('seed','sens.95')
      test.sens.99.spec = test[test$fpr_mean == '0.01',c('seed','tpr_mean')]
      colnames(test.sens.99.spec) = c('seed','sens.99')
      return.df =merge(test.sens.99.spec,test.sens.95.spec, by = 'seed')
      return.df$Cancer = st[1]
      
      return(return.df)
      
      
    }
    
  }
  
}

#
wkdir='/Path/to/prediction/output/'
setwd(wkdir)

perf_masterlist = list()
results_list_all = list()
combined_list = list()
st = c('breast_cancer','control')

sample_info = readRDS('brca.sample.information.RDS')
matrix.raw = readRDS('normalised.matrix.RDS')

#performance summary
model = c('rf')
for (mod in 1:length(model )) {
  setwd(paste0(wkdir,'/',model[mod]))
  tmpdir = getwd()
 
  feature_size = paste0(seq(100,200,25),'.hyper')

  for (feat in 1:length(feature_size)) {
    setwd(paste0(tmpdir,'/',feature_size[feat],'/performance/'))
    print(model[mod])
    print(feature_size[feat])
    indv_performance = list()
    average_performnace = list()
    complete_list = list()
    combined_list =list()
    seed.range = c(1:200)
  
    for (seedno in seed.range) {
      file = '.*read.5_seed.seedno_fold.*'
      targ_file = gsub('seedno',seedno, file)
      performance_files = list.files(pattern = targ_file)
      if (length(performance_files) > 0) {
        fold_perf_list = lapply(performance_files, function(x) {
          tmp = readRDS(x)
          tmp$seed = seedno
          tmp = tmp[order(-tmp$methylation_score),]
          controls = tmp[tmp$reported == st[2],]
          controlno = nrow(controls)
          spec.95.score = controls[controlno - round(0.95*controlno),'methylation_score']
          
          if (controlno - round(0.99*controlno) == 0) {
            spec.99.score = controls[controlno - round(0.99*controlno) + 1,'methylation_score']
          } else {
            spec.99.score = controls[controlno - round(0.99*controlno),'methylation_score']
          }
          
          
          cancers = tmp[tmp$reported == st[1],]
          spec.95.score = nrow(cancers[cancers$methylation_score > spec.95.score,])/nrow(cancers)
          spec.99.score =  nrow(cancers[cancers$methylation_score > spec.99.score,])/nrow(cancers)
          
          tmp$sens.spec.95 = spec.95.score
          tmp$sens.spec.99 = spec.99.score
          return(tmp)
        } )
        for (i in 1:length(fold_perf_list)) {
          fold_perf_list[[i]]$fold = gsub('.*seed.','',gsub('_fold.*','',performance_files[i]))
        }
        indv_performance[[seedno]] = do.call('rbind',lapply(fold_perf_list, function(x) x[1,]))
        average_performnace[[seedno]] = data.frame(model = model[mod], features = feature_size[feat], average_auroc = mean(indv_performance[[seedno]]$auroc),seed = seedno)
        combined_list[[seedno]] = do.call('rbind',fold_perf_list)
        complete_list =c(complete_list, fold_perf_list)
      }
    }
    
    indv_performance = indv_performance[which(sapply(indv_performance,length)>0)]
    combined_list = combined_list[which(sapply(combined_list,length)>0)]
    average_performnace =average_performnace[which(unlist(sapply(average_performnace,ncol))>0)]
    predList = list()
    featScores = list()
    
    
    results_list_all[[length(results_list_all) + 1]] = do.call('rbind',average_performnace)
    tmp_perf = data.frame(mean_auc = mean(do.call('rbind', lapply(average_performnace, function(x) x$average_auroc[1]))), model = model[mod], features = feature_size[feat])
    perf_masterlist[[length(perf_masterlist) + 1]] = tmp_perf
    roc_plot_function.only(complete_list, file_name = 'performance.smooth', dir = getwd(), st = st, results_list_all[[length(results_list_all)]], sample_info, core.number = detectCores()-2,foldno = 10) 

  }
  
}

#AUC across different feature sizes
overall.performance = do.call('rbind',perf_masterlist)
overall.performance = overall.performance[order(-overall.performance[,1]),]
head(overall.performance)

overall.performance$count= as.numeric(gsub('.hyper','',overall.performance$features))
pdf(paste0(wkdir,'feature.size.auc.performance.pdf'),height = 5, width = 5)
plot1 = ggplot(overall.performance, aes(x = count, y= mean_auc)) +  #Diagnosis_Time
  geom_point() + geom_line()+
  theme_bw()+
  theme(text = element_text(size=15),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank()) +
  xlab('Feature Number') + ylab('Mean AUC')
print(plot1)
dev.off()

####performance across subgroups####
sample.information = readRDS('brca.sample.information.RDS')
targ_stage = list('Stage I' = c('I'), 'Stage II-IV' = c('II','III','IV'), 'Stage II' = c('II'),'Stage III-IV' = c('III','IV'))
targ_mammo = c('Pre-dx Mammo < 1 year')
hr_groups = c('HR Positive','HR Negative','HER2 Positive','HER2 Negative') #,'HR Positive, HER2 Negative','HR Negative, HER2 Negative','HR Positive, HER2 Positive','HR Negative, HER2 Positive'
targ_age_list = list('Dx Age Pre-50' =c(1:50),'Dx Age Post-50' = c(50:100))
stage_list = list()
stage_mean_list = list()

auc_calc2 = function(plot_list) {
  tmp_list = lapply(plot_list, function(x) x[order(-(x$methylation_score)),] )
  auc_all = lapply(tmp_list, function(dat)  {
    prauc <- function(dat) {
      x <- dat@x.values[[1]]
      y <- dat@y.values[[1]]
      idx <- which(is.nan(y))
      if (any(idx)) {
        x <- x[-idx]
        y <- y[-idx]
      }
      return(pracma::trapz(x, y))
    }
    dat$STATUS = dat$reported
    dat$PRED_CLASS = dat$predictions
    if(nrow(dat[dat$STATUS == st[1],]) > 1 & nrow(dat[dat$STATUS == st[2],]) > 1) {
      pred <- ROCR::prediction(dat$methylation_score, dat$STATUS ==
                                 st[1])
      c1 <- st[1]
      tp <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS ==
                  c1)
      tn <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS !=
                  c1)
      fp <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS !=
                  c1)
      fn <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS ==
                  c1)
      curRoc <- ROCR::performance(pred, "tpr", "fpr")
      curPr <- ROCR::performance(pred, "prec", "rec")
      auroc <- performance(pred, "auc")@y.values[[1]]
      aupr <- prauc(curPr)
      return_df = data.frame(auroc = auroc, aupr = aupr)
      return(return_df)
    } else {
      return(NULL)
    }
    
    
  })

  return(auc_all)
}
subgroups.order = c(rev(names(targ_stage)), rev(hr_groups), rev(targ_mammo),rev(names(targ_age_list)),'All')
subgroups = c(names(targ_age_list),names(targ_stage), hr_groups, targ_mammo,'All')
results_list = complete_list
combined.auc.summary= NULL
sens.ci.info = NULL

for (sg in subgroups) {
  print(sg)
  predList = list()
  featScores = list()
  
  for (i in 1:length(results_list)) {
    predictions = results_list[[i]]
    #subgrouping
    if (length(predictions) > 5) {
      
      if (sg %in% hr_groups) {
        #if hr
        if (sg == 'HR Positive, HER2 Negative') {
          targ_diag_time_samples = sample.information[sample.information$ER %in% c('positive','low positive') & sample.information$HER2_neu == 'negative' | sample.information$PR %in% c('positive','low positive')  & sample.information$HER2_neu == 'negative',]
        } else if (sg == 'HR Negative, HER2 Negative') {
          targ_diag_time_samples = sample.information[sample.information$ER == 'negative' & sample.information$HER2_neu == 'negative' & sample.information$PR == 'negative',]
        } else if (sg == 'HR Positive, HER2 Positive') {
          targ_diag_time_samples = sample.information[sample.information$ER %in% c('positive','low positive') & sample.information$HER2_neu == 'positive' | sample.information$PR  %in% c('positive','low positive') & sample.information$HER2_neu == 'positive',]
        } else if (sg == 'HR Negative, HER2 Positive') {
          targ_diag_time_samples = sample.information[sample.information$ER %in% c('negative') & sample.information$HER2_neu == 'positive' & sample.information$PR %in% c('negative'),]
        } else if (sg == 'HR Negative') {
          targ_diag_time_samples = sample.information[sample.information$ER %in% c('negative') & sample.information$PR %in% c('negative') ,]
        }else if (sg == 'HR Positive') {
          targ_diag_time_samples = sample.information[sample.information$ER %in% c('positive','low positive') | sample.information$PR %in% c('positive','low positive') ,]
        }  else if (sg == 'HER2 Positive') {
          targ_diag_time_samples = sample.information[ sample.information$HER2_neu == 'positive' ,]
        }  else if (sg == 'HER2 Negative') {
          targ_diag_time_samples = sample.information[ sample.information$HER2_neu == 'negative' ,]
        }
        predictions = predictions[predictions$GRP_Id %in% targ_diag_time_samples$GRP_Id | predictions$reported == 'control',]
      } else if (sg %in% targ_mammo) {
        if (sg == 'Pre-dx Mammo < 1 year') {
          targ.sample = sample.information[as.numeric(sample.information$mammogram.time.median) < 365, 'GRP_Id']
        } 

        predictions = predictions[predictions$GRP_Id %in% targ.sample  | predictions$reported == 'control',] 
        
      } else if (sg %in% names(targ_stage)) {
        targ_stages = sample.information[sample.information$Stage %in% targ_stage[[sg]],]
      
        predictions = predictions[predictions$GRP_Id %in% targ_stages$GRP_Id | predictions$reported == 'control',]
      } else if (sg == '3.7years') {
        targ.sample = sample.information[sample.information$followup_time > 365*3,]
        predictions = predictions[predictions$GRP_Id %in% targ.sample$GRP_Id ,] 
        
      } else if (sg %in% names(targ_age_list)) {
        sample.information$dx.age = ifelse(sample.information$group != 'control', sample.information$Age + abs(sample.information$diff_in_days)/365, sample.information$Age)
        if (sg == 'Dx Age Pre-50') {
          targ.ages = sample.information[sample.information$dx.age <= 50,] 
          
        } else if (sg == 'Dx Age Post-50') {
          targ.ages = sample.information[sample.information$dx.age > 50   ,]
          
        }
        predictions = predictions[predictions$GRP_Id %in% targ.ages$GRP_Id ,] 
        
      }
      
      predictions$reported = factor(predictions$reported, levels = st)
      predictions$predictions = factor(predictions$predictions, levels = st)
      
      #plotting performance
      pred <- predictions;
      # predictions table
      tmp <- pred[,c("GRP_Id","reported","methylation_score","predictions","seed")]
      predList[[i]] = tmp
    }
    

    
  }
  plot.list = predList
  
  return = auc_calc.seed(plot.list)
  return$average_auroc = return$auroc
  return = ddply(return, c('seed'),numcolwise(mean))
  ci.auc = confidence_interval(return$average_auroc)
  tmp.summary.auc = data.frame(mean.auc = mean(return$auroc),mean.ci.u = ci.auc[2],mean.ci.l = ci.auc[1], subgroup = sg)
  combined.auc.summary = rbind(combined.auc.summary,tmp.summary.auc)
  
  sens = sens.calc(plot.list, file_name = gsub(' ','_',sg), dir = getwd(), st = st, results_list_all[[length(results_list_all)]], sample_info, core.number = detectCores()-2,foldno = 10) 
  sens.ci.info = rbind(sens.ci.info, sens)

} 
combined.auc.summary$subgroup = gsub('\\.',' ',combined.auc.summary$subgroup)
sens.ci.info$subgroup =  gsub('_',' ',sens.ci.info$subgroup)

combined.overall = merge(combined.auc.summary, sens.ci.info, by = 'subgroup')
combined.auc.summary$subgroup = factor(combined.auc.summary$subgroup, levels =subgroups.order)
saveRDS(combined.overall, 'auc.performance.summary.RDS')

combined.auc.summary.melt = reshape2::melt(combined.overall[!combined.overall$subgroup %in% c('Stage II','Stage III-IV'),c('subgroup','mean.auc','mean.ci.u','mean.ci.l')],id.vars = 'subgroup')
pdf(paste0("all.subtype.auc.summary.pdf"),height = 3, width = 4)
ggplot(combined.auc.summary.melt[combined.auc.summary.melt$variable != 'mean.auc',], aes(x = subgroup, y = value, group= subgroup)) + geom_line(col = '#334E58',size = 1) + geom_point(data = combined.auc.summary.melt[combined.auc.summary.melt$variable == 'mean.auc',], aes(x = subgroup, y = value),size = 2, col = '#119DA4') + coord_flip() + theme_minimal()+
  theme(text = element_text(size=10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y  = element_blank()) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +ylab('AUROC') + xlab('Subgroup') 
dev.off()

combined.sens.summary.melt = reshape2::melt(combined.overall[combined.overall$subgroup %in% combined.auc.summary.melt$subgroup,c('subgroup','sens.99.mean','sens.99.ciu','sens.99.cil')],id.vars = 'subgroup')
pdf(paste0("all.subtype.sens.summary.pdf"),height = 3, width = 4)
ggplot(combined.sens.summary.melt[combined.sens.summary.melt$variable != 'sens.99.mean',], aes(x = subgroup, y = value, group= subgroup)) + geom_line(col = '#334E58',size = 1) + 
  geom_point(data = combined.sens.summary.melt[combined.sens.summary.melt$variable == 'sens.99.mean',], aes(x = subgroup, y = value),size = 2, col = '#119DA4') + 
  coord_flip() + theme_minimal()+
  theme(text = element_text(size=10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y  = element_blank()) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +ylab('Sensitivity at 99% Specificity') + xlab('Subgroup') 
dev.off()


combined.conc.summary.melt = reshape2::melt(combined.overall[combined.overall$subgroup %in% combined.auc.summary.melt$subgroup,c('subgroup','conc.mean','conc.ciu','con.cil')],id.vars = 'subgroup')
pdf(paste0("all.subtype.conc.summary.pdf"),height = 3, width = 4)
ggplot(combined.conc.summary.melt[combined.conc.summary.melt$variable != 'conc.mean',], aes(x = subgroup, y = value, group= subgroup)) + geom_line(col = '#334E58',size = 1) + 
  geom_point(data = combined.conc.summary.melt[combined.conc.summary.melt$variable == 'conc.mean',], aes(x = subgroup, y = value),size = 2, col = '#119DA4') + 
  coord_flip() + theme_minimal()+
  theme(text = element_text(size=10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y  = element_blank()) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +ylab('Concordance Index') + xlab('Subgroup') 
dev.off()



