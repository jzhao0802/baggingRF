################################################################################
# Lichao: 
# Previously the function was defined inside msOnTest_sep_v2. This is not 
# recommended especially it's a little complex. 
################################################################################

FindPrecisionAtGivenRecall <- function(recall, rec_prec)
{
  idx=which(abs(rec_prec[, 1]-recall)==min(abs(rec_prec[, 1]-recall), na.rm=T))[1]
  ################################################################################
  # Lichao: 
  # Safe guard added. 
  ################################################################################
  if (length(idx) == 0)
    stop("Error! The target recall level doesn't exist in the recal-precision table. ")
  prec_sel <- rec_prec[idx, 2]
  return(prec_sel)
}


msOnTest_sep_v2 <- function(pred, response, recall_tar){
  #pred <- apply(pred, 1, mean, na.rm=T)
  predobj <- prediction(pred, response)
  #add plot
  perf <- performance(predobj, 'ppv', 'sens') # added by jie for recall-precision plot.
  recall <- perf@x.values[[1]]
  precision <- perf@y.values[[1]]
  auc <- performance(predobj, 'auc')@y.values[[1]]
  rec_prec <- data.frame(recall=recall, precision=precision)
  rec_prec_omitMiss <- rec_prec[complete.cases(rec_prec),]
  aupr <- trapz(rec_prec_omitMiss$recall, rec_prec_omitMiss$precision)
  bucket <- cut(recall, breaks=seq(0, 1, 0.01), include.lowest=T,right=F)
  rec_prec_byBucket <- aggregate(rec_prec, by=list(bucket), function(i) mean(i, na.rm=T))
  #write.csv(rec_prec_byBucket, paste('Curve_dong.csv', sep=''), 
  #         row.names=F, quote=T)
  # plot
#     pdf(file=paste(plots_path, 'recall-precision curve on test on 200K.pdf', sep=''))
#     plot(recall, precision, type='l', main=paste('recall-precision curve(5 simulations)'))
#     plot(perf)
  # dev.off()
  
  ##in simulation
  temp4 <- unlist(lapply(recall_tar, FindPrecisionAtGivenRecall, rec_prec=rec_prec))    
  
  ##end
  ms <- c(auc, aupr, temp4)
  names(ms) <- c('auc',"aupr", paste("PPV(recall=", recall_tar,')', sep=''))
  
  return(list(ms=ms, curve=rec_prec_byBucket, rec_prec=rec_prec))
    
}

msOnTest_sep_v3 <- function(pred, response, recall_tar, simu){
    #pred <- apply(pred, 1, mean, na.rm=T)
    predobj <- prediction(pred, response)
    #add plot
    perf <- performance(predobj, 'ppv', 'sens') # added by jie for recall-precision plot.
    recall <- perf@x.values[[1]]
    precision <- perf@y.values[[1]]
    auc <- performance(predobj, 'auc')@y.values[[1]]
    rec_prec <- data.frame(recall=recall, precision=precision)
    rec_prec_omitMiss <- rec_prec[complete.cases(rec_prec),]
    aupr <- trapz(rec_prec_omitMiss$recall, rec_prec_omitMiss$precision)
    bucket <- cut(recall, breaks=seq(0, 1, 0.01), include.lowest=T,right=F)
    rec_prec_byBucket <- aggregate(rec_prec, by=list(bucket), function(i)mean(i, na.rm=T))
    plot(recall, precision, type='l', main=paste0('recall-precision curve simulation', simu))
    #plot(perf)
    #     dev.off()
    
    ##in simulation
    temp4 <- unlist(lapply(recall_tar, function(X){
        #idx <- sample(rep(which(abs(rec_prec[, 1]-X)==min(abs(rec_prec[, 1]-X), na.rm=T)), 2), 1)
        idx=which(abs(rec_prec[, 1]-X)==min(abs(rec_prec[, 1]-X), na.rm=T))[1]
        prec_sel <- rec_prec[idx, 2]
        return(prec_sel)
    }))    
    
    ##end
    ms <- c(auc, aupr, temp4)
    names(ms) <- c('auc',"aupr", paste("PPV(recall=", recall_tar,')', sep=''))
    
    return(list(ms=ms, curve=rec_prec_byBucket, rec_prec=rec_prec))
    
}


split_simulations <- function(n.simu, haeFile, nonhaeFile, haeDir, nonhaeDir, outDir, iters){
  
  dat_hae_1111_rf_nov26 <- 
      read.csv(paste0(haeDir, "dat_hae_1111_rf_nov26_flag.csv")
               , header=T, sep=",", check.names=F)
  
  
  hae <- dat_hae_1111_rf_nov26 %>% 
    filter(FLAG_==0) %>% 
    select(-c(FLAG_, REGION)) 
  hae$GENDERM <- ifelse(hae$GENDER=='M', 1, 0)
  hae$GENDER <- NULL
  
  ################################################################################
  # Lichao: 
  # This line was  
  # HAE_ptid <- readRDS(paste0(haeDir, haeFile, '.RDS'))
  # but there was not such a RDS file at all! Now I have created such an RDS file 
  # in order to execute the code. 
  ################################################################################
  HAE_ptid <- readRDS(paste0(haeDir, haeFile, '.RDS'))
  ################################################################################
  # Lichao: 
  # The following line was  
  # hae <- hae[hae$PATIENT_ID %in% HAE_ptid]
  # It wouldn't run since hae is a data frame. 
  ################################################################################
  hae <- hae[hae$PATIENT_ID %in% HAE_ptid, ]
  # nonhae
  if(nonhaeFile=='nonhae_Dong200K'){
      load(paste(haeDir, "dat_nonhae_1111_rf.RData", sep=""))
      
      nonhae <- dat_nonhae_1111_rf %>% 
        mutate(LOOKBACK_DAYS=lookback_days) %>% 
        select(-c(lookback_days, REGION))
      nonhae <- nonhae %>% 
        mutate(hae_patient_id=HAE_PATIENT)  %>% 
        select(-HAE_PATIENT) 
      nonhae$hae_patient_id <- as.numeric(nonhae$hae_patient_id)
      nonhae$GENDERM <- ifelse(nonhae$GENDER=='M', 1, 0)
      nonhae$GENDER <- NULL
  }else if(grepl('nonhae_200K_A\\dE\\d', nonhaeFile, ignore.case = T)){
      dat_nonhae <- 
        read.table(paste0(nonhaeDir, nonhaeFile, ".csv"), 
                   sep=',', stringsAsFactors = F, head=T)
      nonhae <- dat_nonhae %>% 
        mutate(LOOKBACK_DAYS=lookback_days) %>% 
        select(-c(lookback_days))
  }else if(grepl('300K|nonhae_A_200K_A\\dE\\d', nonhaeFile, ignore.case = T)){
      dat_nonhae <- 
        read.table(paste0(nonhaeDir, nonhaeFile, ".csv"), 
                   sep=',', stringsAsFactors = F, head=T)
      dat_nonhae <- dat_nonhae %>% 
        mutate(LOOKBACK_DAYS=lookback_days) %>% 
        select(-c(lookback_days))
      
  }else{
      stop('the wrong nonhae input!\n')
  }
  
  ################################################################################
  # Lichao: 
  # Previously the following line was only
  # outDir <- paste0(outDir, nonhaeFile, '&', haeFile, '\\')
  # which was not compatible to code in other functions that reads data from the dir
  ################################################################################
  
  ################################################################################
  # Jie: 
  # code back to the following line. because models with different iterations share
  # the same model data.
  # outDir <- paste0(outDir, nonhaeFile, '&', haeFile, '\\')
  ################################################################################
  
  outDir <- paste0(outDir, nonhaeFile, '&', haeFile, "/")
  if(!dir.exists(outDir)) 
    dir.create(outDir, showWarnings = T, recursive = TRUE)
  
  set.seed(20)
  
  for (simu in 1:n.simu)
  {
    if(grepl('300K|nonhae_A_200K_A\\dE\\d', nonhaeFile, ignore.case = T)){
      ################################################################################
      # Lichao: 
      # tr_idx <- createFolds(hae$PATIENT_ID, 5, returnTrain=T)[[simu]]
      # where the number '5' is a constant, which should be replaced by n.simu.
      ################################################################################
      tr_idx <- createFolds(hae$PATIENT_ID, n.simu, returnTrain=T)[[simu]]
      dat_hae_trn <- hae[tr_idx, ]
      dat_hae_tst <- hae[-tr_idx, ]
      
      tr_idx_nonhae <- createFolds(dat_nonhae$patient_id, n.simu, returnTrain=T)[[simu]]
      dat_nonhae_trn <- dat_nonhae[tr_idx_nonhae, ]
      dat_nonhae_tst <- dat_nonhae[-tr_idx_nonhae, ]
      
    }else{
      tr_idx <- createFolds(hae$PATIENT_ID, n.simu, returnTrain=T)[[simu]]
      dat_hae_trn <- hae[tr_idx, ]
      dat_hae_tst <- hae[-tr_idx, ]
      
      dat_nonhae_trn <- nonhae[nonhae$hae_patient_id %in% dat_hae_trn$PATIENT_ID,]
      dat_nonhae_tst <- nonhae[nonhae$hae_patient_id %in% dat_hae_tst$PATIENT_ID,]
    }
    
    dat_nonhae_trn$HAE <- 0
    dat_nonhae_tst$HAE <- 0
    
#     outDirThisSimu <- paste0(outDir, "simu", simu, "/")
    if(!dir.exists(outDir)) 
      dir.create(outDir, showWarnings = T, recursive = TRUE)
    
    saveRDS(dat_hae_trn, file=paste0(outDir, "dat_hae_trn_simu", simu, ".RDS"))
    saveRDS(dat_hae_tst, file=paste0(outDir, "dat_hae_tst_simu", simu, ".RDS"))
    saveRDS(dat_nonhae_trn, file=paste0(outDir, "dat_nonhae_trn_simu", simu, ".RDS"))
    saveRDS(dat_nonhae_tst, file=paste0(outDir, "dat_nonhae_tst_simu", simu, ".RDS"))
#     save(dat_hae_trn, dat_hae_tst, dat_nonhae_trn , dat_nonhae_tst
#          , file=paste(outDirThisSimu, "dat_hae_trn_tst_split_simu", simu, ".RData", sep=""))
  }
}

run_split <- function(n.simu, haeFile, nonhaeFile, haeDir, nonhaeDir, outDir, iters){
  ################################################################################
  # Lichao: 
  # There's no need to do it 5 times. Once is enough. 
  ################################################################################
  
  split_simulations(n.simu=n.simu, haeFile=haeFile, nonhaeFile=nonhaeFile, haeDir=haeDir, 
                    nonhaeDir=nonhaeDir, outDir=outDir, iters=iters)
}



# run the bagging forest model

run_bagging_rf_par_forTrainigFit <- 
  function(simu, dir, lasso_rf_iters, nonhaeFile, haeFile)
{
  ################################################################################
  # Lichao: 
  # Previously there was this line: 
  # dir <- paste0(dir, nonhaeFile, '&', haeFile, '\\')
  # and this shouldn't be there since the same was done in the function that
  # calls this one. Otherwise the directory won't be correct. 
  ################################################################################
 
      ################################################################################
      # Jie: 
      # read in the model data using dir instead of outDir to correspond to the split step
      # & there is no need to define the dataDir in line 239
      ################################################################################
      
  outDir  <- paste0(dir, 'iters=', lasso_rf_iters, '\\simu', simu, '\\')
  if(!dir.exists(outDir)) 
    dir.create(outDir, showWarnings = T, recursive = TRUE)
    
  dat_hae_trn <- readRDS(file=paste0(dir, "dat_hae_trn_simu", simu, ".RDS"))
  dat_hae_tst <- readRDS(file=paste0(dir, "dat_hae_tst_simu", simu, ".RDS"))
  dat_nonhae_trn <- readRDS(file=paste0(dir, "dat_nonhae_trn_simu", simu, ".RDS"))
  dat_nonhae_tst <- readRDS(file=paste0(dir, "dat_nonhae_tst_simu", simu, ".RDS"))
  
    
    
    
  #########################################################################
  ### Model training (all training data)
  #########################################################################
  
  Sys.time()->start
  
  # (1) Underbagging LASSO / Random Forest
  trn_ans_lasso_rf <- 
    undbag_lasso_rf(rf_formula=NA, dat_pos=dat_hae_trn, dat_neg=dat_nonhae_trn, 
                    iters=lasso_rf_iters, mtry=25, ntree=300)
  # trn_undbag_lasso_fit = trn_ans_lasso_rf$undbag_lasso_fit
  trn_undbag_rf_fit = trn_ans_lasso_rf$undbag_rf_fit
  
  print(Sys.time()-start)
  
  #
  ################################################################################
  # Lichao: 
  # No longer named '*Mar31*'
  ################################################################################
  saveRDS(trn_undbag_rf_fit, file=paste0(outDir, "trn_rf_fit.RDS"))
  
  # save(dat_hae_trn, dat_nonhae_trn, trn_undbag_rf_fit, file=paste(outDir, "trn_rf_fit_Mar31.RData", sep=""))
  
  
  #########################################################################
  ### Apply all-training model to testing data
  #########################################################################

  x_hae_trn <- dat_hae_trn[,-match(c('PATIENT_ID', 'HAE'), names(dat_hae_trn))]
  y_hae_trn <- dat_hae_trn[,match('HAE', names(dat_hae_trn))]
  x_hae_tst <- dat_hae_tst[,-match(c('PATIENT_ID', 'HAE'), names(dat_hae_trn))]
  y_hae_tst <- dat_hae_tst[,match('HAE', names(dat_hae_trn))]
  
  dat_nonhae_trn$HAE <- 0
  dat_nonhae_tst$HAE <- 0
  
  x_nonhae_trn <- dat_nonhae_trn[
    , -match(grep('patient_id|hae', names(dat_nonhae_trn), valu=T, perl=T, ignore.case=T), 
             names(dat_nonhae_trn))
    ]
  y_nonhae_trn <- dat_nonhae_trn[,match('HAE', names(dat_nonhae_trn))]
  x_nonhae_tst <- dat_nonhae_tst[
    , -match(grep('patient_id|hae', names(dat_nonhae_tst), valu=T, perl=T, ignore.case = T), 
             names(dat_nonhae_tst))
    ]
  y_nonhae_tst = dat_nonhae_tst[,match('HAE', names(dat_nonhae_tst))]
  
  x_nonhae_tst <- x_nonhae_tst[, match(names(x_hae_tst), names(x_nonhae_tst))]
  x_nonhae_trn <- x_nonhae_trn[, match(names(x_hae_trn), names(x_nonhae_trn))]
  x_tst <- rbind(x_hae_tst, x_nonhae_tst)
  y_tst <- c(y_hae_tst, y_nonhae_tst)
  dat_tst <- data.frame(y_tst, x_tst)
  names(dat_tst)[1]='HAE'
  
  tst_label <- y_tst
  
  tst_prob_rf <- rep(0, length(y_tst))
  
  Sys.time()->start
  
  # Bagging LASSO, Bagging Random Forest
  for (i in 1:lasso_rf_iters){
    tst_prob_rf <- tst_prob_rf + 
      predict(trn_undbag_rf_fit[[i]], dat_tst[, -match('HAE', names(dat_tst))], type = "prob")[,2]/lasso_rf_iters
  }
  
  print(Sys.time()-start)
  # Time difference of 12.26336 mins
  
  ################################################################################
  # Lichao: 
  # Previously responses and predictions are saved saparately, which leads to
  # later loading training + testing data all over again just for combining
  # the responses and predictions. 
  # Now they're combined in this function. 
  ################################################################################
  result <- cbind(tst_label, tst_prob_rf)
  colnames(result) <- c("resp", "prob")
  saveRDS(result, file=paste0(outDir, "tst_rf_prob.RDS"))
  # save(tst_label, tst_prob_rf, file=paste(outDir, "tst_rf_prob.RData", sep=""))
    
}


run_bagging_rf <- function(n.simu, wk_dir, dir, lasso_rf_iters, nonhaeFile, haeFile)
{
  dir <- paste0(dir, nonhaeFile, '&', haeFile, '\\')
  
  trace_path <- paste0(dir, 'iters=', lasso_rf_iters, '\\')
  if(!dir.exists(trace_path)) 
    dir.create(trace_path, showWarnings = T, recursive = TRUE)
  
  traceFile <- paste0(trace_path, 'traceFile.csv')
  cat(file=traceFile, append=T, 'parallele on n.simu simulation start!\n')
  #     wk_dir = "D:\\jzhao\\Shire_followup\\02_Code\\HAE_R_codes_Dec15\\"
  
  sfInit(parallel=TRUE, cpus=n.simu, type='SOCK')
  #sfSource("F:\\Jie\\Shire\\03_code\\subFunc_v3.R")
  
  cat(file=traceFile, append=TRUE, 'n.simu simulations parallel sfExport running!\n')
  sfExport('dir', 'wk_dir', "nonhaeFile", "haeFile")
  sfExport('undbag_lasso_rf')
  
  sfSource(paste0(wk_dir, "scripts/loadpackage.R"))
  # Auxiliary functions
  sfSource(paste0(wk_dir, "functions/auxfunctions.R"))
  # 
  sfSource(paste0(wk_dir, "functions/funs_baggingRF.R"))
  
  sfClusterEval(library(ggplot2))
  sfClusterEval(library(ROCR))
  sfClusterEval(library(PRROC))
  # sfClusterEval(library(FSelector))
  sfClusterEval(library(randomForest))
  sfClusterEval(library(caret))
  sfClusterEval(library(e1071))
  sfClusterEval(library(reshape2))
  sfClusterEval(library(sqldf))
  sfClusterEval(library(glmnet))
  sfClusterEval(library(caTools))
  # sfClusterEval(library(gbm))
  # sfClusterEval(library(xlsx))
  sfClusterEval(library(dplyr))   
  ################################################################################
  # Lichao: 
  # Arguments "nonhaeFile" and "haeFile" were not passed to the function previously. 
  ################################################################################
  temp <- sfClusterApplyLB(1:n.simu, run_bagging_rf_par_forTrainigFit, dir, 
                           lasso_rf_iters, nonhaeFile, haeFile)
  #save(pred_ts_allSim, file=paste0(modelDir, '//pred_allSim.RData'))
  sfStop()
  # run_bagging_rf_par_forTrainigFit(1, dir, lasso_rf_iters, nonhaeFile, haeFile)
  #     cat(unlist(lapply(pred_ts_allSim, length)), '\n')
}

################################################################################
# Lichao: 
# Previously the following function was defined inside of lapply in get_perf_allSimu, 
# which is not recommended. 
################################################################################

ConcatenatePreds <- function(simu, dir, iters)
{
  ################################################################################
  # Lichao: 
  # Previously there was this line: 
  # load(paste0(dir, "dat_hae_trn_tst_split_simu", i, ".RData"))
  # which doesn't correspond to the correct directory defined in training. 
  ################################################################################
  
  ################################################################################
  # Lichao: 
  # Previously there were also lines reading training and testing data. 
  # 1. Loading training data is definitely not needed. This is why we don't want
  #    to save multiple datasets in one RData file. 
  # 2. Loading testing data is also not necessary because it was for generating
  #    the responses, which now has been saved with the predictions in run_bagging_rf_par_forTrainigFit
  ################################################################################
  
  ################################################################################
  # Lichao: 
  # Previously the following line was
  # saveRDS(resp, file = paste0(dir, "resp_simu", simu, ".RData"))
  # which despite using saveRDS, the target file has an extention of .RData. 
  # Also I don't think it's necessary to save a separate file of resp_simu? 
  ################################################################################
  # saveRDS(resp, file = paste0(dir, "resp_simu", simu, ".RDS"))
  tst_prob_rf <- readRDS(paste0(dir, "iters=", iters, "\\simu", simu, "\\tst_rf_prob.RDS"))
  
  return(tst_prob_rf)
}

get_perf_allSimu <- function(dir, iters, n.simu, recall_tar, haeFile, nonhaeFile){
  
  dir <- paste0(dir, nonhaeFile, '&', haeFile, '/')
  
  temp <- lapply(1:n.simu, ConcatenatePreds, dir=dir, iters=iters)
  
  resp_pred <- ldply(temp, rbind)
  ################################################################################
  # Lichao: 
  # Previously the following line was
  # saveRDS(resp_pred, paste0(dir, 'iters=', iters, '\\resp_pred.RData'))
  # which despite using saveRDS, the target file has an extention of .RData. 
  ################################################################################
  saveRDS(resp_pred, paste0(dir, 'iters=', iters, '\\resp_pred.RDS'))
  temp1 <- msOnTest_sep_v2(resp_pred[, 2], resp_pred[, 1], recall_tar)
  perf <- temp1$ms
  write.csv(perf, paste0(dir, 'iters=', iters, '\\performance_onAllSimu.csv'))
  return(perf)
}


get_perf_allSimu_forPRcurve <- function(dataDir, iters, n.simu, recall_tar){
    
    dataDir1 <- outdir <- paste0(dataDir, 'iters=', iters, '\\')
        
    # saveRDS(resp_pred, paste0(outdir, 'iters=', iters, '\\resp_pred.RData'))
    resp_pred <- readRDS(paste0(dataDir1, 'resp_pred.RData'))
    temp1 <- msOnTest_sep_v2(resp_pred[, 2], resp_pred[, 1], recall_tar)
    perf <- temp1$ms
    recPrec <- temp1$rec_prec
    write.csv(recPrec, paste0(outDir, "recall_precision_on200K.csv"))
    write.csv(temp1$curve, paste0(outDir, "recall_precision_byBucket_on200K.csv"), 
              row.names=F, quote=T)
    pdf(file=paste(outDir, 'recall-precision curve on test 200K.pdf', sep=''))
    
    plot(recall=recPrec[, 1], precision=recPrec[, 2], type='l', main='PR curve on 200K test data(5 simulations)')
    dev.off()
    # write.csv(perf, paste0(outdir, 'iters=', iters, '\\performance_onAllSimu.csv'))
    return(perf)
}


get_perf_3M_par_forPRcurve <- function(simu, dir, lasso_rf_iters, recall_tar, fileNm_3M, path_3M){
    
    dataDir <- outDir <- paste0(dir, 'iters=', lasso_rf_iters, '\\simu', simu, '\\')

    load(paste(dataDir, "tst_rf_prob_haeTs&3M.RData", sep=""))
    resp <- c(rep(1, length(tst_prob_rf_hae)), rep(0, length(tst_prob_rf)))
    
    pred <- c(tst_prob_rf_hae, tst_prob_rf)
    
    perf_result <- msOnTest_sep_v3(pred, resp, recall_tar=seq(0.5, 0.05, -0.05), dataDir, simu)
    recPrec <- perf_result$rec_prec
    write.csv(recPrec, paste0(plots_path, "recall_precision_sim", simu, '.csv'))
    write.csv(perf_result$curve, paste0(plots_path, "recall_precision_byBucket_sim", simu, '.csv'), 
              row.names=F, quote=T)
    
    #     write.csv(perf_result$ms, paste0(plots_path, 'perf_on3M.csv'))
    result <- c(simu=simu, perf_result$ms)
    return(result)
}
run_perf_3M_forPRcurve <- function(dir, wk_dir, lasso_rf_iters, n.simu, recall_tar, fileNm_3M, path_3M){
    trace_path <- paste0(dir, 'iters=', lasso_rf_iters, '\\')
    if(!dir.exists(trace_path)) dir.create(trace_path, showWarnings = T, recursive = TRUE)
    #     plot
    pdf(file=paste(trace_path, 'recall-precision curve on test 3M.pdf', sep=''))
    par(mfrow=c(3,2))
    par(pty='m')
    par(cex.lab=1.2, cex.axis=0.9)
    
    for(simu in 1:n.simu){
        temp=get_perf_3M_par_forPRcurve(simu
                                        , trace_path
                                        , lasso_rf_iters
                                        , recall_tar
                                        , fileNm_3M
                                        , path_3M
        )
        #         return('finished\n')
    }
    
    dev.off()
}




get_perf_3M_par <- function(simu, dir, lasso_rf_iters, recall_tar, fileNm_3M, 
                            path_3M, nonhaeFile, haeFile){
  ################################################################################
  # Lichao: 
  # Previously there was this line: 
  # dir <- paste0(dir, nonhaeFile, '&', haeFile, '\\')
  # and this shouldn't be there since the same was done in the function that
  # calls this one. Otherwise the directory won't be correct. 
  ################################################################################
    
  outDir <- dataDir <- paste0(dir, 'iters=', lasso_rf_iters, '\\simu', simu, '\\')
  # if(!dir.exists(plots_path)){dir.create(plots_path, showWarnings = T, recursive = T, model='0777')}
  x_3M <- read.table(paste0(path_3M, fileNm_3M, ".csv"), sep=',', 
                     stringsAsFactors = F, head=T)
  if('lookback_days' %in% names(x_3M)){
    x_3M <- x_3M %>% 
      mutate(LOOKBACK_DAYS=lookback_days) %>% 
      select(-patient_id) %>% 
      select(-lookback_days)
  }
  

  
  ################################################################################
  # Lichao: 
  # Instead of loading all 4 datasets and removing three of them, only one is
  # loaded now. 
  ################################################################################
  dat_hae_tst <- readRDS(file=paste0(dataDir, "dat_hae_tst_simu", simu, ".RDS"))
  # load(paste0(dataDir, 'dat_hae_trn_tst_split_simu', simu, '.RData'))
  # rm(dat_hae_trn, dat_nonhae_trn, dat_nonhae_tst)
  gc()
  
  tst_prob_rf <- rep(0, nrow(x_3M))
  tst_prob_rf_hae <- rep(0, nrow(dat_hae_tst)) #218
  
  Sys.time()->start
  
  # Bagging LASSO, Bagging Random Forest
  trn_undbag_rf_fit <- readRDS(file=paste0(dataDir, 'trn_rf_fit.RDS'))
  # load(paste0(dataDir, 'trn_rf_fit_Mar31.RData'))
  for (i in 1:lasso_rf_iters){
    # 	tst_prob_lasso = tst_prob_lasso + predict(trn_undbag_lasso_fit[[i]], as.matrix(x_tst), s="lambda.min", type="response")/lasso_rf_iters
    cat('i=', i, '!\n')
    tst_prob_rf <- tst_prob_rf + 
      predict(trn_undbag_rf_fit[[i]], x_3M, type = "prob")[,2]/lasso_rf_iters
    tst_prob_rf_hae <- tst_prob_rf_hae + 
      predict(trn_undbag_rf_fit[[i]], dat_hae_tst[,-match(c('HAE', 'PATIENT_ID'), names(dat_hae_tst))], type = "prob")[,2]/lasso_rf_iters
  }
  
  print(Sys.time()-start)
  # Time difference of 12.26336 mins
  
  ################################################################################
  # Lichao: 
  # Previously there was this following line: 
  # save( tst_prob_rf_hae, tst_prob_rf, file=paste(outDir, "tst_rf_prob_haeTs&3M.RData", sep=""))
  # Now the two predicted vectors are combined before being saved (with the responses). 
  ################################################################################
  # save( tst_prob_rf_hae, tst_prob_rf, file=paste(outDir, "tst_rf_prob_haeTs&3M.RData", sep=""))
  
  resp <- c(rep(1, length(tst_prob_rf_hae)), rep(0, length(tst_prob_rf)))
  
  pred <- c(tst_prob_rf_hae, tst_prob_rf)
  resp_pred <- cbind(resp, pred)
  colnames(resp_pred) <- c("resp", "prob")
  saveRDS(resp_pred, file=paste0(outDir, "tst_rf_prob_haeTs&3M.RDS"))
  
  ################################################################################
  # Lichao: 
  # Previously the following line was: 
  # perf_result <- msOnTest_sep_v2(pred, resp, recall_tar=seq(0.5, 0.05, -0.05), oudDir)
  # But function msOnTest_sep_v2 doesn't has a parameter 'outDir'. This is also 
  # inconsistent with other places where the same function is called. 
  ################################################################################
  perf_result <- msOnTest_sep_v2(pred, resp, recall_tar=seq(0.5, 0.05, -0.05))
  write.csv(perf_result$ms, paste0(outDir, 'perf_on3M.csv'))
  result <- c(simu=simu, perf_result$ms)
  return(result)
}

run_perf_3M <- function(dir, wk_dir, lasso_rf_iters, n.simu, recall_tar, 
                        fileNm_3M, path_3M, haeFile, nonhaeFile){
  
  dir <- paste0(dir, nonhaeFile, '&', haeFile, '\\')
  
  trace_path <- paste0(dir, 'iters=', lasso_rf_iters, '\\')
  if(!dir.exists(trace_path)) 
    dir.create(trace_path, showWarnings = T, recursive = TRUE)
  
  traceFile <- paste0(trace_path, 'traceFile_pred3M.csv')
  cat(file=traceFile, append=T, 'parallele on n.simu simulation start!\n')
  
  sfInit(parallel=TRUE, cpus=n.simu, type='SOCK')
  #sfSource("F:\\Jie\\Shire\\03_code\\subFunc_v3.R")
  
  cat(file=traceFile, append=TRUE, 'n.simu simulations parallel sfExport running!\n')
  sfExport('dir', 'wk_dir', "nonhaeFile", "haeFile", "lasso_rf_iters",
           "recall_tar", "fileNm_3M", "path_3M")
  sfExport('get_perf_3M_par', 'msOnTest_sep_v2')
  
  sfSource(paste0(wk_dir, "scripts/loadpackage.R"))
  # Auxiliary functions
  sfSource(paste0(wk_dir, "functions/auxfunctions.R"))
  # 
  sfSource(paste0(wk_dir, "functions/funs_baggingRF.R"))
  
  sfClusterEval(library(ggplot2))
  sfClusterEval(library(ROCR))
  sfClusterEval(library(PRROC))
  # sfClusterEval(library(FSelector))
  sfClusterEval(library(randomForest))
  sfClusterEval(library(caret))
  sfClusterEval(library(e1071))
  sfClusterEval(library(reshape2))
  sfClusterEval(library(sqldf))
  sfClusterEval(library(glmnet))
  sfClusterEval(library(caTools))
  # sfClusterEval(library(gbm))
  # sfClusterEval(library(xlsx))
  sfClusterEval(library(dplyr))    
  ################################################################################
  # Lichao: 
  # Arguments "nonhaeFile" and "haeFile" were not passed to the function previously. 
  ################################################################################
  temp <- sfClusterApplyLB(1:n.simu, get_perf_3M_par, dir, lasso_rf_iters, 
                           recall_tar, fileNm_3M, path_3M, nonhaeFile, haeFile)
  sfStop()
  
    # get_perf_3M_par(1, dir, lasso_rf_iters, recall_tar, fileNm_3M, path_3M, nonhaeFile, haeFile)
#     print("get_perf_3M_par finished.")
  ms_allSimu <- ldply(temp, quickdf)
  
  ################################################################################
  # Lichao: 
  # Previously the following line was: 
  # write.csv(ms_allSimu, pate0(trace_path, 'perf_on3M_allSimu.csv'))
  # It should be 'paste0' instead of 'pate0'
  # Seriously, we need to be way more careful when coding. 
  ################################################################################
  write.csv(ms_allSimu, paste0(trace_path, 'perf_on3M_allSimu.csv'))
  return(ms_allSimu)
    
}


# get_perf_2.5M_par <- function(simu, dir, lasso_rf_iters, recall_tar, fileNm_2.5M, path_2.5M){
#     
#     outDir <- dataDir  <- paste0(dir, 'iters=', lasso_rf_iters, '\\simu', simu, '\\')
#     if(!dir.exists(outDir)){dir.create(outDir, showWarnings = T, recursive = T, model='0777')}
# 
#     x_3M <- lapply(1:51, function(i) {
#         load(paste0(path_2.5M, fileNm_2.5M, '_', i, '.RData'))
#         return(adj_ppv_samp)
#     })
#     x_3M <- ldply(x_3M, rbind)
#     
#     load(paste0(dataDir, 'trn_rf_fit_Mar31.RData'))
#     vars_rf <- rownames(trn_undbag_rf_fit[[1]]$importance)
#     if('LOOKBACK_DAYS' %in% vars_rf & 'lookback_days' %in% names(x_3M)){
#         x_3M <- x_3M %>% mutate(LOOKBACK_DAYS=lookback_days) %>% select(-PATIENT_ID) %>% select(-lookback_days)
#         
#     }
#     load(paste0(dataDir, 'dat_hae_trn_tst_split_simu', simu, '.RData'))
#     rm(dat_hae_trn, dat_nonhae_trn, dat_nonhae_tst)
#     gc()
#     
#     tst_prob_rf = rep(0, nrow(x_3M))
#     tst_prob_rf_hae = rep(0, nrow(dat_hae_tst)) #218
#     
#     Sys.time()->start
#     
#     # Bagging LASSO, Bagging Random Forest
#     for (i in 1:lasso_rf_iters){
#         # 	tst_prob_lasso = tst_prob_lasso + predict(trn_undbag_lasso_fit[[i]], as.matrix(x_tst), s="lambda.min", type="response")/lasso_rf_iters
#         cat('i=', i, '!\n')
#         tst_prob_rf = tst_prob_rf + predict(trn_undbag_rf_fit[[i]], x_3M, type = "prob")[,2]/lasso_rf_iters
#         tst_prob_rf_hae = tst_prob_rf_hae + predict(trn_undbag_rf_fit[[i]], dat_hae_tst[,-match(c('HAE', 'PATIENT_ID'), names(dat_hae_tst))], type = "prob")[,2]/lasso_rf_iters
#     }
#     
#     print(Sys.time()-start)
#     # Time difference of 12.26336 mins
#     
#     save( tst_prob_rf_hae, tst_prob_rf, file=paste(outDir, "tst_rf_prob_haeTs&2.5M.RData", sep=""))
#     
#     resp <- c(rep(1, length(tst_prob_rf_hae)), rep(0, length(tst_prob_rf)))
#     
#     pred <- c(tst_prob_rf_hae, tst_prob_rf)
#     
#     perf_result <- msOnTest_sep_v2(pred, resp, recall_tar=seq(0.5, 0.05, -0.05))
#     write.csv(perf_result$ms, paste0(outDir, 'perf_on2.5M.csv'))
#     result <- c(simu=simu, perf_result$ms)
#     return(result)
# }
# 
# run_perf_2.5M <- function(dir, wk_dir, lasso_rf_iters, n.simu, recall_tar, fileNm_2.5M, path_2.5M){
#     trace_path <- paste0(dir, 'iters=', lasso_rf_iters, '\\')
#     # if(!dir.exists(trace_path)) dir.create(trace_path, showWarnings = T, recursive = TRUE)
#     
#     traceFile <- paste0(trace_path, 'traceFile_pred2.5M.csv')
#     cat(file=traceFile, append=T, 'parallele on n.simu simulation start!\n')
#     #     wk_dir = "D:\\jzhao\\Shire_followup\\02_Code\\HAE_R_codes_Dec15\\"
#     
#     sfInit(parallel=TRUE, cpus=n.simu, type='SOCK')
#     #sfSource("F:\\Jie\\Shire\\03_code\\subFunc_v3.R")
#     
#     cat(file=traceFile, append=TRUE, 'n.simu simulations parallel sfExport running!\n')
#     sfExport(  'outDir', 'wk_dir')
#     sfExport('get_perf_2.5M_par', 'msOnTest_sep_v2'
#     )
#     
#     sfSource(paste0(wk_dir, "scripts/loadpackage.R"))
#     # Auxiliary functions
#     sfSource(paste0(wk_dir, "functions/auxfunctions.R"))
#     # 
#     sfSource(paste0(wk_dir, "functions/funs_baggingRF.R"))
#     
#     sfClusterEval(library(ggplot2))
#     sfClusterEval(library(ROCR))
#     sfClusterEval(library(PRROC))
#     # sfClusterEval(library(FSelector))
#     sfClusterEval(library(randomForest))
#     sfClusterEval(library(caret))
#     sfClusterEval(library(e1071))
#     sfClusterEval(library(reshape2))
#     sfClusterEval(library(sqldf))
#     sfClusterEval(library(glmnet))
#     sfClusterEval(library(caTools))
#     # sfClusterEval(library(gbm))
#     # sfClusterEval(library(xlsx))
#     sfClusterEval(library(dplyr))    
#     temp <- sfClusterApplyLB(1:n.simu, get_perf_2.5M_par
#                              , dir
#                              , lasso_rf_iters 
#                              , recall_tar
#                              , fileNm_2.5M
#                              , path_2.5M
#     )
#     sfStop()
#     
#     ms_allSimu <- ldply(temp, quickdf)
#     
#     return(ms_allSimu)
#     
# }
