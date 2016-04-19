# main function for running bagging forest and performance
rm(list=ls())


source("scripts/loadpackage.R")
# Load R packages
install_load("ggplot2", 'ROCR', 'PRROC', 'randomForest', 'caret', 'e1071', 'reshape2', 'sqldf', 'glmnet', 'caTools'
             , "gbm", "dplyr", "snowfall", 'plyr')
# Auxiliary functions
source("functions/auxfunctions.R")
# 
source("functions/funs_baggingRF.R")


# setup the R work directory
wk_dir <- "./"
setwd(wk_dir)

# 
n.simu = 5
recall_tar <- seq(0.5, 0.05, -0.05)
iters <- 20

################################################################################
# Previously the file names are given as follows: 
# nonhaeFile='nonhae_200K_II2.a'
# haeFile="II2.a"
# However there are no such files in the given path: 
# haeDir="F:\\Jie\\Shire_follow_up\\01_data\\"
# nonhaeDir="F:\\Jie\\Shire_follow_up\\01_data\\newdata_200K_3M\\"
################################################################################
nonhaeFile <- 'nonhae_200K_v2'
haeFile <- "HAE973_ptid"
fileNm_3M <- "neg_3M_clean2335697"

dir <- "F:\\Jie\\Shire_follow_up\\"
outDir <- paste0(dir, "03_Results\\")
haeDir <- paste0(dir, "01_data\\")
nonhaeDir <- paste0(haeDir, "newdata_200K_3M\\")
path_3M <- nonhaeDir


run_split(n.simu=n.simu, nonhaeFile=nonhaeFile, haeFile=haeFile, haeDir=haeDir, 
          nonhaeDir=nonhaeDir, split_fun=split_fun)

t0 <- proc.time()
run_bagging_rf(n.simu=n.simu, wk_dir=wk_dir, dir=outDir, lasso_rf_iters=iters, 
               nonhaeFile=nonhaeFile, haeFile=haeFile)
cat((proc.time()-t0)[3]/60, '!\n')



perf_new200_1233_iter200 <- 
  get_perf_allSimu(dir=outdir, iters=iters, n.simu=n.simu, 
                   recall_tar=recall_tar, nonhaeFile=nonhaeFile, haeFile=haeFile)


run_perf_3M(dir=outDir, wk_dir=wk_dir, lasso_rf_iters=iters, n.simu=n.simu, 
            recall_tar=recall_tar, fileNm_3M=fileNm_3M, path_3M=path_3M, 
            haeFile=haeFile, nonhaeFile=nonhaeFile)
