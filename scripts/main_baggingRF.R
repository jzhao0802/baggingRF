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
main.wk_dir <- "./"
setwd(main.wk_dir)

# 
main.n.simu = 5
main.recall_tar <- seq(0.5, 0.05, -0.05)
main.iters <- 20

################################################################################
# Lichao: 
# Previously the file names are given as follows: 
# nonhaeFile='nonhae_200K_II2.a'
# haeFile="II2.a"
# However there are no such files in the given path: 
# haeDir="F:\\Jie\\Shire_follow_up\\01_data\\"
# nonhaeDir="F:\\Jie\\Shire_follow_up\\01_data\\newdata_200K_3M\\"
################################################################################
main.nonhaeFile <- 'nonhae_200K_v2'
main.haeFile <- "HAE973_ptid"
main.fileNm_3M <- "neg_3M_clean2335697"

main.dir <- "F:\\Jie\\Shire_follow_up\\"
main.haeDir <- paste0(main.dir, "01_data\\")
main.nonhaeDir <- paste0(main.haeDir, "newdata_200K_3M\\")
main.path_3M <- main.nonhaeDir

main.timeStamp <- as.character(Sys.time())
main.timeStamp <- gsub(":", ".", main.timeStamp)  # replace ":" by "."
main.outDir <- paste(main.wk_dir, "Results/", main.timeStamp, "/", sep = '')
dir.create(main.outDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")


run_split(n.simu=main.n.simu, nonhaeFile=main.nonhaeFile, haeFile=main.haeFile, 
          haeDir=main.haeDir, nonhaeDir=main.nonhaeDir, outDir=main.outDir,
          iters=main.iters)
print("run_split finished.")


main.t0 <- proc.time()
run_bagging_rf(n.simu=main.n.simu, wk_dir=main.wk_dir, dir=main.outDir, 
               lasso_rf_iters=main.iters, nonhaeFile=main.nonhaeFile, haeFile=main.haeFile)
cat((proc.time()-main.t0)[3]/60, '!\n')
print("run_bagging_rf finished. ")

# main.outDir <- "F:/Lichao/work/Projects/HAE/code/R/GoFromBaggingForest/Results/2016-04-19 18.07.52/"

main.perf_new200_1233_iter200 <- 
  get_perf_allSimu(dir=main.outDir, iters=main.iters, n.simu=main.n.simu, 
                   recall_tar=main.recall_tar, nonhaeFile=main.nonhaeFile, 
                   haeFile=main.haeFile)
print("get_perf_allSimu finished.")


run_perf_3M(dir=main.outDir, wk_dir=main.wk_dir, lasso_rf_iters=main.iters, n.simu=main.n.simu, 
            recall_tar=main.recall_tar, fileNm_3M=main.fileNm_3M, path_3M=main.path_3M, 
            haeFile=main.haeFile, nonhaeFile=main.nonhaeFile)
print("run_perf_3M finished. ")
