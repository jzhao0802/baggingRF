
# Apr15

# main function for running bagging forest and performance
rm(list=ls())
# R work directory
setwd("F:\\Jie\\Shire_follow_up\\02_Code\\GoFromBaggingForest\\")
wk_dir = "./"

# Setup R work directory
setwd(wk_dir)
# Load R packages
source("scripts/loadpackage.R")
# Auxiliary functions
source("functions/auxfunctions.R")
# 
source("functions/funs_baggingRF.R")

timeStamp <- as.character(Sys.time())
timeStamp <- gsub(":", ".", timeStamp)  # replace ":" by "."
resultDir <- paste("./Results/", timeStamp, "/", sep = '')
dir.create(resultDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")


resultDirDataSpecific <- paste0(resultDir, "for_new200K&973\\")

run_split(n.simu=5, haeFile=973
          , nonhaeFile='for_new200K'
          , dataDir="F:\\Jie\\Shire_follow_up\\01_data\\"
          , outDir=resultDir
          , split_fun=3
)


run_bagging_rf(n.simu = 5, 
               wk_dir = wk_dir, 
               outDir =resultDirDataSpecific, 
               lasso_rf_iters = 20)



perf_new300_iter20 <- 
    get_perf_allSimu(outdir=resultDirDataSpecific
                     , iters=20
                     , n.simu=5
                     , recall_tar=seq(0.5, 0.05, -0.05)
    )

t0 <- proc.time()
perf_new200K_973_on3M_20 <- run_perf_3M(outDir=resultDirDataSpecific
                                         , wk_dir=wk_dir
                                         , lasso_rf_iters=20
                                         , n.simu=5
                                         , recall_tar=seq(0.5, 0.05, -0.05)
                                         , fileNm_3M='neg_3M_clean2335697'
                                         , path_3M='F:\\Jie\\Shire_follow_up\\01_data\\newdata_200K_3M\\'
)
cat((proc.time()-t0)[3]/60, 'min!\n')
write.csv(perf_new200K_973_on3M_20,
          paste0(resultDirDataSpecific, "iters=20\\perf_on3M.csv"))

