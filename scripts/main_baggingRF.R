# main function for running bagging forest and performance
rm(list=ls())
# R work directory
wk_dir = "F:\\Jie\\Shire_follow_up\\02_Code\\HAE_R_codes_Dec15\\"

# Setup R work directory


setwd(wk_dir)

source("loadpackage.R")
# Load R packages
install_load("ggplot2", 'ROCR', 'PRROC', 'randomForest', 'caret', 'e1071', 'reshape2', 'sqldf', 'glmnet', 'caTools'
             , "gbm", "dplyr", "snowfall", 'plyr')
# Auxiliary functions
source("auxfunctions.R")
# 
source("funs_baggingRF.R")


n.simu=5
split_fun=3
recall_tar <- seq(0.5, 0.05, -0.05)
nonhaeFile='nonhae_200K_II2.a'
haeFile="II2.a"
dir="F:\\Jie\\Shire_follow_up\\03_Results\\"
haeDir="F:\\Jie\\Shire_follow_up\\01_data\\"
nonhaeDir="F:\\Jie\\Shire_follow_up\\01_data\\newdata_200K_3M\\"
iters=20
fileNm_3M="neg_3M_clean2335697"
path_3M="F:\\Jie\\Shire_follow_up\\01_data\\newdata_200K_3M\\"    


run_split(n.simu=n.simu, nonhaeFile=nonhaeFile, haeFile=haeFile
          , haeDir=haeDir
          , nonhaeDir=nonhaeDir
          , split_fun=split_fun
          )

t0 <- proc.time()
run_bagging_rf(n.simu=n.simu
               , wk_dir=wk_dir
               , dir=dir
               , lasso_rf_iters=iters
               , nonhaeFile=nonhaeFile
               , haeFile=haeFile
)
cat((proc.time()-t0)[3]/60, '!\n')



perf_new200_1233_iter200 <- get_perf_allSimu(dir=dir
                                             , iters=iters
                                             , n.simu=n.simu
                                             , recall_tar=recall_tar
                                             , nonhaeFile=nonhaeFile
                                             , haeFile=haeFile
)


run_perf_3M(dir=dir
            , wk_dir=wk_dir
            , lasso_rf_iters=iters
            , n.simu=n.simu
            , recall_tar=recall_tar
            , fileNm_3M=fileNm_3M
            , path_3M=path_3M
            , haeFile=haeFile
            , nonhaeFile=nonhaeFile
            )












