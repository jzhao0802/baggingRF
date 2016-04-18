# # load several packages
# library(ggplot2)
# library(ROCR)
# library(PRROC)
# # library(FSelector)
# library(randomForest)
# library(caret)
# library(e1071)
# library(reshape2)
# library(sqldf)
# library(glmnet)
# library(caTools)
# library(gbm)
# # library(xlsx)
# library(dplyr)
# library(snowfall)
# library(plyr)

install_load <- function (package1, ...)  {   
    # convert arguments to vector
    packages <- c(package1, ...)
    # start loop to determine if each package is installed
    for(package in packages){
        # if package is installed locally, load
        if(package %in% rownames(installed.packages()))
            do.call('library', list(package))
        # if package is not installed locally, download, then load
        else {
            install.packages(package)
            do.call("library", list(package))
        }
    }
}


install_load("ggplot2", 'ROCR', 'PRROC', 'randomForest', 'caret', 'e1071', 'reshape2', 'sqldf', 'glmnet', 'caTools'
             , "gbm", "dplyr", "snowfall", 'plyr')

