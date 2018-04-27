normalize_data <- function(data,
                           normalize_factor = 65535){
  
  data.intensity_colnames <- grepl("Intensity", colnames(data)) & !grepl("Location", colnames(data))
  data[,data.intensity_colnames] <- data[,data.intensity_colnames]*normalize_factor
  return(list(data = data))
}

library(gridExtra)
library(reshape2)
library(ggplot2)
library(deamer)
library(foreach)
library(doParallel)
library(flowCore)
# library(doSNOW)

# package.list <- list("ggplot2", "gridExtra", "ggthemes", "grid", "reshape2", 
#                      "deamer", "foreach", "doParallel", "flowcore")

