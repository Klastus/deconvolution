#### pSTAT1 distribution for deconvolution ####
## preparing the distribution of pSTAT1 accumulation at time 0 in nuclei,
## which can be next used for deconvolution of other distribution 
## (e.g. in WT cells)

## microscopic data only ##

#### data preparation- loading and normalization ####
## Important! first load normalization_data() function.
# path - the path to the PT36 shrinked nuclei (or other) .csv

pS1.decon.distrib <- function(path){

  PT36.nuc <- normalize_data(read.table(path, 
                                        header=TRUE, sep=","))$data
  
  #### choosing apropriate columns and renaming them accordingly ####
  ## columns of main interests: ##
  columns <- c("Intensity_IntegratedIntensity_Alexa488",
               "Intensity_MeanIntensity_Alexa488",
               "Intensity_IntegratedIntensity_Alexa555",
               "Intensity_MeanIntensity_Alexa555")
  
  variable.subset <- function(data, columns, new.columns){
    data.2 <- data[, colnames(data) %in% columns]
    colnames(data.2) <- new.columns
    return(data.2)
  }
  
  PT36.nuc_2 <- variable.subset(PT36.nuc, 
                                c(columns, "Intensity_IntegratedIntensity_DAPI"), 
                                c("Integrated_Alexa488_nucleus",
                                  "Integrated_Alexa555_nucleus",
                                  "Integrated_DAPI_nucleus",
                                  "Mean_Alexa488_nucleus",
                                  "Mean_Alexa555_nucleus"))
  
  PT36 <- cbind(PT36.nuc_2, 
                PT36.nuc[, 87:96])
  
  #### subsetting only values for pSTAT1 antibodies, not STAT1 #### 
  PT36.pS1 <- PT36[! PT36$antibody.1.1==0, ]
  
  #### plotting ####

    # ggplot(PT36.pS1)+
    #   geom_boxplot(aes(x=factor(time.1.1), y= Mean_Alexa488_nucleus))+
    #   theme_jetka()+
    #   ylim(0,60)
    # 
    # ggplot(PT36.pS1)+
    #   geom_histogram(aes(x=Mean_Alexa488_nucleus), bins=100)+
    #   theme_jetka()+
    #   xlim(0,60)+
    #   facet_grid(. ~ time.1.1)

  # conclusion from above: distributions of the pS1 in nucleus is the same
  # for all times
  
  return(PT36.pS1$Mean_Alexa488_nucleus)
}

# experiment <- "PT36"
# paths <- paste("Y:/PiotrT/deconvolution/",
#                experiment, sep='')
# normaliz <- "/input/"
# 
# test2 <- pS1.decon.distrib(path=paste(paths, normaliz, "ShrinkedNuclei.csv", sep=''))
