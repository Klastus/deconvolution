#### Script for deconvoluting two distributions 
## function of preparing deamer object ##
# library(foreach)
# library(doParallel)
deconvolution <- function(observ, noise, grid.len=100, from, to, user.scalar=1){
  # parameters: 
    # - noise: vector of pure noise values
    # - observ: vector of noisy observations
  # - grid.len: number of points for which the distribution will be calculated
  # - from, to: from/to which value the deconvoluted distribution will be calculated
  # -user. scalar: sometimes a scaling of values can be needed. By default, there is no scaling 
  
  if(user.scalar == 0){
    scalar <- mean(observ)/log(30)  
  } else {
    scalar <- user.scalar
  }
  est <- deamerSE(y=(observ/scalar), errors = (noise/scalar), grid.length=grid.len, from=from/scalar, to=to/scalar)
  deconvoluted <- cbind(melt(as.data.frame(est$f/scalar))[2], 
                        as.data.frame(est$supp*scalar))
  colnames(deconvoluted)<- c('density', 'model_fluo')
  return(deconvoluted)
}

decon.to.values.par <- function(noise, observ, x_s, no_cores){
  
  ## function to convert distribution of deconvoluted values to 
  ## exact fluorescence values ##
  # parameters: 
  # - noise: vector of pure noise values
  # - observ: vector of noisy observations
  # - x_s: deamer object of deconvoluted values (in principle, a distribution)
  # - no_cores: number of cores that will be used in parallel computations
  
  decon.working1 <- function(xi, noise, x_s){
    a <- foreach(epsi=noise$x) %do%{
      signal <- (xi-epsi)
      Pepsi <- approx(x=noise$x, y=noise$y, xout=epsi)
      Psignal <- approx(x=x_s$model_fluo, y=x_s$dens, xout=signal)
      vector <- c(xi, signal, epsi, Pepsi$y, Psignal$y, Pepsi$y*Psignal$y)
    }
    big.data <- as.data.frame(do.call(rbind, a))
    colnames(big.data) <- c("xi", "signal", "epsi", 
                          "Pepsi", "Psignal", "Pepsi_Psignal")
    return(big.data)
  }
  eps <- density(noise)
  x_m <- density(observ)
  len.o <- length(observ)
  len.e <- length(eps$x)
  
  # calculating all probabilities of x_s
  
  registerDoParallel(no_cores)
  big.list <- foreach(xi=observ, .packages = "foreach") %dopar% {
    decon.working1(xi, eps, x_s)
  }
  stopImplicitCluster()
  # finding the highest probability of x_s
  dane <- list()
  for(i in (1:len.o)){
    no.na <- big.list[[i]][complete.cases(big.list[[i]]$Pepsi_Psignal), ]
    a <- max(no.na$Pepsi_Psignal)
    if(a==-Inf){
      print(i)
    }
    b <- big.list[[i]][which(big.list[[i]]$Pepsi_Psignal==a), ]
    dane[[i]] <- c(b$xi, b$signal)
  }
  dane.df <- as.data.frame(do.call(rbind, dane))
  colnames(dane.df) <- c("fluorescence", "t.signal")
  dane.df$roznica <- dane.df$fluorescence-dane.df$t.signal
  return(dane.df)
}

decon.treatment <- function( data.decon, col.name, treatment, exper.noise.log,
                              grid.len = 100, from = 0, to=10, no_cores = 20) {
  #function to deconvolute values for each separate group ov values from df
  # ('treatments')
  ## args:
  # - data.decon; df with convoluted values
  # - colname; column with convoluted values
  # - treatment; column which indicates to which group an observation belongs to
  # - exper.noise.log; a vector of log(values) of pure errors
  # - grid.len, from, to, no_cores; all args come from 'deconvolution' function
  
  new.col.name.log <- paste("decon.", col.name, ".log", sep='')
  new.col.name <- paste("decon.", col.name, sep='')
  data.decon[[new.col.name.log]] <- -100
  unique.treatment <- unique(data.decon[[treatment]])

  for(i in unique.treatment){
    
    exper.observ.log <- log(data.decon[data.decon[[treatment]]==i, ][[col.name]])
    deconvoluted <- deconvolution(observ=exper.observ.log, noise=exper.noise.log, 
                                  grid.len = grid.len, 
                                  from = from, to=to)
    deconvoluted.values <- decon.to.values.par(noise=exper.noise.log,
                                               observ = exper.observ.log, 
                                               x_s=deconvoluted, no_cores = no_cores)
    
    data.decon[data.decon[[treatment]]==i, ][[new.col.name.log]] <- deconvoluted.values$t.signal
  }
  data.decon[[new.col.name]] <- exp(data.decon[[new.col.name.log]])
  return(data.decon)
}






















#### testing the speed of calculations: ####
# blad <- noise[1:100]
# proba <- observ[1:100]
# start_time <- Sys.time()
# test.values <-decon.to.values.par(noise=blad,
#                                   observ = proba,
#                                   x_s=deconvoluted, no_cores = 20)
# end_time <- Sys.time()
# end_time - start_time
# # 
# ggplot()+
# geom_line(aes(x=fun.research$x, y=fun.research$y), color= 'blue')+
# geom_line(aes(x=deconvoluted$model_fluo, y=deconvoluted$density), color= 'green')+
# geom_density(aes(x=test.values$t.signal), color="orange")+
# ggtitle('deconvoluted vs cytometry data')+
# theme_jetka()+
# scale_colour_discrete(name="sample type", breaks=c('red', 'blue', 'green'),
#                       labels=c('\ncytometry\nsecondary\ncontrol\n',
#                                "cytometry\nprotein\n",
#                                "deconvoluted\n"))+
# xlab("fluorescence")+
# ylab('density')+
# xlim(0,2000)
# # 
# # 
# # ggplot(data=test.values)+
# #   geom_histogram(aes(x=fluorescence), bins=100, fill='darkred', alpha=0.5)+
# #   geom_histogram(aes(x=t.signal), bins=100, fill='blue', alpha=0.5)+
# #   theme_jetka()+
# #   xlim(0,2000)
