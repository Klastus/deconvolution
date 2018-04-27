set.seed(23456)
n = 1000; M = 500
x <- rchisq(n, 100)
b=0.5
e <- rlaplace(n, 30, b)
y <- x + e
# y <- y[y>0]
eps <- rlaplace(M, 30, b)

est <- deamerSE(y, eps, from=0)
est
log.y <- log(y)
log.eps <- log(eps)
est.log <- deamerSE(log.y, log.eps, from = 0)

est.log <- deconvolution(log.y, log.eps, from=0, to=6)
values.log <- decon.to.values.par(noise = log.eps, observ = log.y, x_s = est.log, no_cores = 20)
values.log$log.correction <- exp(values.log$fluorescence) - (exp(values.log$fluorescence)/exp(values.log$t.signal))


est <- deconvolution(y, eps, from=0, to=180)
values <- decon.to.values.par(noise = eps, observ = y, x_s = est, no_cores = 20)

alfa = 0.5
ggplot()+
  geom_histogram(aes(log(y)), bins=50, fill="green", alpha=alfa)+
  geom_histogram(aes(log(eps)), bins=50, fill="brown", alpha=alfa)+
  geom_histogram(aes(values.log$t.signal), bins=50, fill="orange", alpha=alfa)+
  geom_histogram(aes(log(x)), bins=50, fill="red", alpha=alfa)
  # geom_histogram(aes(log(values.log$log.correction)), bins=50, fill="yellow")

ggplot()+
  geom_histogram(aes(y), bins=50, fill="green", alpha=alfa)+
  geom_histogram(aes(eps), bins=50, fill="brown", alpha=alfa)+
  # geom_histogram(aes(exp(values.log$t.signal)), bins=50, fill="orange", alpha=alfa)+
  geom_histogram(aes(x), bins=50, fill="red", alpha=alfa)+
  geom_histogram(aes(values.log$log.correction), bins=50, fill="yellow", alpha=alfa)+
  geom_histogram(aes(values$t.signal), bins=50, fill="blue", alpha=alfa)+
  xlim(0,180)

merged <- data.frame(fluo=values$fluorescence, 
                     decon.correction=values.log$log.correction, 
                     decon.normal=values$t.signal)

ggplot(merged)+
  geom_point(aes(x=decon.normal, y=decon.correction))











set.seed(23456)
n = 1000; M = 500
x <- rchisq(n, 100)
b=3
e <- rlaplace(n, 30, b)
y <- x + e
# y <- y[y>0]
eps <- rlaplace(M, 30, b)

est <- deamerSE(y, eps, from=0)
est
log.y <- log(y)
log.eps <- log(eps)
est.log <- deamerSE(log.y, log.eps, from = 0)

est.log <- deconvolution(log.y, log.eps, from=0, to=6)
values.log <- decon.to.values.par(noise = log.eps, observ = log.y, x_s = est.log, no_cores = 20)
values.log$log.correction <- exp(values.log$fluorescence) - (exp(values.log$fluorescence)/exp(values.log$t.signal))


est <- deconvolution(y, eps, from=0, to=180)
values <- decon.to.values.par(noise = eps, observ = y, x_s = est, no_cores = 20)

ggplot()+
  geom_histogram(aes(log(y)), bins=50, fill="green", alpha=alfa)+
  geom_histogram(aes(log(eps)), bins=50, fill="brown", alpha=alfa)+
  geom_histogram(aes(values.log$t.signal), bins=50, fill="orange", alpha=alfa)+
  geom_histogram(aes(log(x)), bins=50, fill="red", alpha=alfa)
# geom_histogram(aes(log(values.log$log.correction)), bins=50, fill="yellow")

ggplot()+
  geom_histogram(aes(y), bins=50, fill="green", alpha=alfa)+
  geom_histogram(aes(eps), bins=50, fill="brown", alpha=alfa)+
  # geom_histogram(aes(exp(values.log$t.signal)), bins=50, fill="orange", alpha=alfa)+
  geom_histogram(aes(x), bins=50, fill="red", alpha=alfa)+
  geom_histogram(aes(values.log$log.correction), bins=50, fill="yellow", alpha=alfa)+
  geom_histogram(aes(values$t.signal), bins=50, fill="blue", alpha=alfa)+
  xlim(0,180)

merged <- data.frame(fluo=values$fluorescence, 
                     decon.correction=values.log$log.correction, 
                     decon.normal=values$t.signal)

ggplot(merged)+
  geom_point(aes(x=decon.normal, y=decon.correction))



















set.seed(23456)
n = 1000; M = 500
x <- rnorm(n, 150, sd = 50)
e <- rnorm(n, 50, sd = 20)
y <- x + e
y[y<0]

eps <- rnorm(n=M, mean=50, sd = 15)

log.y <- log(y)
log.eps <- log(eps)

est.log <- deconvolution(log.y, log.eps, from=0, to=9)
values.log <- decon.to.values.par(noise = log.eps, observ = log.y, x_s = est.log, no_cores = 20)
values.log$log.correction <- exp(values.log$fluorescence) - (exp(values.log$fluorescence)/exp(values.log$t.signal))


est <- deconvolution(y, eps, from=0, to=370)
values <- decon.to.values.par(noise = eps, observ = y, x_s = est, no_cores = 20)

ggplot()+
  geom_histogram(aes(log(y)), bins=50, fill="green", alpha=alfa)+
  geom_histogram(aes(log(eps)), bins=50, fill="brown", alpha=alfa)+
  geom_histogram(aes(values.log$t.signal), bins=50, fill="orange", alpha=alfa)+
  geom_histogram(aes(log(x)), bins=50, fill="red", alpha=alfa)
# geom_histogram(aes(log(values.log$log.correction)), bins=50, fill="yellow")

ggplot()+
  geom_histogram(aes(y), bins=50, fill="green", alpha=alfa)+
  geom_histogram(aes(eps), bins=50, fill="brown", alpha=alfa)+
  # geom_histogram(aes(exp(values.log$t.signal)), bins=50, fill="orange", alpha=alfa)+
  geom_histogram(aes(x), bins=50, fill="red", alpha=alfa)+
  geom_histogram(aes(values.log$log.correction), bins=50, fill="yellow", alpha=alfa)+
  geom_histogram(aes(values$t.signal), bins=50, fill="blue", alpha=alfa)+
  xlim(0,400)

merged <- data.frame(fluo=values$fluorescence, 
                     decon.correction=values.log$log.correction, 
                     decon.normal=values$t.signal,
                     x=x)

ggplot(merged)+
  geom_point(aes(x=x, y=decon.correction))
