#Reading in the data

mastitis <- read.csv("~/Documents/SECOND/bayesian/mastitis.dat")
rabbit <- read.delim("~/Documents/SECOND/bayesian/rabbit.txt")

#Required Library
library(MCMCpack) ; library(mcmc) ; library(mcmcplots) ; library(R2jags)
library(rjags) ; library(coda) ;  library(ggmcmc);  library(runjags) ; library(snow)
library(superdiag) ; library(tidyr) ; library(tidyverse) ; library(parfm) 

mastitis$logMidpoint = NULL ; mastitis$logMidpoint = log(mastitis$Midpoint)


#classical nonlinear model for question 1ajkgy
mm = parfm(Surv(Midpoint, Status) ~  Heifer + Quarter , cluster = "Cowid" , data=mastitis, dist = "lognormal",frailty = "lognormal" )
mm
ci.parfm(mm, level=0.05)["Heifer",] ; ci.parfm(mm, level=0.05)["Quarter",]


#classical nonlinear model for question 1b
Lens = rabbit$Lens
Age = rabbit$Age
p = nls(Lens2 ~ (theta1 * exp((-1*theta2) / (theta3+Age2))), start=list(theta1=270, theta2= 120, theta3=30))
summary(p)

q = nls(Lens2 ~  exp(theta1 + exp((-1*theta2) / (theta3+Age2))), start=list(theta1=199, theta2= 0.19, theta3=79 ), control = list(maxiter = 500), trace = TRUE)
nls(rabbit$loglens ~  theta1 + exp((-1*theta2) / (theta3+Age2)), start=list(theta1=4.254, theta2= 76.58, theta3=0.358 ), control = list(maxiter = 500), trace = TRUE)


require(ssym)

data("Erabbits", package="ssym")
plot(Erabbits$age, Erabbits$wlens, type="p", cex=0.3, lwd=3,
     ylab="Dry weight of eye lens (in milligrams)",
     xlab="Age of the animal (in days)")

#Now the Bayesian approach
#Question 1
set.seed(1234)
cat(
  "model{
    for (i in 1:400) {
    epsalum[i] = (logMidpoint[i] - beta0 - beta1*heifer[i] - beta2*Q1[i] - beta3*Q2[i] - beta4*Q3[i] - b[cowid[i]])/sigma.e
    f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
    s_o[i] = 1 - phi(epsalum[i])
    ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])

#the poisson zero tricks
   constant[i] = 10000 - ll[i]
   samp[i] ~ dpois(constant[i])
    }

#random effects 
for (j in 1:100) {
  b[j] ~ dnorm(0,tau.b)
}

#priors  (fixed effect)
beta0 ~ dnorm(0, 1.0E-3)
beta1 ~ dnorm(0, 1.0E-3)
beta2 ~ dnorm(0, 1.0E-3)
beta3 ~ dnorm(0, 1.0E-3)
beta4 ~ dnorm(0, 1.0E-3)

#priors  (within subject variance)
sigma.e ~ dunif(0,100)

#priors  (between subject variance)
tau.b = pow(sigma.b, -2)
sigma.b ~ dunif(0,100)
}
  ",file="mastitis.jags")

#the data to supply
mdata = list(
  logMidpoint = mastitis$logMidpoint, 
  heifer = mastitis$Heifer,
  Q1 = as.numeric(mastitis$Quarter==1),
  Q2 = as.numeric(mastitis$Quarter==2),
  Q3 = as.numeric(mastitis$Quarter==3),
  Midpoint = mastitis$Midpoint,
  status = mastitis$Status,
  samp = rep(0,400),
  cowid = mastitis$Cowid
  
)

#the initial values
inits = list(
  list(beta0 = 4.9,
       beta1 = -0.1,
       beta2 = 0.07,
       beta3 = 0.15,
       beta4 = -0.05,
       sigma.e = 0.5,
       sigma.b = 0.82), 
  list(beta0 = 4.0,
       beta1 = -0.24,
       beta2 = 0.09,
       beta3 = 0.25,
       beta4 = -0.055,
       sigma.e = 0.7,
       sigma.b = 0.90),
  list(beta0 = 3.9,
       beta1 = -0.21,
       beta2 = 0.02,
       beta3 = 0.35,
       beta4 = -0.09,
       sigma.e = 0.59,
       sigma.b = 0.42)
)

myjags <- jags(data=mdata, inits = inits, parameters.to.save = c("beta0","beta1","beta2","beta3","beta4","sigma.e","sigma.b"),n.iter = 50000, n.burnin = 25000,n.chains = 3,model.file = "mastitis.jags")


myjags_mcmc = as.mcmc(myjags)
myjags_mcmcggs <- ggs(myjags_mcmc)
ggs_traceplot(myjags_mcmcggs)
ggs_autocorrelation(myjags_mcmcggs)
ggs_running(myjags_mcmcggs)
gelman.plot(myjags_mcmc,autoburnin = F)
gelman.diag(myjags2_mcmc,autoburnin = F)

library(xtable)
effectiveSize(myjags_mcmc)
effectiveSize(myjags2_mcmc)

gelman.diag(myjags2_mcmc)

#hierachical centering

cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] - beta1*heifer[i] - beta2*Q1[i] - beta3*Q2[i] - beta4*Q3[i] - b[cowid[i]])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #random effects 
  for (j in 1:100) {
  b[j] ~ dnorm(beta0,tau.b)
  }
  
  #priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-3)
  beta1 ~ dnorm(0, 1.0E-3)
  beta2 ~ dnorm(0, 1.0E-3)
  beta3 ~ dnorm(0, 1.0E-3)
  beta4 ~ dnorm(0, 1.0E-3)
  
  #priors  (within subject variance)
  sigma.e ~ dunif(0,100)
  
  #priors  (between subject variance)
  tau.b = pow(sigma.b, -2)
  sigma.b ~ dunif(0,100)
  }
  ",file="mastitis_hcentering.jags")

mdata2 = list(
  logMidpoint = mastitis$logMidpoint, 
  heifer = mastitis$Heifer,
  Q1 = as.numeric(mastitis$Quarter==1),
  Q2 = as.numeric(mastitis$Quarter==2),
  Q3 = as.numeric(mastitis$Quarter==3),
  Midpoint = mastitis$Midpoint,
  status = mastitis$Status,
  samp = rep(0,400),
  cowid = mastitis$Cowid
)

myjags2 <- jags(data=mdata2, inits = inits, 
                parameters.to.save = c("beta0","beta1","beta2","beta3","beta4",
                                       "sigma.e","sigma.b"),n.iter = 50000, 
                n.burnin = 25000,n.chains = 3,n.thin=1,
                model.file = "mastitis_hcentering.jags")

# converting to mcmc so we can use functions in coda
# traceplot(myjags2)
# 
myjags2_mcmc = as.mcmc(myjags2)
myjags2_mcmcggs <- ggs(myjags2_mcmc)
ggs_traceplot(myjags2_mcmcggs)
ggs_autocorrelation(myjags2_mcmcggs)
ggs_running(myjags2_mcmcggs)
gelman.plot(myjags2_mcmc)


summary(myjags2_mcmc)
acfplot(myjags2_mcmc)

### saving only the random effect

myjags2_b1 <- jags(data=mdata, inits = inits, 
                  parameters.to.save = c("b"),n.iter = 50000, 
                  n.burnin = 25000,n.chains = 3,n.thin=1,
                  model.file = "mastitis_hcentering.jags")


myjags2_b2 <- jags(data=mdata2, inits = inits, 
                parameters.to.save = c("b"),n.iter = 50000, 
                n.burnin = 25000,n.chains = 3,n.thin=1,
                model.file = "mastitis.jags")

#checking the histogram/normality of the random effects
myjags2_b2_mcmc <- as.mcmc(myjags2_b2)
myjags2_b2_mcmc_summary <- summary(myjags2_b2_mcmc)
head(myjags2_b_mcmc_summary$statistics)
hist(myjags2_b2_mcmc_summary$statistics[1:100,1])

#for convergence
#traceplot(myjags2_b)


### using t random effect
cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] - beta1*heifer[i] - beta2*Q1[i] - beta3*Q2[i] - beta4*Q3[i] - b[cowid[i]])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #random effects 
  for (j in 1:100) {
  b[j] ~ dt(beta0,tau.b,nu.b)
  }
  
  #priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-3)
  beta1 ~ dnorm(0, 1.0E-3)
  beta2 ~ dnorm(0, 1.0E-3)
  beta3 ~ dnorm(0, 1.0E-3)
  beta4 ~ dnorm(0, 1.0E-3)
  
  #priors  (within subject variance)
  sigma.e ~ dunif(0,100)
  
  #priors  (between subject variance)
  tau.b = pow(sigma.b, -2)
  sigma.b ~ dunif(0,100)

  #prior for the df
  nu.b <- 1/nu.inverse.b
  nu.inverse.b ~ dunif(0,1)
  }
  ",file="mastitis_hcenteringt.jags")

mdata2t = list(
  logMidpoint = mastitis$logMidpoint, 
  heifer = mastitis$Heifer,
  Q1 = as.numeric(mastitis$Quarter==1),
  Q2 = as.numeric(mastitis$Quarter==2),
  Q3 = as.numeric(mastitis$Quarter==3),
  Midpoint = mastitis$Midpoint,
  status = mastitis$Status,
  samp = rep(0,400),
  cowid = mastitis$Cowid
)

myjags2t <- jags(data=mdata2t, inits = inits, 
                parameters.to.save = c("beta0","beta1","beta2","beta3","beta4",
                                       "sigma.e","sigma.b","b"),n.iter = 50000, 
                n.burnin = 25000,n.chains = 3,n.thin=1,
                model.file = "mastitis_hcenteringt.jags")

#For each monitored node,  effectiveSize()  returns the ???effective sample size??? ??? 
#an estimate of the number of independent samples from the posterior that would 
#hold the same information.


#tppc
cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] - beta1*heifer[i] - beta2*Q1[i] - beta3*Q2[i] - beta4*Q3[i] - b[cowid[i]])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #random effects 
  for (j in 1:100) {
  b[j] ~ dt(beta0,tau.b,nu.b)
b.new[j] ~ dnorm(beta0,tau.b)
}

#maximum
b_max <- max(b)
b.new_max <- max(b.new)
test1 <- step(b.new_max - b_max)
#minimum
b_min <- min(b)
b.new_min <- min(b.new)
test2 <- step(b.new_min - b_min)
#skewness
for (j in 1:100) {
  b_skew[j] <- pow((b[j] - mean(b[]))/sigma.b,3)
  b_skew_new[j] <- pow((b.new[j] - mean(b.new[]))/sigma.b,3)
  #kurtosis
  b_kurt[j] <- pow((b[j] - mean(b[]))/sigma.b,4) - 3
  b_kurt_new[j] <- pow((b.new[j] - mean(b.new[]))/sigma.b,4) - 3
}
skewb <- mean(b_skew[])
skewb_new <- mean(b_skew_new[])
test3 <- step(skewb_new - skewb)
kurtb <- mean(b_kurt[])
kurtb_new <- mean(b_kurt_new[])
test4 <- step(kurtb_new - kurtb)
#ksmethod
sortb <- sort(b[])
sortb_new <- sort(b.new[])
for (j in 1:100) {
  pb[j] <- phi(sortb[j])
  pbnew[j] <- phi(sortb_new[j])
  maxpb[j] <- max(pb[j] - ((j-1)/100) , (j/100) - pb[j])
  maxpbnew[j] <- max(pbnew[j] - ((j-1)/100) , (j/100) - pbnew[j])
}
ks_b <- max(maxpb) 
ks_bnew <- max(maxpbnew)
test5 <- step(ks_bnew - ks_b)
#sinharay and stern
absb <- abs(b_max - sortb[50]) - abs(b_min - sortb[50])
absb_new <- abs(b.new_max - sortb_new[50]) - abs(b.new_min - sortb_new[50])
test6 <- step(absb_new - absb)

## alltest
ppc_test[1] <- test1 
ppc_test[2] <- test2
ppc_test[3] <- test3
ppc_test[4] <- test4
ppc_test[5] <- test5
ppc_test[6] <- test6

#priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-3)
beta1 ~ dnorm(0, 1.0E-3)
beta2 ~ dnorm(0, 1.0E-3)
beta3 ~ dnorm(0, 1.0E-3)
beta4 ~ dnorm(0, 1.0E-3)

#priors  (within subject variance)
sigma.e ~ dunif(0,100)

#priors  (between subject variance)
tau.b = pow(sigma.b, -2)
sigma.b ~ dunif(0,100)
}
",file="mastitis_hpcc.jags")

myjags2ppc <- jags(data=mdata2, inits = inits, 
                   parameters.to.save = c("beta0","beta1","beta2", "beta3","beta4","b","sigma.e","sigma.b","ppc_test"),n.iter = 50000, 
                   n.burnin = 25000,n.chains = 3,n.thin=1,
                   model.file = "mastitis_hpcc.jags")







myjags2tmcmc <- as.mcmc(myjags2t)
myjags2tsummary <- summary(myjags2tmcmc)

par(mfrow=c(1,2))
hist(myjags2_b2_mcmc_summary$statistics[1:100,1],xlab="b",main="Histogram of Posterior Means \nof Normal Random Intercepts",probability=T,col="cyan4")
qqnorm(myjags2_b2_mcmc_summary$statistics[1:100,1],col="cyan4")
qqline(myjags2_b2_mcmc_summary$statistics[1:100,1])

hist(myjags2tsummary$statistics[1:100,1],xlab="b",main="Histogram of Posterior Means \nof t Random Intercepts",probability=T,col="cyan4")



####Posterior Predictive Check

cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] - beta1*heifer[i] - beta2*Q1[i] - beta3*Q2[i] - beta4*Q3[i] - b[cowid[i]])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #random effects 
  for (j in 1:100) {
  b[j] ~ dnorm(beta0,tau.b)
  b.new[j] ~ dnorm(beta0,tau.b)
  }
  
  #maximum
  b_max <- max(b)
  b.new_max <- max(b.new)
  test1 <- step(b.new_max - b_max)
    #minimum
    b_min <- min(b)
    b.new_min <- min(b.new)
    test2 <- step(b.new_min - b_min)
        #skewness
        for (j in 1:100) {
        b_skew[j] <- pow((b[j] - mean(b[]))/sigma.b,3)
        b_skew_new[j] <- pow((b.new[j] - mean(b.new[]))/sigma.b,3)
        #kurtosis
        b_kurt[j] <- pow((b[j] - mean(b[]))/sigma.b,4) - 3
        b_kurt_new[j] <- pow((b.new[j] - mean(b.new[]))/sigma.b,4) - 3
        }
       skewb <- mean(b_skew[])
       skewb_new <- mean(b_skew_new[])
       test3 <- step(skewb_new - skewb)
       kurtb <- mean(b_kurt[])
       kurtb_new <- mean(b_kurt_new[])
       test4 <- step(kurtb_new - kurtb)
            #ksmethod
            sortb <- sort(b[])
            sortb_new <- sort(b.new[])
            for (j in 1:100) {
            pb[j] <- phi(sortb[j])
            pbnew[j] <- phi(sortb_new[j])
            maxpb[j] <- max(pb[j] - ((j-1)/100) , (j/100) - pb[j])
            maxpbnew[j] <- max(pbnew[j] - ((j-1)/100) , (j/100) - pbnew[j])
            }
            ks_b <- max(maxpb) 
            ks_bnew <- max(maxpbnew)
            test5 <- step(ks_bnew - ks_b)
                #sinharay and stern
                absb <- abs(b_max - sortb[50]) - abs(b_min - sortb[50])
                absb_new <- abs(b.new_max - sortb_new[50]) - abs(b.new_min - sortb_new[50])
                test6 <- step(absb_new - absb)

## alltest
ppc_test[1] <- test1 
ppc_test[2] <- test2
ppc_test[3] <- test3
ppc_test[4] <- test4
ppc_test[5] <- test5
ppc_test[6] <- test6


  #priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-3)
  beta1 ~ dnorm(0, 1.0E-3)
  beta2 ~ dnorm(0, 1.0E-3)
  beta3 ~ dnorm(0, 1.0E-3)
  beta4 ~ dnorm(0, 1.0E-3)
  
  #priors  (within subject variance)
  sigma.e ~ dunif(0,100)
  
  #priors  (between subject variance)
  tau.b = pow(sigma.b, -2)
  sigma.b ~ dunif(0,100)
  }
  ",file="mastitis_hpcc.jags")

myjags2ppc <- jags(data=mdata2, inits = inits, 
                parameters.to.save = c("beta0","beta1","beta2", "beta3","beta4","b","sigma.e","sigma.b","ppc_test"),n.iter = 50000, 
                n.burnin = 25000,n.chains = 3,n.thin=1,
                model.file = "mastitis_hpcc.jags")

# myjags2ppc$BUGSoutput
asmcmc <- as.mcmc(myjags2ppc); summary(asmcmc)
# 
#ppc.check(myjags2ppc ,b_min ,b.new_min)

#####varying prior specification.
cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] - beta1*heifer[i] - beta2*Q1[i] - beta3*Q2[i] - beta4*Q3[i] - b[cowid[i]])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #random effects 
  for (j in 1:100) {
  b[j] ~ dnorm(beta0,tau.b)
  }
  
  #priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-9)
  beta1 ~ dnorm(0, 1.0E-9)
  beta2 ~ dnorm(0, 1.0E-9)
  beta3 ~ dnorm(0, 1.0E-9)
  beta4 ~ dnorm(0, 1.0E-9)
  
  #priors  (within subject variance)
  sigma.e ~ dunif(0,100)
  
  #priors  (between subject variance)
  tau.b = pow(sigma.b, -2)
  sigma.b ~ dunif(0,100)
  }
  ",file="mastitis_hcprior2.jags")


myjags2prior <- jags(data=mdata2, inits = inits, 
                parameters.to.save = c("beta0","beta1","beta2","beta3","beta4","b","sigma.e","sigma.b"),n.iter = 50000, 
                n.burnin = 25000,n.chains = 3,n.thin=1,
                model.file = "mastitis_hcprior2.jags")

#####varying sigmas prior specification
cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] - beta1*heifer[i] - beta2*Q1[i] - beta3*Q2[i] - beta4*Q3[i] - b[cowid[i]])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #random effects 
  for (j in 1:100) {
  b[j] ~ dnorm(beta0,tau.b)
  }
  
  #priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-3)
  beta1 ~ dnorm(0, 1.0E-3)
  beta2 ~ dnorm(0, 1.0E-3)
  beta3 ~ dnorm(0, 1.0E-3)
  beta4 ~ dnorm(0, 1.0E-3)
  
  #priors  (within subject variance)
  sigma.e ~ dunif(0,1.0E6)
  
  #priors  (between subject variance)
  tau.b = pow(sigma.b, -2)
  sigma.b ~ dunif(0,1.E6)
  }
  ",file="mastitis_hcprior3.jags")

myjags3prior <- jags(data=mdata2, inits = inits, 
                     parameters.to.save = c("beta0","beta1","beta2","beta3","beta4","b","sigma.e","sigma.b"),n.iter = 50000, 
                     n.burnin = 25000,n.chains = 3,n.thin=1,
                     model.file = "mastitis_hcprior3.jags")

# converting to mcmc so we can use functions in coda
traceplot(myjags3prior)
myjags3prior_mcmc = as.mcmc(myjags3prior)
summary(myjags3prior_mcmc)
acfplot(myjags3prior_mcmc)


##Removing unimportant covariate
cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] - beta2*Q1[i] - beta3*Q2[i] - beta4*Q3[i] - b[cowid[i]])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #random effects 
  for (j in 1:100) {
  b[j] ~ dnorm(beta0,tau.b)
  }
  
  #priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-3)
  beta2 ~ dnorm(0, 1.0E-3)
  beta3 ~ dnorm(0, 1.0E-3)
  beta4 ~ dnorm(0, 1.0E-3)
  
  #priors  (within subject variance)
  sigma.e ~ dunif(0,100)
  
  #priors  (between subject variance)
  tau.b = pow(sigma.b, -2)
  sigma.b ~ dunif(0,100)
  }
  ",file="mastitis_hcenteringhf.jags")

mdata2hf = list(
  logMidpoint = mastitis$logMidpoint, 
  Q1 = as.numeric(mastitis$Quarter==1),
  Q2 = as.numeric(mastitis$Quarter==2),
  Q3 = as.numeric(mastitis$Quarter==3),
  Midpoint = mastitis$Midpoint,
  status = mastitis$Status,
  samp = rep(0,400),
  cowid = mastitis$Cowid
)

#the initial values
initshf = list(
  list(beta0 = 4.9,
       beta2 = 0.07,
       beta3 = 0.15,
       beta4 = -0.05,
       sigma.e = 0.5,
       sigma.b = 0.82), 
  list(beta0 = 4.0,
       beta2 = 0.09,
       beta3 = 0.25,
       beta4 = -0.055,
       sigma.e = 0.7,
       sigma.b = 0.90),
  list(beta0 = 3.9,
       beta2 = 0.02,
       beta3 = 0.35,
       beta4 = -0.09,
       sigma.e = 0.59,
       sigma.b = 0.42)
)

myjags2hf <- jags(data=mdata2hf, inits = initshf, 
                parameters.to.save = c("beta0","beta2","beta3","beta4",
                                       "sigma.e","sigma.b", "b"),n.iter = 50000, 
                n.burnin = 25000,n.chains = 3,n.thin=1,
                model.file = "mastitis_hcenteringhf.jags")

# converting to mcmc so we can use functions in coda
traceplot(myjags2hf)
myjags3prior_mcmchf = as.mcmc(myjags2hf)
summary(myjags3prior_mcmchf)
acfplot(myjags3prior_mcmchf)

# without quater

cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] - beta1*heifer[i]- b[cowid[i]])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #random effects 
  for (j in 1:100) {
  b[j] ~ dnorm(beta0,tau.b)
  }
  
  #priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-3)
  beta1 ~ dnorm(0, 1.0E-3)
  
  #priors  (within subject variance)
  sigma.e ~ dunif(0,100)
  
  #priors  (between subject variance)
  tau.b = pow(sigma.b, -2)
  sigma.b ~ dunif(0,100)
  }
  ",file="mastitis_hcenteringqt.jags")

mdata2qt = list(
  logMidpoint = mastitis$logMidpoint, 
  heifer = mastitis$Heifer,
  Midpoint = mastitis$Midpoint,
  status = mastitis$Status,
  samp = rep(0,400),
  cowid = mastitis$Cowid
)

#the initial values
initsqt = list(
  list(beta0 = 4.9,
       beta1 = -0.1,
       sigma.e = 0.5,
       sigma.b = 0.82), 
  list(beta0 = 4.0,
       beta1 = -0.24,
       sigma.e = 0.7,
       sigma.b = 0.90),
  list(beta0 = 3.9,
       beta1 = -0.21,
       sigma.e = 0.59,
       sigma.b = 0.42)
)

myjags2qt <- jags(data=mdata2qt, inits = initsqt, 
                  parameters.to.save = c("beta0","beta1",
                                         "sigma.e","sigma.b","b"),n.iter = 50000, 
                  n.burnin = 25000,n.chains = 3,n.thin=1,
                  model.file = "mastitis_hcenteringqt.jags")

## No covariates
cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] - b[cowid[i]])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #random effects 
  for (j in 1:100) {
  b[j] ~ dnorm(beta0,tau.b)
  }
  
  #priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-3)
  
  #priors  (within subject variance)
  sigma.e ~ dunif(0,100)
  
  #priors  (between subject variance)
  tau.b = pow(sigma.b, -2)
  sigma.b ~ dunif(0,100)
  }
  ",file="mastitis_hcenteringno.jags")

mdata2no = list(
  logMidpoint = mastitis$logMidpoint, 
  Midpoint = mastitis$Midpoint,
  status = mastitis$Status,
  samp = rep(0,400),
  cowid = mastitis$Cowid
)

#the initial values
initsno = list(
  list(beta0 = 4.9,
       sigma.e = 0.5,
       sigma.b = 0.82), 
  list(beta0 = 4.0,
       sigma.e = 0.7,
       sigma.b = 0.90),
  list(beta0 = 3.9,
       sigma.e = 0.59,
       sigma.b = 0.42)
)

myjags2no <- jags(data=mdata2no, inits = initsno, 
                  parameters.to.save = c("beta0",
                                         "sigma.e","sigma.b","b"),n.iter = 50000, 
                  n.burnin = 25000,n.chains = 3,n.thin=1,
                  model.file = "mastitis_hcenteringno.jags")

#Dzero = D + 2nC ; D = Dzero - 2nc ; c= 10000 ; n= 100
(DICqt1 =myjags2qt$BUGSoutput$DIC - (2*400*10000))
(DIChf1 =myjags2hf$BUGSoutput$DIC - (2*400*10000))
(DICno1 =myjags2no$BUGSoutput$DIC - (2*400*10000))
(DICal1 =myjags2$BUGSoutput$DIC - (2*400*10000))

#Without random effect (Independence)

cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] -beta0 - beta1*heifer[i] - beta2*Q1[i] - beta3*Q2[i] - beta4*Q3[i])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-3)
  beta1 ~ dnorm(0, 1.0E-3)
  beta2 ~ dnorm(0, 1.0E-3)
  beta3 ~ dnorm(0, 1.0E-3)
  beta4 ~ dnorm(0, 1.0E-3)
  
  #priors  (within subject variance)
  sigma.e ~ dunif(0,100)
  
  }
  ",file="mastitis_hcenteringnorf.jags")

mdata2norf = list(
  logMidpoint = mastitis$logMidpoint, 
  heifer = mastitis$Heifer,
  Q1 = as.numeric(mastitis$Quarter==1),
  Q2 = as.numeric(mastitis$Quarter==2),
  Q3 = as.numeric(mastitis$Quarter==3),
  Midpoint = mastitis$Midpoint,
  status = mastitis$Status,
  samp = rep(0,400)
)

initsnorf = list(list(beta0=4.2,
           beta1=-0.43,
           beta2=0.05,
           beta3=0.10,
           beta4=-0.01,
           sigma.e=0.45
  ),
  list(beta0=3.1,
       beta1=-0.13,
       beta2=0.06,
       beta3=0.23,
       beta4=-0.15,
       sigma.e=0.31
  ),
  list(beta0=7.5,
       beta1=-0.45,
       beta2=0.10,
       beta3=0.05,
       beta4=-0.012,
       sigma.e=0.76)
)

myjags2norf <- jags(data=mdata2norf, inits = initsnorf, 
                parameters.to.save = c("beta0","beta1","beta2","beta3","beta4",
                                       "sigma.e"),n.iter = 50000, 
                n.burnin = 25000,n.chains = 3,n.thin=1,
                model.file = "mastitis_hcenteringnorf.jags")

(DICnorf1 =myjags2norf$BUGSoutput$DIC - (2*400*10000))
traceplot(myjags2norf)
acfplot(myjags2norf)

#meian survival time for each subgroup

#hierachical centering

cat(
  "model{
  for (i in 1:400) {
  
  epsalum[i] =  (logMidpoint[i] - beta1*heifer[i] - beta2*Q1[i] - beta3*Q2[i] - beta4*Q3[i] - b[cowid[i]])/sigma.e
  #
  f_o[i] = pow(2*3.142,-0.5) * exp(-0.5 * pow(epsalum[i],2))
  s_o[i] = 1 - phi(epsalum[i])
  ll[i] = status[i] * log(f_o[i]/(Midpoint[i]*sigma.e)) + (1-status[i]) * log(s_o[i])
  
  #the poisson zero tricks
  constant[i] = 10000 - ll[i]
  samp[i] ~ dpois(constant[i])
  }
  
  #random effects 
  for (j in 1:100) {
  b[j] ~ dnorm(beta0,tau.b)
  }
  
# h=0 q=1
med1 = exp(beta0 + beta2)
# h=0 q=2
med2 = exp(beta0 + beta3)
# h=0 q=3
med3 = exp(beta0 + beta4)
#h=0 q=4
med4 = exp(beta0)
#h=1 q=1
med5 = exp(beta0 + beta1 + beta2)
#h=1 q=2
med6 = exp(beta0 + beta1 + beta3)
#h=1 q=3
med7 = exp(beta0 + beta1 + beta4)
#h=1 q=4
med8 = exp(beta0 + beta1)
 
  #priors  (fixed effect)
  beta0 ~ dnorm(0, 1.0E-3)
  beta1 ~ dnorm(0, 1.0E-3)
  beta2 ~ dnorm(0, 1.0E-3)
  beta3 ~ dnorm(0, 1.0E-3)
  beta4 ~ dnorm(0, 1.0E-3)
  
  #priors  (within subject variance)
  sigma.e ~ dunif(0,100)
  
  #priors  (between subject variance)
  tau.b = pow(sigma.b, -2)
  sigma.b ~ dunif(0,100)
  }
  ",file="mastitis_hcentering.jags")

mdata2 = list(
  logMidpoint = mastitis$logMidpoint, 
  heifer = mastitis$Heifer,
  Q1 = as.numeric(mastitis$Quarter==1),
  Q2 = as.numeric(mastitis$Quarter==2),
  Q3 = as.numeric(mastitis$Quarter==3),
  Midpoint = mastitis$Midpoint,
  status = mastitis$Status,
  samp = rep(0,400),
  cowid = mastitis$Cowid
)

myjags2surv <- jags(data=mdata2, inits = inits, 
                parameters.to.save = c("med1", "med2", "med3", "med4", "med5","med6","med7","med8"),n.iter = 50000, 
                n.burnin = 25000,n.chains = 3,n.thin=1,
                model.file = "mastitis_hcentering.jags")n

summary(as.mcmc(myjags2surv))

myjags2surv_mcmc = as.mcmc(myjags2surv)
myjags2surv_mcmcggs <- ggs(myjags2surv_mcmc)
ggs_traceplot(myjags2surv_mcmcggs)
ggs_autocorrelation(myjags2_mcmcggs)
ggs_running(myjags2_mcmcggs)
gelman.plot(myjags2_mcmc)


#Question 2 

cat(
  "model{
  for (i in 1:71) {
  Lens2[i] ~ dnorm(mui[i] , tau.ep)
  mui[i] = theta1 * exp((-1*theta2) / (theta3 + Age2[i]))
  esp[i] = Lens2[i] - mui[i]
  }
  
  #priors  (fixed effect)
  theta1 ~ dnorm(0, 1.0E-3)
  theta2 ~ dnorm(0, 1.0E-3)
  theta3 ~ dnorm(0, 1.0E-3)
  
  #priors  for variance
  tau.ep = pow(sigma.ep, -2)
  sigma.ep ~ dunif(0,100)
  }
  ",file="rabbit.jags")

mdata3 = list(
  Age2 = rabbit$Age, 
  Lens2 = rabbit$Lens
)

#the initial values
inits2 = list(
  list(theta1 = 270,
       theta2 = 120,
       theta3 = 30,
       sigma.ep = 60), 
  list(theta1 = 260,
       theta2 = 150,
       theta3 = 50,
       sigma.ep = 82),
  list(theta1 = 250,
       theta2 = 100,
       theta3 = 60,
       sigma.ep = 70)
)

myjags3 <- jags(data=mdata3, inits = inits2, 
                parameters.to.save = c("theta1","theta2","theta3","sigma.ep"),n.iter = 100000, 
                n.burnin = 50000,n.chains = 3,n.thin=10,
                model.file = "rabbit.jags")

myjags3_mcmc = as.mcmc(myjags3) ; acfplot(myjags_mcmc)
traceplot(myjags3)


rabbit$loglens = log(rabbit$Lens)

cat(
  "model{
  for (i in 1:71) {
loglens[i] ~ dnorm(mui[i] , tau.ep)
  mui[i] = theta1 + exp((-1*theta2) / (theta3 + Age2[i]))
  esp[i] = loglens[i] - mui[i]
  }
  
  #priors  (fixed effect)
  theta1 ~ dnorm(0, 1.0E-3)
  theta2 ~ dnorm(0, 1.0E-3)
  theta3 ~ dnorm(0, 1.0E-3)
  
  #priors  for variance
  tau.ep = pow(sigma.ep, -2)
  sigma.ep ~ dunif(0,100)
  }
  ",file="rabbit2.jags")

mdata4 = list(
  Age2 = rabbit$Age, 
  loglens = rabbit$loglens
)

#the initial values
inits3 = list(
  list(theta1 = 0.4,
       theta2 = 120,
       theta3 = 30,
       sigma.ep = 60), 
  list(theta1 = 260,
       theta2 = 150,
       theta3 = 50,
       sigma.ep = 82),
  list(theta1 = 250,
       theta2 = 100,
       theta3 = 60,
       sigma.ep = 70)
)

myjags4 <- jags(data=mdata4, inits = inits3, 
                parameters.to.save = c("theta1","theta2","theta3","sigma.ep", "esp"),n.iter = 50000, 
                n.burnin = 25000,n.chains = 3,n.thin=1,
                model.file = "rabbit2.jags")
myjags4$BUGSoutput

myjags3_mcmc = as.mcmc(myjags3) ; acfplot(myjags_mcmc)
traceplot(myjags3)
