# Biostats project: icddr,b Shigella analysis
# 3/18/23
# Mac LaPrete

require(lme4)
require(latex2exp)
library(scales)

setwd("/Users/JayLaPrete/Desktop/BIOL 6500 Biostats/project/")
source("dhaka_data.R")
source("matlab_data.R")

################################################################################
## Recreate paper figures (pdf) ################################################
################################################################################
PLOT <- FALSE                                                                   # TRUE if want to draw plots
strains_dhaka <- get_strains_dhaka()
strains_matlab <- get_strains_matlab()
if(PLOT){
  plot_strains_dhaka(strains_dhaka)
  plot_strains_matlab(strains_matlab)
}

resistance_dhaka <- get_resistance_dhaka()
resistance_matlab <- get_resistance_matlab()
# plot both locations together
plot_resistance_both <- function(resistance_dhaka,resistance_matlab){
  titles <- c("Ciprofloxacin All","Ciprofloxacin S. flexneri","Ciprofloxacin S. sonnei",
              "Azithromycin All","Azithromycin S. flexneri","Azithromycin S. sonnei",
              "Mecillinam All","Mecillinam S. flexneri","Mecillinam S. sonnei",
              "Ceftriaxone All","Ceftriaxone S. flexneri","Ceftriaxone S. sonnei",
              "Multidrug All")
  pdf("resistance_both.pdf", width = 12, height = 20)
  par(mfrow=c(5,3))
  for(i in 2:dim(resistance_dhaka)[2]){
    plot(resistance_dhaka$times, resistance_dhaka[,i],
         col="darkgreen", type = "l", ylim = c(0,1),
         ylab = "Portion of Resistance", xlab = "Year",
         main = titles[i-1])
    lines(resistance_matlab$times, resistance_matlab[,i],
         col="red", type = "l")
    legend("topleft",c("Dhaka","Matlab"),col = c("darkgreen","red"), lty = c(1,1))
  }
  dev.off()
}
if(PLOT){
  plot_resistance_both(resistance_dhaka,resistance_matlab)
}

# combine to one data set
rd <- resistance_dhaka
rd$loc <- rep("d",length(resistance_dhaka$times))
rm <- resistance_matlab
rm$loc <- rep("m",length(resistance_matlab$times))
resistance_both <- rbind(rd,rm)
resistance_both <- within(resistance_both, {loc <- factor(loc)})

# estimated count data set
nsamples <- 1000
res_count <- cbind(subset(resistance_both, select=-c(times,loc))*nsamples,      # gets count data from percents (R)
                   loc=resistance_both[,"loc"],times=resistance_both[,"times"])
sus_count <- cbind(nsamples - subset(res_count, select=-c(times,loc)),          # gets count of susceptible from R and n
                   loc=resistance_both[,"loc"],times=resistance_both[,"times"])
all_count <- data.frame()
for(i in 1:dim(res_count)[1]){                                                  # turns counts into 0s and 1s
  all_count <- rbind(all_count,
                     cbind(rep(res_count$times[i],res_count$cipro_all[i]),
                           rep(1,res_count$cipro_all[i])),
                     cbind(rep(sus_count$times[i],sus_count$cipro_all[i]),
                           rep(0,sus_count$cipro_all[i])))
}
colnames(all_count) <- c("t","R")

# the big ugly plot of all the 1s and 0s jittered for visibility
plot_all_dots <- function(){
  counts_plot <- data.frame(rep(0,dim(all_count)[2],dim(all_count)[1]))
  counts_plot$t <- jitter(all_count$t, 2.5)
  counts_plot$R <- jitter(all_count$R, 0.5)
  plot(counts_plot$t, counts_plot$R,
       xlab = "Years", ylab = TeX("Proportion resistant ($P_L$)"),
       col = alpha("black",0.5), cex = 0.5)
}


################################################################################
## Estimating r_R and r_S ######################################################
################################################################################
ts <- resistance_both$times                                                     # store times to plot
uts <- c(unique(ts),2020,2025)                                                  # times for plotting

########################
### glm proportions ####
########################
get_P_0 <- function(beta0) 1/(1+exp(-beta0))                                    # get parameters from betas
get_rdiffs <- function(beta1) beta1

Ps <- resistance_both$cipro_all                                                 # store P values to plot

PtoP_L <- function(P) log(P/(1-P))                                              # transform data to log form
plot(ts, PtoP_L(Ps),
     xlab = "Years", ylab = TeX("log Proportion resistant ($P_L$)"))            # plot transformed data to check linearity
m_lognorm <- glm(PtoP_L(cipro_all)~times, data=resistance_both, family=gaussian)# run glm with normal
summary(m_lognorm)

glm_b0 <- m_lognorm$coefficients[1]
glm_b1 <- m_lognorm$coefficients[2]
lines(uts,glm_b0+glm_b1*uts, col="blue")                                        # plot linear regression on linearized data

P_LtoP <- function(t,b0,b1) 1/(1+exp(-(b0+b1*t)))                               # reverse data transform
plot(ts, Ps,
     xlab = "Years", ylab = "Proportion resistant (P)",
     ylim = c(0,1), xlim = c(min(uts),max(uts)))                                # plot the data
lines(uts,P_LtoP(uts,glm_b0,glm_b1), col="blue")                                # plot fit line

psR2_lognorm <- with(summary(m_lognorm), 1 - deviance/null.deviance)            # R^2? pseudo R^2?

get_P_0(glm_b0)                                                                 # P0 parameter estimate
get_rdiffs(glm_b1)                                                              # r_R - r_S parameter estimate

# R^2
1 - (m_lognorm$deviance/m_lognorm$null.deviance)
# adjusted R^2
1 - ((m_lognorm$deviance/m_lognorm$df.residual)/(m_lognorm$null.deviance/m_lognorm$df.null))


############
### MLE ####
############
q <- function(b0,b1,ts) 1/(1+exp(-(b0+b1*ts)))                                  # sample mean of simulation
Cs <- res_count$cipro_all                                                       # store P values to plot

negloglik <- function(params){                                                  # negative log likelihood function
  b0 <- params[1]
  b1 <- params[2]
  -sum(dbinom(Cs,
              size = nsamples,
              prob = q(b0,b1,ts),
              log=TRUE))
}

initvals <- c(glm_b0,glm_b1)
optimalparams <- optim(par=initvals, negloglik)
mle_b0 <- optimalparams$par[1]
mle_b1 <- optimalparams$par[2]

plot(ts, Ps,
     xlab = "Years", ylab = "Proportion resistant (P)",
     ylim = c(0,1), xlim=c(2000,2025))                                          # plot the data
lines(uts,q(mle_b0,mle_b1,uts),col="darkgreen")                                 # plot fit line
# stop and save fig 1

lines(uts,P_LtoP(uts,glm_b0,glm_b1), col="blue")                                # add glm line
legend("topleft",c("glm","MLE"), col = c("blue","darkgreen"), lty = c(1,1))

mle_b0
mle_b1
get_P_0(mle_b0)                                                                 # P0 parameter estimate
get_rdiffs(mle_b1)                                                              # r_R - r_S parameter estimate

nsamples
negloglik(c(mle_b0,mle_b1))

# AIC
p <- -sum(dbinom(Cs,size = nsamples,prob = q(mle_b0,mle_b1,ts),log=TRUE))
2*p + 4

p<-1
for(i in 1:8){
  p <- p*choose(1000,Cs[i])*q(mle_b0,mle_b1,ts[i])^Cs[i]*(1-q(mle_b0,mle_b1,ts[i]))^(1000-Cs[i])
}
-2*log(p) + 4

# Cs
# p<-1
# for(i in 1:8000){
#   p <- p*q(mle_b0,mle_b1,all_count$t[i])^(all_count$R[i])*(1-q(mle_b0,mle_b1,all_count$t[i]))^(1-(all_count$R[i]))
# }
# -2*log(p) + 4

# R^2
deviance <- 0
for(i in 1:8000){
  deviance <- deviance + (all_count$R[i] - q(mle_b0,mle_b1,all_count$t[i]))^2
}
null.deviance <- 0
the_mean <- sum(q(mle_b0,mle_b1,all_count$t))/8000
for(i in 1:8000){
  null.deviance <- null.deviance + (all_count$R[i] - the_mean)^2
}
df.r <- 7999
df.n <- 7998
1 - (deviance/null.deviance)
# adjusted R^2
1 - ((deviance/df.r)/(null.deviance/df.n))


##########################
### glm estimate data ####
##########################
m_glm_est <- glm(R~t, data=all_count, family=binomial)                          # run glm with normal
summary(m_glm_est)

glm_est_b0 <- m_glm_est$coefficients[1]
glm_est_b1 <- m_glm_est$coefficients[2]
Ps_est <- 1/(1+exp(-(glm_est_b0 + glm_est_b1*uts)))

plot_all_dots()                                                                 # plot transformed data to check linearity
lines(uts,Ps_est, col="orange", lwd=2)                                          # plot fit line

points(ts, Ps, pch=19, cex=1, col="red")                                        # plot the data
lines(uts,q(mle_b0,mle_b1,uts),col="darkgreen", lwd=2, lty=2)                   # plot MLE fit line
lines(uts,P_LtoP(uts,glm_b0,glm_b1), col="blue", lwd=2)                         # add glm line
legend("topleft",c("glm","MLE","glm est"),
       col = c("blue","darkgreen","orange"), lty = c(1,2,2), lwd = c(1,1,1))

get_P_0(glm_est_b0)                                                             # P0 parameter estimate
get_rdiffs(glm_est_b1)                                                          # r_R - r_S parameter estimate

# R^2
1 - (m_glm_est$deviance/m_glm_est$null.deviance)
# adjusted R^2
1 - ((m_glm_est$deviance/m_glm_est$df.residual)/(m_glm_est$null.deviance/m_glm_est$df.null))


#############
### MCMC ####
#############
sapply(c('nimble','dplyr','ggplot2','MCMCvis'),require,character.only = TRUE)

model <- nimbleCode({                                                           # specify likelihood
  for(i in 1:ndat){
    predicted.y[i] <- 1/(1+exp(-(intercept + slope*x[i])))                      # equation for regression line
    y[i] ~ dbinom(predicted.y[i], nsamp)                                        # likelihood for each data point (each row)
  }
  
  intercept ~ dnorm(0, sd = 1000) # dbeta(1,1)
  slope ~ dbeta(1,1)
  nsamp ~ dbinom(0.5,1000) # dbeta(1,1)
})

#p_data <- list(y=Ps,x=ts)                                                       # proportion data to list form
p_data <- list(y=all_count$R,x=all_count$t)                                                       # proportion data to list form
ndat <- length(p_data$y)

parameters.to.save <- c("intercept", "slope", "nsamp")

initial.values <- function(){                                                   # random initial values for each parameter
  list(intercept = runif(1,-1000,0),
       slope = runif(1,0,1),
       nsamp = runif(1,500,1500))
}
initial.values()

n.iter <- 7000                                                                  # MCMC tech specs
n.burnin <- 2000
n.chains <- 3

p_mcmc.output <- nimbleMCMC(code = model,                                       # run MCMC on p_data
                            data = p_data,
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            nburnin = n.burnin,
                            nchains = n.chains)
str(p_mcmc.output)                                                              # displays structure of mcmc output
dim(p_mcmc.output$chain1)                                                       # 5000x3 (3 parameters we are keeping track of)
head(p_mcmc.output$chain1)                                                      # look at chain 1
mean(p_mcmc.output$chain1[,"slope"])                                            # expected value of posterior for slope in chain 1
mean(p_mcmc.output$chain2[,"slope"])
mean(p_mcmc.output$chain3[,"slope"])

par(mfrow=c(1,2))
hist(p_mcmc.output$chain1[,'slope'],
     xlab = "Slope", ylab = "Frequency", main = "Posterior")                    # histogram of slope
hist(p_mcmc.output$chain1[,'intercept'],
     xlab = "Intercept", ylab = "Frequency", main = "Posterior")                # histogram of intercept
par(mfrow=c(1,1))

MCMCsummary(object = p_mcmc.output, round = 2)                                  # mean of posterior

library(HDInterval)
hdi(p_mcmc.output$chain1[,'slope'], credMass = 0.95)                            # credible intervals
hdi(p_mcmc.output$chain1[,'intercept'], credMass = 0.95)

MCMCtrace(object = p_mcmc.output,
          pdf = FALSE,
          ind = TRUE,                                                           # separate density lines per chain
          params = "slope")
MCMCtrace(object = p_mcmc.output,
          pdf = FALSE,
          ind = TRUE,                                                           # separate density lines per chain
          params = "intercept")

slp <- MCMCsummary(object = p_mcmc.output, round = 2)[3,1]
intcpt <- MCMCsummary(object = p_mcmc.output, round = 2)[1,1]
pred_line <- 1/(1+exp(-(intcpt + slp*uts)))
plot(ts, Ps,
     xlab = "Years", ylab = "Proportion resistant (P)",
     ylim = c(0,1), xlim=c(2000,2025))                                          # plot the data
lines(uts, pred_line, col="purple",lwd=2)                                       # plot MCMC line
lines(uts,q(mle_b0,mle_b1,uts),col="darkgreen", lwd=2, lty=1)                   # plot MLE fit line
lines(uts,P_LtoP(uts,glm_b0,glm_b1), col="blue", lwd=2)                         # add glm line
lines(uts,Ps_est, col="orange", lwd=2, lty=2)                                   # plot fit line
legend("topleft",c("glm","MLE","glm est","MCMC"),
       col = c("blue","darkgreen","orange","purple"), lty = c(1,2,2,1))







