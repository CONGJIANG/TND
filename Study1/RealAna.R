## Preparation
#Before running this VE analyse code,  the TND dataset should be well prepared: 
#Outcome is binary (SARS-Cov-2 infection or not), Covid19 Vaccination status is binary, e.g., unvaccinated (V=0) vs fully vaccinated + 14 days (V=1)
#Covariates that we will use to adjust is selected, such as ages group (e.g., < 18, 18~60, or >=60), gender, race, 
#and the calendar month of a test-positive subject's first positive COVID test and a test-negative subject's last COVID test, etc. 
require("dplyr") # to use select function
require("sandwich")
# read TND data, for example, ... is the Path you store the TND data, and TNDdat.tex is the name (change it!)
# The data need to be cleaned: the structure is like Y is the outcome, V is vaccination status, and other columns are covariates we will use to adjust. 
TNDdat <-read.table(".../TNDdat.txt",header=F)
head(TNDdat)
# isolate the names of baseline covariates
baselinevars <- names(dplyr::select(TNDdat, !c(V,Y)))
baselinevars
######## fit models
# adjust the exposure variable (primary interest) + covariates
# fit mu model
mu.formula <- as.formula(paste("Y~ V +", 
                               paste(baselinevars, 
                                     collapse = "+")))

mu.fit <- glm(mu.formula,family="binomial", data=TNDdat)
coef(mu.fit)
mu1 <- predict(mu.fit,newdata=as.data.frame(cbind(V=1,select(TNDdat, !c(V,Y)))),type="response")
mu0 <- predict(mu.fit,newdata=as.data.frame(cbind(V=0,select(TNDdat, !c(V,Y)))),type="response")
# fit m model
m.formula <- as.formula(paste("Y~", 
                              paste(baselinevars, 
                                    collapse = "+")))
m.fit <- glm(mu.formula,family="binomial", data=TNDdat)
coef(m.fit)
m <- predict(m.fit,type="response")
m0 <- 1 - m
# fit propensity models using control data
ps.formula <- as.formula(paste("V ~",
                               paste(baselinevars,
                                     collapse = "+")))
ps.fit <- glm(ps.formula,family="binomial", data=TNDdat, subset=(TNDdat$Y==0))
g1_cont<-predict(ps.fit,type="response",newdata= TNDdat);
g0_cont <- 1 - g1_cont

# we put these elements together as the following function:
element <- function(TNDdat){
  baselinevars <- names(dplyr::select(TNDdat, !c(V,Y)))
  ######## fit models
  # adjust the exposure variable (primary interest) + covariates
  # fit mu model
  mu.formula <- as.formula(paste("Y~ V +", 
                                 paste(baselinevars, 
                                       collapse = "+")))
  
  mu.fit <- glm(mu.formula,family="binomial", data=TNDdat)
  mu1 <- predict(mu.fit,newdata=as.data.frame(cbind(V=1,select(TNDdat, !c(V,Y)))),type="response")
  mu0 <- predict(mu.fit,newdata=as.data.frame(cbind(V=0,select(TNDdat, !c(V,Y)))),type="response")
  # fit m model
  m.formula <- as.formula(paste("Y~", 
                                paste(baselinevars, 
                                      collapse = "+")))
  m.fit <- glm(mu.formula,family="binomial", data=TNDdat)
  m0 <- 1 - predict(m.fit,type="response")
  # debiasing weights
  #w1 <-m0/(1 - mu1)
  #w0 <-m0/(1 - mu0)
  # fit propensity models using control data
  ps.formula <- as.formula(paste("V ~",
                                 paste(baselinevars,
                                       collapse = "+")))
  ps.fit <- glm(ps.formula,family="binomial", data=TNDdat, subset=(TNDdat$Y==0))
  g1_cont<-predict(ps.fit,type="response",newdata= TNDdat);
  g0_cont <- 1 - g1_cont
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 =g1_cont, g0 = g0_cont, w1 = m0/(1 - mu1), w0 = m0/(1 - mu0)))
}
######## Methods for VE
res <- element(TNDdat)
# IPW estimate and SE
mod_IPW <- mean(TNDdat$Y*TNDdat$V/res$g1)/mean(TNDdat$Y*(1-TNDdat$V)/res$g0)
mod_IPW
TNDdat$ipw <- ifelse(TNDdat$V == 1, 1/res$g1, 1/res$g0)
modY.ipw <- glm(Y ~ V, family=binomial(link = "log"), weights = ipw, data=TNDdat)
# if does not work try:
# modY.ipw <- glm(Y ~ V, family=binomial(), weights = ipw, data=TNDdat)
(est.ipw <- exp(modY.ipw$coefficients[2]))
(se.ipw <- vcovHC(modY.ipw)[2,2])

# proposed eif estimator
A1 <- (1 - TNDdat$Y)*(TNDdat$V - res$g1)/(res$g1* res$m0)
A0 <- (1 - TNDdat$Y)*((1-TNDdat$V) - res$g0)/(res$g0* res$m0)
psi1 <- mean(TNDdat$Y*TNDdat$V/res$g1 - res$w1* res$mu1*A1)
psi0 <- mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$w0* res$mu0*A0)
(mod_eif <- psi1/psi0)
var1 <- mean((TNDdat$Y*TNDdat$V/res$g1 - res$w1* res$mu1*A1)^2)
var0 <- mean((TNDdat$Y*(1-TNDdat$V)/res$g0 - res$w0* res$mu0*A0)^2)
se1 <- sqrt(var1/(psi1^2) + var0/(psi0^2))
CI_l1 <- mod_eif *exp(- 1.96*se1/sqrt(nrow(TNDdat)))
CI_u1 <- mod_eif *exp(1.96*se1/sqrt(nrow(TNDdat)))
mod_eif; c(CI_l1, CI_u1)

# proposed eif estimator 2: SAME as the previous method, i.e.,eif
A.1 <- (1 - TNDdat$Y)*(TNDdat$V - res$g1)/(res$g1* (1 -res$mu1))
A.0 <- (1 - TNDdat$Y)*((1-TNDdat$V) - res$g0)/(res$g0* (1 - res$mu0))
psi.1 <- mean(TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)
psi.0 <- mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)
(mod_eif2 <- psi.1/psi.0)
var.1 <- mean((TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)^2)
var.0 <- mean((TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)^2)
se2 <- sqrt(var.1/(psi.1^2) + var.0/(psi.0^2))
CI_l2 <- mod_eif2 *exp(- 1.96*se2/sqrt(nrow(TNDdat)))
CI_u2 <- mod_eif2 *exp( 1.96*se2/sqrt(nrow(TNDdat)))
mod_eif2; c(CI_l2, CI_u2)


# outcome regression
mod_OR <- mean(res$mu1 * res$w1)/mean(res$mu0 * res$w0)
mod_OR
nbs <- 50
bsest<-rep(NA,nbs)
for(i in 1:nbs){
  resamps<-sample(1:nrow(TNDdat),size=nrow(TNDdat),replace=T)
  datk<-TNDdat[resamps,]
  res <- element(datk)
  bsest[i] <- mean(res$mu1 * res$w1)/mean(res$mu0 * res$w0)
}
bs_var <- var(bsest)
CI = quantile(bsest,c(0.025,0.975))
mod_OR; CI
