require("dplyr")
# read TND data, for example, ... is the Path you store the TND data, and TNDdat.tex is the name (change it)
# The data need to be cleaned: the structure is like Y is outcome, V is vaccination status, and other columns are covariates that we will use to adjust. 
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
# IPW
mod_IPW <- mean(TNDdat$Y*TNDdat$V/res$g1)/mean(TNDdat$Y*(1-TNDdat$V)/res$g0)
mod_IPW
# proposed eif estimator
A1 <- (1 - TNDdat$Y)*(TNDdat$V - res$g1)/(res$g1* res$m0)
A0 <- (1 - TNDdat$Y)*((1-TNDdat$V) - res$g0)/(res$g0* res$m0)
mod_eif <- mean(TNDdat$Y*TNDdat$V/res$g1 - res$w1* res$mu1*A1)/mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$w0* res$mu0*A0)
mod_eif

# proposed eif estimator 2
A.1 <- (1 - TNDdat$Y)*(TNDdat$V - res$g1)/(res$g1* (1 -res$mu1))
A.0 <- (1 - TNDdat$Y)*((1-TNDdat$V) - res$g0)/(res$g0* (1 - res$mu0))
mod_eif2 <- mean(TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)/mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)
mod_eif2
