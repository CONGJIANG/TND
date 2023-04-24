#Common data-generating function

#popsize set at 1500000, can increase in large ssize is needed
#ssize set at 500, may have trouble when not enough patient available

#Defaults (scenario 1) true marg RR = 
#OR_C<-3 #effect of vaccine on covid
#OR_WI<-1 #no effect of vaccine on W from other infection
#OR_WC<-5 #effect of vaccine on covid symptoms
#OR_H<-1.5 #effect of vaccine on hospitalization among those with symptoms
datagen<-function(seed=sample(1:1000000,size=1),ssize=500,popsize=1500000,OR_C=3,OR_WI=1,OR_WC=5,OR_H=1.5,em=0,cfV0=F,cfV1=F,return_full=F){
  set.seed(seed)
  
  #generate data
  C<-runif(n=popsize, -3, 3)
  U1<-rbinom(n=popsize,size=1,prob=0.5) #affects both
  U2<-rbinom(n=popsize,size=1,prob=0.5) #affects covid
  
  if(cfV0==T) V=rep(0,popsize); if(cfV1==T) V=rep(1,popsize);
  if(cfV0==F&cfV1==F){
    V<-rbinom(prob=plogis(0.5 + 0.3*C + abs(C) - sin(pi*C) ),size=1,n=popsize) #prevalence is around 0.61
  }
  #Infection (with something) has some common risk factors U1 and C
  Infec<-rbinom(prob=plogis(0.5*C-5+0.5*U1),size=1,n=popsize) #current prevalence around 0.007
  
  #Infected with COVID
  Infec_COVID<- rbinom(prob=plogis( -log(OR_C)*V -5 + 2*C - 0.15*exp(C)+log(3)*U2*(1.5-V)-2*U1), size=1,n=popsize) #0.009
  #symptoms based on infection
  #can come from either or both infections, if present
  W=rep(0,popsize)
  W[Infec==1]<-rbinom(prob=plogis(2+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
  W[Infec_COVID==1]<-rbinom(prob=plogis(-3+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
  #mean(W[Infec==1|Infec_COVID==1]) #25%
  #mean(W[Infec_COVID==1]) #39%
  #mean(W[Infec==1]) #12%
  #mean(W[Infec_COVID==1&V==1]) #22%
  
  #hospitalization, only possible if symptoms present
  H=rep(0,popsize)
  H[W==1]<-rbinom(prob=plogis(1.5+0.5*C[W==1]+log(OR_H)*V[W==1]-0.5*U1[W==1]),size=1,n=sum(W==1))
  #mean(H[W==1]) #83% with severe symptoms go to hospital
  
  #selection on outcome for testing (does not condition on infectious status, just being in the hospital)
  R<-sample(which(H==1),ssize) #sample randomly from those in hospital to control the study size
  
  if(return_full==F){
    dat<-as.data.frame(cbind(Y=Infec_COVID,V=V,C=C)[R,])
  } else{dat<-as.data.frame(cbind(Infec_COVID=Infec_COVID,Infec=Infec,H=H,W=W,V=V,C=C))}
  return(dat)
}


dat<-datagen(ssize=1000, em=0)
sum(dat$Y)/1000

######################################################################
###Run Simulation Study
# IPW estimator
library(sandwich)
# IPW correct
mod_IPW_c <- function(dat){
  TNDmod_g_col<-glm(V ~  C + abs(C) +  sin(pi*C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  ipw <- ifelse(dat$V == 1, 1/g1_cont, 1/(1 - g1_cont))
  modY.ipw <- glm(Y ~ V, family=binomial(link = "logit"), weights = ipw, data=dat)
  est <- exp(modY.ipw$coefficients[2])
  se <- sqrt(vcovHC(modY.ipw)[2,2])
  return(list(est = est, CI = c(est - 1.96*se/sqrt(nrow(dat)), est + 1.96*se/sqrt(nrow(dat))) ))
}
# IPW wrong
mod_IPW_w <- function(dat){
  TNDmod_g_col<-glm(V ~ 1 ,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  ipw <- ifelse(dat$V == 1, 1/g1_cont, 1/(1 - g1_cont))
  modY.ipw <- glm(Y ~ V, family=binomial(link = "logit"), weights = ipw, data=dat)
  est <- exp(modY.ipw$coefficients[2])
  se <- sqrt(vcovHC(modY.ipw)[2,2])
  return(list(est = est, CI = c(est - 1.96*se/sqrt(nrow(dat)), est + 1.96*se/sqrt(nrow(dat))) ))
}

# outcome regression correct
mod_OR_c <- function(dat){
  TNDmod<-glm(Y~ V + C + exp(C),family=binomial(),data=dat)
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  est <- mean(Q1)/mean(Q2)
  return(list(est = est))
}
#  outcome regression wrong
mod_OR_w <- function(dat){
  TNDmod<-glm(Y~ 1,family=binomial(),data=dat)
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  est <- mean(Q1)/mean(Q2)
  return(list(est = est))
}
######################################################################
# EIF estimator 1 (Equation 10):
# S1) both are correct
modEIF1a<-function(dat){
  TNDmod_g_col<-glm(V~  C + abs(C) +  sin(pi*C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  TNDmod1<-glm(Y ~ C + exp(C),family=binomial(),data=dat)
  preY=predict(TNDmod1,type="response")
  
  TNDmod<-glm(Y~ V + C + exp(C),family=binomial(),data=dat)
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - preY))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - preY))
  psi1 <- mean(dat$Y*dat$V/g1_cont - Q1*A1)
  psi0 <- mean(dat$Y*(1-dat$V)/g0_cont - Q0*A0)
  var1 <- mean((dat$Y*dat$V/g1_cont - Q1*A1)^2)
  var0 <- mean((dat$Y*(1-dat$V)/g0_cont - Q0*A0)^2)
  se <- sqrt(var1/(psi1^2) + var0/(psi0^2))
  CI_l <- (psi1/psi0)* exp(- 1.96*se/sqrt(nrow(dat)))
  CI_u <- (psi1/psi0)* exp(1.96*se/sqrt(nrow(dat)))
  return(list(est = psi1/psi0, CI = c(CI_l, CI_u) ) )
}
# S2) Out is correct, but PS is wrong
modEIF1b<-function(dat){
  TNDmod_g_col<-glm( V ~ 1 ,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod1<-glm(Y ~ C + exp(C),family=binomial(),data=dat)
  preY=predict(TNDmod1,type="response")
  
  TNDmod<-(glm(Y~ V + C + exp(C),family=binomial(),data=dat)) 
  
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - preY))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - preY))
  
  psi1 <- mean(dat$Y*dat$V/g1_cont - Q1*A1)
  psi0 <- mean(dat$Y*(1-dat$V)/g0_cont - Q0*A0)
  var1 <- mean((dat$Y*dat$V/g1_cont - Q1*A1)^2)
  var0 <- mean((dat$Y*(1-dat$V)/g0_cont - Q0*A0)^2)
  se <- sqrt(var1/(psi1^2) + var0/(psi0^2))
  CI_l <- (psi1/psi0)* exp(- 1.96*se/sqrt(nrow(dat)))
  CI_u <- (psi1/psi0)* exp(1.96*se/sqrt(nrow(dat)))
  return(list(est = psi1/psi0, CI = c( CI_l, CI_u) ))
}
# S3) PS is correct, but Out is wrong
modEIF1c<-function(dat){
  TNDmod_g_col<-glm(V~  C + abs(C) +  sin(pi*C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod1<-glm(Y ~ 1,family=binomial(),data=dat)
  preY=predict(TNDmod1,type="response")
  
  TNDmod<-(glm(Y~ V,family=binomial(),data=dat)) 
  
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - preY))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - preY))
  
  psi1 <- mean(dat$Y*dat$V/g1_cont - Q1*A1)
  psi0 <- mean(dat$Y*(1-dat$V)/g0_cont - Q0*A0)
  var1 <- mean((dat$Y*dat$V/g1_cont - Q1*A1)^2)
  var0 <- mean((dat$Y*(1-dat$V)/g0_cont - Q0*A0)^2)
  se <- sqrt(var1/(psi1^2) + var0/(psi0^2))
  CI_l <- (psi1/psi0)* exp(- 1.96*se/sqrt(nrow(dat)))
  CI_u <- (psi1/psi0)* exp(1.96*se/sqrt(nrow(dat)))
  return(list(est = psi1/psi0, CI = c( CI_l, CI_u)))
}
# S4) Both are wrong
modEIF1d<-function(dat){
  TNDmod_g_col<-glm( V ~ 1 ,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod1<-glm(Y ~ 1,family=binomial(),data=dat)
  preY=predict(TNDmod1,type="response")
  
  TNDmod<-(glm(Y~ V,family=binomial(),data=dat)) 
  
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - preY))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - preY))
  
  psi1 <- mean(dat$Y*dat$V/g1_cont - Q1*A1)
  psi0 <- mean(dat$Y*(1-dat$V)/g0_cont - Q0*A0)
  var1 <- mean((dat$Y*dat$V/g1_cont - Q1*A1)^2)
  var0 <- mean((dat$Y*(1-dat$V)/g0_cont - Q0*A0)^2)
  se <- sqrt(var1/(psi1^2) + var0/(psi0^2))
  CI_l <- (psi1/psi0)* exp(- 1.96*se/sqrt(nrow(dat)))
  CI_u <- (psi1/psi0)* exp(1.96*se/sqrt(nrow(dat)))
  return(list(est = psi1/psi0, CI = c( CI_l, CI_u)) )
}


######################################################################
### EIF 2 (Equation 11)
# EIF estimator 2:S1) both are correct
modEIF2a<-function(dat){
  TNDmod_g_col<-glm(V~  C + abs(C) +  sin(pi*C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V + C + exp(C),family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  psi1 <- mean(dat$Y*dat$V/g1_cont - mu1*A1)
  psi0 <- mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  var1 <- mean((dat$Y*dat$V/g1_cont - mu1*A1)^2)
  var0 <- mean((dat$Y*(1-dat$V)/g0_cont - mu0*A0)^2)
  se <- sqrt(var1/(psi1^2) + var0/(psi0^2))
  CI_l <- (psi1/psi0)* exp(- 1.96*se/sqrt(nrow(dat)))
  CI_u <- (psi1/psi0)* exp(1.96*se/sqrt(nrow(dat)))
  return(list(est = psi1/psi0, CI = c( CI_l, CI_u)))
}

modEIF2b<-function(dat){
  TNDmod_g_col<-glm( V ~ 1 ,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V + C + exp(C),family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  psi1 <- mean(dat$Y*dat$V/g1_cont - mu1*A1)
  psi0 <- mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  var1 <- mean((dat$Y*dat$V/g1_cont - mu1*A1)^2)
  var0 <- mean((dat$Y*(1-dat$V)/g0_cont - mu0*A0)^2)
  se <- sqrt(var1/(psi1^2) + var0/(psi0^2))
  CI_l <- (psi1/psi0)* exp(- 1.96*se/sqrt(nrow(dat)))
  CI_u <- (psi1/psi0)* exp(1.96*se/sqrt(nrow(dat)))
  return(list(est = psi1/psi0, CI = c( CI_l, CI_u)))
}

modEIF2c<-function(dat){
  TNDmod_g_col<-glm(V~  C + abs(C) +  sin(pi*C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V,family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  psi1 <- mean(dat$Y*dat$V/g1_cont - mu1*A1)
  psi0 <- mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  var1 <- mean((dat$Y*dat$V/g1_cont - mu1*A1)^2)
  var0 <- mean((dat$Y*(1-dat$V)/g0_cont - mu0*A0)^2)
  se <- sqrt(var1/(psi1^2) + var0/(psi0^2))
  CI_l <- (psi1/psi0)* exp(- 1.96*se/sqrt(nrow(dat)))
  CI_u <- (psi1/psi0)* exp(1.96*se/sqrt(nrow(dat)))
  return(list(est = psi1/psi0, CI = c( CI_l, CI_u)))
}

modEIF2d<-function(dat){
  TNDmod_g_col<-glm( V ~ 1 ,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V,family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  psi1 <- mean(dat$Y*dat$V/g1_cont - mu1*A1)
  psi0 <- mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  var1 <- mean((dat$Y*dat$V/g1_cont - mu1*A1)^2)
  var0 <- mean((dat$Y*(1-dat$V)/g0_cont - mu0*A0)^2)
  se <- sqrt(var1/(psi1^2) + var0/(psi0^2))
  CI_l <- (psi1/psi0)* exp(- 1.96*se/sqrt(nrow(dat)))
  CI_u <- (psi1/psi0)* exp(1.96*se/sqrt(nrow(dat)))
  return(list(est = psi1/psi0, CI = c( CI_l, CI_u)))
}

CI1=CI2=CI3=CI4=CI5=CI6=CI7=CI8=CI9=CI10=c(0,0)

for (i in 1:1000){
  dat<-datagen(ssize=1000, em=0)
  #######################################################
  # Marginal RR 
  # IPW ps correct
  est1<-mod_IPW_c(dat)$est
  CI1 <- mod_IPW_c(dat)$CI
  # IPW ps wrong
  est2<-mod_IPW_w(dat)$est
  CI2 <- mod_IPW_w(dat)$CI
  #######################################################
  # EIF estimator 1
  #Case 1: All correct
  est3 <- modEIF1a(dat)$est
  CI3 <-  modEIF1a(dat)$CI
  #######################################################
  #Case 2: Outcome model is right, but PS is not
  est4 <- modEIF1b(dat)$est
  CI4 <-  modEIF1b(dat)$CI
  #######################################################
  # Case 3: (PS is correct, but the outcome model is not correct)
  est5 <- modEIF1c(dat)$est
  CI5 <-  modEIF1c(dat)$CI
  #######################################################
  # Both are not correct
  est6 <- modEIF1d(dat)$est
  CI6 <-  modEIF1d(dat)$CI
  #######################################################
  # EIF estimator 2
  #Case 1: All correct
  est7 <- modEIF2a(dat)$est
  CI7 <-  modEIF2a(dat)$CI
  #######################################################
  #Case 2: Outcome model is right, but PS is not
  est8 <- modEIF2b(dat)$est
  CI8 <-  modEIF2b(dat)$CI
  #######################################################
  # Case 3: (PS is correct, but the outcome model is not correct)
  est9 <- modEIF2c(dat)$est
  CI9 <-  modEIF2c(dat)$CI
  #######################################################
  # Both are not correct
  est10 <- modEIF2d(dat)$est
  CI10 <-  modEIF2d(dat)$CI
  #######################################################
  # OR estimator 
  est11 <- mod_OR_c(dat)
  est12 <- mod_OR_w(dat)
  write(c(i,est1,CI1, est2,CI2, est3,CI3,est4,CI4, est5,CI5, est6,CI6, est7,CI7,est8,CI8, est9,CI9, est10,CI10, est11, est12),file="Study2results1000.txt",ncolumns=40,append=T)
}


res1<-read.table("Study2results1000.txt",header=F)
head(res1)
dim(res1)
#truth mRR
psi = 0.04664484
psi = 0.047
colnames(res1)[seq(2, 30,3)] <- c("IPWCorrect", "IPWrong", "BothCorrect1", "OutCorrect1", "PSCorrect1", "BothWrong1", "BothCorrect2", "OutCorrect2", "PSCorrect2", "BothWrong2")
(apply(res1[,seq(2, 30,3)], 2, mean) - psi)
apply(res1[,seq(2, 30,3)], 2, sd)/sqrt(1000)
library(vioplot)
vioplot( res1$IPWCorrect, res1$IPWrong, res1$BothCorrect1,res1$OutCorrect1 , res1$PSCorrect1 , res1$BothWrong1, names=c("IPWCorrect", "IPWrong","BothCorrect", "OutCorrect", "PSCorrect", "BothWrong1"),  
         xlab= "n = 2000", ylab = expression(paste( psi, " estimates")) )
abline(h = psi, col = "red")


#CIs
mean(psi<=res1$V4 & psi>=res1$V3) # 50
mean(psi<=res1$V7 & psi>=res1$V6) # 0
mean(psi<=res1$V10 & psi>=res1$V9) #88
mean(psi<=res1$V13 & psi>=res1$V12) #50
mean(psi<=res1$V16 & psi>=res1$V15) # 73
mean(psi<=res1$V19 & psi>=res1$V18) #0

