#Common data-generating function

#popsize set at 1500000, can increase in large ssize is needed
#ssize set at 500, may have trouble when not enough patient available

#Defaults (scenario 1) true marg RR = 
#OR_C<-3 #effect of vaccine on covid
#OR_WI<-1 #no effect of vaccine on W from other infection
#OR_WC<-5 #effect of vaccine on covid symptoms
#OR_H<-1.5 #effect of vaccine on hospitalization among those with symptoms
#em=0
#true marg RR = 0.04 ; true cond RR = 0.04 ; #ONLY GOOD TO TWO DECIMALS

#Scenario 2: Effect modification of V by C on Infec_COVID
#em=1
#true marg RR = 0.25; true cond RR = 0.23;

datagen<-function(seed=sample(1:1000000,size=1),ssize=500,popsize=1500000,OR_C=3,OR_WI=1,OR_WC=5,OR_H=1.5,em=0,cfV0=F,cfV1=F,return_full=F){
  set.seed(seed)
  
  #generate data
  C<-rnorm(n=popsize)
  U1<-rbinom(n=popsize,size=1,prob=0.5) #affects both
  U2<-rbinom(n=popsize,size=1,prob=0.5) #affects covid
  
  if(cfV0==T) V=rep(0,popsize); if(cfV1==T) V=rep(1,popsize);
  if(cfV0==F&cfV1==F){
    V<-rbinom(prob=plogis(0.5+0.3*C+sin(C)),size=1,n=popsize) #prevalence is around 0.61
  }
  #Infection (with something) has some common risk factors U1 and C
  Infec<-rbinom(prob=plogis(0.5*C-5+0.5*U1),size=1,n=popsize) #current prevalence around 0.007
  
  #Infected with COVID
  Infec_COVID<- rbinom(prob=plogis( 1-log(OR_C)*V -6 + 2*C - 0.15*exp(C)+log(3)*U2*(1.5-V)-2*U1), size=1,n=popsize) #0.009
  #symptoms based on infection
  #can come from either or both infections, if present
  W=rep(0,popsize)
  W[Infec==1]<-rbinom(prob=plogis(-4+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
  W[Infec_COVID==1]<-rbinom(prob=plogis(-2+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
  #mean(W[Infec==1|Infec_COVID==1]) #9%
  #mean(W[Infec_COVID==1]) #27%
  #mean(W[Infec==1]) #2%
  #mean(W[Infec_COVID==1&V==1]) #8%
  
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


######################################################################
###Run Simulation Study
modEIF1<-function(dat){
  TNDmod_g_col<-glm(V~C + sin(C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V + C + exp(C),family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  est=mean(dat$Y*dat$V/g1_cont - mu1*A1)/mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  return(est)
}

modEIF2<-function(dat){
  TNDmod_g_col<-glm(V~C,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V + C + exp(C),family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  est=mean(dat$Y*dat$V/g1_cont - mu1*A1)/mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  return(est)
}

modEIF3<-function(dat){
  TNDmod_g_col<-glm(V~C + sin(C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V + C,family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  est=mean(dat$Y*dat$V/g1_cont - mu1*A1)/mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  return(est)
}

modEIF4<-function(dat){
  TNDmod_g_col<-glm(V~C,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V + C,family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  est=mean(dat$Y*dat$V/g1_cont - mu1*A1)/mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  return(est)
}


CI0=CI1=CI2=CI3=CI4=CI5=CI6=CI7=c(0,0)
#seeds<-read.table("seeds.txt", header=F)$V1
seeds <- seq(1, 1000, 1)
for (i in 1:1000){
  dat<-datagen(ssize=500, em=0)
  #######################################################
  # Marginal RR 
  #Case 1: All correct
  est1<-modEIF1(dat)
  #######################################################
  #Case 2: Outcome model is right, but PS is not
  est2<-modEIF2(dat)
  #######################################################
  # Case 3: (PS is correct, but the outcome model is not correct)
  est3<-modEIF3(dat)
  #######################################################
  # Both are not correct
  est4<-modEIF4(dat)
  write(c(i,est1,est2,est3,est4),file="/Users/congjiang/Downloads/Study2results2v2.txt",ncolumns=40,append=T)
}


res1<-read.table("/Users/congjiang/Downloads/Study2results2.txt",header=F)
head(res1)
dim(res1)
#em=0
#true cRR = 0.04166 ; true mRR = 0.04363 ;
colnames(res1)[seq(2, 5)] <- c("BothCorrect", "OutCorrect", "PSCorrect", "BothWrong")
psi = 0.0469;
apply(res1, 2, mean) - psi
- mRR
1-RR1[c(3:8)]
apply(res1, 2, sd)/sqrt(nrow(res1))




# Library
library(ggplot2)

# create a dataset
data <- data.frame(
  name=c( rep("A",500), rep("B",500), rep("B",500), rep("C",20), rep('D', 100)  ),
  value=c( rnorm(500, 10, 5), rnorm(500, 13, 1), rnorm(500, 18, 1), rnorm(20, 25, 4), rnorm(100, 12, 1) )
)

# Most basic violin chart
p <- ggplot(res1, aes(x=name, y=value, fill=name)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin()
p

library(tidyr)
library(ggplot2)
library(dplyr)
data_wide <- iris[ , 1:4]
data_wide %>% 
  gather(key="MesureType", value="Val") %>%
  ggplot( aes(x=MesureType, y=Val, fill=MesureType)) +
  geom_violin()

data_wide <- res1[ , 2:5]
data_wide %>% 
  gather(key="MesureType", value="Val") %>%
  ggplot( aes(x=MesureType, y=Val, fill=MesureType)) +
  geom_violin()
