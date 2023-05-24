#Common data-generating function

#popsize set at 1500000, can increase in large ssize is needed
#ssize set at 500, may have trouble when not enough patient available

datagen<-function(seed=sample(1:1000000,size=1),ssize=5000,popsize=15000000,OR_C=2,OR_WI=1,OR_WC=4,OR_H=1.5,em=1,cfV0=F,cfV1=F,return_full=F){
  set.seed(seed)
  #generate data
  C<- runif(n=popsize, -1,1)
  U1<-rbinom(n=popsize,size=1,prob=0.5) #affects both
  U2<-rbinom(n=popsize,size=1,prob=0.5) #affects covid
  
  if(cfV0==T) V=rep(0,popsize); if(cfV1==T) V=rep(1,popsize);
  if(cfV0==F&cfV1==F){
    p = plogis(1 -  C - C^2 - 2 * sin(C))
    V<-rbinom(prob=p,size=1,n=popsize) #prevalence is around 0.71
  }
  #Infection (with something) has some common risk factors U1 and C
  Infec<-rbinom(prob=plogis(0.5*C-3+0.5*U1),size=1,n=popsize) #current prevalence around 0.007
  
  #Infected with COVID
  pcovid <- plogis( -log(OR_C)*V -2+ C - C^3 - abs(C)+em*V*C+log(3)*U2*(1.5-V)-2*U1)
  #symptoms based on infection
  #can come from either or both infections, if present
  W=rep(0,popsize)
  W[Infec==1]<-rbinom(prob=plogis(0.5+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
  W[Infec_COVID==1]<-rbinom(prob=plogis(-2+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
  #mean(W[Infec==1|Infec_COVID==1]) #9%
  #mean(W[Infec_COVID==1]) #27%
  #mean(W[Infec==1]) #2%
  #mean(W[Infec_COVID==1&V==1]) #8%
  
  #hospitalization, only possible if symptoms present
  H=rep(0,popsize)
  H[W==1]<-rbinom(prob=plogis(1.5+0.5*C[W==1]+log(OR_H)*V[W==1]-0.5*U1[W==1]),size=1,n=sum(W==1))
  #mean(H[W==1]) #83% with severe symptoms go to hospital
  #sum(H==1)
  #selection on outcome for testing (does not condition on infectious status, just being in the hospital)
  R<-sample(which(H==1),ssize) #sample randomly from those in hospital to control the study size
  
  if(return_full==F){
    dat<-as.data.frame(cbind(Y=Infec_COVID,V=V,C=C)[R,])
  } else{dat<-as.data.frame(cbind(Infec_COVID=Infec_COVID,Infec=Infec,H=H,W=W,V=V,C=C))}
  return(dat)
}



#true marginal RR
orvect<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T)
  datfull0<-datagen(cfV0=T,return_full=T)
  orvect[j]<-mean(datfull1$H*datfull1$Infec_COVID)/mean(datfull0$H*datfull0$Infec_COVID)
}
(psi = mean(orvect)) #2000, em = 1, 0.04230757; em = 2 0.07174366, em = 3, 0.1597932;
hist(orvect)

#true marginal RR for Infec_COVID
orvect1<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T)
  datfull0<-datagen(cfV0=T,return_full=T)
  orvect1[j]<-mean(datfull1$Infec_COVID)/mean(datfull0$Infec_COVID)
}
mean(orvect1) #0.16
hist(orvect1)

#true marginal RR for Infec_COVID * W (severe COVID disease)
orvect2<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T)
  datfull0<-datagen(cfV0=T,return_full=T)
  orvect2[j]<-mean(datfull1$W*datfull1$Infec_COVID)/mean(datfull0$W*datfull0$Infec_COVID)
}
mean(orvect2) #0.04
hist(orvect2)
