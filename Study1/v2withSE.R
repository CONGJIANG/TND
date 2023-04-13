
# Method 1: IPW estimator

mod_IPW<-function(TNDdf_train){
  gmodTND_cont<-mod_g_control(TNDdf_train)
  gTND_cont<-gmodTND_cont$g1_cont
  est=mean(TNDdf_train$Y*TNDdf_train$V/gTND_cont)/mean(TNDdf_train$Y*(1-TNDdf_train$V)/(1-gTND_cont))
  ipw <- ifelse(TNDdf_train$V == 1, 1/gTND_cont, 1/(gTND_cont))
  modY.ipw <- glm(Y ~ V, family=binomial(), weights = ipw, data=TNDdf_train)
  est.ipw <- exp(modY.ipw$coefficients[2])
  se.ipw <- vcovHC(modY.ipw)[2,2]
  return(list(est = est, CI = c(est - 1.96*se.ipw/sqrt(nrow(TNDdf_train)), est + 1.96*se.ipw/sqrt(nrow(TNDdf_train))) ))
}
mod_IPW(TNDdf_train)
# Method 2: Outcome Regression
# option for debiasing weights, which is optional and defaults to "option1" if not specified.
mod_OutReg <- function(TNDdf_train, deWet = c("PS", "Out")){
  mu <- mod_mu(TNDdf_train)
  if (deWet == "PS") {
    wet <- de_wetPS(TNDdf_train, TNDdf_train)
    w1 <- wet$dewet_v1; w0 <- wet$dewet_v0; 
  } else if (deWet == "Out") {
    wet <- de_wetOut(TNDdf_train, TNDdf_train)
    w1 <- wet$dewet_v1; w0 <- wet$dewet_v0; 
  } else {
    # Throw an error if an invalid option is selected
    stop("Invalid option selected.")
  }
  Q1 <- mu$mu1*w1
  Q0 <- mu$mu0*w0
  return(list(est=mean(Q1)/mean(Q0), Q1 = Q1, Q0 = Q0))
}

mod_OutReg(TNDdf_train, deWet = "PS")$est
mod_OutReg(TNDdf_train, deWet = "Out")$est

# Method 3: The proposed estimator that is based on the EIF
mod_eifsub <- function(TNDdf_train1, TNDdf_train2, deWet = c("PS", "Out")){
  # models 1) mu model
  mu <- mod_mu(TNDdf_train1)
  if (deWet == "PS") {
    wet <- de_wetPS(TNDdf_train1, dfpred = TNDdf_train2)
    w1 <- wet$dewet_v1; w0 <- wet$dewet_v0; 
  } else if (deWet == "Out") {
    wet <- de_wetOut(TNDdf_train1, dfpred = TNDdf_train2)
    w1 <- wet$dewet_v1; w0 <- wet$dewet_v0; 
  } else {
    # Throw an error if an invalid option is selected
    stop("Invalid option selected.")
  }
  Q1 <- predict(mu$mod,new_data=as.data.frame(cbind(V=1,C=TNDdf_train2$C)),type="response")*w1
  Q0 <- predict(mu$mod,new_data=as.data.frame(cbind(V=0,C=TNDdf_train2$C)),type="response")*w0
  #2) m model
  m <- mod_m(TNDdf_train1)
  m=predict(m$mod,type="response", new_data=as.data.frame(cbind(C=TNDdf_train2$C,V=rep(1, nrow(TNDdf_train2)),Y=TNDdf_train2$Y)))
  #3) PS model
  g1Cont <- mod_g_control(TNDdf_train1)
  g1_v1<-predict(g1Cont$mod,type="response",new_data=as.data.frame(cbind(C=TNDdf_train2$C,V=rep(1, nrow(TNDdf_train2)),Y=TNDdf_train2$Y)))
  g1_v0 <- 1 - g1_v1
  A1 <- (1 - TNDdf_train2$Y)*(TNDdf_train2$V - g1_v1)/(g1_v1 * (1 - m))
  A0 <- (1 - TNDdf_train2$Y)*((1-TNDdf_train2$V) - g1_v0)/(g1_v0* (1 - m))
  est=mean(TNDdf_train2$Y*TNDdf_train2$V/g1_v1 - Q1*A1)/mean(TNDdf_train2$Y*(1-TNDdf_train2$V)/g1_v0 - Q0*A0)
  var1 <- mean((TNDdf_train2$Y*TNDdf_train2$V/g1_v1 - Q1*A1)^2)/(mean(TNDdf_train2$Y*TNDdf_train2$V/g1_v1 - Q1*A1))^2
  var0 <- mean((TNDdf_train2$Y*(1-TNDdf_train2$V)/g1_v0 - Q0*A0)^2)/(mean(TNDdf_train2$Y*(1-TNDdf_train2$V)/g1_v0 - Q0*A0))^2
  selnRR <- sqrt(var1 + var0)
  return(list(est = est, se = selnRR))
}

mod_EIF <- function(TNDdf_train, deWet = c("PS", "Out")){
  # Data splitting
  TNDdf_split <- rsample::initial_split(TNDdf_train, prop = 0.5)
  TNDdf_train1 <- rsample::training(TNDdf_split)
  TNDdf_train2  <- rsample::testing(TNDdf_split)
  res1<- mod_eifsub(TNDdf_train1, TNDdf_train2, deWet)
  res2<- mod_eifsub(TNDdf_train2, TNDdf_train1, deWet)
  est <- nrow(TNDdf_train1)*res1$est/nrow(TNDdf_train) + nrow(TNDdf_train2)*res2$est/nrow(TNDdf_train)
  se <- nrow(TNDdf_train1)*res1$se/nrow(TNDdf_train) + nrow(TNDdf_train2)*res2$se/nrow(TNDdf_train)
  upper <- est*exp(+ 1.96*se/sqrt(nrow(TNDdf_train)))
  lower <- est*exp(- 1.96*se/sqrt(nrow(TNDdf_train)))
  return(list(est = est, se = se, CI =c(lower, upper)) )
}


mod_EIF(TNDdf_train, deWet = "PS")
mod_EIF(TNDdf_train, deWet = "Out")

# Method 4: The proposed estimator that is based on the EIF
mod_eifsub2 <- function(TNDdf_train1, TNDdf_train2, deWet = c("PS", "Out")){
  mu <- mod_mu(TNDdf_train1)
  mu1=predict(mu$mod,new_data=as.data.frame(cbind(V=1,C=TNDdf_train2$C)),type="response")
  mu0=predict(mu$mod,new_data=as.data.frame(cbind(V=0,C=TNDdf_train2$C)),type="response")
  g1Cont <- mod_g_control(TNDdf_train1)
  g1_v1<-predict(g1Cont$mod,type="response",new_data=as.data.frame(cbind(C=TNDdf_train2$C,V=rep(1, nrow(TNDdf_train2)),Y=TNDdf_train2$Y)));
  g1_v0 <- 1 - g1_v1
  
  A1 <- (1 - TNDdf_train2$Y)*(TNDdf_train2$V - g1_v1)/(g1_v1 * (1 - mu1))
  A0 <- (1 - TNDdf_train2$Y)*((1-TNDdf_train2$V) - g1_v0)/(g1_v0* (1 - mu0))
  est=mean(TNDdf_train2$Y*TNDdf_train2$V/g1_v1 - mu1*A1)/mean(TNDdf_train2$Y*(1-TNDdf_train2$V)/g1_v0 - mu0*A0)
  var1 <- mean((TNDdf_train2$Y*TNDdf_train2$V/g1_v1 - mu1*A1)^2)/(mean(TNDdf_train2$Y*TNDdf_train2$V/g1_v1 - mu1*A1))^2
  var0 <- mean((TNDdf_train2$Y*(1-TNDdf_train2$V)/g1_v0 - mu0*A0)^2)/(mean(TNDdf_train2$Y*(1-TNDdf_train2$V)/g1_v0 - mu0*A0))^2
  selnRR <- sqrt(var1 + var0)
  return(list(est = est, se = selnRR))
}

mod_EIF2 <- function(TNDdf_train){
  # Data splitting
  TNDdf_split <- rsample::initial_split(TNDdf_train, prop = 0.5)
  TNDdf_train1 <- rsample::training(TNDdf_split)
  TNDdf_train2  <- rsample::testing(TNDdf_split)
  # models
  res1<- mod_eifsub2(TNDdf_train1, TNDdf_train2)
  res2<- mod_eifsub2(TNDdf_train2, TNDdf_train1)
  est <- nrow(TNDdf_train1)*res1$est/nrow(TNDdf_train) + nrow(TNDdf_train2)*res2$est/nrow(TNDdf_train)
  se <- nrow(TNDdf_train1)*res1$se/nrow(TNDdf_train) + nrow(TNDdf_train2)*res2$se/nrow(TNDdf_train)
  upper <- est*exp( + 1.96*se/sqrt(nrow(TNDdf_train)) )
  lower <- est*exp( - 1.96*se/sqrt(nrow(TNDdf_train)) )
  return(list(est = est, se = se, CI =c(lower, upper)) )
}

mod_EIF2(TNDdf_train)

CI1=CI2=CI3=CI4=CI5=c(0,0)
#seeds<-read.table("seeds.txt", header=F)$V1
seeds <- seq(1, 1000, 1)
for (i in 1:1000){
  TNDdf_train <-datagen(ssize=2000, em=3)
  #######################################################
  #IPW with controls to fit g
  est.1<-mod_IPW(TNDdf_train)
  est1<- est.1$est
  CI1<- est.1$CI
  #######################################################
  # Outcome regression with propensity ratio
  est2<-mod_OutReg(TNDdf_train, deWet = "PS")$est
  # Outcome regression with outcome ratio
  est3<-mod_OutReg(TNDdf_train, deWet = "Out")$est
  #######################################################
  ## Method 3: The proposed estimator that is based on the EIF
  est.4 <- mod_EIF(TNDdf_train, deWet = "PS")
  est4<- est.4$est
  CI2<- est.4$CI
  est.5 <- mod_EIF(TNDdf_train, deWet = "Out")
  est5<- est.5$est
  CI3<- est.5$CI
  est.6 <- mod_EIF2(TNDdf_train)
  est6<- est.6$est
  CI4<- est.6$CI
  write(c(i,est1,CI1,est2,est3,est4,CI2,est5,CI3,est6,CI4),file="TNDStudy1hal2000_em3.txt",ncolumns=10,append=T)
}

res1<-read.table("TNDStudy1hal2000.txt",header=F)
res1 <- na.omit(res1)
head(res1)
dim(res1)
#em=1 true mRR = 0.04363 ;
colnames(res1) <- c("i", "IPW", "out-PS","out-Out", "EIF-PS", "EIF-Out", "EIF2")
mRR.em3 = 0.159513;
mRR.em1 = 0.04363
apply((res1- mRR.em1), 2, median)
apply(res1- mRR, 2, sd)
