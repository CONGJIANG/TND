
# Method 1: IPW estimator

mod_IPW<-function(TNDdf_train){
  gmodTND_cont<-mod_g_control(TNDdf_train)
  gTND_cont<-gmodTND_cont$g1_cont
  est=mean(TNDdf_train$Y*TNDdf_train$V/gTND_cont)/mean(TNDdf_train$Y*(1-TNDdf_train$V)/(1-gTND_cont))
  return(est)
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
mod_EIF <- function(TNDdf_train, deWet = c("PS", "Out")){
  # Data splitting
  TNDdf_split <- rsample::initial_split(TNDdf_train, prop = 0.5)
  TNDdf_train1 <- rsample::training(TNDdf_split)
  TNDdf_train2  <- rsample::testing(TNDdf_split)
  # models 1) mu model
  mu <- mod_mu(TNDdf_train1)
  if (deWet == "PS") {
    wet <- de_wetPS(TNDdf_train2, dfpred = TNDdf_train)
    w1 <- wet$dewet_v1; w0 <- wet$dewet_v0; 
  } else if (deWet == "Out") {
    wet <- de_wetOut(TNDdf_train1, dfpred = TNDdf_train)
    w1 <- wet$dewet_v1; w0 <- wet$dewet_v0; 
  } else {
    # Throw an error if an invalid option is selected
    stop("Invalid option selected.")
  }
  Q1 <- predict(mu$mod,newdata=as.data.frame(cbind(V=1,C=TNDdf_train$C)),type="response")*w1
  Q0 <- predict(mu$mod,newdata=as.data.frame(cbind(V=0,C=TNDdf_train$C)),type="response")*w0
  #2) m model
  m <- mod_m(TNDdf_train1)
  m=predict(m$mod,type="response", newdata=as.data.frame(cbind(C=TNDdf_train$C,V=rep(1, nrow(TNDdf_train)),Y=TNDdf_train$Y)))
  #3) PS model
  g1Cont <- mod_g_control(TNDdf_train2)
  g1_v1<-predict(g1Cont$mod,type="response",newdata=as.data.frame(cbind(C=TNDdf_train$C,V=rep(1, nrow(TNDdf_train)),Y=TNDdf_train$Y)))
  g1_v0 <- 1 - g1_v1
  A1 <- (1 - TNDdf_train$Y)*(TNDdf_train$V - g1_v1)/(g1_v1 * (1 - m))
  A0 <- (1 - TNDdf_train$Y)*((1-TNDdf_train$V) - g1_v0)/(g1_v0* (1 - m))
  est=mean(TNDdf_train$Y*TNDdf_train$V/g1_v1 - Q1*A1)/mean(TNDdf_train$Y*(1-TNDdf_train$V)/g1_v0 - Q0*A0)
  return(est)
}


mod_EIF(TNDdf_train, deWet = "PS")
mod_EIF(TNDdf_train, deWet = "Out")

# Method 4: The proposed estimator that is based on the EIF
mod_EIF2 <- function(TNDdf_train){
  # Data splitting
  TNDdf_split <- rsample::initial_split(TNDdf_train, prop = 0.5)
  TNDdf_train1 <- rsample::training(TNDdf_split)
  TNDdf_train2  <- rsample::testing(TNDdf_split)
  # models
  mu <- mod_mu(TNDdf_train1)
  mu1=predict(mu$mod,newdata=as.data.frame(cbind(V=1,C=TNDdf_train$C)),type="response")
  mu0=predict(mu$mod,newdata=as.data.frame(cbind(V=0,C=TNDdf_train$C)),type="response")
  g1Cont <- mod_g_control(TNDdf_train2)
  g1_v1<-predict(g1Cont$mod,type="response",newdata=as.data.frame(cbind(C=TNDdf_train$C,V=rep(1, nrow(TNDdf_train)),Y=TNDdf_train$Y)));
  g1_v0 <- 1 - g1_v1
  
  A1 <- (1 - TNDdf_train$Y)*(TNDdf_train$V - g1_v1)/(g1_v1 * (1 - mu1))
  A0 <- (1 - TNDdf_train$Y)*((1-TNDdf_train$V) - g1_v0)/(g1_v0* (1 - mu0))
  est=mean(TNDdf_train$Y*TNDdf_train$V/g1_v1 - mu1*A1)/mean(TNDdf_train$Y*(1-TNDdf_train$V)/g1_v0 - mu0*A0)
  return(est)
}
mod_EIF2(TNDdf_train)

CI1=CI2=CI3=CI4=CI5=c(0,0)
#seeds<-read.table("seeds.txt", header=F)$V1
seeds <- seq(1, 1000, 1)
for (i in 1:1000){
  TNDdf_train <-datagen(ssize=5000, em=3)
  #######################################################
  #IPW with controls to fit g
  est1<-mod_IPW(TNDdf_train)
  #bsobj<-bootstrap(dat,nbs=500,method1=mod_IPW,method2=mod_IPW_all)
  #var2<-bsobj$bs_var1
  #CI2<-bsobj$CI1
  #######################################################
  # Outcome regression with propensity ratio
  est2<-mod_OutReg(TNDdf_train, deWet = "PS")$est
  # Outcome regression with outcome ratio
  est3<-mod_OutReg(TNDdf_train, deWet = "Out")$est
  #######################################################
  ## Method 3: The proposed estimator that is based on the EIF
  est4 <- mod_EIF(TNDdf_train, deWet = "PS")
  est5 <- mod_EIF(TNDdf_train, deWet = "Out")
  est6 <- mod_EIF2(TNDdf_train)
  write(c(i,est1,est2,est3,est4,est5,est6),file="TNDStudy1earth5000_em3.txt",ncolumns=10,append=T)
}

res1<-read.table("TNDStudy1Hal2000.txt",header=F)
res1 <- na.omit(res1)
head(res1)
dim(res1)
#em=1 true mRR = 0.04363 ;
colnames(res1) <- c("i", "IPW", "out-PS","out-Out", "EIF-PS", "EIF-Out", "EIF2")
#mRR.em3 = 0.159513;
mRR.em1 = 0.04363
colMeans(res1- mRR.em1)
apply(res1- mRR, 2, sd)
