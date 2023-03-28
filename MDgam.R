dat<-datagen(ssize=1000, em=1)
# Create training (70%) and test (30%) sets for the rsample::attrition data.
# Use set.seed for reproducibility
set.seed(123)
TNDdf_split <- rsample::initial_split(dat, prop = 0.5)
TNDdf_train <- rsample::training(TNDdf_split)
TNDdf_test  <- rsample::testing(TNDdf_split)
library(earth)
library(caret)
library(dplyr)
library(ggplot2)
library(resample)
library(h2o)
library(tidymodels)
#################################################################
# Outcome model
mod_mu<-function(TNDdf_train){
  formula <- Y ~ V + s(C)
  Out_mu <- gam(formula ,
               data = TNDdf_train,
               family = binomial)
  mu1=predict(Out_mu,newdata=as.data.frame(cbind(V=1,C=TNDdf_train$C)),type="response")
  mu0=predict(Out_mu,newdata=as.data.frame(cbind(V=0,C=TNDdf_train$C)),type="response")
  return(list(mod=Out_mu,mu1=mu1, mu0=mu0))
}


mod_m<-function(TNDdf_train){
  TNDdf_train = subset(TNDdf_train, select = -V)
  formula <- Y ~ s(C) 
  Out_m <- gam(formula ,
                data = TNDdf_train,
                family = binomial)
  m=predict(Out_m,type="response")
  return(list(mod=Out_m,m=m))
}

mod_mu(TNDdf_train)
mod_m(TNDdf_train)
#################################################################
# PS model with whole sampled data
mod_g<-function(TNDdf_train){
  # training
  formula <- V ~ s(C) 
  TNDmod_g <- gam(formula ,
               data = TNDdf_train,
               family = binomial)
  g1=predict(TNDmod_g,type="response")
  g1_v1=predict(TNDmod_g,type="response", newdata=as.data.frame(cbind(C=TNDdf_train$C,V=rep(1, nrow(TNDdf_train)))) )
  return(list(mod=TNDmod_g,g1=g1, g1_v1 = g1_v1, g1_v0 = (1-g1_v1)))
}

# PS model with only control data
mod_g_control<-function(TNDdf_train){
  TND_train_ctr <- subset(TNDdf_train, Y==0)
  # training
  formula <- V ~ s(C) 
  TNDmod_g <- gam(formula ,
                  data = TND_train_ctr,
                  family = binomial)
  g1_cont<-predict(TNDmod_g,type="response",newdata=as.data.frame(cbind(C=TNDdf_train$C,V=TNDdf_train$V, Y=TNDdf_train$Y)))
  g1_v1<-predict(TNDmod_g,type="response",newdata=as.data.frame(cbind(C=TNDdf_train$C,V=rep(1, nrow(TNDdf_train)),Y=TNDdf_train$Y)));
  return(list(mod=TNDmod_g,g1_cont=g1_cont, g1Cont_v1 = g1_v1, g1Cont_v0 = (1-g1_v1)))
}

mod_g(TNDdf_train)
mod_g_control(TNDdf_train)

###################################################################
###################################################################
de_wetPS <- function(TNDdf_train, dfpred){
  g1 <- mod_g(TNDdf_train)
  g1Cont <- mod_g_control(TNDdf_train)
  
  g1_v1=predict(g1$mod,type="response", newdata=as.data.frame(cbind(C=dfpred$C,V=rep(1, nrow(dfpred)))) )
  g1_v0 = 1-g1_v1
  
  g1Cont_v1<-predict(g1Cont$mod,type="response",newdata=as.data.frame(cbind(C=dfpred$C,V=rep(1, nrow(dfpred)),Y=dfpred$Y)));
  g1Cont_v0 = 1-g1Cont_v1
  
  
  w1 <- g1_v1/g1Cont_v1
  w0 <- g1_v0/g1Cont_v0
  return(list(dewet_v1 = w1, dewet_v0 = w0))
}

de_wetOut <- function(TNDdf_train, dfpred){
  mu <- mod_mu(TNDdf_train)
  m_mod <- mod_m(TNDdf_train)
  
  m=predict(m_mod$mod,type="response", newdata=as.data.frame(cbind(C=dfpred$C,V=rep(1, nrow(dfpred)),Y=dfpred$Y)))
  mu1=predict(mu$mod,newdata=as.data.frame(cbind(V=1,C=dfpred$C)),type="response")
  mu0=predict(mu$mod,newdata=as.data.frame(cbind(V=0,C=dfpred$C)),type="response")
  w1 <- (1 - m)/(1 - mu1)
  w0 <- (1 - m)/(1 - mu0)
  return(list(dewet_v1 = w1, dewet_v0 = w0))
}


###################################################################
###################################################################
wet <- de_wetPS(TNDdf_train, TNDdf_train)
w1 <- wet$dewet_v1; w0 <- wet$dewet_v0
