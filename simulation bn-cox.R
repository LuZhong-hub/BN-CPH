############simulation#################
rm(list = ls())
library(coxed):library(Rmpfr)
library(dplyr);library(plyr);library(survival);library(bnlearn);library(parallel);library(bnlearn);library(survey)
library(pROC);library(ResourceSelection);library(rmda);library(pec); library(survcomp)
library(missMethods)
library(DMwR)
library("VIM")
library(sparsebn)
library(missForest)
library(mice)
library(Publish)
library(igraph)
library(dplyr)
library(data.table)
library(tictoc)
library(survival)
library(fastDummies)
library(pec)
library(caret)
library(ggpubr)




tabel <-  data  %>% group_by(gender,Age,Smoking,drinking,family_history,NSCLC,Radiotherapy,CHEMO,targeted_therapy,COPD,pneumonia,pleural_effusion,hypertension,STAGE,URI,lung_abscess,pulmonary_embolism,Pulmonary_heart_disease,ILD,respiratory_failure,pneumothorax,Redbloodcell,Monocyte,Eosinophil,Fibrinogen,Directbilirubin) %>% dplyr::summarise(n = n())
tabel$prop <- tabel$n / nrow(data)

tabel$case <- 1:nrow(tabel)
lists <- split(tabel[, var], f = tabel$case)


data$time <- data$time/365.25
data <- data%>%filter(time != 0)



FML <- as.formula(paste0("Surv(time, status) ~",paste0(var_step, collapse = "+")))
logreg <-  survreg(FML,
                   data = data, dist = "lognormal",
                   robust = TRUE)

coefs <- matrix(logreg$coefficients) # intercept is parameter mu
sigma <- as.vector(logreg$scale)


cens_times <- Surv(data$time, data$status == 0)
plot(cens_times)
summary(cens_times)
distr <- survreg(cens_times ~ 1, dist = "weibull")
cens <- rweibull(n = nrow(data), shape = 1 / distr$scale, scale = exp(distr$coefficients))
summary(cens)

# data simulation procedure

data_simulator <- function(n = 2500) {
  
  ind <- sample(x = 1:nrow(tabel), size = n, prob = tabel$prop, replace = TRUE)
  data_list <- vector(mode = "list", length = length(ind))
  
  for (i in 1:length(ind)) {
    
    value <- ind[i]
    line <- lists[[value]]
    data_list[[i]] <- line
  }
  
  data_sim <- as.data.frame(rbindlist(data_list))
  
  mat_data <- model.matrix(~., data_sim[, var_step])
  

  time <- exp(mat_data %*% coefs + sigma*rnorm(n = n, 0, 1))#summary(time)
  cens <- rweibull(n = n, shape = 1 / distr$scale, scale = exp(distr$coefficients)) 
  cens[cens == 0] <- (1/365.25)
  status <- as.numeric(time <= cens)
  time <- pmin(time, cens)
  data_sim <- cbind(data_sim, time, status)
  
  return(data_sim)
}
data1 <- data_simulator(n = 5000)

# step1

big_list <- replicate(n = 500, {
  list(data_simulator(n = 6000))
})

MNAR_gen_imp2 <- function(data_c, mmode="mar_wrapper", mis_var_prop=seq(0.1,0.4,0.1),
                          impute_list, cores=5,var_step=var_step) {
  if (cores==1) {
    results <- list()
    #train test
    index <- createDataPartition(data_c$status, p = 1/2, list = FALSE)
    train <- data_c[index, ]
    test  <- data_c[-index,]
    var <- setdiff(names(data_c),c("status","time"))
    # paste0(paste0("'",var_step,"'"),collapse = ",")
    var_step <- c('gender','Age','Smoking','drinking','NSCLC','COPD','pneumonia','STAGE','ILD','respiratory_failure','Fibrinogen')
    for (i in seq_along(mis_var_prop)) {
      prop <- mis_var_prop[i]
      
      test_m_res <- prodNA(test[,var], noNA = prop)
      data <- rbind(train[,var], test_m_res )
      data[,!names(data) %in% c("status","time")] <- lapply(data[,!names(data) %in% c("status","time")], 
                                                            as.factor)
      result_list <- list()
      tests <- data[-c(1:nrow(train)),]
      for (j in seq_along(impute_list)) {
        method=impute_list[j]
        print(method)
        result_list[[j]] <- do.call(method, list(tests[,var_step]))
      }
      data1 <- rbind(train,cbind(test_m_res,test[,c("status","time")]))
      data1[,!names(data1) %in% c("status","time")] <- lapply(data1[,!names(data1) %in% c("status","time")], 
                                                              as.factor)
      train <- data1[c(1:nrow(train)),]
      test1 <- data1[-c(1:nrow(train)),]
      result_list[[j+1]] <- test1
      results[[i]] <- list(train=train,test=result_list,prop=prop)
    }
  } else if (cores>1) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    results <- foreach(prop=mis_var_prop, .export=c(impute_list,mmode), .packages=c('missForest', 'randomForest', 'mice','VIM','caret')) %dopar% {
      source('~/zhonglu/code/da/mnar_generate.R')
      source('~/zhonglu/code/da/miss type produce.R')
      index <- createDataPartition(data_c$status, p = 1/2, list = FALSE)
      train <- data_c[index, ]
      test  <- data_c[-index,]
      var <- setdiff(names(data_c),c("status","time"))
      var_step <- c('gender','Age','Smoking','drinking','NSCLC','COPD','pneumonia','STAGE','ILD','respiratory_failure','Fibrinogen')
      test_m_res <- do.call(mmode, list(test[,var],prop,var_step))
      
      
      data <- rbind(train[,var],test_m_res)
      data[,!names(data) %in% c("status","time")] <- lapply(data[,!names(data) %in% c("status","time")], 
                                                            as.factor)
      result_list <- list()
      tests <- data[-c(1:nrow(train)),]
      for (j in seq_along(impute_list)) {
        method=impute_list[j]
        print(method)
        result_list[[j]] <- do.call(method, list(tests[,var_step]))
      }
      data1 <- rbind(train,cbind(test_m_res,test[,c("status","time")]))
      data1[,!names(data1) %in% c("status","time")] <- lapply(data1[,!names(data1) %in% c("status","time")], 
                                                              as.factor)
      train <- data1[c(1:nrow(train)),]
      test1 <- data1[-c(1:nrow(train)),]
      result_list[[j+1]] <- test1
      result_list <- list(train=train,test=result_list,prop=prop)
      result_list
    }
    stopCluster(cl)
  } else {cat('Improper argument: cores!')}
  return(results=results)
}

MNAR_gen_imp_big_lc<- function(datas,mmode) {
  
  datas[,!names(datas) %in% c("status","time")] <- lapply(datas[,!names(datas) %in% c("status","time")], 
                                                          as.numeric)
  aa <- MNAR_gen_imp2(datas, mmode=mmode, mis_var_prop=seq(0.1,0.4,0.1), 
                      impute_list=c('RF_wrapper', 'knn_wrapper','mice_wrapper'), cores=10,var_step=var_step)
  
  return(aa)
}
# step2
big_aa_mnar <- lapply(big_list,MNAR_gen_imp_big_lc,mmode="mnar_wrapper")



mnar_wrapper <- function (data_c, prop = .3,var_step=var_step) {
  result2 <- produce_NA(data = data_c[,var_step], mechanism="MNAR", perc.missing = prop)
  data_c[,var_step] <- result2$data.incomp
  return (data_c)
}

#######Evaluation of simulation results#################################################
csbn23 <- function(train,test,c1,var_step, t = 3) {
  H0<-basehaz(c1, centered=T)
  cases = expand.grid(lapply(train[,var_step],levels))  
  cases$prob_surv3y <- predictSurvProb(object = c1, newdata = cases,
                                       times = 3)
  survCPT <- c()
  survCPT <-c(t(cbind(cases$prob_surv3y,1-cases$prob_surv3y)))
  lev <- c()
  for(x in var_step) {
    a <- length(levels(train[,x]))
    lev <- c(lev,a)
  }
  dim(survCPT) = c(2,lev)
  train$status <-ifelse(train$status==1 & train$time<= t,1,0)
  train$status <- as.factor(train$status)
  dimnms <- lapply(train[,c('status',var_step)],levels)
  dimnames(survCPT) = dimnms
  # BN
  var <- setdiff(names(train),c('status','time'))
  bnn = tabu(train[,var], score="aic")
  fit=bn.fit(bnn,train[,var],replace.unidentifiable=TRUE)
  CPTs = lapply(fit, function(n) { n$prob })
  CPTs[["status"]] = survCPT
  st=modelstring(fit)
  st2=paste(st, "[status|",paste0(var_step, collapse = ":"),"]",sep = "")
  mybn=as.bn(st2)
  cfit = custom.fit(mybn, dist = CPTs)
  # prediction
  train$status <- as.numeric(as.character(train$status))
  vars_evidence <- var
  myprocess <- function(x){
    vars_evidence_nomi <- vars_evidence[which(is.na(x[vars_evidence])==F)]
    set.seed(1)
    risk <- cpquery(cfit, event = eval(parse(text = "(status=='1')")),
                    evidence = as.list(x[vars_evidence_nomi]), method = "lw")
    return(risk)
  }
  test <- test[(test$status==1)|(test$time>=3),]
  test$status <- ifelse(test$status==1 & test$time<=t,1,0)
  pre_te <- apply(test[,vars_evidence],MARGIN = 1,myprocess)
  ares.coxs2 <- roc(test$status,pre_te,plot=TRUE,legacy.axes=T,print.thres=T,print.auc=T,main="b")
  auc_cli <-ares.coxs2$auc
  # cindex_cli <- concordance.index(pre_te,surv.time = test$time, surv.event = test$status, method = "noether")$"c.index"
  return(auc_cli)
}


MNAR_evalue <- function(data,var_step, cores=5,t=3) {
  if (cores==1) {
    results_cox <- vector(mode = "list", length = length(data))
    for (i in 1:length(data)) {
      train <- data[[i]]$train
      test <- data[[i]]$test
      test4 <- test[[length(test)]]
      var <- setdiff(names(train),c('status','time'))
      FML <- as.formula(paste0("Surv(time, status) ~",paste0(var_step, collapse = "+")))
      cox_all <- coxph(FML,
                       data = train,
                       method = "breslow", x = TRUE, y = TRUE)
      auc_cox <- list()
      for (j in 1:(length(test)-1)){
        
        test_new <- cbind(test[[j]],test4[,c("status","time")])
        test_new <- test_new[(test_new$status==1)|(test_new$time>=t),]
        test_new$status <- ifelse(test_new$status==1 & test_new$time<=t,1,0)
        lp <- predict(cox_all, newdata = test_new[,var_step] , type = "lp")
        probs_cox_all <- predictSurvProb(object = cox_all, newdata = test_new[,var_step] ,
                                         times = 3)
        ares.coxs2 <- roc(test_new$status,1-probs_cox_all,plot=TRUE,legacy.axes=T,print.thres=T,print.auc=T)
        auc_cox[[j]] <- ares.coxs2$auc
        # cindex_cox[[j]] <- concordance.index(x = lp,
        #                                      surv.time = test_new$time,
        #                                      surv.event = test_new$status)$c.index
      }
      auc_cox[[j+1]] <- csbn(train,test4,cox_all,var_step,t=3)
      results_cox[[i]] <- list(auc=auc_cox,missrate=data[[i]]$prop)
    }
  } else if (cores>1) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    results_cox <- foreach(i=1:length(data), .packages=c('bnlearn','pec', 'survcomp', 'survival','pROC')) %dopar% {
      source('~/zhonglu/code/da/mnar_generate.R')
      train <- data[[i]]$train
      test <- data[[i]]$test
      test4 <- test[[length(test)]]
      var <- setdiff(names(train),c('status','time'))
      FML <- as.formula(paste0("Surv(time, status) ~",paste0(var_step, collapse = "+")))
      cox_all <- coxph(FML,
                       data = train,
                       method = "breslow", x = TRUE, y = TRUE)
      auc_cox <- list()
      for (j in 1:(length(test)-1)){
        test_new <- cbind(test[[j]],test4[,c("status","time")])
        test_new <- test_new[(test_new$status==1)|(test_new$time>=t),]
        test_new$status <- ifelse(test_new$status==1 & test_new$time<=t,1,0)
        probs_cox_all <- predictSurvProb(object = cox_all, newdata = test_new[,var_step] ,
                                         times = 3)
        lp <- predict(cox_all, newdata = test_new[,var_step] , type = "lp")
        ares.coxs2 <- roc(test_new$status,1-probs_cox_all,plot=TRUE,legacy.axes=T,print.thres=T,print.auc=T)
        auc_cox[[j]] <- ares.coxs2$auc
        # cindex_cox[[j]] <- concordance.index(x = lp,
        #                                      surv.time = test_new$time,
        #                                      surv.event = test_new$status)$c.index
      }
      auc_cox[[j+1]] <- csbn23(train,test4,cox_all,var_step,t=3)
      results <- list(auc=auc_cox,missrate=data[[i]]$prop)
      results
    }
    stopCluster(cl)
  } else {cat('Improper argument: cores!')}
  return(results_cox)
}

bb1 <- MNAR_evalue(big_aa_mnar[[1]],var_step, cores=5)


MNAR_evalue_big<- function(datas) {
  
  aa <- MNAR_evalue(datas, var_step, cores=10)
  
  return(aa)
}
# step3
#非随机缺失
big_bb1 <- lapply(big_aa_mnar,MNAR_evalue_big)


#######结果整理#################################################

save(SUMMX, file = "~/zhonglu/data/da/simulations_cox.RData")


result <- big_bb1
SUMMX <- data.frame()
modeln <- c('CPH with missForest', 'CPH with KNN','CPH with MICE',"CSBN")
for (i in 1:length(result)){
  SUMM <- data.frame()
  for (j in 1:length(result[[i]])){
    bb <- data.frame()
    for (k in seq(length(result[[i]][[j]]$auc))){
      cc <- result[[i]][[j]]$auc
      cidx<- cbind( cc[[k]],modeln[k],result[[i]][[j]]$missrate,i)
      bb <- rbind(bb, cidx)
    }
    
    SUMM <- rbind(SUMM, bb)
  }
  SUMMX <- rbind(SUMMX, SUMM)
}
names(SUMMX) <- c("AUC","model","missrate","Group")
SUMM <- SUMMX

library(ggplot2)
SUMM$model <- as.factor(SUMM$model)
SUMM$missrate <- as.factor(SUMM$missrate)
SUMM$AUC <- as.numeric(SUMM$AUC)
SUMM$missrate <- revalue(SUMM$missrate,c("0.1"='10%',"0.2"='20%',"0.3"='30%',"0.4"="40%"))

ggplot(SUMM,aes(missrate,AUC))+
  stat_boxplot(aes(fill=model),geom="errorbar",width=0.1,size=0.3,position=position_dodge(0.6),color="black")+
  geom_boxplot(aes(fill=model),
               position=position_dodge(0.6),
               size=0.3,
               width=1,
               color="black",
               #outlier.color = "black",
               outlier.fill = "red",
               outlier.shape = 19,
               outlier.size = 1.5,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  theme_bw()+xlab("Missing Prop") + ylab("AUC")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14))
