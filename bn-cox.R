library(missMethods)
library(DMwR)
library("VIM")
library(sparsebn)
library(coxed);library(Rmpfr);library(missForest);library(glmnet);library(survminer);library(bnviewer)
library(dplyr);library(plyr);library(survival);library(bnlearn);library(parallel);library(bnlearn);library(survey)
library(pROC);library(ResourceSelection);library(rmda);library(pec); library(survcomp)
library(randomForest)
library(missForest)
library(mice)
library(Publish)
library(igraph)

#The data set is divided into training set and verification set
#and the training set is used to construct and learn the network
data$id <- 1:nrow(data)
summary(data$time)
train <- na.omit(data)
trainid <- train$id
test <- data[-trainid,c("status","time",var)]
data1 <- train[,c("status","time",var)]
##############univariate cox##################
UniCox_adjust <- function(x){
  data1 <- data1[complete.cases(data1[,x]),]
  y<- Surv(time=data1$time,event=data1$status==1)
  FML <- as.formula(paste0 ("y~",x))
  GCox<- coxph(FML,data=data1[complete.cases(data1[,x]),])
  GSum <- summary(GCox)
  HR <- sprintf("%.3f", GSum$coefficients[,2])
  P <- sprintf("%.3f", GSum$coefficients[,5])
  Wald <- sprintf("%.3f", GSum$coefficients[,4])
  Beta <- sprintf("%.3f", GSum$coefficients[,1])
  SE <- sprintf("%.3f", GSum$coefficients[,3])
  HR.confint.lower <- sprintf("%.3f", GSum$conf.int[,"lower .95"])
  HR.confint.upper <- sprintf("%.3f", GSum$conf.int[,"upper .95"])
  # HR <- paste0(HR, " (",
  #              HR.confint.lower, "-", HR.confint.upper, ")")
  # Mulcox <- data.frame("Characteristics" = rownames(GSum$coefficients),
  #                      Beta,SE,HR,HR.confint.lower,HR.confint.upper,Wald,P)
  Mulcox <- data.frame("Characteristics" = rownames(GSum$coefficients),
                       Beta,HR,HR.confint.lower,HR.confint.upper,P)
  var_step <- c(x)
  Mulcox$Var <- x
  return(Mulcox)
}
covariates1 <- var
time_3 <- system.time(UniVar_adjust <- lapply(covariates1, UniCox_adjust))
time_4 <- system.time(UniVar_adjust <- ldply(UniVar_adjust,data.frame))
UniVar <- UniVar_adjust[ifelse(UniVar_adjust$P < 0.05,TRUE,FALSE),]
UniVar <- UniVar[as.numeric(UniVar$HR.confint.lower)>1|as.numeric(UniVar$HR.confint.upper)<1,]
var <- unique(UniVar$Var)

###################################################
# LASSO

x=data.matrix(train[,setdiff(var,vardrop)])
x=data.matrix(train[,var])
y=as.double(train[,2])
status=as.double(train[,1])
surv=Surv(y,status)
cv.fit =cv.glmnet(x, surv, family="cox", nfolds=10)
fit <- glmnet(x, surv,family = "cox") 
tiff("F:/writing/supplement/lasso1.tiff", res = 300, width=17, height = 12, units = "cm")

plot(fit, xvar = "lambda",label = TRUE, las=1)

dev.off()
tiff("F:/writing/supplement/lasso2.tiff", res = 300, width=17, height = 12, units = "cm")
plot(cv.fit, las =1)
dev.off()
var_step = rownames(as.data.frame(which(coef(cv.fit, s = "lambda.1se")[,1]!=0)))
####################################

#################################
# Conditional probability table------
FML <- as.formula(paste0("Surv(time, status) ~",paste0(var_step, collapse = "+")))
c1 <- coxph(FML, data=train, x = TRUE)
summary(c1)
H0<-basehaz(c1, centered=T)
cases = expand.grid(lapply(train[,var_step],levels))
pre1 <- predict(c1,type="lp",newdata=cases)
t <- 3
cuttime <- 365*t
b <- which(abs(H0$time-cuttime)==min(abs(H0$time-cuttime)))
s0<-exp(-H0[b,1])
cases$prob_surv3y<-s0^exp(pre1)
cases$prob_surv180days<-s0^exp(pre1)


# Structure learning: Building networks from variables filtered by single factors
#-----------------bootstrap--------
boot.model.averaged <- function(data, R = 200, threshold = 0.7, blacklist=NULL,whitelist=NULL,
                                algo = c("pc.stable","fast.iamb","iamb","gs","hc","tabu",
                                         "mmhc","rsmax2")){
  if (!require('bnlearn')) {
    stop('The package bnlearn was not installed')
  }
  netlist <- list(NULL)
  length(netlist) <- length(algo)
  options(warn = -1)
  for (i in 1:length(algo)) {
    netlist[[i]] <- boot.strength(data, R = R, algorithm = algo[i],
                                  algorithm.args = list(blacklist=blacklist, whitelist=whitelist))
  }
  options(warn = 1)
  for (i in 1:length(algo)) {
    netlist[[i]] <- averaged.network(netlist[[i]], threshold = threshold)
  }
  arcs.model <- custom.strength(netlist, nodes = colnames(data),cpdag = FALSE)
  res <- averaged.network(arcs.model, threshold = threshold)
  return(res)
}

iter.model.average <- function(data, R = 200, threshold = 0.7, blacklist=NULL,whitelist=NULL,
                               algo = c("pc.stable","iamb","fast.iamb","gs","hc","tabu","mmhc","rsmax2"), 
                               iter.max=Inf){
  if (!require('dplyr')) {
    stop('The package dplyr was not installed')
  }
  iter <- 0
  improve <- 1
  while (improve > 0 & iter < iter.max) {
    iter <- iter + 1
    if (iter==1) {
      bn.iter <- boot.model.averaged(data, R = R, threshold = threshold, blacklist=blacklist,
                                     whitelist=whitelist, algo = algo)
      new.arc <- data.frame(arcs(bn.iter))
      arc.imp <- new.arc
    } else {
      bn.iter <- boot.model.averaged(data, R = R, threshold = threshold, blacklist=blacklist,
                                     whitelist=new.arc, algo = algo)
      new.arc <- data.frame(arcs(bn.iter))
      arc.imp <- suppressMessages(dplyr::anti_join(new.arc, arc.imp))
    }
    improve <- nrow(arc.imp)
  }
  return(bn.iter)
}

#--------------------------------------

#--------blacklist----------
lists <- lapply(zlvar, function(r) t(sapply(basename, function(x) c(r, x), USE.NAMES = FALSE)))
lists <- do.call("rbind", lists)
bl2 <-  matrix(lists, nodes, ncol = 2)
black_list <- rbind(black_list,bl2)
#-------------------------

#------white list-------------
white_list <- read.csv(file="F:/writing/blacklist.csv")
#------------------------
bnn = tabu(train[,var], score="aic",blacklist = black_list,whitelist = white_list)
fit=bn.fit(bnn,train[,var],replace.unidentifiable=TRUE)
#model average
abn <- boot.strength(train[,var], R = 200, algorithm = "tabu",
                     algorithm.args = list(blacklist=black_list, whitelist=white_list))
fita <- averaged.network(abn, threshold = 0.7)
#add status to the obtained BN--------------- 
st=modelstring(fit)
paste0(var_step, collapse = ":")
st2=paste(st, "[status|gender:Age:Smoking:drinking:NSCLC:targeted_therapy:COPD:pneumonia:STAGE:ILD:respiratory_failure:Fibrinogen]",sep = "")
mybn=as.bn(st2)
#generate a CPT for this STATUS node:
survCPT <- c()
survCPT <-c(t(cbind(cases$prob_surv180days,1-cases$prob_surv180days)))

lev <- c()
for(x in var_final) {
  a <- length(levels(train[,var_step][,x]))
  lev <- c(lev,a)
}
dim(survCPT) = c(2,lev)

t <- 3
train$status <-ifelse(train$status==1 & train$time<=as.vector(t %*% 365),1,0)
train$status <- as.factor(train$status)
dimnms <- lapply(train[,c('status',var_final)],levels)
dimnames(survCPT) = dimnms

# Add it to a list of the other CPTs
CPTs = lapply(fit, function(n) { n$prob })
CPTs[["status"]] = survCPT

cfit = custom.fit(mybn, dist = CPTs)
#-------------------------
type_test <- as.data.frame(unlist(sapply(train, class)))
names(type_test ) <- "TYPE"
num_names <- rownames(subset(type_test, TYPE=="numeric"))
num_names <- setdiff(num_names,c("status"))
for (i in num_names){
  train[,i] <- as.factor(train[,i])
}

vars_evidence <- names(train[,var_step])
pred <- c()

#lw
#-----
str2 <- "(status=='1')"
myprocess <- function(x){
  vars_evidence_nomi <- vars_evidence[which(is.na(x[vars_evidence])==F)]
  str2 <- "(status=='1')"
  set.seed(1)
  risk <- cpquery(cfit, event = eval(parse(text = str2)),
                  evidence = as.list(x[vars_evidence_nomi]), method = "lw")
  return(risk)
}
time_1 <- system.time(pre_tr <- apply(train[,var_step],MARGIN = 1,myprocess))
#-------
############################################


########Evaluation###################
#AUC
ares.coxs2 <- roc(train$status,pre_tr,plot=TRUE,legacy.axes=T,print.thres=T,print.auc=T)
ares.coxs2$auc
ci(ares.coxs2)

cuttime <- 365*t
train$status1 <- ifelse(train$time>cuttime,0,train$status)
train$time1 <- ifelse(train$time>=cuttime,cuttime,train$time)
cindex_cli <- concordance.index(aabb,surv.time = train$time, surv.event = train$status, method = "noether")$"c.index"

#Calibration curve
calibration_plot <- function(data, pred){
  pihat <- pred
  breaks=c(0,quantile(pihat, probs=seq(0.1,0.9,0.1)),1)
  pihatcat <- cut(pihat, breaks, labels=FALSE)
  meanprobs <- array(0, dim=c(10,2))
  expevents <- array(0, dim=c(10,2))
  obsevents <- array(0, dim=c(10,2))
  stdprobs <- array(0, dim=c(10,2))
  obsprobs <- array(0, dim=c(10,2))
  for (i in 1:10) {
    meanprobs[i,1] <- mean(pihat[pihatcat==i])
    stdprobs[i,1] <- sd(pihat[pihatcat==i])
    expevents[i,1] <- sum(pihatcat==i)*meanprobs[i,1]
    obsevents[i,1] <- sum(data$status[pihatcat==i])
    if(sum(pihatcat==i)==0){
      obsprobs[i,1] <- 0
    }
    else {
      obsprobs[i,1] <- sum(data$status[pihatcat==i])/sum(pihatcat==i)
    }
    
    meanprobs[i,2] <- mean(1-pihat[pihatcat==i])
    stdprobs[i,2] <- sd(1-pihat[pihatcat==i])
    expevents[i,2] <- sum(pihatcat==i)*meanprobs[i,2]
    obsevents[i,2] <- sum(1-data$status[pihatcat==i])
    if(sum(pihatcat==i)==0){
      obsprobs[i,2] <- 1
    }
    else {
      obsprobs[i,2] <- 1-(sum(data$status[pihatcat==i])/sum(pihatcat==i))
    }
  }
  xx <- meanprobs[,c(1)]
  yy <- obsprobs[,c(1)]
  plot(xx, yy, type = "n", xlab = "Predicted probability", ylab = "Oberserved probability",xlim = c(0,1),ylim = c(0,1),main = "calibration")
  points(xx, yy, type = "p", pch = 1, col = "black", lty = 1)
  abline(0,1,lty=3,lwd=1,col=c(rgb(0,0,0,maxColorValue=255)))
  lm9<-lm(yy~xx)##
  summary(lm9)
  Corner_text <- function(text, location="topleft"){
    legend(location,legend=text, bty ="n", pch=NA) 
  }
  return(summary(lm9))
}
calibration_plot(train,aabb)

#DCA --------------
a3 <- as.data.frame(cbind(test$status,pre_te))
a1 <- as.data.frame(cbind(test$status,pre_coxte))
model1 <- decision_curve(V1~pre_coxte,
                         data = a1,
                         fitted.risk = TRUE, 
                         thresholds = seq(0, 1, by = .05),
                         bootstraps = 25) 
model3 <- decision_curve(V1~pre_te,
                         data = a3,
                         fitted.risk = TRUE, 
                         thresholds = seq(0, 1, by = .05),
                         bootstraps = 25) 
plot_decision_curve( list(model3, model1), 
                     main="Decision Curve",
                     curve.names = c("CSBN model","CPH"),
                     cost.benefit.axis=FALSE,
                     col = c("blue", "red"), 
                     lty = c(1,2), 
                     lwd = c(3,2, 2, 1),
                     confidence.intervals=FALSE,
                     # xlim = c(0, 1), ylim = c(-.2,2),
                     # lengend.size
                     legend.position = "topright"
                     #standardize = FALSE
)
dev.off()

#------------------

#############################



#validation dataset
test <- test[which((test$status==1)|(test$time>=1095)),]
test$status <- ifelse(test$status==1 & test$time<=as.vector(t %*% 365),1,0)
time_2 <- system.time(pre_test <- apply(test[,var_step],MARGIN = 1,myprocess))


##########3-year OS Prediction Performance#############
#AUC
ares.coxs2 <- roc(test$status,aabb,plot=TRUE,legacy.axes=T,print.thres=T,print.auc=T)
ares.coxs2$auc
ci(ares.coxs2)

test$status1 <- ifelse(test$time>cuttime,0,test$status)
test$time1 <- ifelse(test$time>=cuttime,cuttime,test$time)
cindex_cli <- concordance.index(aabb,surv.time = test$time, surv.event = test$status, method = "noether")$"c.index"

calibration_plot(test,pre_test)

