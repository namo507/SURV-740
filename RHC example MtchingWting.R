#RHC Example for PS matching and weighting
# install packages
# install.packages("tableone")
# install.packages("Matching")
# install.packages("glmmTMB")
# load packages
library(tableone)
library(Matching)
library(glmmTMB)
library(survey)
library(ipw)
library(sandwich)
#read in data
load("/Users/yli/Dropbox/COURSES/2024 Fall courses/SURV740/rhc.sav")

#treatment variable is swang1: Primary disease category
#x variables that we will use
#cat1: primary disease category #age #sex #meanbp1: mean blood pressure ...
#create a data set with just these variables, for simplicity
ARF<-as.numeric(rhc$cat1=='ARF')
CHF<-as.numeric(rhc$cat1=='CHF')
Cirr<-as.numeric(rhc$cat1=='Cirrhosis')
colcan<-as.numeric(rhc$cat1=='Colon Cancer')
Coma<-as.numeric(rhc$cat1=='Coma')
COPD<-as.numeric(rhc$cat1=='COPD')
lungcan<-as.numeric(rhc$cat1=='Lung Cancer')
MOSF<-as.numeric(rhc$cat1=='MOSF w/Malignancy')
sepsis<-as.numeric(rhc$cat1=='MOSF w/Sepsis')

female<-as.numeric(rhc$sex=='Female')
died<-as.numeric(rhc$death=='Yes')
age<-as.numeric(rhc$age)
treatment<-as.numeric(rhc$swang1=='RHC')
meanbp1<-rhc$meanbp1

#new dataset
mydata<-cbind(ARF,CHF,Cirr,colcan,Coma,lungcan,MOSF,sepsis, age,female,meanbp1,treatment,died)
mydata<-data.frame(mydata)

#covariates we will use 
xvars<-c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis", "age","female","meanbp1")

#look at a table 1
table1.raw<- CreateTableOne(vars=xvars, strata="treatment", data=mydata, test=FALSE)
## include standardized mean difference (SMD)
print(table1.raw,smd=TRUE)

############################################
# do greedy matching on Mahalanobis distance
############################################
greedymatch<-Match(Tr=treatment, M=1, X=mydata[xvars], replace=FALSE)
matched<-mydata[unlist(greedymatch[c("index.treated","index.control")]), ]
#get table 1 for matched data with standardized mean differences
table1.mtch.Mahal<-CreateTableOne(vars=xvars, strata ="treatment", 
                                  data=matched, test = FALSE)
print(table1.mtch.Mahal, smd = TRUE)

#outcome analysis
y_trt<-matched$died[matched$treatment==1]
y_con<-matched$died[matched$treatment==0]

#pairwise difference
diffy<-y_trt-y_con
#paired t-test
t.test(diffy)

###########################
# propensity score matching
###########################
#fit a propensity score model - logistic regression
psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+age+female+meanbp1,
      family=binomial(),data=mydata)
#show coefficients etc
summary(psmodel)

# #create propensity score
# pscore<-psmodel$fitted.values
# #do greedy matching on logit(PS) using Match with a caliper
# logit <- function(p) {log(p)-log(1-p)}
# psmatch<-Match(Tr=mydata$treatment, M=1, X=logit(pscore), replace=FALSE, caliper=.2)

psmatch<-Match(Tr=mydata$treatment, M=1, 
               X=psmodel$linear.predictors, replace=FALSE, caliper=.2)

matched<-mydata[unlist(psmatch[c("index.treated","index.control")]), ]

#get standardized differences
table1.mtch.pslp<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = FALSE)
print(table1.mtch.pslp, smd = TRUE)

# outcome analysis
y_trt<-matched$died[matched$treatment==1]
y_con<-matched$died[matched$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)
#t = 3.0396, df = 1931, p-value = 0.0046

###########################
# propensity score matching
###########################
psmodel<-glm(treatment~age+female+meanbp1+ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis,
      data=mydata, family = binomial(link ="logit"))

## value of propensity score for each subject
ps <-predict(psmodel, type = "response")

#create IP weights
weight<-ifelse(treatment==1,1/ps,1/(1-ps))

#apply weights to data
weighteddata<-svydesign(ids = ~ 1, data =mydata, weights = ~ weight)

#weighted table 
table1.wt.ps <-svyCreateTableOne(vars = xvars, strata="treatment",
                                  data=weighteddata, test=FALSE)
## Show table with SMD
print(table1.wt.ps, smd = TRUE)

#---get causal risk difference (RD)
glm.obj<-glm(died~treatment, weights=weight, 
             family= quasibinomial (link = "identity"))
#summary(glm.obj)
betaiptw<-coef(glm.obj)
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))
causalrd<-(betaiptw[2])
lcl<-(betaiptw[2]-1.96*SE[2])
ucl<-(betaiptw[2]+1.96*SE[2])
c(lcl,causalrd,ucl)

#---get causal relative risk(RR)
glm.obj<-glm(died~treatment, weights=weight, 
             family=quasibinomial ( link = "log"))

betaiptw<-coef(glm.obj)

#to properly account for weighting, use asymptotic (sandwich) variance
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))

#get point estimate and CI for relative risk (need to exponentiate)
causalrr<-exp(betaiptw[2])
lcl<-exp(betaiptw[2]-1.96*SE[2])
ucl<-exp(betaiptw[2]+1.96*SE[2])
c(lcl,causalrr,ucl)

#---truncate weights at 10
truncweight<-replace(weight, weight>10,10)
#get causal risk difference
glm.obj<-glm(died~treatment,weights=truncweight, family=
               quasibinomial(link="identity"))

#summary(glm.obj)
betaiptw<-coef(glm.obj)
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))
causalrd<-(betaiptw[2])
lcl<-(betaiptw[2]-1.96*SE[2])
ucl<-(betaiptw[2]+1.96*SE[2])
c(lcl,causalrd,ucl)
#0.055 (0.028, 0.082)

#############################
#alternative: use ipw package
#############################
#first fit propensity score model to get weights
weightmodel<-ipwpoint(exposure= treatment, family = "binomial", link ="logit",
  denominator=~age+female+meanbp1+ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis, 
  data=mydata)

#numeric summary of weights
summary(weightmodel$ipw.weights)

ipwplot(weights = weightmodel$ipw.weights, logscale = FALSE, 
        main ="weights", xlim = c(0, 22))

mydata$wt<-weightmodel$ipw.weights

#--Estimate risk difference
rd <- (svyglm(died ~ treatment, 
              design = svydesign(~1, weights = ~wt, data=mydata)))
coef(rd); confint(rd)
#0.052 (0.023, 0.080)

# fit propensity score model to get weights, but truncated
weightmodel<-ipwpoint(exposure= treatment, family = "binomial", link ="logit",
           denominator=~age+female+meanbp1+ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis, 
           data=mydata, trunc=.01)

#numeric summary of weights
summary(weightmodel$weights.trun)

#plot of weights
ipwplot(weights = weightmodel$weights.trun, logscale = FALSE, main =
          "weights", xlim = c(0, 22))

mydata$wt<-weightmodel$weights.trun
#--estimate risk difference with truncation
rd.trk <- (svyglm(died ~ treatment, 
                  design = svydesign(~ 1, weights = ~wt,data =mydata)))
coef(rd.trk); confint(rd.trk)

