################################
####Simulation of Population####
################################
set.seed(49)
nPOP<-1000000
tOR<-1.1

####Creating population####
pop <- data.frame(X1 = rbinom(nPOP, 1, .5), X2 = rbinom(nPOP, 1, .75),
  X3 = rbinom(nPOP, 1, .25), X4 = rbinom(nPOP, 1, .4), X5 = rbinom(nPOP, 1, .1),
  X6 = rbinom(nPOP, 1, .8), X7 = rbinom(nPOP, 1, .3), X8 = rbinom(nPOP, 1, .5))
pop = transform(pop,        #add A
  A = rbinom(nPOP, 1, 1/(1+exp(-(3*X1+X2+2*X3-2*X4-3*X5+3*X6-4*X7-4*X8)))))
pop = transform(pop,        #add Y
  Y = rbinom(nPOP, 1, 1/(1+exp(-(tOR*A+X1-1.5*X2-3*X3-X4-3*X5-3*X6+X7-3*X8)))))
sum(pop$Y)/nPOP             #probability of outcome

####################################################################################
####Simulation and Analyses of Nested Case-Control Samples Drawn from EHR Cohort####
####################################################################################
set.seed(75)
n<-100000
nS<-500                      #number of simulations
M<-3                         #number of controls matched to each case
nColY<-dim(pop)[2]           #number of columns and location of Y
nColX1<-which(colnames(pop)=="X1") #location of matching variable X1
OR<-matrix(NA, nrow=nS, ncol=3)
colnames(OR)<-c("Database","Unmatched", "Matched")
n0<-matrix(NA, nrow=nS, ncol=1)

for(p in 1:nS){
####Construct EHR cohort sample####
EHRcohort<-pop[(sample((1:nPOP),n, replace=FALSE)),]
nCC<-sum(EHRcohort$Y)        #number of cases
nCoC<-n-nCC                  #number of controls in EHR cohort
nCoM<-M*nCC                  #number of controls in sample
n0[p]<-nCC*(M+1)             #number of observations in study sample

####Subset the EHRvcohort into case and control groups by X1####
casesEHRcohort1<-EHRcohort[EHRcohort[,nColY]==1 & EHRcohort[,nColX1]==1,]
casesEHRcohort0<-EHRcohort[EHRcohort[,nColY]==1 & EHRcohort[,nColX1]==0,]
contEHRcohort1<-EHRcohort[EHRcohort[,nColY]==0 & EHRcohort[,nColX1]==1,]
contEHRcohort0<-EHRcohort[EHRcohort[,nColY]==0 & EHRcohort[,nColX1]==0,]
nC1<-dim(casesEHRcohort1)[1]
nC0<-dim(casesEHRcohort0)[1]
nCoC1<-dim(contEHRcohort1)[1]
nCoC0<-dim(contEHRcohort0)[1]

####Construct a matched case-control dataset####
contsample<-data.frame(matrix(NA, nrow=nCoM, ncol=nColY))
colnames(contsample)<-colnames(EHRcohort)
rstart<-1;rend<-M
for(i in 1:nC1){
  contsample[rstart:rend,]<-contEHRcohort1[(sample(1:nCoC1, M, replace=FALSE)),]
  rstart<-rend+1;rend<-rend+M}
for(i in 1:nC0){
  contsample[rstart:rend,]<-contEHRcohort0[(sample(1:nCoC0, M, replace=FALSE)),]
  rstart<-rend+1;rend<-rend+M}
ccDatamatched<-rbind(casesEHRcohort1,casesEHRcohort0,contsample)

####Construct an unmatched case-control dataset####
contEHRcohort<-rbind(contEHRcohort1,contEHRcohort0)
contsampleum<-contEHRcohort[(sample((1:nCoC),nCoM, replace=FALSE)),]
ccDataunmatched<-rbind(casesEHRcohort1,casesEHRcohort0,contsampleum)

####Regressions####
fit1<-glm(Y~A+X1+X2+X3+X4+X5+X6+X7+X8, family='binomial', data=EHRcohort)
fit2<-glm(Y~A+X1+X2+X3+X4+X5+X6+X7+X8, family='binomial', data=ccDataunmatched)
fit3<-glm(Y~A+X1+X2+X3+X4+X5+X6+X7+X8, family='binomial', data=ccDatamatched)
OR[p,1]<-coef(fit1)[2];OR[p,2]<-coef(fit2)[2];OR[p,3]<-coef(fit3)[2]
}

#########################
####Results Summaries####
#########################
mean(n0)
AverageBias<-apply(OR, 2, function (x) mean(abs(x-tOR)/tOR))
AverageBias
