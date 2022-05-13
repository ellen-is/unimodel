


library(Rcpp)
library(tidyverse)
library(gplots)
source("R/plotfunctions.R") 

starttime=Sys.time()


params <- list()
params <- within(params, {

  SIMTIME=300
  iperiod=5.0; 
  load("data/betamat.Rdata")
  load("data/Npop.Rdata")
  years=read.csv("data/yearnumbers1.csv")
  nages=length(Npop)
  
  gam_p=1/2;
  beta=betamat/(iperiod+1/gam_p)
  gamma=1/iperiod;  sigma=1/(3.2);  gam_a=1/(iperiod+1/gam_p); gam_h=1/3; gam_q=1/14 ##recovery rates
  testrate=1/2; testrate_a=rep(0,nages); testrate_a2=0; testrate_a3=0; ##testing parameters
  eps=0.3; eps_q=0.5 ##
  h=0.002 #hospitalisation rate
  f=0.75 #fraction asymptomatic
  mrate=0.038 #mortality rate
  backgroundrate = 1e-4 #background infection rate
  waning = 0
  ve_inf = 0.5
  ve_trans = 0.1
    
  ## initial conditions, list within list
  init <- within(list(), {
    S=Npop
    E=rep(0,nages)
    E[sample(1:nages,10)]=1
    A=rep(0,nages)
    A[sample(1:nages,10)]=1
    I=rep(0,nages)
    I[sample(1:nages,10)]=1
    H=rep(0,nages)
    R=rep(0,nages)
    R[sample(1:nages,10)]=1
    P=rep(0,nages)
    D=rep(0,nages)
    Hinc=rep(0,nages)
    Q=rep(0,nages)
    QA=rep(0,nages)
    EV=rep(0,nages)
    N=S+E+A+I+H+R+P+Q+QA+EV
  })
})
attach(params)
set.seed(params$seed)
sourceCpp('R/studentmodel_reactive2.cpp')

nsim=10
par(mar=c(5,5,1,1),mfrow=c(1,2))

result.rep=c()
for(nn in 1:nsim)
{
  set.seed(nn)
  result <- c19uni(params)
  result$nsim <- nn
  result$time[params$SIMTIME]=(params$SIMTIME-1)
  result.rep=rbind(result.rep,result)
}


Infecteds=getinfecteds(result.rep,years)
allinf=getasymp(result.rep,years)

Infecteds %>%
  group_by(time) %>%
  summarise_all(mean) -> meanout

maxbaseline=max(allinf$total)*1.2
plotcases(allinf,Infecteds,maxbaseline,'baseline')
plotyears(meanout)


end.time = Sys.time()

print(end.time-starttime)
