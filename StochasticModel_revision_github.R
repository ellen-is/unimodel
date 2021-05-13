library(Rcpp)
library(tidyverse)
library(gplots)
library(plot3D)
library(socialmixr)


####this code creates figure 2 from High COVID-19 transmission potential associated with re-opening universities can be mitigated with layered interventions from Brooks-Pollock et al. 

#change to local directory
create_mixingmatrix = function(cs=24,meanstudy=20,meanextra=4)
{
  #cs = living circle size
  load(file=paste("ContactMatrix",cs,".Rdata",sep=""))
  load(file=paste("Npop",cs,".Rdata",sep=""))
  numgroups = dim(Contactmatrix)[1]
  extra1=rnorm(1,mean=meanextra,sd=1)
  while(extra1<0){extra1=rnorm(1,mean=meanextra,sd=1)}
  extracontacts = matrix(extra1*Npop/sum(Npop),nrow=numgroups,ncol=numgroups,byrow=T)
  betamatrix=Contactmatrix + extracontacts
  studycontacts=rnorm(1,mean=meanstudy,sd=4)
  while(studycontacts<0){studycontacts=rnorm(1,mean=meanstudy,sd=4)}
  
  diag(betamatrix) = diag(betamatrix) + studycontacts
  return(betamatrix)
}
normaliseNGMbaseline = function(NGM,rzero)
{
  #normalised to living circle size = 24
  mev=max(Re(eigen(NGM)$value))
  multiplier = rzero/mev
  return(multiplier)
}
normaliseNGM = function(NGM,multiplier)
{
  #normalised to living circle size = 24
  mev=max(Re(eigen(NGM)$value))
  betamat = multiplier*NGM
  return(betamat)
}

params <- list()
params <- within(params, {
  
  ## set rng state
  seed <- 0

  SIMTIME=300
  iperiod=5.0; 
  load("Npopsize.Rdata")
  years=names(Npop)
  nages=length(Npop)
  meanextra=4
  meanstudy=20
  cs=24
  NGM = create_mixingmatrix() 
  rzero=2.26
  fac1=1
  multiplier=normaliseNGMbaseline(NGM,rzero)
  betamat = multiplier*NGM
  gam_p=1/2;
  beta=betamat/(iperiod+1/gam_p)
  
  gamma=1/iperiod; 
  sigma=1/(3.2); 
  gam_a=1/(iperiod+1/gam_p); 
  gam_h=1/3
 
  gam_q=1/14
  testrate=1/2;
  #testrate=0;
  testrate_a=rep(0,nages) #reactive testing of asymptomatic cases
  testrate_a2=0;
  testrate_a3=0;
  eps=0.5
  eps_q=0.5
  h=0.002
  f=0.75 #fraction asymptomatic
  mrate=0.038
  backgroundrate = 1e-4
  #backgroundrate = (rzero/mev)/100
  
  ## initial conditions, list within list
  ## use within() to modify empty list, as above
  
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
    N=S+E+A+I+H+R+P+Q+QA
  })
})
studentcolsa=hcl.colors(10, "YlOrRd", rev = FALSE)
studentcols=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.1)
set.seed(params$seed)

sourceCpp('studentmodel_reactive2.cpp')

nages=params$nages
years=params$years
startofterm=as.Date("28/09/2020",format="%d/%m/%Y")
Xmas=as.Date("21/12/2020",format="%d/%m/%Y")-startofterm

getdoubling = function(result.rep,statenums)
{
  s1=result.rep[,c('time','nsim')]
  for(i in statenums)
  {
    s1 = cbind(s1,result.rep[,(i*nages+2):((i+1)*nages+1)])
  }
  s1$totalcount = rowSums(s1[,3:length(s1[1,])])
  
  growthrate=rep(0,max(s1$nsim))
  PT=rep(0,max(s1$nsim))
  for(n1 in 1:max(s1$nsim))
  {
    yall=log(s1$totalcount[s1$nsim==n1])
    xmin=1
    peaktime=xmin+20
    y=yall[xmin:peaktime]
    x=xmin:peaktime
    m1=lm(y~x)
    growthrate[n1]=m1$coefficients[2]
    PT[n1]=which.max(s1$totalcount[s1$nsim==n1])
  }
  timestats=c(mean(growthrate),min(growthrate),max(growthrate),
              mean(PT),min(PT),max(PT))
  return(timestats)
}

getstate = function(result.rep,statenums)
{
  s1=result.rep[,c('time','nsim')]
  for(i in statenums)
  {
    s1 = cbind(s1,result.rep[,(i*nages+2):((i+1)*nages+1)])
  }
  s1$totalcount = rowSums(s1[,3:length(s1[1,])])
  s1$year1 = rowSums(s1[c(F,F,rep(years=="1",length(statenums)))])
  s1$year2 = rowSums(s1[c(F,F,rep(years=="2",length(statenums)))])
  s1$year3 = rowSums(s1[c(F,F,rep(years=="3",length(statenums)))])
  s1$year4 = rowSums(s1[c(F,F,rep(years=="4",length(statenums)))])
  s1$yearR = rowSums(s1[c(F,F,rep(years=="R",length(statenums)))])
  s1$yearT = rowSums(s1[c(F,F,rep(years=="T",length(statenums)))])
  
  s1 %>%
    group_by(time) %>%
    summarise(total=mean(totalcount),min=min(totalcount),max=max(totalcount),
              my1=mean(year1),my2=mean(year2),my3=mean(year3),my4=mean(year4),
              myr=mean(yearR),myt=mean(yearT)) -> meanout
  
  return(meanout)
}
getstateyears = function(result.rep,years,statenums)
{
  s1=result.rep[,c('time','nsim')]
  for(i in statenums)
  {
    s1y=result.rep[,(i*nages+2):((i+1)*nages+1)]
    s1 = cbind(s1,s1y)
    
  }
  s1$total = rowSums(s1[,3:length(s1[1,])])
  s1$year1 = rowSums(s1[c(F,F,rep(years=="1",length(statenums)))])
  s1$year2 = rowSums(s1[c(F,F,rep(years=="2",length(statenums)))])
  s1$year3 = rowSums(s1[c(F,F,rep(years=="3",length(statenums)))])
  s1$year4 = rowSums(s1[c(F,F,rep(years=="4",length(statenums)))])
  s1$yearR = rowSums(s1[c(F,F,rep(years=="R",length(statenums)))])
  s1$yearT = rowSums(s1[c(F,F,rep(years=="T",length(statenums)))])
  return(s1)
}



plotyears = function(meanout,plotlegend=F)
{
  studentcolsa=hcl.colors(10, "YlOrRd", rev = FALSE)
  studentcols=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.1)
  startofterm=as.Date("28/09/2020",format="%d/%m/%Y")
  Xmas=as.Date("21/12/2020",format="%d/%m/%Y")-startofterm
  days=1:length(meanout$my1)
  
  studentcolsa=hcl.colors(10, "YlOrRd", rev = FALSE)
  studentcols=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.1)
  
  plot(days+startofterm,meanout$my1,type="l",lwd=2,col=studentcolsa[1],
       xlab="",cex.axis=1.1,cex.lab=1.3,ylim=c(0,max(meanout$my1)),ylab="# symptomatic")
  
  lines(days+startofterm,meanout$my2,type="l",lwd=2,col=studentcolsa[2])
  lines(days+startofterm,meanout$my3,type="l",lwd=2,col=studentcolsa[3])
  lines(days+startofterm,meanout$my4,type="l",lwd=2,col=studentcolsa[4])
  lines(days+startofterm,meanout$myt,type="l",lwd=2,col=studentcolsa[5])
  lines(days+startofterm,meanout$myr,type="l",lwd=2,col=studentcolsa[6])
  
  if(plotlegend==T)legend('topright',c(paste('year',1:4),"PGT","PGR"),col=studentcolsa,pch=19,cex=1,box.col=NA,bg=NA)
  
  abline(v=as.Date("21/12/2020",format="%d/%m/%Y"),lty=2)
}


plotcases = function(allinf,Infecteds,maxbaseline,stratname,plotlegend=F)
{
  studentcolsa=hcl.colors(10, "YlOrRd", rev = FALSE)
  studentcols=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.1)
  startofterm=as.Date("28/09/2020",format="%d/%m/%Y")
  Xmas=as.Date("21/12/2020",format="%d/%m/%Y")-startofterm
  n=1
  days=1:length(allinf$total[allinf$nsim==n])
  plot(startofterm+days,allinf$total[allinf$nsim==n],type="l",lwd=2,col=studentcols[2],
       xlab="",cex.axis=1.1,cex.lab=1.3,
       ylim=c(1,1.0*maxbaseline),ylab="# infected",log="")
  
  
  lines(startofterm+days,Infecteds$total[Infecteds$nsim==n],type="l",lwd=2,col=studentcols[1])
  
  for(n in 1:nsim)
  {
    lines(startofterm+days,Infecteds$total[Infecteds$nsim==n],type="l",lwd=2,col=studentcols[1])
    
    lines(startofterm+days,allinf$total[allinf$nsim==n],type="l",lwd=2,col=rgb(0,0.5,0,0.1))
  }
  
  #lines(startofterm+days,meanout$total,type="l",lwd=2,col="black")
  abline(v=as.Date("21/12/2020",format="%d/%m/%Y"))
  legend('topright',stratname,box.col = NA)
  #mtext(stratname,side=4)
}
plotallcases = function(allinf,maxbaseline,mycolour)
{
  startofterm=as.Date("28/09/2020",format="%d/%m/%Y")
  Xmas=as.Date("21/12/2020",format="%d/%m/%Y")-startofterm
  n=1
  days=1:length(allinf$total[allinf$nsim==n])
  plot(startofterm+days,allinf$total[allinf$nsim==n],type="l",lwd=2,col=mycolour[1],
       xlab="",cex.axis=1.1,cex.lab=1.3,
       ylim=c(1,1.0*maxbaseline),ylab="# infected",log="")
  
  for(n in 1:nsim)
  {
    lines(startofterm+days,allinf$total[allinf$nsim==n],type="l",lwd=2,col=mycolour[1])
  }
  
}

plotallcases_add = function(allinf,mycolour)
{
  startofterm=as.Date("28/09/2020",format="%d/%m/%Y")
  n=1
  days=1:length(allinf$total[allinf$nsim==n])
  for(n in 1:nsim)
  {
    lines(startofterm+days,allinf$total[allinf$nsim==n],type="l",lwd=2,col=mycolour[1])
  }
}



runplotmodel1 = function(params,stratname,maxbaseline=800)
{
  #jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  mycols=jet.col(6,alpha=0.3)
  result.rep=c()
  for(nn in 1:nsim)
  {
    set.seed(nn)
    NGM = create_mixingmatrix(cs=params$cs,meanstudy = params$meanstudy,meanextra = params$meanstudy)
    NGM=normaliseNGM(NGM,params$multiplier)
    params$beta=params$fac1*NGM/(params$iperiod+(1/params$gam_p))
    result <- c19uni(params)
    result$nsim <- nn
    result$time[params$SIMTIME]=(params$SIMTIME-1)
    result.rep=rbind(result.rep,result)
    
  }
  #0=S, 1=E, 2=A, 3=I, 4=R, 5=H, 6=P, 7=D, 8=Q, 9=QA, 10=N
  
  Infecteds=getstateyears(result.rep,years,3)
  Asymp=getstateyears(result.rep,years,c(2,6))
  allinf=getstateyears(result.rep,years,c(1,2,3,6))
  symp=getstate(result.rep,c(3))
  asymp1=getstate(result.rep,c(2,6))
  #plotcases(Asymp,Infecteds,maxbaseline,stratname)
  maxbaseline=max(allinf$total)
  plotallcases(allinf,maxbaseline,mycols[1])
  abline(v=as.Date("21/12/2020",format="%d/%m/%Y"),lty=2)
  legend('topright',stratname,box.col = NA,inset=c(0,0.0),bg=NA,pch=19,col=mycols[1],pt.cex=1.5)
  legend('topleft','a',box.col = NA,inset=c(-0.1,-0.1),cex=2,bg=NA)
  plotyears(symp,plotlegend = T)
  legend('topleft','b',box.col=NA,inset=c(-0.1,-0.1),cex=2,bg=NA)
  iw=which.max(symp$total)
  asymp1=getstate(result.rep,c(2,6))
  modeltimes=getdoubling(result.rep,c(2,3,6))

  out1=c(as.numeric(symp[84,2:4]),as.numeric(asymp1[84,2:4]),as.numeric(symp[iw,1:4]),
         summary(asymp1$total[50:300]/symp$total[50:300])[c(1,3,4,6)],
      modeltimes)
  return(out1)
}
runplotmodel3 = function(params,params2,params3,stratname,maxbaseline=800,mycols1)
{
  #jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),alpha=T)
  mycols=jet.col(6,alpha=0.3)
  result.rep=c()
  result.rep2=c()
  result.rep3=c()
  for(nn in 1:nsim)
  {
    set.seed(nn)
    NGM = create_mixingmatrix(cs=params$cs,meanstudy = params$meanstudy,meanextra = params$meanstudy)
    NGM=normaliseNGM(NGM,params$multiplier)
    params$beta=params$fac1*NGM/(params$iperiod+(1/params$gam_p))
    
    result <- c19uni(params)
    result$nsim <- nn
    result$time[params$SIMTIME]=(params$SIMTIME-1)
    result.rep=rbind(result.rep,result)
    
    
    NGM = create_mixingmatrix(cs=params2$cs,meanstudy = params2$meanstudy,meanextra = params2$meanstudy)
    NGM=normaliseNGM(NGM,params2$multiplier)
    params2$beta=params$fac1*NGM/(params2$iperiod+(1/params2$gam_p))
    
    result2 <- c19uni(params2)
    result2$nsim <- nn
    result2$time[params$SIMTIME]=(params2$SIMTIME-1)
    result.rep2=rbind(result.rep2,result2)
    
    
    NGM = create_mixingmatrix(cs=params3$cs,meanstudy = params3$meanstudy,meanextra = params3$meanstudy)
    NGM=normaliseNGM(NGM,params3$multiplier)
    params3$beta=params3$fac1*NGM/(params3$iperiod+(1/params3$gam_p))
    
    result3 <- c19uni(params3)
    result3$nsim <- nn
    result3$time[params$SIMTIME]=(params3$SIMTIME-1)
    result.rep3=rbind(result.rep3,result3)
    
  }
  #0=S, 1=E, 2=A, 3=I, 4=R, 5=H, 6=P, 7=D, 8=Q, 9=QA, 10=N
  
  Infecteds=getstateyears(result.rep,years,3)
  Asymp=getstateyears(result.rep,years,c(2,6))
  allinf=getstateyears(result.rep,years,c(1,2,3,6))
  allinf2=getstateyears(result.rep2,years,c(1,2,3,6))
  allinf3=getstateyears(result.rep3,years,c(1,2,3,6))
  totalR1=getstate(result.rep,c(1:9))
  totalR2=getstate(result.rep2,c(1:9))
  totalR3=getstate(result.rep3,c(1:9))
  symp=getstate(result.rep2,c(3))
  
  #plotcases(Asymp,Infecteds,maxbaseline,stratname)
  maxbaseline=max(allinf$total)
  plotallcases(allinf,maxbaseline,mycols[1])
  plotallcases_add(allinf2,mycols1[2])
  plotallcases_add(allinf3,mycols1[3])
  abline(v=as.Date("21/12/2020",format="%d/%m/%Y"),lty=2)
  legend('topright',c('baseline',stratname),box.col = NA,inset=c(0,0),col=c(mycols[1],mycols1[2:3]),pch=19,bg=NA,pt.cex=1.5)
  #plotyears(symp)
  
  symp1=getstate(result.rep,c(3))
  iw=which.max(symp1$total)
  asymp1=getstate(result.rep,c(2,6))
  modeltimes=getdoubling(result.rep,c(2,3,6))
  out1=c(as.numeric(symp1[84,2:4]),as.numeric(asymp1[84,2:4]),as.numeric(symp1[iw,1:4]),summary(asymp1$total[50:300]/symp1$total[50:300])[c(1,3,4,6)],
         modeltimes,100*totalR1[84,2:4]/sum(params$Npop),100*totalR1[300,2:4]/sum(params$Npop))
  
  symp2=getstate(result.rep2,c(3))
  iw=which.max(symp2$total)
  asymp2=getstate(result.rep2,c(2,6))
  modeltimes=getdoubling(result.rep2,c(2,3,6))
  out2=c(as.numeric(symp2[84,2:4]),as.numeric(asymp2[84,2:4]),as.numeric(symp2[iw,1:4]),summary(asymp2$total[50:300]/symp2$total[50:300])[c(1,3,4,6)],
         modeltimes,100*totalR2[84,2:4]/sum(params$Npop),100*totalR2[300,2:4]/sum(params$Npop))
  
  symp3=getstate(result.rep3,c(3))
  iw=which.max(symp3$total)
  asymp3=getstate(result.rep3,c(2,6))
  modeltimes=getdoubling(result.rep3,c(2,3,6))
  out3=c(as.numeric(symp3[84,2:4]),
         as.numeric(asymp3[84,2:4]),
         as.numeric(symp3[iw,1:4]),
         summary(asymp3$total[50:300]/symp3$total[50:300])[c(1,3,4,6)],
         modeltimes,100*totalR3[84,2:4]/sum(params$Npop),100*totalR3[300,2:4]/sum(params$Npop))
  
  outall=rbind(out1,out2,out3)
  
  return(outall)
}


nsim=10

mylabels=c("CS 25%","CS 50%",
           "f2f teaching 15","f2f teaching 5",
          "Living circles 20","Living circles 14",
           "Mass testing 1 week","Mass testing 2 days")
summaryout = matrix(0,nrow=params$SIMTIME,ncol=3*4)



params$eps = 0.5
par(mar=c(2,5,1,1))
layout(matrix(c(1,2,3,4,5,6,7,7,7,7), 5, 2, byrow = TRUE))
####1 baseline
fac1=1
if(params$eps == 0.1){params$fac1=0.77}
if(params$eps == 0.2){params$fac1=0.9}
if(params$eps == 0.3){params$fac1=1}
if(params$eps == 0.4){params$fac1=1.097}
if(params$eps == 0.5){params$fac1=1.18}
if(params$eps == 0.6){params$fac1=1.26}
if(params$eps == 0.7){params$fac1=1.34}
if(params$eps == 0.8){params$fac1=1.4}
if(params$eps == 0.9){params$fac1=1.46}
if(params$eps == 1){params$fac1=1.52}

mycols1=rainbow(8,alpha=0.3)
params$beta=fac1*params$beta

out1=runplotmodel1(params,"baseline",maxbaseline=maxY)
log(2)/out1[15:17]

#####2. covid security
params2=params
params2$meanstudy=15
params2$meanextra=3

params3=params
params3$meanstudy=10
params3$meanextra=2

OUT=runplotmodel3(params,params2,params3,mylabels[1:2],maxbaseline=maxY,mycols1=mycols1[c(1,1,2)])
legend('topleft','c',box.col = NA,inset=c(-0.1,-0.1),cex=2,bg=NA)

#####3. reduced contact
params2=params
params2$meanstudy=15
params2$meanextra=4

params3=params
params3$meanstudy=5
params3$meanextra=4

OUT=rbind(OUT,runplotmodel3(params,params2,params3,mylabels[3:4],maxbaseline=maxY,mycols1=mycols1[c(1,3,4)]))
legend('topleft','d',box.col = NA,inset=c(-0.1,-0.1),cex=2,bg=NA)

#####4. reduced living circles
params2=params
params2$cs=20

params3=params
params3$cs=14

OUT=rbind(OUT,runplotmodel3(params,params2,params3,mylabels[5:6],maxbaseline=maxY,mycols=mycols1[c(1,5,6)]))
legend('topleft','e',box.col = NA,inset=c(-0.1,-0.1),cex=2,bg=NA)

#####5. testing
params2=params
params2['testrate_a2']=1/7
params3=params
params3['testrate_a2']=1/2
OUT=rbind(OUT,runplotmodel3(params,params2,params3,mylabels[7:8],maxbaseline=maxY,mycols=mycols1[c(1,7,8)]))
legend('topleft','f',box.col = NA,inset=c(-0.1,-0.1),cex=2,bg=NA)

colnames(OUT)[1:3]=c("mean_sympxmas","min","max")
colnames(OUT)[4:6]=c("mean_Axmas","min","max")
colnames(OUT)[7:10]=c("peaknumber_time","mean","min","max")
colnames(OUT)[11:14]=c("s/a_min","s/a_med","s/a_mean","s/a_max")
colnames(OUT)[15:17]=c("growthrate","min","max")
colnames(OUT)[18:20]=c("peaktime","min","max")
colnames(OUT)[21:23]=c("R_xmas","min","max")
colnames(OUT)[24:26]=c("R_end","min","max")


rank1=read.csv("ranking.csv")

C1=contact_matrix(polymod, countries = "United Kingdom", age.limits = c(0,5,18, 30, 50, 70),quiet=T)
betamatrix=C1$matrix
X=rowSums(betamatrix)
eps=seq(0.1,1,0.1)
rstudents=eps
r0=2.7
f1=0.6
f2=0.25
for(i in 1:10)
{
  
  rs=r0/(f1+(1-f1)*eps[i])
  rstudents[i]=rs*(f2+(1-f2)*eps[i])*(X[3]^2/(mean(X)^2))
}
rank2=rank1
rank1[2,]=rank2[1,]#I got them the wrong way round
rank1[1,]=rank2[2,]
rank1[5,]=rank2[6,]
rank1[6,]=rank2[5,]


par(mar=c(5.5,9.5,5.5,3.5))
mycols2=rainbow(8,alpha=1)
order1=c(4,2,3,8,1,7,6,5)

plot(eps,as.numeric(rank1[4,2:11]),type="o",ylim=c(0,9),col=mycols2[4],yaxt="n",ylab="",xlab="")
mtext(expression(paste("asymptomatic infectiousness, ",epsilon)),side=1,line=2.5,cex=1.1)
for(i in 1:8)
{
  points(eps,as.numeric(rank1[i,2:11]),type="o",col=mycols2[i],lwd=2)
  
}

axis(2,at=1:8,labels=c("f2f teaching 5","COVID security 50%","f2f teaching 15","Mass testing 2 days","COVID security 25%","Mass testing 1 week","Living circles 14","Living circles 20"),las=2)
axis(4,at=1:8,labels=1:8,las=2)
axis(3,at=eps,labels=signif(rstudents,2))
mtext("Reproduction number",side=3,line=2.5,cex=1.1)

mtext("g",side=2,at=9.5,las=2,cex=1.5,line=2)

corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
par(xpd = TRUE) #Draw outside plot area
text(x = 1.1, y = 5, "Rank", srt = 270,cex=1.5)

