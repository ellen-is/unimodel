
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



getinfecteds = function(result.rep,years)
{
  Infecteds = result.rep[,c((3*nages+2):((3+1)*nages+1),1,length(result.rep[1,]))]
  Infecteds$year1 = rowSums(Infecteds[which(years=="1")])
  Infecteds$year2 = rowSums(Infecteds[which(years=="2")])
  Infecteds$year3 = rowSums(Infecteds[which(years=="3")])
  Infecteds$year4 = rowSums(Infecteds[which(years=="4")])
  Infecteds$yearR = rowSums(Infecteds[which(years=="R")])
  Infecteds$yearT = rowSums(Infecteds[which(years=="T")])
  Infecteds$total = rowSums(Infecteds[1:nages])
  return(Infecteds)
}
getasymp = function(result.rep,years)
{
  
  allinf = result.rep[,c((1*nages+2):((1+1)*nages+1),#E compartment
                         (2*nages+2):((2+1)*nages+1),#A compartment
                         (6*nages+2):((6+1)*nages+1), #P compartment
                         (8*nages+2):((8+1)*nages+1), #Q compartment
                         (9*nages+2):((9+1)*nages+1), #QA compartment
                         (11*nages+2):((11+1)*nages+1), #EV compartment
                         1,length(result.rep[1,]))] 
  
  allinf$total = rowSums(allinf[1:(length(allinf[1,])-2)])
  
  return(allinf)
}
getallinf = function(result.rep,years)
{
  #the result.rep columns are time, S, E, A, I, R, H, P, D, Q, QA, N, EV, nsim, scenario
  allinf = result.rep[,c((1*nages+2):((1+1)*nages+1), #E compartment
                         (2*nages+2):((2+1)*nages+1), #A compartment
                         (3*nages+2):((3+1)*nages+1), #I compartment
                         (6*nages+2):((6+1)*nages+1), #P compartment
                         (8*nages+2):((8+1)*nages+1), #Q compartment
                         (9*nages+2):((9+1)*nages+1), #QA compartment
                         (11*nages+2):((11+1)*nages+1), #Ev compartment
                         1,length(result.rep[1,]))] #time and nsim
  allinf$total = rowSums(allinf[1:(length(allinf[1,])-2)])
  
  return(allinf)
}

getrecoverds = function(result.rep,years)
{
  #the result.rep columns are time, S, E, A, I, R, H, P, D, Q, QA, N, nsim, scenario
  recoverds = result.rep[,c((4*nages+2):((4+1)*nages+1), #R compartment
                        1,length(result.rep[1,]))]
  recoverds$total = rowSums(recoverds[1:(length(recoverds[1,])-2)])
  
  return(allinf)
}


getiso = function(result.rep,years)
{
  #the result.rep columns are time, S, E, A, I, R, H, P, D, Q, QA, N, nsim, scenario
  isolated = result.rep[,c((8*nages+2):((8+1)*nages+1), #Q compartment
                         (9*nages+2):((9+1)*nages+1), #QA compartment
                         1,length(result.rep[1,]))]
  isolated$total = rowSums(isolated[1:(length(isolated[1,])-2)])
  
  return(isolated)
}

plotyears = function(meanout)
{
  studentcolsa=hcl.colors(10, "YlOrRd", rev = FALSE)
  studentcols=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.1)
  startofterm=as.Date("28/09/2020",format="%d/%m/%Y")
  Xmas=as.Date("21/12/2020",format="%d/%m/%Y")-startofterm
  days=1:length(meanout$year1)
  
  studentcolsa=hcl.colors(10, "YlOrRd", rev = FALSE)
  studentcols=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.1)
  
  plot(days+startofterm,meanout$year1,type="l",lwd=2,col=studentcolsa[1],
       xlab="",cex.axis=1.1,cex.lab=1.3,ylab="# symptomatic")
  
  lines(days+startofterm,meanout$year2,type="l",lwd=2,col=studentcolsa[2])
  lines(days+startofterm,meanout$year3,type="l",lwd=2,col=studentcolsa[3])
  lines(days+startofterm,meanout$year4,type="l",lwd=2,col=studentcolsa[4])
  lines(days+startofterm,meanout$yearT,type="l",lwd=2,col=studentcolsa[5])
  lines(days+startofterm,meanout$yearR,type="l",lwd=2,col=studentcolsa[6])
  legend('topright',c(paste('year',1:4),"PGT","PGR"),col=studentcolsa,lwd=3,cex=1.0,box.col=NA)
  
  abline(v=as.Date("21/12/2020",format="%d/%m/%Y"),lty=2)
}


plotcases = function(allinf,Infecteds,maxbaseline,stratname)
{
  studentcolsa=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.6)
  studentcols=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.1)
  startofterm=as.Date("20/09/2021",format="%d/%m/%Y")
  Xmas=as.Date("21/12/2020",format="%d/%m/%Y")-startofterm
  n=1
  days = 1:91
  allinf
  #days=1:length(allinf$total[allinf$nsim==n])
  plot(startofterm+days,allinf$total[allinf$nsim==n],type="l",lwd=2,col=rgb(0,0.5,0,0.1),
  #plot(startofterm+days,allinf$total[allinf$nsim==n]+Infecteds$total[Infecteds$nsim==n],type="l",lwd=1,col= "tan1",
      xlab="Date",cex.axis=1.1,cex.lab=1.3,
       #ylim=c(1,1.0*maxbaseline)log="")
       ylab = "Number of infected students", ylim=c(1,15000),log="")
       
  
  for(n in 1:nsim)
  {
    #lines(startofterm+days,Infecteds$total[Infecteds$nsim==n],type="l",lwd=2,col=studentcols[1])
    lines(startofterm+days,Infecteds$total[Infecteds$nsim==n],type="l",lwd=1,col="tan1")
    
    #lines(startofterm+days,allinf$total[allinf$nsim==n],type="l",lwd=2,col=rgb(0,0.5,0,0.1))
    lines(startofterm+days,allinf$total[allinf$nsim==n],type="l",lwd=1,col="royalblue")
    
  }
  #abline(a= 3, b = 0) #was using to try to calcualte the outbreak length
  #abline(v=as.Date("21/12/2020",format="%d/%m/%Y"),lty=2)
  #legend('topright',c("symptomatic","asymptomatic"),col=c(studentcolsa[1],rgb(0,0.5,0,0.5)),pch=19,box.col = NA)
}

plotcases_percent = function(allinf,Infecteds,maxbaseline,stratname, N)
{
  studentcolsa=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.6)
  studentcols=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.1)
  startofterm=as.Date("20/09/2021",format="%d/%m/%Y")
  Xmas=as.Date("21/12/2020",format="%d/%m/%Y")-startofterm
  n=1
  days = 1:91
  allinf
  #days=1:length(allinf$total[allinf$nsim==n])
  plot(startofterm+days,(allinf$total[allinf$nsim==n]/sum(N))*100,type="l",lwd=2,col=rgb(0,0.5,0,0.1),
       xlab="",cex.axis=1.1,cex.lab=1.3,
       #ylim=c(1,1.0*maxbaseline),ylab="",log="")
       ylim=c(0,5),ylab="",log="")
  
  
  for(n in 1:nsim)
  {
    lines(startofterm+days,(Infecteds$total[Infecteds$nsim==n]/sum(N))*100,type="l",lwd=2,col=studentcols[1])
    
    lines(startofterm+days,(allinf$total[allinf$nsim==n]/sum(N))*100,type="l",lwd=2,col=rgb(0,0.5,0,0.1))
  }
  #abline(v=as.Date("21/12/2020",format="%d/%m/%Y"),lty=2)
  #legend('topright',c("symptomatic","asymptomatic"),col=c(studentcolsa[1],rgb(0,0.5,0,0.5)),pch=19,box.col = NA)
}

plotcasesSA = function(allinf,Infecteds,maxbaseline,stratname)
{
  studentcolsa=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.6)
  studentcols=hcl.colors(10, "YlOrRd", rev = FALSE,alpha=0.1)
  startofterm=as.Date("28/09/2020",format="%d/%m/%Y")
  Xmas=as.Date("21/12/2020",format="%d/%m/%Y")-startofterm
  n=1
  days=1:length(allinf$total[allinf$nsim==n])
  plot(startofterm+days,allinf$total[allinf$nsim==n],type="l",lwd=2,col=rgb(0,0.5,0,0.1),
       xlab="",cex.axis=1.1,cex.lab=1.3,
       ylim=c(1,1.0*maxbaseline),ylab="# infected",log="")
  
  for(n in 1:nsim)
  {
    lines(startofterm+days,Infecteds$total[Infecteds$nsim==n],type="l",lwd=2,col=studentcols[1])
    
    lines(startofterm+days,allinf$total[allinf$nsim==n],type="l",lwd=2,col=rgb(0,0.5,0,0.1))
  }
  abline(v=as.Date("21/12/2020",format="%d/%m/%Y"),lty=2)
  legend('topright',c("symptomatic","asymptomatic"),col=c(studentcolsa[1],rgb(0,0.5,0,0.5)),pch=19,box.col = NA)
}
  
  

  
  
  
 