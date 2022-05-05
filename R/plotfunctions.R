
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
  
  allinf = result.rep[,c((2*nages+2):((2+1)*nages+1),
                         (6*nages+2):((6+1)*nages+1),1,length(result.rep[1,]))]
  
  allinf$total = rowSums(allinf[1:(length(allinf[1,])-2)])
  
  return(allinf)
}
getallinf = function(result.rep,years)
{
  
  allinf = result.rep[,c((1*nages+2):((1+1)*nages+1),
                         (2*nages+2):((2+1)*nages+1),
                         (3*nages+2):((3+1)*nages+1),
                         (6*nages+2):((6+1)*nages+1),
                         (8*nages+2):((8+1)*nages+1),
                         (9*nages+2):((9+1)*nages+1),
                         1,length(result.rep[1,]))]
  allinf$total = rowSums(allinf[1:(length(allinf[1,])-2)])
  
  return(allinf)
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

