library(Rcpp)
library(tidyverse)
library(gplots)
library(ggpubr)
source("R/plotfunctions_withgetstate.R") 

load("data/betamat.RData")
load("data/Npop.Rdata")
years=read.csv("data/yearnumbers1.csv")
nages=length(Npop)
############################################################################
###Setting different levels of asymptomatic testing and vaccination uptake ##
############################################################################
#These are currently set as for Figure 2 in the submitted paper
#The way I've made the lists makes sure there are all combinations of the three parameter sets

#Vaccine uptake informs the proportion of students who are in the R compartment at the beginning of the simulation

vaccine_uptake_list <- c(rep(rep(c(0.3, 0.5, 0.7, 0.9), 4),1))

#Asymptomatic testing uptake 
a_test_uptake_list <-  c((rep(c(rep(0.9,4),rep(0.5,4),rep(0.2,4), rep(0,4)),1)))


#Range of percentage of student population initially infected
#Low- 0.05%
#Medium- 0.1%
#High- 0.2% 
#infected <- c(rep(0.005,16), rep(0.01,16), rep(0.02,16))
#Used the average for the paper results

infected <- c(rep(0.01,16))




############################################################################
###Code currently commented out but can be used to save data to make other figures ##
############################################################################

#make data frames for saving data so I can make other plots (e.g Figure 3), not currently in use while I am doing exploratory analysis
#Ntot <- data.frame()
#cumulative_incidence <- data.frame() #this will record the cumulative incidence of each scenario
#Isolated_allsims <- data.frame() #this will calcualte the number in isolation in all scenarios
#Infected_allsims <- data.frame()
#allinf_total <- data.frame()
#Cumulative_incidence_total <- data.frame()
#Recoverds_total <- data.frame()
#Isolated_total <- data.frame()
#this is the start of a for loop which can be used to do sensitivity analysis with the different parameter sets


#############################################################################################
###Running the model under different levels of vaccination and asymptomatic testing uptake##
############################################################################################

#Setting layout for figure

par(mar=c(2, 3, 2, 2),mfrow=c(4,4))



for (i in 1:length(vaccine_uptake_list)) {
  
  starttime=Sys.time()
  params <- list()
  params <- within(params, {
    
    SIMTIME=91
    iperiod=5.0; 
     load("data/betamat.RData")
     load("data/Npop.Rdata")
     years=read.csv("data/yearnumbers1.csv")
     nages=length(Npop)
    
    gam_p=1/2;
    vaccineuptake <-vaccine_uptake_list[i]
    a_test_uptake <- a_test_uptake_list[i]
    ve_inf <- 0.45  #vaccine effectiveness against infection
    #ve_trans <- 0.4 #proportion reduction in the probability of symptoms.
    ve_trans <- 0.16 #proportion reduction in the probability of symptoms.
    beta=(betamat/(iperiod+1/gam_p))*1.5*1.6 #Setting this for delta, hence multiplying by 1.5 and 1.6 (the increase in transmisbility from WT to alpha and then alpha to delta respectively)
    gamma=1/iperiod;  sigma=1/(3.2);  gam_a=1/(iperiod+1/gam_p); gam_h=1/3; gam_q=1/10 ##recovery rates, not we are now assuming 10 days are spent in isolation (this is what we assumed in the vaccination paper but we assumed 14 in the original paper)
    daysto_asy_test <- 1/4
    LFTsensitivity <- 0.58 #(for asymptomatic people according to this Cochrane report https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.CD013705.pub2/full)
    testrate_a_single_value <- a_test_uptake * daysto_asy_test * LFTsensitivity #higher sensitivity reduces the time to test
    testrate=1/2; testrate_a=rep(testrate_a_single_value,nages); testrate_a2=0; testrate_a3=0; ##testing parameters, we aren't using testrate_a2 or testrate_a3 anymore
    eps=0.5; eps_q=0.5 ##
    h = 0 #hospitalisation switched off
    
    #I think f is already being scaled in the c code so I have commented out the ways this was done in the vaccination paper and left it as it was in the original paper
    f = 0.75
    #vaccine_symptoms <- 0.84 #2nd dose pfzier 84-95% We had this in the 
    #f=0.75/(1+vaccine_uptake*vaccine_symptoms) #((0.25 * vaccine_symptoms)*vaccine_uptake) #fraction asymptomatic - increase this because of vaccination
    #f = 0.75*(1-vaccineuptake) + 0.75*vaccineuptake/ve_trans
    #f = 0.75*(1-vaccine_uptake) + 0.75*vaccine_uptake/vaccine_symptoms
    #mrate=0.038 * 0.11 #mortality rate adjusted according to vaccine efficacy
    mrate = 0 #mortality rate switched off
    backgroundrate = 1e-3 #background infection rate/seasonal importation from the community
    waning <- 1/243
    
    

    ### initial conditions, list within list ##
    
    #These have been set up so that they are allocated across all school/year groups as evenly as possible 
    
    ## S = Number of students minus number in all other compartments
    ## E = 50% of all infections
    ## A = 40% of all infections
    ## P = 5% of all infections
    ## I = 5% of all infections
    ## Q = 0
    ## R = the percentage of students initially vaccinated at the start of the simulation
    ## EV = 0
    
    
    init <- within(list(), {
      S=Npop
      N= sum(Npop)
      E=rep(0,nages)
      A=rep(0,nages)
      I=rep(0,nages)
      P=rep(0,nages)
      
      if((infected[i]*N/2) <= 153) {
        E[sample(1:nages,round(((infected[i]*N))*0.5))]=1 
        A[sample(1:nages,round(((infected[i]*N))*0.4))]=1
        I[sample(1:nages,round(((infected[i]*N))*0.05))]=1
        P[sample(1:nages,round(((infected[i]*N))*0.05))] = 1
      } else {
        
        E[sample(1:nages,153)]=round(((infected[i]*N)*0.5)/153)
        A[sample(1:nages,153)]=round(((infected[i]*N)*0.4)/153)
        I[sample(1:nages,153)]=round(((infected[i]*N)*0.05)/153)
        P[sample(1:nages,153)]=round(((infected[i]*N)*0.05)/153)
      }
      
      H=rep(0,nages)
      R=rep(0,nages)
      
      #this is adjusted based on the vaccine uptake and divided evenly between groups
        if(vaccineuptake*N <= 153) {
          
      R[sample(1:nages,round(vaccineuptake*N))]=1 
      
      
      } else {
       R[sample(1:nages,153)]=round(vaccineuptake*N/153)
       
      }
      
      Q=rep(0,nages)
      D=rep(0,nages)
      Hinc=rep(0,nages)
      QA=rep(0,nages)
      EV = rep(0,nages) 
      
      #Adjust population
      S = S-(E+A+I+P+H+R+Q+D+Hinc+QA+EV)
     
       for(i in 1:length(S)) {
        
        if(S[i] < 0) {
          
          R[i] <- R[i] + S[i]#this makes sure the population size remains constant (although I guess means less are in R if it was E or I that is the problem?) 
          S[i] <- 0 #this makes sure there are no negatives
          
        }
        
        
         
      }
      
      
      N <- S+E+A+I+P+H+R+Q+D+Hinc+QA+EV
      
      #For counting N in each scenario
      #Ntot <<- data.frame(scenario = paste(infected[i], a_test_uptake_list[i], vaccine_uptake_list[i]), Ntot = N)
      
    })
  })
  #attach(params)
  set.seed(params$seed)
  sourceCpp('R/studentmodel_reactive3_Emily9May22.cpp')
  
  nsim=100
  
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
  
  #Lots of the code edited out here is for when I want to make other figures in the paper
  #Cumulative_incidence=getstate(result.rep,c(1,2,3,4,6,8,9)) #all infected compartments E, A, I, P, Q ,QA and R on the last time step
  
  # Cumulative_incidence %>%
  #   group_by(time) %>%
  #   summarise_all(mean) -> meanout
  # 
  #maxbaseline=max(allinf$total)*1.2
  
  #calculating the cumulative incidence for this run
  
  #get all the recovereds
  #Recovered <-  select(result.rep, starts_with("R"))
  
  #Recovered$total <- rowSums(Recovered)
  
 # Recovered <- cbind(result.rep$nsim, result.rep$time, Recovered)
  #find R on the last time step
  #Recovered_last_step <- filter(Recovered, Recovered$`result.rep$time` %in% 90)
  
  #add to all the infecteds on the last time step
  #asymptomatic_last_time_step <- filter(allinf, allinf$time %in% 90)
  #symptomatic_last_time_step <- filter(Infecteds, Infecteds$time %in% 90)
  
  #cumulative incidence
 # cumu_incid <- data.frame(total =(Recovered_last_step$total + asymptomatic_last_time_step$total + symptomatic_last_time_step$total),
                          # sim =seq(1,nsim),
                          # scenario = paste(infected[i], a_test_uptake_list[i], vaccine_uptake_list[i]))
  #cumu_incid$prop <- (Recovered_last_step$total + asymptomatic_last_time_step$total + symptomatic_last_time_step$total)/sum(Ntot$Ntot)
  #cumu_incid$Ntot <- sum(Ntot$Ntot)
  
  
  
  #cumu_incid <- data.frame()
  #sim_nos <- unique(allinf$nsim)
  #for (j in 1:nsim) {
   # cumu_incid[j,1] <- j
    #cumu_incid$cumulativeincidence[j] <- sum(Infecteds$total[Infecteds$nsim==j])+ sum(allinf$total[allinf$nsim==j])
    #cumu_incid$cumulativeincidence[j] <- sum(Infecteds$total[Infecteds$nsim==j])+ sum(allinf$total[allinf$nsim==j])
    #cumu_incid$initially_infected <- infected[i]
    #cumu_incid$asymptomatic_test_uptake <-  a_test_uptake_list[i]
    #cumu_incid$vaccine_uptake <- vaccine_uptake_list[i]
    
    # }
  
  #cumulative_incidence <- rbind(cumulative_incidence, cumu_incid)
  
  ##get all the isolated
  
 #Isolated <-  select(result.rep, starts_with(c("Q", "QA")))
  
  #Isolated$total <- rowSums(Isolated)
  
 # Isolated <- data.frame(scenario = paste(infected[i], a_test_uptake_list[i], vaccine_uptake_list[i]),
                      #   sim = result.rep$nsim,
                      #   time = result.rep$time,
                      #  total = Isolated$total)
  
#  Isolated_allsims <- rbind(Isolated_allsims, Isolated)
  
  ##get all the infected
  
  #Infected <-  select(result.rep, starts_with(c("E", "A", "I", "H", "P","Q")))
  
  #Infected$total <- rowSums(Infected)
  
  #Infected <- data.frame(scenario = paste(backgroundrate_list[i], a_test_uptake_list[i]),
                     #    backgroundrate = backgroundrate_list[i],
                     #    atestuptake = a_test_uptake_list[i],
                      #   sim = result.rep$nsim,
                      #   time = result.rep$time,
                      #   total = Infected$total)
 # print(a_test_uptake_list[i])
  
 # Infected_allsims <- rbind(Infected_allsims, Infected)
  
  
  

  plotcases(allinf,Infecteds,maxbaseline,'baseline')
 
 # # allinf$atestuptake <- a_test_uptake_list[i]
 #  #allinf$backgroundrate <- backgroundrate_list[i]
 #  
 #  Cumulative_incidence$atestuptake <- a_test_uptake_list[i]
 #  Cumulative_incidence$backgroundrate <- backgroundrate_list[i]
 #  Cumulative_incidence$initially_infected <- infected[i]  
 #  Cumulative_incidence$vaccine_uptake <- vaccine_uptake_list[i]
 #  
 #  #just have the last time step
 #  
 #  Cumulative_incidence <- filter(Cumulative_incidence, Cumulative_incidence$time == 90)
 #  
 #  #allinf_total <- rbind(allinf_total, allinf)
 #  Cumulative_incidence_total <- rbind(Cumulative_incidence_total, Cumulative_incidence)
  end.time = Sys.time()
 #  
  #print(end.time-starttime)
  print(vaccine_uptake_list[i])
  
  #mtext("Vaccine Uptake:Low to High", side = 1, line = 1, outer = TRUE, cex = 1.3)
  #mtext("Numbers of initially infected students: Low to High", side = 2, line =0, outer=TRUE, cex = 1.3, las =0)

  
  
  
}

