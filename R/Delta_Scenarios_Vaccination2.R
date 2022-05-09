library(Rcpp)
library(tidyverse)
library(gplots)
library(ggpubr)
source("R/plotfunctions_withgetstate.R") 

starttime=Sys.time()




load("data/betamat.RData")
load("data/Npop.Rdata")
years=read.csv("data/yearnumbers1.csv")
nages=length(Npop)


#Low natural immunity levels (5% in recovered compartment)
#Closed student population (background infection rate = 0)

#Range of vaccine uptake
#Low- 50%
#Medium- 70%
#High- 90% 

#The way I've made the lists makes usre there are all combinations of the three parameter sets
vaccine_uptake_list <- c(rep(rep(c(0.3, 0.5, 0.7, 0.9), 4),1))
a_test_uptake_list <-  c((rep(c(rep(0.9,4),rep(0.5,4),rep(0.2,4), rep(0,4)),1)))

#Range of percentage of student population initially infected
#Low- 0.05%
#Medium- 0.1%
#High- 0.2% 
#infected <- c(rep(0.005,16), rep(0.01,16), rep(0.02,16))
infected <- c(rep(0.01,16))


#Range of student population with natural immunity
#Low- 5%
#Medium- 10%
#High- 20%
#natural_immunity <- c(rep(c(rep(0.05,3),rep(0.1,3), rep(0.2,3)),3))
natural_immunity <- c(rep(rep(0.1,length(vaccine_uptake_list)),1))
#natural_immunity <- c(rep(0.05,9), rep(0.1,9), rep(0.2,9))
#natural_immunity_reinfection <- 0.97 #this was a parameter that was used to assume the reinfection rate for those who have natural immunity
#natural_immunity_reinfection <- 0.97 #if 1 then it's not having an impact
#Repeat with and without background infection on

backgroundrate_list <- c(rep(1e-3, length(vaccine_uptake_list)))
#backgroundrate_list <- c(1e-3, 0, 1e-3,0, 1e-3,0, 1e-3, 0)

par(mar=c(2, 3, 2, 2),mfrow=c(4,4))


Ntot <- data.frame()
#cumulative_incidence <- data.frame() #this will record the cumulative incidence of each scenario
#Isolated_allsims <- data.frame() #this will calcualte the number in isolation in all scenarios

#Infected_allsims <- data.frame()

#allinf_total <- data.frame()
Cumulative_incidence_total <- data.frame()
#Recoverds_total <- data.frame()
#Isolated_total <- data.frame()
#this is the start of a for loop which can be used to do sensitivity analysis with the different parameter sets
for (i in 1:length(vaccine_uptake_list)) {
  
  
  params <- list()
  params <- within(params, {
    
    SIMTIME=91
    iperiod=5.0; 
    load("data/betamat.RData")
    load("data/Npop.Rdata")
    years=read.csv("data/yearnumbers1.csv")
    nages=length(Npop)
    
    gam_p=1/2;
    
    
    #Including the impact of vaccination on transmission
    #*(1-(uptake of vaccine*(1-efficacy of vaccine in reducing transmission))
    #*
    
     vaccine_uptake <- vaccine_uptake_list[i]
    # vaccine_transmission <- 0.45
    # vaccine_susceptibility <- 0.8  #Warwick estimates for alpha or delta
    # v_s <- 1-vaccine_susceptibility
    # v_t <- 1-vaccine_transmission
    ve_inf <- 0.45
    ve_trans <- 0.8
    
    
    #beta=(betamat/(iperiod+1/gam_p))*(1-(vaccine_uptake*vaccine_transmission*vaccine_susceptibility))
    #beta= ((betamat/(iperiod+1/gam_p))*(v_s*v_t*(vaccine_uptake^2)+ vaccine_uptake*(1-vaccine_uptake)*(v_s+v_t) + ((1-vaccine_uptake)^2)))*1.5*1.6
    beta=betamat/(iperiod+1/gam_p)
    
    #beta=(betamat/(iperiod+1/gam_p))*(1-(vaccine_uptake*(1-vaccine_transmission)))*(1-(vaccine_uptake*(1-vaccine_susceptibility)))
    #beta=betamat/(iperiod+1/gam_p)
    gamma=1/iperiod;  sigma=1/(3.2);  gam_a=1/(iperiod+1/gam_p); gam_h=1/3; gam_q=1/10 ##recovery rates
    
    #uptake of asymptomatic testing
    a_test_uptake <- a_test_uptake_list[i]
    daysto_asy_test <- 1/4
    LFTsensitivity <- 0.58 #(for asymptomatic people according to this Cochrane report https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.CD013705.pub2/full)
    
    testrate_a_single_value <- a_test_uptake * daysto_asy_test * LFTsensitivity #higher sensitivity reduces the time to test
    
    testrate=1/2; testrate_a=rep(testrate_a_single_value,nages); testrate_a2=0; testrate_a3=0; ##testing parameters
    

    
    eps=0.5; eps_q=0.5 ##
    #h=0.002 * 0.04 #hospitalisation rate adjusted according to vaccine efficacy
    h = 0 #hospitalisation switched off
    vaccine_symptoms <- 0.84 #2nd dose pfzier 84-95% 
    #f=0.75/(1+vaccine_uptake*vaccine_symptoms) #((0.25 * vaccine_symptoms)*vaccine_uptake) #fraction asymptomatic - increase this because of vaccination
    f = 0.75*(1-vaccine_uptake) + 0.75*vaccine_uptake/vaccine_symptoms
    #mrate=0.038 * 0.11 #mortality rate adjusted according to vaccine efficacy
    mrate = 0 #mortality rate swtiched off
    backgroundrate = backgroundrate_list[i] #background infection rate
    waning <- 1/243
    
    

    ## initial conditions, list within list
    
    ## S = Number of students - number in all other compartments
    ## E = 50% of all infections
    ## A = 40% of all infections
    ## P = 5% of all infections
    ## I = 5% of all infections
    ## Q = 0
    ## R = 5%, 10% or 20% of the student population
    
    
    init <- within(list(), {
      S=Npop
      N= sum(Npop)
      E=rep(0,nages)
      
      
      #E[sample(1:nages,round(((infected[i]*N))/3))]=1 #at the moment the proportion infected is divided between E A and P
      A=rep(0,nages)
      #A[sample(1:nages,round(((infected[i]*N))/3))]=1
      
      
      I=rep(0,nages)
      #I[sample(1:nages,round(((infected[i]*N))/3))]=1
      
      P=rep(0,nages)
      
      if((infected[i]*N/2) <= 153) {
        E[sample(1:nages,round(((infected[i]*N))*0.5))]=1 #problem line
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
      
      #this is adjusted based on the natural_immunity levels and divided evenly between groups
      #if(natural_immunity[i]*natural_immunity_reinfection*N <= 153) {
        if(natural_immunity[i]*N <= 153) {
          
      #R[sample(1:nages,round(natural_immunity[i]*natural_immunity_reinfection*N))]=1 
      R[sample(1:nages,round(natural_immunity[i]*N))]=1 
      
     
      
      } else {
       #R[sample(1:nages,153)]=round(natural_immunity[i]*natural_immunity_reinfection*N/153)
       R[sample(1:nages,153)]=round(natural_immunity[i]*N/153)
       
      }
      
      #testing natural immunity
      #R[sample(1:nages,153)]=9
      
      
      Q=rep(0,nages)
      D=rep(0,nages)
      Hinc=rep(0,nages)
     
      QA=rep(0,nages)
      #Adjust population
      S = S-(E+A+I+H+R+P+Q+QA)
     
       for(i in 1:length(S)) {
        
        if(S[i] < 0) {
          
          R[i] <- R[i] + S[i]#this makes sure the population size remains constant (although I guess means less are in R if it was E or I that is the problem?) 
          S[i] <- 0 #this makes sure there are no negatives
          
        }
        
      EV = rep(0,nages)   
         
      }
      
      
      N <- S+E+A+I+H+R+P+Q+QA+EV
      Ntot <<- data.frame(scenario = paste(infected[i], a_test_uptake_list[i], vaccine_uptake_list[i]), Ntot = N)
      
    })
  })
  #attach(params)
  set.seed(params$seed)
  sourceCpp('R/studentmodel_reactive2_Emily9May22.cpp')
  
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
  #Cumulative_incidence=getstate(result.rep,c(1,2,3,4,6,8,9)) #all infected compartments E, A, I, P, Q ,QA and R on the last time step
  
  # Cumulative_incidence %>%
  #   group_by(time) %>%
  #   summarise_all(mean) -> meanout
  # 
  maxbaseline=max(allinf$total)*1.2 # why?
  
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


 




write.csv(Cumulative_incidence_total, "Cumulative_incidence_total_9thNov21_endtime.csv")
#ggsave("/Users/em15116/OneDrive\ -\ University\ of\ Bristol/Covid19/Universities/Uni_Vaccination/Fig5_90percent_vaccinated.png", last_plot(), width=10, height=10, dpi = 500)
#plotting figure 5 - infecteds 


#combining data




#Infected_allsimssave <- Infected_allsims  

#Infected_allsims$atestuptake <- as.factor(Infected_allsims$atestuptake)


#Infected_allsims$date <- c(rep(seq(as.Date("20/09/2021"), by = "day", length.out = 91),800))


 cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  #make into percentage
  Cumulative_incidence_total$atestuptake <- Cumulative_incidence_total$atestuptake * 100
  #allinf_total$atestuptake <-  allinf_total$atestuptake * 100
  
  #allinf_total$atestuptake <- as.factor(allinf_total$atestuptake)
  Cumulative_incidence_total$atestuptake <- as.factor(Cumulative_incidence_total$atestuptake)
  #Infecteds_total$asymp_and_symp <- Infecteds_total$total + allinf_total$total
  #Cumulative_incidence_total$date <- c(rep(seq(as.Date("20/09/2021"), by = "day", length.out = 91),48))
  
  Cumulative_incidence_total <- as_tibble(Cumulative_incidence_total)
  
  #withbackgroundrate <- Isolated_total %>%
   # filter(., backgroundrate %in% 0.001) %>%
    #group_by(atestuptake,`date`) %>%
    #summarise(median = median(`total`),
           #   lll = min(`total`),
          #    hhh = max(`total`)) %>%
 #   ggplot(aes(fill=`atestuptake`))+
    #scale_fill_manual(values=cbPalette)+
 #   geom_ribbon(aes(x = `date`, ymin =`min`, ymax = `max`))+
    #geom_line (aes(x = `date`, y = `median`, colour= `atestuptake`),  show.legend =FALSE) +
   # scale_fill_manual(values=cbPalette,  aesthetics = c("colour", "fill"))+
   # ggtitle("(a)")+
  #  xlab("Date")+
  #  ylab("Daily number of students infected")+
  #  ylim(0,6000)+
    #xlim("20/09/2021", "0020-12-19")+
    #xlim(0,91)+
   # theme_bw()+
    #theme(plot.title = element_text(size = 12, family = "Helvetica", colour = "black"), 
      #    axis.title = element_text(size = 18, family = "Helvetica", colour = "black", face = "bold"),
     #     axis.text = element_text(size = 12, family = "Helvetica", colour = "black"))+
   # guides(fill=guide_legend(title="Asymptomatic testing uptake (%)"))
  
  
  #withoutbackgroundrate <- Infecteds_total %>%
  #  filter(., backgroundrate %in% 0) %>%
    #group_by(atestuptake, date) %>%
    #summarise(median = median(`asymp_and_symp`),
           #   lll = min(`asymp_and_symp`),
            #  hhh = max(`asymp_and_symp`)) %>%
   # ggplot(aes(fill=`atestuptake`))+
    #scale_fill_manual(values=cbPalette)+
   # geom_ribbon(aes(x = `date`, ymin =`min`, ymax = `max`)) +
    #geom_line ( aes(x = `time`, y = `total`, colour= `atestuptake`),  show.legend =FALSE) +
   # scale_fill_manual(values=cbPalette,  aesthetics = c("colour", "fill"))+
   # ggtitle("(b)")+
   # xlab("Date")+
    #ylab("")+
    #ylim(0,6000)+
   # theme_bw()+
    #theme(plot.title = element_text(size = 12, family = "Helvetica", colour = "black"), 
      #    axis.title = element_text(size = 18, family = "Helvetica", colour = "black", face = "bold"),
      #    axis.text = element_text(size = 12, family = "Helvetica", colour = "black"))+
   # guides(fill=guide_legend(title="Asymptomatic testing uptake (%)"))
  


  
  #ggarrange(withbackgroundrate, withoutbackgroundrate, ncol=2, nrow= 1, common.legend = TRUE, legend = "bottom", align="v")
  
 # ggsave("/Users/em15116/OneDrive\ -\ University\ of\ Bristol/Covid19/Universities/Uni_Vaccination/Fig5_9thNov21.png", last_plot(), width=10, height=10, dpi = 500)
 # write.csv(Isolated_total, "Isolated_total_15thNov21.csv")
#Plotting summary of cumulative results for each scenario

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)


#using the function to calculate the mean and se of the simulations
#summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                    #  conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
#write.csv(result.rep, "result.rep_univacpaper_21stOct2021_v1.csv")
#write.csv(cumulative_incidence, "cumulativeincidence_21stOct2021_v1.csv")

#write.csv(Isolated, "Isolated_21OCt2021_v1.csv")

#plot isolated

#Isolated_allsims <- as_tibble(Isolated_allsims)
#Isolated_allsims %>% group_by(`scenario`, `time`) %>%
 # summarise(median = median(`total`),
        #    lll = quantile(`total`, 0.025),
         #   hhh = quantile(`total`, 0.975)) -> Isolated_summary
#

## BE CAREFUL as these don't automatically update

#SIMTIME <- 91
#Isolated_summary$initially_infected <- c( rep(0.005,20*SIMTIME), rep(0.01,20*SIMTIME), rep(0.02,20*SIMTIME))
#
#Isolated_summary$asymptomatic_test_uptake <- c(rep(rep(c(rep("0%",5*SIMTIME),rep("20%",5*SIMTIME),rep("50%",5*SIMTIME),rep("90%",5*SIMTIME))),3))

#Isolated_summary$vaccine_uptake <- c(rep(rep(c(rep("30%", SIMTIME),  rep("50%", SIMTIME), rep("60%", SIMTIME), rep("70%", SIMTIME), rep("80%", SIMTIME)), 4),3))

#Isolated_summary$asymptomatic_test_uptake <- factor(Isolated_summary$asymptomatic_test_uptake)

colnames(Cumulative_incidence_total)[11] <- "Asymptomatic test uptake"
colnames(Cumulative_incidence_total)[14] <- "Vaccine uptake"

#cbPalette_2 <- c("#999999", "#999999", "#E69F00","#E69F00", "#56B4E9", "#56B4E9", "#009E73", "#009E73", "#F0E442", "#F0E442") #, "#0072B2", "#D55E00", "#CC79A7")

#make it have dates

#Isolated_summary$date <- c(rep(seq(as.Date("20/09/2021"), by = "day", length.out = 91),60))
  
# Isolated_total$`Vaccine uptake` <- as.factor(Isolated_total$`Vaccine uptake`)
# Isolated_total$`Asymptomatic test uptake` <- as.factor(Isolated_total$`Asymptomatic test uptake`)
# 
# initialinf5e04_isolate <- Isolated_total %>%
#   filter(., initially_infected %in% "0.005") %>%
#   ggplot() +
#   geom_ribbon(aes(x = `date`, ymin = (`min`/28116)*100, ymax = (`max`/28116)*100, fill = `Vaccine uptake`, lty = `Asymptomatic test uptake`), alpha = 0.3)+
#   geom_line(aes(x= `date`, y= (`total`/28116)*100, colour= `Vaccine uptake`, lty = `Asymptomatic test uptake`))+
#   ylim(0,max((Isolated_total$max/28116)*100))+
#   theme_minimal()+
#   labs(y = "", x = " Day of term")+
#   ggtitle("Initially infected: 0.5%")+
#   theme(legend.position="bottom")+
#   theme(plot.title = element_text(size = 12, family = "Helvetica", colour = "black"),
#         axis.text.y = element_text(size = 12, family = "Helvetica", colour = "black"),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         legend.title =element_text(size = 12, family = "Helvetica", colour = "black"))
  
# initialinf1e03_isolate <- Isolated_total %>%
#   filter(., initially_infected %in% "0.01") %>%
#   ggplot() +
#   geom_ribbon(aes(x = `date`, ymin = (`min`/28116)*100, ymax = (`max`/28116)*100,fill = `Vaccine uptake`, lty = `Asymptomatic test uptake`), alpha = 0.3)+
#   geom_line(aes(x= `date`, y= (`total`/28116)*100, colour= `Vaccine uptake`, lty = `Asymptomatic test uptake`))+
#   ylim(0,max((Isolated_total$max/28116)*100))+
#   theme_minimal()+
#   labs(y = "% of students in isolation", x = " Day of term")+
#   ggtitle("Initially infected: 1%")+
#   theme(legend.position="bottom")+
#   theme(plot.title = element_text(size = 12, family = "Helvetica", colour = "black"),
#         axis.title.y = element_text(size = 18, family = "Helvetica", colour = "black", face = "bold"),
#         axis.text.y = element_text(size = 12, family = "Helvetica", colour = "black"),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())


# initialinf2e03_isolate <- Isolated_total %>%
#   filter(., initially_infected %in% "0.02") %>%
#   ggplot() +
#   geom_ribbon(aes(x = `date`, ymin = (`min`/28116)*100, ymax = (`max`/28116)*100,fill = `Vaccine uptake`, lty = `Asymptomatic test uptake`), alpha = 0.3)+
#   geom_line(aes(x= `date`, y= (`total`/28116)*100, colour= `Vaccine uptake`, lty = `Asymptomatic test uptake`))+
#   ylim(0,max((Isolated_total$max/28116)*100))+
#   theme_minimal()+
#   labs(y = "", x = " Day of term")+
#   ggtitle("Initially infected: 2%") + 
#   theme(legend.position="bottom")+
#   theme(plot.title = element_text(size = 12, family = "Helvetica", colour = "black"), 
#         axis.title = element_text(size = 18, family = "Helvetica", colour = "black", face = "bold"),
#         axis.text = element_text(size = 12, family = "Helvetica", colour = "black"))
#   
# 
# #plot them all on the same plot
# 
# ggarrange(initialinf5e04_isolate, initialinf1e03_isolate, initialinf2e03_isolate, ncol=1, nrow= 3, common.legend = TRUE, legend = "bottom", align="v")
# 
# ggsave("/Users/em15116/OneDrive\ -\ University\ of\ Bristol/Covid19/Universities/Uni_Vaccination/Fig3_isolated_15Nov21.png", last_plot(), width=10, height=10, dpi = 500)
# 
# 




# summary(filter(Isolated_total, scenario %in% "0.02 0 0.8"))
# summary(filter(Isolated_summary, scenario %in% "0.005 0 0.8"))
# 
# summary(filter(Isolated_summary, scenario %in% "0.02 0.9 0.8"))
# summary(filter(Isolated_summary, scenario %in% "0.005 0.9 0.8"))

#ggsave("/Users/em15116/OneDrive\ -\ University\ of\ Bristol/Covid19/Universities/Uni_Vaccination/Fig4_isolationpercent_v1_2percent_initially_infected.png", last_plot(), width=10, height=10, dpi = 500)

#plot all on the same plot later as you have to use package ggarrange which messes up other packages you need first
  


#cumulative_incidence$scenario <- paste(cumulative_incidence$initially_infected, cumulative_incidence$asymptomatic_test_uptake, cumulative_incidence$vaccine_uptake, sep = "")

Cumulative_incidence_total$percent <- (Cumulative_incidence_total$total/sum(Npop))*100

#min and max percents
Cumulative_incidence_total$minpercent <- (Cumulative_incidence_total$min/sum(Npop))*100

Cumulative_incidence_total$maxpercent <- (Cumulative_incidence_total$max/sum(Npop))*100

#tgc <- summarySE(cumulative_incidence, measurevar="percent", groupvars=c("scenario"))

# ###THIS IS ALL WRONG AND I NEED TO SORT IT!!! 
# cumulative_incidence2 <- cumulative_incidence[,2:ncol(cumulative_incidence)]
# #the 2nd lowest sim value, this groups by scenario, picks the two lowest then picks the highest of those 
# PILow <- cumulative_incidence2 %>% group_by(`scenario`) %>% slice_min(., order_by = percent, n = 2) %>% slice_max(., order_by = percent, n = 1)
# 
# colnames(PILow)[5] <- "PILow"
# 
# PILow <- cbind(PILow[,2], PILow[,5])
# 
# #this groups by scenario, for each group picks the three highest and then of those picks the lowest (98th Prediction interval)
# #why does this lead to 62 instead of 60? because for two scenarios two runs led to the same 98th PI 
# PIhigh <- cumulative_incidence2 %>% group_by(`scenario`) %>% slice_max(., order_by = percent, n = 3) %>% slice_min(., order_by = percent, n = 1)
# 
# colnames(PIhigh)[5] <- "PIHigh"
# 
# PIhigh <- cbind(PIhigh[,2], PIhigh[,5])
# 
# #join the prediction interval
# tgc <- right_join(tgc, PILow, by = "scenario")
# 
# tgc <- right_join(tgc, PIhigh, by = "scenario")
# 
# #remove duplicate rows 
# tgc <- distinct(tgc)
# 
# 
# #CAUTION
# #this isn't an automated way of doing this so need to be careful
# #making it so I can plot based on these things
# tgc$initially_infected <- c( rep(0.005,20), rep(0.01,20), rep(0.02,20))
# 
# tgc$asymptomatic_test_uptake <- c(rep(rep(c(rep("0%",5),rep("20%",5),rep("50%",5),rep("90%",5))),3))
# 
# tgc$vaccine_uptake <- c(rep(rep(c("30%", "50%", "60%", "70%", "80%"), 4),3))
# 
# tgc$asymptomatic_test_uptake <- factor(tgc$asymptomatic_test_uptake)
# 
# #make scenario name vaccine uptake just because of labelling
# tgc$scenario <- factor(tgc$vaccine_uptake)





#colour blind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Cumulative_incidence_total$`Vaccine uptake` <- Cumulative_incidence_total$`Vaccine uptake`*100
Cumulative_incidence_total$`Vaccine uptake` <- as.factor(Cumulative_incidence_total$`Vaccine uptake`)

#Cumulative_incidence_total$scenario <- gl(48,91)#change for different combinations

# Cumulative_incidence_total %>% group_by(scenario) %>%
#   summarise(yeartotal = sum(total),
#             yearmin = sum(min),
#             yearmax = sum(max),
#             atestuptake = `Asymptomatic test uptake`,
#             backgroundrate= `backgroundrate`,
#             initially_infected = `initially_infected`,
#             Vaccine_uptake = `Vaccine uptake`) -> CI_summary

#CI_summary <- distinct(CI_summary)

#colnames(CI_summary)[5] <- "Asymptomatic testing uptake"
#colnames(CI_summary)[8] <- "Vaccine uptake"

#CI_summary$percent <- (CI_summary$yeartotal/sum(Npop))*100

#make plot for initially infected 0.5%
# Error bars represent standard error of the mean
initialinf5e04 <- filter(Cumulative_incidence_total, Cumulative_incidence_total$initially_infected %in% 5e-03) %>%
ggplot(., aes(x=`Vaccine uptake`, y=`percent`, fill = `Asymptomatic test uptake`)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=`minpercent`, ymax=`maxpercent`),#change to be the prediction interval
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  scale_fill_manual(values=cbPalette)+
  xlab("Vaccine uptake (%)")+
  ylab("")+
  theme_bw()+
  scale_y_continuous(
    limits =   c(0,80),
    expand = expansion(mult = c(0,0.05)))+
  guides(fill=guide_legend(title="Asymptomatic testing uptake"))+
  ggtitle("Initially infected: 0.5%")+
  theme(plot.title = element_text(size = 12, family = "Helvetica", colour = "black"),
        axis.text.y = element_text(size = 12, family = "Helvetica", colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title =element_text(size = 12, family = "Helvetica", colour = "black"))

#make plot for initially infected 1%
initialinf1e03 <- filter(Cumulative_incidence_total, Cumulative_incidence_total$initially_infected %in% 1e-02) %>%
  ggplot(., aes(x=`Vaccine uptake`, y=`percent`, fill = `Asymptomatic test uptake`)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=`minpercent`, ymax=`maxpercent`),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  scale_fill_manual(values=cbPalette)+
  xlab("Vaccine uptake (%)")+
  ylab("% of students infected")+
  theme_bw()+
  scale_y_continuous(
    limits =   c(0,80),
    expand = expansion(mult = c(0,0.05)))+
guides(fill=guide_legend(title="Asymptomatic testing uptake"))+
ggtitle("Initially infected: 1%")+
  theme(plot.title = element_text(size = 12, family = "Helvetica", colour = "black"),
        axis.title.y = element_text(size = 18, family = "Helvetica", colour = "black", face = "bold"),
        axis.text.y = element_text(size = 12, family = "Helvetica", colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#make plot for initially infected 2%
initialinf2e03 <- filter(Cumulative_incidence_total, Cumulative_incidence_total$initially_infected %in% 2e-02) %>%
  ggplot(., aes(x=`Vaccine uptake`, y=`percent`, fill = `Asymptomatic test uptake`)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=`minpercent`, ymax=`maxpercent`),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  scale_fill_manual(values=cbPalette)+
  xlab("Vaccine uptake (%)")+
  ylab("")+
  theme_bw()+
  scale_y_continuous(
    limits =   c(0,80),
    expand = expansion(mult = c(0,0.05)))+
guides(fill=guide_legend(title="Asymptomatic testing uptake"))+
ggtitle("Initially infected: 2%")+
  theme(plot.title = element_text(size = 12, family = "Helvetica", colour = "black"), 
        axis.title = element_text(size = 18, family = "Helvetica", colour = "black", face = "bold"),
        axis.text = element_text(size = 12, family = "Helvetica", colour = "black"))

#plot them all on the same plot
#library(ggpubr)
ggarrange(initialinf5e04, initialinf1e03, initialinf2e03, ncol=1, nrow= 3, common.legend = TRUE, legend = "bottom", align="v")
ggsave("/Users/em15116/OneDrive\ -\ University\ of\ Bristol/Covid19/Universities/Uni_Vaccination/Fig_cumulative incidence_all_v1.png", last_plot(), width=10, height=10, dpi = 500)

#plot all the isolated on the same plor

ggsave("/Users/em15116/OneDrive\ -\ University\ of\ Bristol/Covid19/Universities/Uni_Vaccination/Fig4_isolationpercent_all_initially_infected_v2.png", last_plot(), width=10, height=10, dpi = 500)








####Old code #######
#calculate isolation 

Isolated <- select(result.rep, starts_with("I"))
Isolated$nsim <- result.rep$nsim
Isolated$time <- result.rep$time

pivot_longer(Isolated, cols = starts_with("I"), names_to = "Group")

Isolated$total <- rowSums(Isolated[,1:153])



ggplot(Isolated) +
  geom_line(aes(x = time, y = total))

plot(Isolated$time,  Isolated$total, type = "l")

Isolated_PIlow <- Isolated %>% group_by(`time`) %>% slice_min(., order_by = total, n = 2) %>% slice_max(., order_by = total, n = 1)

colnames(PILow)[5] <- "PILow"

PILow <- cbind(PILow[,2], PILow[,5])

#this groups by scenario, for each group picks the three highest and then of those picks the lowest (98th Prediction interval)
#why does this lead to 62 instead of 60? because for two scenarios two runs led to the same 98th PI 
Isolated_PIhigh <- Isolated %>% group_by(`time`) %>% slice_max(., order_by = total, n = 3) %>% slice_min(., order_by = total, n = 1)

colnames(PIhigh)[5] <- "PIHigh"

PIhigh <- cbind(PIhigh[,2], PIhigh[,5])




#cumulative incidence
cumu_incid <- data.frame(total =(Recovered_last_step$total + asymptomatic_last_time_step$total + symptomatic_last_time_step$total),
                         sim =seq(1,nsim),
                         scenario = paste(infected[i], a_test_uptake_list[i], vaccine_uptake_list[i]))
cumu_incid$prop <- (Recovered_last_step$total + asymptomatic_last_time_step$total + symptomatic_last_time_step$total)/sum(Ntot$Ntot)
cumu_incid$Ntot <- sum(Ntot$Ntot)


