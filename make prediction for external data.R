#this program calculate risk score for external data provided by users using the models fitted using the SPRINT data. 
##users only need to modify the 8 lines below before running the code#####
event.name="pd_mci_protocol_death"  ##What composite event you want to predict? options are "pd_mci_amnestic_death" or "pd_death" or "pd_mci_protocol_death" 
predictor.set="full set"  # What set of predictors you want to use? Options are "full set" or  "reduced set 1",  "reduced set 2"
new.data.file="example data.csv"   #this is the input file you want to make prediction on. Check the file "variables needed for prediction.csv" to see what variables are needed 
export.file.name="output.csv"      #this is the output file you want to save the results as.
#if the outcome variables, pd_mci_amnestic_death,pd_mci_amnestic_death_yrs, pd_death and pd_death_yrs
# are available in the input fie, you can calibrate the baseline survival S0 and calculate the risk. Do you want to do so? Otherwise the baseline survival estimated using the SPRINT data will be used.
calibrate="NO"                    #calibrate and calculate risk?
t.assess=4.13                     #at what time point (years) you want to assess risk?
###############End input setup##########################################################

#########################################################################
###############Users don't touch the code below!!!!!#####################
#########################################################################
library(splines)
library(dplyr)
require(tidyr)
library(sas7bdat)
library(fastDummies)
formulas=read.csv("./formulas.csv")

#The models contain spline terms of age and MOCA. The following lines create spline terms from spline term coefficients
beta.spline=read.csv( "./coefficients for spline terms.csv" )
beta.spline=beta.spline%>%separate(segments, sep="--",into=c("lower","upper"), remove = FALSE)%>%
  mutate(lower=as.numeric(lower),upper=as.numeric(upper))

get.splines <- function(b, x) {
  lower.bound=min(b$lower)
  upper.bound=max(b$upper)
  if(x>=lower.bound & x<=upper.bound){
    s=b %>%filter((x>=lower) & (x<=upper))%>%mutate(terms=beta*x^power)%>%
      group_by(spline.term)%>%summarise(terms=sum(terms))
    output=t(s$terms)
  }else if(x<lower.bound){ #linear trend when x is out of bounds
    delta=0.000001
    x1=lower.bound
    x2=x1+delta
    s1=b %>%filter((x1>=lower) & (x1<=upper))%>%mutate(terms=beta*x1^power)%>%
      group_by(spline.term)%>%summarise(terms=sum(terms))
    s2=b %>%filter((x2>=lower) & (x2<=upper))%>%mutate(terms=beta*x2^power)%>%
      group_by(spline.term)%>%summarise(terms=sum(terms))
    slope=(s2$terms-s1$terms)/delta
    output=t(s1$terms+slope*(x-x1))
  }else if(x>upper.bound){
    delta=0.000001
    x1=upper.bound
    x2=x1-delta
    s1=b %>%filter((x1>=lower) & (x1<=upper))%>%mutate(terms=beta*x1^power)%>%
      group_by(spline.term)%>%summarise(terms=sum(terms))
    s2=b %>%filter((x2>=lower) & (x2<=upper))%>%mutate(terms=beta*x2^power)%>%
      group_by(spline.term)%>%summarise(terms=sum(terms))
    slope=(s2$terms-s1$terms)/(-delta)
    output=t(s1$terms+slope*(x-x1))
  }
  output
}

new.data=read.csv(new.data.file)

#create the spline terms using the coefficients
b1=beta.spline%>%filter(var=="BLage")%>%select(-var)
b2=beta.spline%>%filter(var=="BLMoCA")%>%select(-var)
new.data=new.data%>%rowwise()%>%mutate(age.spline1=get.splines(b1,x=BLage)[1],
                                       age.spline2=get.splines(b1,x=BLage)[2],
                                       age.spline3=get.splines(b1,x=BLage)[3],
                                       MoCA.spline1=get.splines(b2,x=BLMoCA)[1],
                                       MoCA.spline2=get.splines(b2,x=BLMoCA)[2],
                                       MoCA.spline3=get.splines(b2,x=BLMoCA)[3] )

#####end create spline terms########

new.data$BLeducation=as.factor(new.data$BLeducation)
dummy.edu=dummy_cols(new.data$BLeducation)
dummy.edu=dummy.edu[,3:5]
names(dummy.edu)=paste("education",2:4,sep="")
new.data$BLcaresource=as.factor(new.data$BLcaresource)
dummy.BLcaresource=dummy_cols(new.data$BLcaresource)
dummy.BLcaresource=dummy.BLcaresource[,3:4]
names(dummy.BLcaresource)=paste("BLcaresource",2:3,sep="")
new.data=data.frame(new.data,dummy.edu,dummy.BLcaresource )

#make prediction of the risk score (X*beta) using the formula
formula.m.z=formulas%>%filter(outcome==event.name & predictor==predictor.set)%>%select(formula.m.z) 
formula.risk.score=formulas%>%filter(outcome==event.name & predictor==predictor.set)%>%select(formula.risk.score)
new.data=new.data%>%mutate(m.z=eval(parse(text=formula.m.z)))
new.data=new.data%>% mutate(risk.score=eval(parse(text=formula.risk.score)))
 


#make prediction of risk
if(toupper(calibrate)=="NO"){
  S0_t=read.csv( "./baseline survival.csv" )
  S0_t=S0_t%>%filter(event_name==event.name & predictors==predictor.set)
  S0=S0_t$S0.survival[which(abs(S0_t$S0.time-t.assess)==min(abs(S0_t$S0.time-t.assess)))]
  
}
if(toupper(calibrate)=="YES"){#calibrate and calculate risk?
  #this function calculates the predicted baseline survival given risk score (X*Beta)
  calculate.survival0=function(followup.time,status,xbeta){
    t=followup.time
    
    data<-data.frame(t_event=followup.time, event=status,xbeta)
    
    #------preparation------
    tab <- data.frame(table(data[data$event == 1, "t_event"])) 
    y <- as.numeric(levels(tab[, 1]))[tab[, 1]] #ordered distinct event times
    d <- tab[, 2]                               #number of events
    
    #------my quick implementation to calculate h0------
    
    
    h0 <- rep(NA, length(y))
    for(l in 1:length(y))
    {
      h0[l] <- d[l] / sum(exp(data$xbeta[data$t_event >= y[l]]))
    }
    H0=cumsum(h0)
    S0=exp(-H0)
    #plot(y,S0)
    
    survival=data.frame(time=y,survival=S0)
    survival
  }
  
  
  S0_t=calculate.survival0(followup.time=new.data[,paste(event.name,"_yrs",sep="")],status=new.data[,event.name],xbeta=new.data$risk.score)
  S0=S0_t$survival[which(abs(S0_t$time-t.assess)==min(abs(S0_t$time-t.assess)))]
  S0
}


new.data=new.data%>%mutate(survival=S0^exp(risk.score), risk=1-survival)
write.csv(new.data,export.file.name,row.names = F)
