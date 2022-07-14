#this program calculate risk score for data provided by users
##users only need to modify the 4 lines below before running the code#####
event.name="pd_mci_protocol_death" #or  "pd_mci_amnestic_death" # or "pd_death"
predictor.set="full set"  # or  "reduced set 1",  "reduced set 2"
new.data.file="example data.csv"   #this is the input file we want to make prediction on. Check he file "variables needed for prediction.csv" to see what variables are needed 
export.file.name="output.csv"      #this is the output file
#if the outcome variables, pd_mci_amnestic_death,pd_mci_amnestic_death_yrs, pd_death and pd_death_yrs
# are available in the input fie, you can calibrate the baseline survival S0 and calculate the risk. Do you want to do so?
calibrate="NO"                    #calibrate and calculate risk?
t.assess=4.13                     #at what time point you want to assess risk?
###############End##########################################################

library(splines)
library(dplyr)
require(tidyr)
library(sas7bdat)
library(fastDummies)
library(ggplot2)

new.data=read.csv(export.file.name)
all.predictions=read.csv(paste("../results/no cross validation/prediction for individuals all models iterate.csv",sep=""))
prediction=all.predictions%>%filter(event_name==event.name & maskid!="S54755")
plot(prediction$m.z.for.check, new.data$m.z)
abline(a=0, b=1,col="red") 

plot(prediction$PredictedXBeta, new.data$risk.score)
abline(a=0, b=1,col="red") 
hist(prediction$PredictedXBeta- new.data$risk.score)
plot(prediction$PredictedXBeta,prediction$PredictedXBeta- new.data$risk.score)

# for (var in names(prediction)[5:68]){
#   plot(prediction[,var],prediction$PredictedXBeta- new.data$risk.score)
# }



plot(prediction$PredictedSurvival, new.data$survival)
abline(a=0, b=1,col="red") 

plot((prediction$PredictedSurvival), (new.data$survival))
abline(a=0, b=1,col="red") 


plot(prediction$PredictedRisk, new.data$risk)
abline(a=0, b=1,col="red") 




