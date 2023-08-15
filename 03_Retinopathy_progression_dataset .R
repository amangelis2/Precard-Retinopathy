################################################################################
## File: 03_Retinopathy progression dataset 
## Purpose: Compute transition to advance retinopathy and mark as event, 
##          censored or pre event death and time-to-event
## Author: Tasos Mangelis
## Last update: 22.09.2022
################################################################################

# Libraries ----

library("readxl")
library("tidyverse")
library("dplyr")
library("lubridate")

# Load dataset ----

# Clean and define event and time to event for transitions of patients ----

# from baseline retinopathy assessments R01 and M0 to final assessments R2, R3 and/or M1

Ret_prog_DTI<-retinopathy_DTI_full%>%
  # Create an indicator and sort to find the first time the patients progressed
  mutate(rpg=ifelse((RM_bsl=="R0M0"|RM_bsl=="R1M0") & 
                            (RM=="R1M1"|RM=="R2M0"|RM=="R2M1"|
                               RM=="R3M0"|RM=="R3M1"),"R01->23", "no"),
         # flag as 1-2 for sorting
         rpg01=ifelse(rpg=="no", 2, 1))%>%
  group_by(R50_ID)%>%
  # Identify latest test date for all patients
  mutate(TEST_DATE_final=max(TEST_DATE))%>%
  # Identify test dates that an event appears for the first time
  arrange(rpg01, TEST_DATE)%>%slice(1)%>%ungroup()%>%
  # Event and time to event calculations
  mutate(Eventg=ifelse(rpg=="no" & !is.na(YODeath), 2,
                             ifelse(rpg=="no" & is.na(YODeath), 0, 1)),
         TEST_DATE_bsl=as.POSIXct(TEST_DATE_bsl,format="%Y-%m-%d", origin = "1970-01-01"),
         TEST_DATE=as.POSIXct(TEST_DATE,format="%Y-%m-%d", origin = "1970-01-01"),
         T2Eg= ifelse(rpg=="no" & !is.na(YODeath), (as.numeric(YODeath) - as.numeric(YEAR_TEST_bsl)),
                      ifelse(rpg=="no" & is.na(YODeath), round((difftime(TEST_DATE_final, TEST_DATE_bsl, units = "days"))/365,1),
                             round((difftime(TEST_DATE, TEST_DATE_bsl, units = "days"))/365,1))),
         EventCox=ifelse(Eventg==1&rpg=="R01->23", 2,ifelse(Eventg==2, NA, 1)),
         # Other calculations and modifications
         ACR_G=ifelse(ACR<3,"A1", ifelse(ACR>=30,"A3","A2")),
         Duration_diabetes_groups=ifelse(Diabetes_Dur==0, "0",
                                         ifelse(Diabetes_Dur>0&Diabetes_Dur<10, "1-9",
                                                ifelse(Diabetes_Dur>9&Diabetes_Dur<20, "10-20", "20+"))),
         Age_groups=as.factor(ifelse(Age<30, "0-30", ifelse(Age>29&Age<60, "31-60", "60+"))),
         Afrocaribbean_Ethnicity=as.factor(ifelse(Ethnicity_groups=="BLACK", "afro-caribbean", "other")))
# Change level of ethnicity for survival analytics 
Ret_prog_DTI <- within(Ret_prog_DTI, Afrocaribbean_Ethnicity<-relevel(Afrocaribbean_Ethnicity, ref = 2))

write.csv(Ret_prog_DTI, paste0(Datasave_file_location, "Retinopathy_Type1_progression.csv"))
# End of code ----
