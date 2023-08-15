################################################################################
## File: 02_Retinopathy data Imputation 
## Purpose: Load and impute Type 1 data 
## Author: Tasos Mangelis
## Last update: 22.09.2022
################################################################################

# Libraries ----

library("readxl")
library("tidyverse")
library("dplyr")
library("lubridate")
library("mice")

# Load dataset ----

# Define file location
Dataload_file_location <- "E:\\Precard project\\Data\\Initial_merge\\"
Datasave_file_location <- "E:\\Precard project\\Data\\Retinopathy_DTI_files\\"

# Load data and keep Type I patients only
Retinopathy_Data_merge<-readRDS(paste0(Data_file_location,"Retinopathy_Data_merge.rds"))

retinopathy_DTI<-Retinopathy_Data_merge$retinopathy_merged%>%filter(DM_TYPE_ID=="DT_I")

# impute missing values ----

# Step1: Check missingness xclude variables with more than 40% missing

retinopathy_DTI<-retinopathy_DTI%>%
  select(-c(PCR, Smoking, DM_TYPE_ID, HDL, LDL))
# Step 2: - define variable that need to be imputed, 
#         - variables not to be imputed but used as predictors and 
#         - variables to be completely excluded from the process
#         - impute missing values using pmm method from mice package

prediction<-mice(retinopathy_DTI%>%
                   select(R50_ID, TEST_DATE, Age, Gender, Ethnicity_groups, eGFR, 
                          Weight, ACR, Triglycerids, BMI, SYSTOLIC, DIASTOLIC, 
                          Cholesterol, Hba1c, RM), 
                 m=5, maxit = 15, method = c('pmm'), seed = 500, 
                 predictorMatrix = quickpred(retinopathy_DTI%>%
                                               select(R50_ID, TEST_DATE, Age, Gender, Ethnicity_groups, eGFR, 
                                                      Weight, ACR, Triglycerids, BMI, SYSTOLIC, DIASTOLIC, 
                                                      Cholesterol, Hba1c, RM), 
                                             include = c("Age", "Gender", "Ethnicity_groups", "Cholesterol", "Hba1c", "eGFR", "Weight", "ACR", "Triglycerids", 
                                                         "BMI", "SYSTOLIC", "DIASTOLIC"),
                                             exclude = c("R50_ID", "TEST_DATE", "RM")))

# Summary of the imputation and density plot for checking new distribution against the old in 04a

# Assignment of the iteration for the imputed dataset
retinopathy_DTI_full<-mice::complete(prediction, 5)
retinopathy_DTI_full<-retinopathy_DTI_full%>%
  inner_join(retinopathy_DTI%>%select(R50_ID, TEST_DATE, TEST_DATE_bsl, YEAR_TEST_bsl, YODeath, YOBirth, DIAGNOSIS_YEAR, Diabetes_Dur, RM_bsl), by = c("R50_ID", "TEST_DATE"))%>%
  distinct(R50_ID, TEST_DATE, eGFR, .keep_all = T)

write.csv(retinopathy_DTI_full, paste0(Datasave_file_location, "Retinopathy_Type1_full.csv"))
# End of code ----
