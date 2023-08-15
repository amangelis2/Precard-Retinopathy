################################################################################
## File: 01_Retinopathy Baseline extraction
## Purpose: Import and combine retinopathy raw data, prepare datasets 
## for analysis 
## Author: Tasos Mangelis
## Last update: 22.09.2022
################################################################################

# Libraries ----

library("readxl")
library("readr")
library("tidyverse")
library("dplyr")
library("lubridate")

# Data loading ----

sheet_names <- c("weight", "acr", "trig", "bmi", "demographics", 
                 "pcr", "bp", "tot_chol", "hba1c", "hdl", "ldl", 
                 "ret", "creat_serum", "smoke", "medsfull")

dataset_names <- c("raw_weight", "raw_acr", "raw_trig", "raw_BMI", "raw_demogs", 
                   "raw_pcr", "raw_BP", "raw_chol", "raw_HBA1C", "raw_HDL", 
                   "raw_LDL", "raw_ret", "raw_creat", "raw_smoke", "raw_med")

imported_data <- list()

# Define file location
Dataload_file_location <- "E:\\Precard project\\Data\\Initial_merge\\"
Datasave_file_location <- "E:\\Precard project\\Data\\Retinopathy_DTI_files\\"


for (i in seq_along(sheet_names)) {
  sheet_name <- sheet_names[i]
  dataset_name <- dataset_names[i]
  
  imported_data[[dataset_name]] <- read_excel(file.path(Data_file_location, "Precardc_Full.xlsx"), sheet = sheet_name)
}

# Data cleaning Retinopathy dataset ----

retinopathy_main<-imported_data$raw_ret%>%
  filter(RM!="U")%>%
  mutate(RM=ifelse(RM=="R3S"|RM=="R3A", "R3", RM),
         RM=ifelse(RM=="M1"|RM=="R3"|RM=="R1", NA, RM),
         R=ifelse(R=="M1", NA, R),
         M=ifelse(M=="3S"|M=="3A"|M=="R3"|M=="R1", NA, M))%>%
  inner_join(imported_data$raw_demogs%>%filter(!is.na(Ethnicity_groups)), by = "R50_ID")%>%
  mutate(Age=as.numeric(as_date(YEAR_TEST-YOBirth)),
         TEST_DATE=as.Date(TEST_DATE,format="%Y-%m-%d", origin = '1970-01-01'), 
         Diabetes_Dur=as.numeric(YEAR_TEST-DIAGNOSIS_YEAR))%>%
  group_by(R50_ID)%>% 
  # Calculate the follow up time per patient
  mutate(FU_time = round(as.numeric(difftime(max(as.Date(TEST_DATE)), min(as.Date(TEST_DATE)),
                         units = "days")/365), 2))%>%
  ungroup()%>%
  filter(YEAR_TEST>=2004)

write.csv(retinopathy_main, paste0(Datasave_file_location, "Retinopathy_initial_data.csv"))

# Baseline retinopathy data ----

Baseline_dataset<-retinopathy_main%>%
  group_by(R50_ID)%>%
  arrange(TEST_DATE)%>%slice(1)%>%
  ungroup()%>%
  # Exclude those patients that have R2 or higher at baseline
  filter(R!="R3"| R!="R2"|M!="M1" )%>%
  rename("TEST_DATE_bsl"="TEST_DATE", "YEAR_TEST_bsl"="YEAR_TEST", "RM_bsl"="RM",
         "Diabetes_Dur_bsl"="Diabetes_Dur", "R_bsl"="R", "M_bsl"="M")%>%
  # replace baseline test date with Diagnosis year - 06 - 01 in case Duration<0
  mutate(Diabetes_Dur_bsl=ifelse(Diabetes_Dur_bsl<0, NA, Diabetes_Dur_bsl))


# Add baseline information to main dataset for later calculation of duration, timet to event etc

retinopathy_bsl<-retinopathy_main%>%
  inner_join(Baseline_dataset%>%select(R50_ID, TEST_DATE_bsl, YEAR_TEST_bsl, RM_bsl, R_bsl, M_bsl), by = "R50_ID")

write.csv(retinopathy_bsl, paste0(Datasave_file_location, "REtinopathy_Baseline.csv"))

# BASELINE Datasets by closest test date to baseline RET 2 years ----

retinopathy_merged<-retinopathy_bsl
ret_ID_TD<-retinopathy_bsl%>%
  select(R50_ID, TEST_DATE_bsl)%>%
  mutate(TEST_DATE_bsl=as.Date(TEST_DATE_bsl))

# Exclude raw_demogs (demographic data) from the list of datasets as not needed

datasets_loop<-imported_data
datasets_loop[["raw_demogs"]]<-NULL
datasets_loop[["raw_med"]]<-NULL
datasets_loop[["raw_ret"]]<-NULL

for (i in 1:length(datasets_loop)) {
  print(i)
  df<-datasets_loop[[i]]%>%
    inner_join(ret_ID_TD, by = "R50_ID", all = TRUE)%>%
    mutate(TEST_DATE=as.Date(TEST_DATE), 
           date_diff=round((as.Date(TEST_DATE)-as.Date(TEST_DATE_bsl))/365, 1))%>%
    filter(abs(date_diff)<= 2)%>%
    group_by(R50_ID)%>%
    arrange(TEST_DATE)%>%
    slice(1)%>%
    select(-c(TEST_DATE, TEST_DATE_bsl, date_diff, YEAR_TEST))
  retinopathy_merged<-retinopathy_merged%>%left_join(df, by = "R50_ID")
}

# Calculation of eGFR values according to the CKD-EPI equation ----

# Function of CKD-epi equation calculation
calculate_eGFR <- function(creat, age, gender, ethn) {
  if (ethn != "BLACK") {
    if (gender != "FEMALE") {
      k <- 61.9
      a <- -0.329
    } else {
      k <- 79.6
      a <- -0.411
    }
  } else {
    if (gender != "FEMALE") {
      k <- 61.9
      a <- -0.329
    } else {
      k <- 79.6
      a <- -0.411
    }
  }
  
  egfr <- 141 * ((min(creat / k, 1)) ** a) * ((max(creat / k, 1)) ** (-1.209)) * (0.993 ** age)
  return(egfr)
}

# Initialize a vector to store eGFR values
gfr <- numeric(length(retinopathy_merged$Creatinine))  # Preallocate space for efficiency

# extraction of total CKD_EPI eGFR using calculate_eGFR() function
for (i in 1:length(retinopathy_merged$Creatinine)) {
  creat_val <- retinopathy_merged$Creatinine[i]
  age_val <- retinopathy_merged$Age[i]
  ethn_val <- retinopathy_merged$Ethnicity_groups[i]
  gender_val <- retinopathy_merged$Gender[i]
  gfr[i] <- calculate_eGFR(creat_val, age_val, gender_val, ethn_val)
}

retinopathy_merged$eGFR <- round(gfr, 1)

# Additional changes

retinopathy_merged<-retinopathy_merged%>%
  mutate(ACR=ifelse(ACR>200, 200, ACR), 
         Cholesterol=ifelse(Cholesterol>10, NA, Cholesterol),
         eGFR=ifelse(eGFR>160, 160, eGFR),
         RM=as.factor(RM),
         DM_TYPE_ID=as.factor(DM_TYPE_ID ),
         Ethnicity_groups=as.factor(Ethnicity_groups),
         Gender=as.factor(Gender))

write.csv(retinopathy_merged, paste0(Datasave_file_location, "Retinopathy_full_inclusive_Baseline.csv"))

# Save all data to RDS file to ensure proper class of the variables ----
data_store_list<- list(retinopathy_bsl=retinopathy_bsl, retinopathy_main=retinopathy_main, 
                       retinopathy_merged=retinopathy_merged, Baseline_dataset=Baseline_dataset)
saveRDS(data_store_list, paste0(Datasave_file_location,"Retinopathy_Data_merge.rds"))

# End of code ----
