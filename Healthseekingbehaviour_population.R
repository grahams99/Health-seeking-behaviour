################################################################################
# Author: Sophie Graham
# Date: 20/03/2023
# Version: R 4.2.2
# File name: healthseekingbehaviour_population.R
# Status: in progress
# CPRD version: May 2022
# Data sets used: patient raw files from CPRD Aurum merged with practice files 
# to obtain practice id and last collection date. Also merged with CPRD Aurum 
# linkage elegibility files and ONS to retrieve death date. 
# R scripts needed: global
# Data sets created: "study_population_covid.parquet", 
# "study_population_flu.parquet" and "study_population_neg.parquet" which are
# the analysis specific populations. It should be noted that for the paper 
# only the flu analysis population was used. 
# Description of file: This script brings in the patient files and checks that 
# the inclusion/exclusion criteria were applied. Then merges the 
# patient file with the practice file to bring in pracid and last collection 
# date. Then it links with ONS data and HES data, removes patients without HES 
# linkage and then creates each of the analysis specific populations (covid, 
# flu, negative exposure). 
################################################################################

#get all patient files from the parquet files
setwd(parquet)
patient <- list.files(path = parquet, pattern = "\\Patient")
patient_df <- list()

#reading in the parquet files and fixing data formatting
for(i in 1:length(patient)){
  
  df <- as.data.table(arrow::read_parquet(file = patient[i]))
  df$regstartdate <- as.Date(df$regstartdate, format = "%d/%m/%Y")
  df$regenddate <- as.Date(df$regenddate, format = "%d/%m/%Y")
  df$cprd_ddate <- as.Date(df$cprd_ddate, format = "%d/%m/%Y")
  df$yob <- as.integer(df$yob)
  df$mob <- as.integer(df$mob)
  df$patienttypeid <- as.integer(df$patienttypeid)
  df$usualgpstaffid <- as.integer(df$usualgpstaffid)
  df$acceptable <- as.factor(df$acceptable)
  df$pracid = as.integer(df$pracid) 
  df$patid <- as.character(df$patid)
  
  patient_df[[i]] <- df
  
}

#total patient population
sum(nrow(patient_df[[1]]), nrow(patient_df[[2]]), nrow(patient_df[[3]]))

#remove unnecessary columns and make date of birth column with imputed month and day and make registration start variable 
for(i in 1:length(patient_df)){
  patient_df[[i]][ ,`:=`(usualgpstaffid = NULL, mob = NULL)]
  patient_df[[i]][, dob := as.IDate(paste0(yob, "0701"), format = "%Y%m%d")]
}

#remove duplicates
sum(nrow(patient_df[[1]]), nrow(patient_df[[2]]), nrow(patient_df[[3]]))
for(i in 1:length(patient_df)){
  patient_df[[i]] <- unique(patient_df[[i]])
}

#check that all patients were aged 66 on 1 September 2019
for(i in 1:length(patient_df)){
  patient_df[[i]] <- patient_df[[i]][dob <= as.Date("01/09/1953", format = "%d/%m/%Y")] 
}

#check that patients have a registration start date that occurs on or before 8 December 2019
for(i in 1:length(patient_df)){
  patient_df[[i]] <- patient_df[[i]][regstartdate <= as.Date("08/12/2019", format = "%d/%m/%Y")] 
}

#check that all patients are acceptable and patient type is all 3 (regular)
for(i in 1:length(patient_df)){
  patient_df[[i]] <- patient_df[[i]][acceptable ==1]
  patient_df[[i]] <- patient_df[[i]][patienttypeid ==3]
}

#check that all patients have registration end date before 1 September 2019
for(i in 1:length(patient_df)){
  patient_df[[i]] <- patient_df[[i]][regenddate >= as.Date("01/09/2019", format = "%d/%m/%Y") | is.na(regenddate)] 
}

#remove patients with indetermined gender
for(i in 1:length(patient_df)){
  patient_df[[i]] <- patient_df[[i]][!gender == 3] 
}

#appending all the patient files together
patient_population <- patient_df[[1]]
for(i in 2:length(patient_df)){
  patient_population <- rbind(patient_population, patient_df[[i]])
}

#now pull in the practice files 
setwd(parquet)
practice <- list.files(path = parquet, pattern = "\\Practice")
practice_df <- list()

#reading in the parquet files and fix the data formatting
for(i in 1:length(practice)){
  df <- as.data.table(arrow::read_parquet(file = practice[i]))
  df$lcd <- as.Date(df$lcd, format = "%d/%m/%Y")
  df$region <- as.integer(df$region)
  df$pracid <- as.integer(df$pracid)
  practice_df[[i]] <- df
}

#quality check - lcd not before start of study
for(i in 1:length(practice_df)){
  n1 <- nrow(practice_df[[i]])
  practice_df[[i]][, lcd >= as.Date("01/09/2019", format = "%d/%m/%Y")] 
  n2 <- nrow(practice_df[[i]])
  print(paste0("drop lcd check: ", (n1-n2)))
}

#appending all the practice files together
all_practices <- practice_df[[1]]
for(i in 2:length(practice_df)){
  all_practices <- rbind(all_practices, practice_df[[i]])
}

nrow(all_practices) == (sum(nrow(practice_df[[1]]), nrow(practice_df[[2]]),
                            nrow(practice_df[[3]])))

#remove duplicates 
all_practices <- unique(all_practices)

#merging with the patient files 
setkey(patient_population, pracid)
setkey(all_practices, pracid)

study_pop <- patient_population[all_practices]

population <- study_pop[, list(patid, pracid, gender, yob, regstartdate,
                                regenddate, cprd_ddate, dob, lcd, region)]

#bring in death date from ons
setwd(raw_data_linked)
ons_dt_1 <- readr::read_tsv("death_patient_22_002202_DM.txt", col_type = cols(.default = col_character()))

#format dod as date
ons_dt_1 <- as.data.table(ons_dt_1)
ons_dt_1[, dod := as.Date(dod, format = "%d/%m/%Y")]

#merge death with patient file
population[ons_dt_1, `:=`(dod = i.dod), on = "patid"]

#generate start of LB and end of follow-up
population[, start_lb := max(dob, regstartdate, na.rm = T), by = "patid"]
population[, end_fu := min(lcd, dod, regenddate, as.Date("31/03/2021", format = "%d/%m/%Y"), na.rm = T), by ="patid"]

#check if any individuals have start of LB before end of FU
n <- population[start_lb <= end_fu]

#then flag individuals identified in each population
population[, covidVE := regstartdate <= as.Date("08/12/2019", format = "%d/%m/%Y") & (regenddate >= as.Date("08/12/2020", format = "%d/%m/%Y") | is.na(regenddate))]
population[, influenzaVE := regstartdate <= as.Date("01/09/2018", format = "%d/%m/%Y") & (regenddate >= as.Date("01/09/2019", format = "%d/%m/%Y") | is.na(regenddate))]
population[, negcontrolVE := regstartdate <= as.Date("01/01/2019", format = "%d/%m/%Y") & (regenddate >= as.Date("01/01/2020", format = "%d/%m/%Y") | is.na(regenddate))]

#count of people in each pop
n <- population[covidVE == "TRUE"]
n <- population[influenzaVE == "TRUE"]
n <- population[negcontrolVE == "TRUE"]

N<-population[covidVE == "TRUE" & influenzaVE == "TRUE" & negcontrolVE == "TRUE"]

#identify patients with HES linkage
setwd("xxx")
HES_linkage <- read.table("Aurum_enhanced_eligibility_January_2022.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
HES_linkage <- as.data.table(HES_linkage)

#merge to get hes linkage 
population[HES_linkage, `:=`(hes_apc_e = i.hes_apc_e, ons_death_e = i.ons_death_e, lsoa_e = i.lsoa_e), on = "patid"]

#save this as parquet file
setwd(datafiles)
write_parquet(population, 
              paste0(datafiles, "study_population.parquet"))

uniqueN(population[!is.na(patid), patid]) 

#now need to create analysis specific populations
covid_population <- population[covidVE == TRUE]
flu_population <- population[influenzaVE == TRUE]
neg_population <- population[negcontrolVE == TRUE]

#apply population specific criteria - covid
covid_population <- covid_population[is.na(regenddate) | regenddate >= as.Date("2020/12/08")]
covid_population <- covid_population[is.na(dod) | dod >= as.Date("2020/12/08")]
covid_population <- covid_population[is.na(cprd_ddate) | cprd_ddate >= as.Date("2020/12/08")]
covid_population <- covid_population[regstartdate < as.Date("2019/12/08")]
covid_population <- covid_population[is.na(lcd) | lcd >= as.Date("2020/12/08")]

#flu
flu_population <- flu_population[is.na(regenddate) | regenddate >= as.Date("2019/09/01")]
flu_population <- flu_population[is.na(dod) | dod >= as.Date("2019/09/01")]
flu_population <- flu_population[is.na(cprd_ddate) | cprd_ddate >= as.Date("2019/09/01")]
flu_population <- flu_population[regstartdate < as.Date("2018/09/01")]
flu_population <- flu_population[is.na(lcd) | lcd >= as.Date("2019/09/01")]

#neg
neg_population <- neg_population[is.na(regenddate) | regenddate >= as.Date("2020/01/01")]
neg_population <- neg_population[is.na(dod) | dod >= as.Date("2020/01/01")]
neg_population <- neg_population[is.na(cprd_ddate) | cprd_ddate >= as.Date("2020/01/01")]
neg_population <- neg_population[regstartdate < as.Date("2019/01/01")]
neg_population <- neg_population[is.na(lcd) | lcd >= as.Date("2020/01/01")]

#save these 
setwd(datafiles)
write_parquet(covid_population, 
              paste0(datafiles, "study_population_covid.parquet"))
write_parquet(flu_population, 
              paste0(datafiles, "study_population_flu.parquet"))
write_parquet(neg_population, 
              paste0(datafiles, "study_population_neg.parquet"))