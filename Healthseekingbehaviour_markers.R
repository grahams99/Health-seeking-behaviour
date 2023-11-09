################################################################################
# Author: Sophie Graham
# Date: 19/06/2023
# Version: R 4.2.2
# File name: healthseekingbehaviour_markers.R
# Status: in progress
# CPRD version: May 2022
# Data sets used: CPRD Aurum and HES APC
# R scripts needed: global (not included as has links to secure severs), data 
# prep new (not included as merging code lists with raw data files) and population.
# Data sets created:
# Description of file: this file creates all the markers for all the populations
# and then merges with the patient file - creating one for each analysis 
# population. 
################################################################################

#Read in part 1 medcode markers - AAA, breast cancer, bowel cancer, cervical screen, pneu vacc, PSA, bone scans and NHS health checks  
setwd(datafiles)
markers_dt_1 <- read_parquet("observations_markers_pt1_tmp1.parquet")
markers_dt_2 <- read_parquet("observations_markers_pt1_tmp2.parquet")
markers_dt_3 <- read_parquet("observations_markers_pt1_tmp3.parquet")

#remove marker_num as we now have separate lists for narrow and broad
markers_dt_1[, marker_num := NULL]
markers_dt_2[, marker_num := NULL]
markers_dt_3[, marker_num := NULL]

#then add in marker num again - load code list
setwd(codelists)

#read in code lists from multiple sheets  
markers <- map(set_names(excel_sheets("Markers_medcodes_part1.xlsx"), c("Aurum_AAA_narrow", "Aurum_AAA_broad", "Aurum_breast_screen_narrow", "Aurum_breast_screen_broad", "Aurum_bowel_screen", "Aurum_cervical_screen_narrow", "Aurum_cervical_screen_broad", "Aurum_pneu_vacc", "Aurum_PSA", "Aurum_bone_scan", "Aurum_NHS_checks")),
               read_excel, path = "Markers_medcodes_part1.xlsx"
)

#merge code lists 
for (i in markers){
  allmarkers <- rbindlist(markers, use.names = TRUE, idcol = TRUE, fill=TRUE)
}
rm(i)

# Create a subset of "allmarkers" with only "medcodeid" and "marker_num" columns
subset_allmarkers_dt <- allmarkers[, .(medcodeid, marker_num)]

# Set the key column for efficient merging
setkey(markers_dt_1, medcodeid)
setkey(markers_dt_2, medcodeid)
setkey(markers_dt_3, medcodeid)
setkey(subset_allmarkers_dt, medcodeid)

# Merge the data tables and allow for multiple rows per marker_num
markers_dt_1 <- subset_allmarkers_dt[markers_dt_1, allow.cartesian = TRUE]
markers_dt_2 <- subset_allmarkers_dt[markers_dt_2, allow.cartesian = TRUE]
markers_dt_3 <- subset_allmarkers_dt[markers_dt_3, allow.cartesian = TRUE]
rm(subset_allmarkers_dt)

#write this for later use
data_tables <- list (markers_dt_1, markers_dt_2, markers_dt_3)
write_data <- function(data, datafiles) {
  iteration <- seq(length(data)) + 0.1
  for (i in seq_along(data)) {
    write_parquet(data[[i]], paste0(datafiles, "/observations_markers_pt1_tmp", iteration[i], ".parquet"))
  }
}
write_data(data_tables, datafiles)
rm(write_data, data_tables)

#read in the data
setwd(datafiles)
markers_dt_1 <- read_parquet("observations_markers_pt1_tmp1.1.parquet")
markers_dt_2 <- read_parquet("observations_markers_pt1_tmp2.1.parquet")
markers_dt_3 <- read_parquet("observations_markers_pt1_tmp3.1.parquet")

################ALL LOOK BACK################
#define those that are all lookback - AAA 
obs_dt_1 <- markers_dt_1[marker_num == 1 | marker_num == 15]
obs_dt_2 <- markers_dt_2[marker_num == 1 | marker_num == 15]
obs_dt_3 <- markers_dt_3[marker_num == 1 | marker_num == 15]

#deduplicate keeping only the earliest 
process_dataset <- function(dataset) {
  dataset <- dataset[order(patid)]
  dataset[, duplicates := .N > 1, by = c("marker_num", "patid")]
  dataset[duplicates == TRUE, first_entry := lapply(.SD, min), by = c("patid", "marker_num"), .SDcols = "eventdate"]
  dataset[duplicates == TRUE, keep_row := eventdate == first_entry]
  dataset[duplicates == FALSE, keep_row := TRUE]
  dataset <- dataset[keep_row == TRUE]
  dataset[, c("duplicates", "first_entry", "keep_row") := NULL]
  return(dataset)
}

# Apply the function to each data set
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#add in analysis population from patient file
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
setkey(patient, patid)

#add info from patient file into observational files 
obs_dt_1[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
obs_dt_2[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
obs_dt_3[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
rm(patient)

#then apply condition for each 
apply_conditions <- function(data) {
  data[covidVE == TRUE & marker_num == 1, COVIDVE_AAA_narrow := eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  data[covidVE == TRUE & marker_num == 15, COVIDVE_AAA_broad := eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  data[influenzaVE == TRUE & marker_num == 1, fluVE_AAA_narrow := eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  data[influenzaVE == TRUE & marker_num == 15, fluVE_AAA_broad := eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  data[negcontrolVE == TRUE & marker_num == 1, negcon_AAA_narrow := eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  data[negcontrolVE == TRUE & marker_num == 15, negcon_AAA_broad := eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  data <- data[COVIDVE_AAA_narrow == TRUE | COVIDVE_AAA_broad == TRUE | fluVE_AAA_narrow == TRUE | fluVE_AAA_broad == TRUE | negcon_AAA_narrow == TRUE | negcon_AAA_broad == TRUE]
  data[, `:=` (obsid = NULL, obsdate = NULL, marker_num = NULL, enterdate = NULL, obstypeid = NULL, medcodeid = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL)]
  return(data)
}

# Apply the function to each dataset
obs_dt_1 <- apply_conditions(obs_dt_1)
obs_dt_2 <- apply_conditions(obs_dt_2)
obs_dt_3 <- apply_conditions(obs_dt_3)

#unique columns 
obs_dt_1 <- unique(obs_dt_1)
obs_dt_2 <- unique(obs_dt_2)
obs_dt_3 <- unique(obs_dt_3)
#length(unique(obs_dt_1$patid))

#combine into one
aa_dt_all <-rbind(obs_dt_1, obs_dt_2, obs_dt_3)
#length(unique(aa_dt_all$patid))

#first remove event date and make all columns character 
col_char <- function(data){
  data[, eventdate := NULL]
  data[, COVIDVE_AAA_narrow := as.character(COVIDVE_AAA_narrow)]
  data[, COVIDVE_AAA_broad := as.character(COVIDVE_AAA_broad)]
  data[, fluVE_AAA_narrow := as.character(fluVE_AAA_narrow)]
  data[, fluVE_AAA_broad := as.character(fluVE_AAA_broad)]
  data[, negcon_AAA_narrow := as.character(negcon_AAA_narrow)]
  data[, negcon_AAA_broad := as.character(negcon_AAA_broad)]
}
aa_dt_all <- col_char(aa_dt_all)

#then we want to collapse into one row per patid
collapse_rows <- function(x) {
  if (any(x == "TRUE", na.rm = TRUE)) {
    return(TRUE)
  } else {
    return(NA)
  }
}
aa_dt_all <- aa_dt_all[, lapply(.SD, collapse_rows), by = patid]

#then want to change the format of the data
aa_dt_all[, COVIDVE_AAA_narrow := ifelse(COVIDVE_AAA_narrow == TRUE, 1, 0)]
aa_dt_all[, COVIDVE_AAA_broad := ifelse(COVIDVE_AAA_broad == TRUE, 1, 0)]
aa_dt_all[, fluVE_AAA_narrow  := ifelse(fluVE_AAA_narrow  == TRUE, 1, 0)]
aa_dt_all[, fluVE_AAA_broad := ifelse(fluVE_AAA_broad == TRUE, 1, 0)]
aa_dt_all[, negcon_AAA_narrow := ifelse(negcon_AAA_narrow == TRUE, 1, 0)]
aa_dt_all[, negcon_AAA_broad := ifelse(negcon_AAA_broad == TRUE, 1, 0)]

#NA to zero
aa_dt_all[is.na(aa_dt_all)] <- 0

#write this for later use
setwd(datafiles)
write_parquet(aa_dt_all, paste0(datafiles, "aa_markers.parquet"))
############################################################

################SHORT LOOK BACK################ 
#short look back conditions - breast (narrow and broad), bowel, cervical (narrow and broad), PSA, bone scan and nhs health checks 
obs_dt_1 <- markers_dt_1[marker_num == 2 | marker_num == 3 | marker_num == 4 | marker_num == 6 | marker_num == 7 | marker_num == 8| marker_num == 16 | marker_num == 17]#17,143,155
obs_dt_2 <- markers_dt_2[marker_num == 2 | marker_num == 3 | marker_num == 4 | marker_num == 6 | marker_num == 7 | marker_num == 8| marker_num == 16 | marker_num == 17]#17,152,804
obs_dt_3 <- markers_dt_3[marker_num == 2 | marker_num == 3 | marker_num == 4 | marker_num == 6 | marker_num == 7 | marker_num == 8| marker_num == 16 | marker_num == 17]#3,221,974

#remove events that occur after 8 dec 2020 
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  return(dt)
}
obs_dt_1 <- remove_rows(obs_dt_1)
obs_dt_2 <- remove_rows(obs_dt_2)
obs_dt_3 <- remove_rows(obs_dt_3)
rm(remove_rows)

#add in analysis population from patient file
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
setkey(patient, patid)

#add info from patient file into observational files 
obs_dt_1[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
obs_dt_2[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
obs_dt_3[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
rm(patient)

#save this for later use
setwd(datafiles)
data_tables <- list (obs_dt_1 , obs_dt_2, obs_dt_3)

write_data <- function(data, datafiles) {
  for (i in seq_along(data)) {
    write_parquet(data[[i]], paste0(datafiles, "/observations_markers_short_lookback_tmp", i, ".parquet"))
  }
}

write_data(data_tables, datafiles)
rm(write_data, data_tables) 

#if starting from here - load 
#obs_dt_1 <- read_parquet("observations_markers_short_lookback_tmp1.parquet")
#obs_dt_2 <- read_parquet("observations_markers_short_lookback_tmp2.parquet")
#obs_dt_3 <- read_parquet("observations_markers_short_lookback_tmp3.parquet")

#deduplicate each patient for each marker num keeping only the latest event 
process_dataset <- function(dataset) {
  dataset <- dataset[order(patid)]
  dataset[, duplicates := .N > 1, by = c("marker_num", "patid")]
  dataset[duplicates == TRUE, last_entry := lapply(.SD, max), by = c("patid", "marker_num"), .SDcols = "eventdate"]
  dataset[duplicates == TRUE, keep_row := eventdate == last_entry]
  dataset[duplicates == FALSE, keep_row := TRUE]
  dataset <- dataset[keep_row == TRUE]
  dataset[, duplicates := NULL]
  dataset[, last_entry := NULL]
  dataset[, keep_row := NULL]
  return(dataset)
}
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#covid - original definition
#process_covid_data <- function(dataset) {
#  dataset[covidVE == TRUE & marker_num == 2, COVIDVE_breast_narrow := eventdate >= date65 & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
#  dataset[covidVE == TRUE & marker_num == 16, COVIDVE_breast_broad := eventdate >= date65 & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
#  dataset[covidVE == TRUE & marker_num == 3, COVIDVE_bowel := eventdate >= as.Date("01/09/2017", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
#  dataset[covidVE == TRUE & marker_num == 4, COVIDVE_cervical_narrow := eventdate >= date54 & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[covidVE == TRUE & marker_num == 17, COVIDVE_cervical_broad := eventdate >= date54 & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[covidVE == TRUE & marker_num == 6, COVIDVE_PSA := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
#  dataset[covidVE == TRUE & marker_num == 7, COVIDVE_bone := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
#  dataset[covidVE == TRUE & marker_num == 8, COVIDVE_NHS := eventdate >= date64 & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
#  dataset <- dataset[COVIDVE_breast_narrow == TRUE | COVIDVE_breast_broad == TRUE | COVIDVE_bowel == TRUE | COVIDVE_cervical_narrow == TRUE | COVIDVE_cervical_broad == TRUE | COVIDVE_PSA == TRUE | COVIDVE_bone == TRUE | COVIDVE_NHS == TRUE]
#  dataset[, `:=` (obsid = NULL, obsdate = NULL, enterdate = NULL, eventdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL, date54 = NULL, date65 = NULL, date64 = NULL)]
#  return(dataset)
#}
#covid_dt_1 <- process_covid_data(obs_dt_1)
#covid_dt_2 <- process_covid_data(obs_dt_2)
#covid_dt_3 <- process_covid_data(obs_dt_3)

#then run for COVID - permissive 
process_covid_data <- function(dataset) {
  dataset[, dob := as.Date(dob, format = "%d/%m/%Y")]
  dataset[, age := as.numeric(difftime(as.Date("08/12/2020", format = "%d/%m/%Y"), dob, units = "days")) / 365.25]
  eligible_age_end_breast <- 71
  dataset[, min_end_breast := pmin(eligible_age_end_breast, age)]
  dataset[, last_4_years_end := as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset[, last_4_years_start := last_4_years_end - 1461]
  dataset[min_end_breast == eligible_age_end_breast, last_4_years_start := dob+(eligible_age_end_breast*365.25) - 1461]
  dataset[, last_4_years_end := as.IDate(last_4_years_end, format = "%d/%m/%Y")]
  dataset[, last_4_years_start := as.IDate(last_4_years_start, format = "%d/%m/%Y")]
  eligible_age_end_cervical <- 64
  dataset[, min_end_cervical := pmin(eligible_age_end_cervical, age)]
  dataset[, last_6_years_end := as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset[, last_6_years_start := last_6_years_end - 1461]
  dataset[min_end_cervical == eligible_age_end_cervical, last_6_years_start := dob+(eligible_age_end_cervical*365.25) - 1461]
  dataset[, last_6_years_end := as.IDate(last_6_years_end, format = "%d/%m/%Y")]
  dataset[, last_6_years_start := as.IDate(last_6_years_start, format = "%d/%m/%Y")]
  eligible_age_end <- 74
  dataset[, min_end_bowel := pmin(eligible_age_end, age)]
  dataset[, last_3_years_end := as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset[, last_3_years_start := last_3_years_end - 1095.75]
  dataset[min_end_bowel == eligible_age_end, last_3_years_start := dob+(eligible_age_end*365.25) - 1095.75]
  dataset[, last_3_years_end := as.IDate(last_3_years_end, format = "%d/%m/%Y")]
  dataset[, last_3_years_start := as.IDate(last_3_years_start, format = "%d/%m/%Y")]
  dataset[, min_end_NHS := pmin(eligible_age_end, age)]
  dataset[, last_6_years_end_NHS := as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset[, last_6_years_start_NHS := last_6_years_end_NHS - 2191.5]
  dataset[min_end_NHS == eligible_age_end, last_6_years_start_NHS := dob+(eligible_age_end*365.25) - 2191.5]
  dataset[, last_6_years_end_NHS := as.IDate(last_6_years_end_NHS, format = "%d/%m/%Y")]
  dataset[, last_6_years_start_NHS := as.IDate(last_6_years_start_NHS, format = "%d/%m/%Y")]
  dataset <- dataset[covidVE == TRUE & (marker_num == 2 | marker_num == 16 | marker_num == 3 | marker_num == 4 | marker_num == 17 | marker_num == 6 | marker_num == 7 | marker_num == 8)]
  dataset[marker_num == 2, COVIDVE_breast_narrow := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 16, COVIDVE_breast_broad := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 3, COVIDVE_bowel := eventdate >= last_3_years_start & eventdate <= last_3_years_end, by = .(patid)]
  dataset[marker_num == 4, COVIDVE_cervical_narrow := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 17, COVIDVE_cervical_broad := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 6, COVIDVE_PSA := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset[marker_num == 7, COVIDVE_bone := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset[marker_num == 8, COVIDVE_NHS := eventdate >= last_6_years_start_NHS & eventdate <= last_6_years_end_NHS, by = .(patid)]
  dataset <- dataset[COVIDVE_breast_narrow == TRUE | COVIDVE_breast_broad == TRUE | COVIDVE_bowel == TRUE | COVIDVE_cervical_narrow == TRUE | COVIDVE_cervical_broad == TRUE | COVIDVE_PSA == TRUE | COVIDVE_bone == TRUE | COVIDVE_NHS == TRUE]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, enterdate = NULL, eventdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL, age = NULL, min_end_breast = NULL, min_end_bowel = NULL, min_end_NHS = NULL, last_4_years_end = NULL, last_4_years_start = NULL, last_3_years_end = NULL, last_3_years_start = NULL, last_6_years_end = NULL, last_6_years_start = NULL, min_end_cervical = NULL, last_6_years_end_NHS = NULL, last_6_years_start_NHS = NULL)]
  return(dataset)
}
covid_dt_1 <- process_covid_data(obs_dt_1)
covid_dt_2 <- process_covid_data(obs_dt_2)
covid_dt_3 <- process_covid_data(obs_dt_3)

#unique columns 
covid_dt_1 <- unique(covid_dt_1)
covid_dt_2 <- unique(covid_dt_2)
covid_dt_3 <- unique(covid_dt_3)

#then turn cols to char
col_char <- function(data){
  data[, COVIDVE_breast_narrow := as.character(COVIDVE_breast_narrow)]
  data[, COVIDVE_breast_broad := as.character(COVIDVE_breast_broad)]
  data[, COVIDVE_bowel := as.character(COVIDVE_bowel)]
  data[, COVIDVE_cervical_narrow := as.character(COVIDVE_cervical_narrow)]
  data[, COVIDVE_cervical_broad := as.character(COVIDVE_cervical_broad)]
  data[, COVIDVE_PSA := as.character(COVIDVE_PSA)]
  data[, COVIDVE_bone := as.character(COVIDVE_bone)]
  data[, COVIDVE_NHS := as.character(COVIDVE_NHS)]
}
covid_dt_1 <- col_char(covid_dt_1)
covid_dt_2 <- col_char(covid_dt_2)
covid_dt_3 <- col_char(covid_dt_3)

#then create one row per patid
covid_dt_1 <- covid_dt_1[, lapply(.SD, collapse_rows), by = patid]
#length(unique(covid_dt_1$patid))
covid_dt_2 <- covid_dt_2[, lapply(.SD, collapse_rows), by = patid]
#length(unique(covid_dt_2$patid))
covid_dt_3 <- covid_dt_3[, lapply(.SD, collapse_rows), by = patid]
#length(unique(covid_dt_3$patid))

#change the format of the data
format_data <- function(data){
  data[, COVIDVE_breast_narrow := ifelse(COVIDVE_breast_narrow == TRUE, 1, 0)]
  data[, COVIDVE_breast_broad := ifelse(COVIDVE_breast_broad == TRUE, 1, 0)]
  data[, COVIDVE_bowel  := ifelse(COVIDVE_bowel  == TRUE, 1, 0)]
  data[, COVIDVE_cervical_narrow := ifelse(COVIDVE_cervical_narrow == TRUE, 1, 0)]
  data[, COVIDVE_cervical_broad := ifelse(COVIDVE_cervical_broad == TRUE, 1, 0)]
  data[, COVIDVE_PSA := ifelse(COVIDVE_PSA == TRUE, 1, 0)]
  data[, COVIDVE_bone := ifelse(COVIDVE_bone == TRUE, 1, 0)]
  data[, COVIDVE_NHS := ifelse(COVIDVE_NHS == TRUE, 1, 0)]
  data[is.na(data)] <- 0
  return(data)
}
covid_dt_1 <- format_data(covid_dt_1)
covid_dt_2 <- format_data(covid_dt_2)
covid_dt_3 <- format_data(covid_dt_3)

#rbind for covid
covid_dt_all <- rbind(covid_dt_1, covid_dt_2, covid_dt_3)
rm(covid_dt_1, covid_dt_2, covid_dt_3)

#now do for negative control
obs_dt_1 <- read_parquet("observations_markers_short_lookback_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_short_lookback_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_short_lookback_tmp3.parquet")

#remove events that occur after 1 jan 2020
obs_dt_1 <- obs_dt_1[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
obs_dt_2 <- obs_dt_2[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
obs_dt_3 <- obs_dt_3[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]

#then run process data function 
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#neg - original definition
#process_neg_data <- function(dataset) {
#  dataset[negcontrolVE == TRUE & marker_num == 2, negVE_breast_narrow := eventdate >= date65 & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[negcontrolVE == TRUE & marker_num == 16, negVE_breast_broad := eventdate >= date65 & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[negcontrolVE == TRUE & marker_num == 3, negVE_bowel := eventdate >= as.Date("01/09/2017", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[negcontrolVE == TRUE & marker_num == 4, negVE_cervical_narrow := eventdate >= date54 & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[negcontrolVE == TRUE & marker_num == 17, negVE_cervical_broad := eventdate >= date54 & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[negcontrolVE == TRUE & marker_num == 6, negVE_PSA := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[negcontrolVE == TRUE & marker_num == 7, negVE_bone := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[negcontrolVE == TRUE & marker_num == 8, negVE_NHS := eventdate >= date64 & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset <- dataset[negVE_breast_narrow == TRUE | negVE_breast_broad == TRUE | negVE_bowel == TRUE | negVE_cervical_narrow == TRUE | negVE_cervical_broad == TRUE | negVE_PSA == TRUE | negVE_bone == TRUE | negVE_NHS == TRUE]
#  dataset[, `:=` (obsid = NULL, obsdate = NULL, enterdate = NULL, eventdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL, date54 = NULL, date65 = NULL, date64 = NULL)]
#  return(dataset)
#}
#neg_dt_1 <- process_neg_data(obs_dt_1)
#neg_dt_2 <- process_neg_data(obs_dt_2)
#neg_dt_3 <- process_neg_data(obs_dt_3)

#neg - permissive 
process_neg_data <- function(dataset) {
  dataset[, dob := as.Date(dob, format = "%d/%m/%Y")]
  dataset[, age := as.numeric(difftime(as.Date("01/01/2020", format = "%d/%m/%Y"), dob, units = "days")) / 365.25]
  eligible_age_end_breast <- 71
  dataset[, min_end_breast := pmin(eligible_age_end_breast, age)]
  dataset[, last_4_years_end := as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset[, last_4_years_start := last_4_years_end - 1461]
  dataset[min_end_breast == eligible_age_end_breast, last_4_years_start := dob+(eligible_age_end_breast*365.25) - 1461]
  dataset[, last_4_years_end := as.IDate(last_4_years_end, format = "%d/%m/%Y")]
  dataset[, last_4_years_start := as.IDate(last_4_years_start, format = "%d/%m/%Y")]
  eligible_age_end_cervical <- 64
  dataset[, min_end_cervical := pmin(eligible_age_end_cervical, age)]
  dataset[, last_6_years_end := as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset[, last_6_years_start := last_6_years_end - 1461]
  dataset[min_end_cervical == eligible_age_end_cervical, last_6_years_start := dob+(eligible_age_end_cervical*365.25) - 1461]
  dataset[, last_6_years_end := as.IDate(last_6_years_end, format = "%d/%m/%Y")]
  dataset[, last_6_years_start := as.IDate(last_6_years_start, format = "%d/%m/%Y")]
  eligible_age_end <- 74
  dataset[, min_end_bowel := pmin(eligible_age_end, age)]
  dataset[, last_3_years_end := as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset[, last_3_years_start := last_3_years_end - 1095.75]
  dataset[min_end_bowel == eligible_age_end, last_3_years_start := dob+(eligible_age_end*365.25) - 1095.75]
  dataset[, last_3_years_end := as.IDate(last_3_years_end, format = "%d/%m/%Y")]
  dataset[, last_3_years_start := as.IDate(last_3_years_start, format = "%d/%m/%Y")]
  dataset[, min_end_NHS := pmin(eligible_age_end, age)]
  dataset[, last_6_years_end_NHS := as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset[, last_6_years_start_NHS := last_6_years_end_NHS - 2191.5]
  dataset[min_end_NHS == eligible_age_end, last_6_years_start_NHS := dob+(eligible_age_end*365.25) - 2191.5]
  dataset[, last_6_years_end_NHS := as.IDate(last_6_years_end_NHS, format = "%d/%m/%Y")]
  dataset[, last_6_years_start_NHS := as.IDate(last_6_years_start_NHS, format = "%d/%m/%Y")]
  dataset <- dataset[negcontrolVE == TRUE & (marker_num == 2 | marker_num == 16 | marker_num == 3 | marker_num == 4 | marker_num == 17 | marker_num == 6 | marker_num == 7 | marker_num == 8)]
  dataset[marker_num == 2, negVE_breast_narrow := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 16, negVE_breast_broad := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 3, negVE_bowel := eventdate >= last_3_years_start & eventdate <= last_3_years_end, by = .(patid)]
  dataset[marker_num == 4, negVE_cervical_narrow := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 17, negVE_cervical_broad := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 6, negVE_PSA := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset[marker_num == 7, negVE_bone := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset[marker_num == 8, negVE_NHS := eventdate >= last_6_years_start_NHS & eventdate <= last_6_years_end_NHS, by = .(patid)]
  dataset <- dataset[negVE_breast_narrow == TRUE | negVE_breast_broad == TRUE | negVE_bowel == TRUE | negVE_cervical_narrow == TRUE | negVE_cervical_broad == TRUE | negVE_PSA == TRUE | negVE_bone == TRUE | negVE_NHS == TRUE]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, enterdate = NULL, eventdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL, age = NULL, min_end_breast = NULL, min_end_bowel = NULL, min_end_NHS = NULL, last_4_years_end = NULL, last_4_years_start = NULL, last_3_years_end = NULL, last_3_years_start = NULL, last_6_years_end = NULL, last_6_years_start = NULL, min_end_cervical = NULL, last_6_years_end_NHS = NULL, last_6_years_start_NHS = NULL)]
  return(dataset)
}
neg_dt_1 <- process_neg_data(obs_dt_1)
neg_dt_2 <- process_neg_data(obs_dt_2)
neg_dt_3 <- process_neg_data(obs_dt_3)

#unique
neg_dt_1 <- unique(neg_dt_1)
neg_dt_2 <- unique(neg_dt_2)
neg_dt_3 <- unique(neg_dt_3)

#then turn cols to char
col_char <- function(data){
  data[, negVE_breast_narrow := as.character(negVE_breast_narrow)]
  data[, negVE_breast_broad := as.character(negVE_breast_broad)]
  data[, negVE_bowel := as.character(negVE_bowel)]
  data[, negVE_cervical_narrow := as.character(negVE_cervical_narrow)]
  data[, negVE_cervical_broad := as.character(negVE_cervical_broad)]
  data[, negVE_PSA := as.character(negVE_PSA)]
  data[, negVE_bone := as.character(negVE_bone)]
  data[, negVE_NHS := as.character(negVE_NHS)]
}
neg_dt_1 <- col_char(neg_dt_1)
neg_dt_2 <- col_char(neg_dt_2)
neg_dt_3 <- col_char(neg_dt_3)

#then create one row per patid
neg_dt_1 <- neg_dt_1[, lapply(.SD, collapse_rows), by = patid]
#length(unique(neg_dt_1$patid))
neg_dt_2 <- neg_dt_2[, lapply(.SD, collapse_rows), by = patid]
#length(unique(neg_dt_2$patid))
neg_dt_3 <- neg_dt_3[, lapply(.SD, collapse_rows), by = patid]
#length(unique(neg_dt_3$patid))

#change the format of the data
format_data <- function(data){
  data[, negVE_breast_narrow := ifelse(negVE_breast_narrow == TRUE, 1, 0)]
  data[, negVE_breast_broad := ifelse(negVE_breast_broad == TRUE, 1, 0)]
  data[, negVE_bowel  := ifelse(negVE_bowel  == TRUE, 1, 0)]
  data[, negVE_cervical_narrow := ifelse(negVE_cervical_narrow == TRUE, 1, 0)]
  data[, negVE_cervical_broad := ifelse(negVE_cervical_broad == TRUE, 1, 0)]
  data[, negVE_PSA := ifelse(negVE_PSA == TRUE, 1, 0)]
  data[, negVE_bone := ifelse(negVE_bone == TRUE, 1, 0)]
  data[, negVE_NHS := ifelse(negVE_NHS == TRUE, 1, 0)]
  data[is.na(data)] <- 0
  return(data)
}
neg_dt_1 <- format_data(neg_dt_1)
neg_dt_2 <- format_data(neg_dt_2)
neg_dt_3 <- format_data(neg_dt_3)

#rbind for neg
neg_dt_all <- rbind(neg_dt_1, neg_dt_2, neg_dt_3)
rm(neg_dt_1, neg_dt_2, neg_dt_3)

#now do for flu
obs_dt_1 <- read_parquet("observations_markers_short_lookback_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_short_lookback_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_short_lookback_tmp3.parquet")

#restrict to observations before 1 september 2019
obs_dt_1 <- obs_dt_1[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
obs_dt_2 <- obs_dt_2[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
obs_dt_3 <- obs_dt_3[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]

#then run process data function 
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#flu - original definition
#process_flu_data <- function(dataset) {
#  dataset[influenzaVE == TRUE & marker_num == 2, fluVE_breast_narrow := eventdate >= date65 & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
#  dataset[influenzaVE == TRUE & marker_num == 16, fluVE_breast_broad := eventdate >= date65 & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
#  dataset[influenzaVE == TRUE & marker_num == 3, fluVE_bowel := eventdate >= as.Date("01/09/2017", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
#  dataset[influenzaVE == TRUE & marker_num == 4, fluVE_cervical_narrow := eventdate >= date54 & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[influenzaVE == TRUE & marker_num == 17, fluVE_cervical_broad := eventdate >= date54 & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
#  dataset[influenzaVE == TRUE & marker_num == 6, fluVE_PSA := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
#  dataset[influenzaVE == TRUE & marker_num == 7, fluVE_bone := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
#  dataset[influenzaVE == TRUE & marker_num == 8, fluVE_NHS := eventdate >= date64 & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
#  dataset <- dataset[fluVE_breast_narrow == TRUE | fluVE_breast_broad == TRUE | fluVE_bowel == TRUE | fluVE_cervical_narrow == TRUE | fluVE_cervical_broad == TRUE | fluVE_PSA == TRUE | fluVE_bone == TRUE | fluVE_NHS == TRUE]
#  dataset[, `:=` (obsid = NULL, obsdate = NULL, enterdate = NULL, eventdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL, date54 = NULL, date65 = NULL, date64 = NULL)]
#  return(dataset)
#}
#flu_dt_1 <- process_flu_data(obs_dt_1)
#flu_dt_2 <- process_flu_data(obs_dt_2)
#flu_dt_3 <- process_flu_data(obs_dt_3)

#flu - permissive 
process_flu_data <- function(dataset) {
  dataset[, dob := as.Date(dob, format = "%d/%m/%Y")]
  dataset[, age := as.numeric(difftime(as.Date("01/09/2019", format = "%d/%m/%Y"), dob, units = "days")) / 365.25]
  eligible_age_end_breast <- 71
  dataset[, min_end_breast := pmin(eligible_age_end_breast, age)]
  dataset[, last_4_years_end := as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset[, last_4_years_start := last_4_years_end - 1461]
  dataset[min_end_breast == eligible_age_end_breast, last_4_years_start := dob+(eligible_age_end_breast*365.25) - 1461]
  dataset[, last_4_years_end := as.IDate(last_4_years_end, format = "%d/%m/%Y")]
  dataset[, last_4_years_start := as.IDate(last_4_years_start, format = "%d/%m/%Y")]
  eligible_age_end_cervical <- 64
  dataset[, min_end_cervical := pmin(eligible_age_end_cervical, age)]
  dataset[, last_6_years_end := as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset[, last_6_years_start := last_6_years_end - 1461]
  dataset[min_end_cervical == eligible_age_end_cervical, last_6_years_start := dob+(eligible_age_end_cervical*365.25) - 1461]
  dataset[, last_6_years_end := as.IDate(last_6_years_end, format = "%d/%m/%Y")]
  dataset[, last_6_years_start := as.IDate(last_6_years_start, format = "%d/%m/%Y")]
  eligible_age_end <- 74
  dataset[, min_end_bowel := pmin(eligible_age_end, age)]
  dataset[, last_3_years_end := as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset[, last_3_years_start := last_3_years_end - 1095.75]
  dataset[min_end_bowel == eligible_age_end, last_3_years_start := dob+(eligible_age_end*365.25) - 1095.75]
  dataset[, last_3_years_end := as.IDate(last_3_years_end, format = "%d/%m/%Y")]
  dataset[, last_3_years_start := as.IDate(last_3_years_start, format = "%d/%m/%Y")]
  dataset[, min_end_NHS := pmin(eligible_age_end, age)]
  dataset[, last_6_years_end_NHS := as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset[, last_6_years_start_NHS := last_6_years_end_NHS - 2191.5]
  dataset[min_end_NHS == eligible_age_end, last_6_years_start_NHS := dob+(eligible_age_end*365.25) - 2191.5]
  dataset[, last_6_years_end_NHS := as.IDate(last_6_years_end_NHS, format = "%d/%m/%Y")]
  dataset[, last_6_years_start_NHS := as.IDate(last_6_years_start_NHS, format = "%d/%m/%Y")]
  dataset <- dataset[influenzaVE == TRUE & (marker_num == 2 | marker_num == 16 | marker_num == 3 | marker_num == 4 | marker_num == 17 | marker_num == 6 | marker_num == 7 | marker_num == 8)]
  dataset[marker_num == 2, fluVE_breast_narrow := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 16, fluVE_breast_broad := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 3, fluVE_bowel := eventdate >= last_3_years_start & eventdate <= last_3_years_end, by = .(patid)]
  dataset[marker_num == 4, fluVE_cervical_narrow := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 17, fluVE_cervical_broad := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 6, fluVE_PSA := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset[marker_num == 7, fluVE_bone := eventdate >= as.Date("01/09/2016", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset[marker_num == 8, fluVE_NHS := eventdate >= last_6_years_start_NHS & eventdate <= last_6_years_end_NHS, by = .(patid)]
  dataset <- dataset[fluVE_breast_narrow == TRUE | fluVE_breast_broad == TRUE | fluVE_bowel == TRUE | fluVE_cervical_narrow == TRUE | fluVE_cervical_broad == TRUE | fluVE_PSA == TRUE | fluVE_bone == TRUE | fluVE_NHS == TRUE]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, enterdate = NULL, eventdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL, age = NULL, min_end_breast = NULL, min_end_bowel = NULL, min_end_NHS = NULL, last_4_years_end = NULL, last_4_years_start = NULL, last_3_years_end = NULL, last_3_years_start = NULL, last_6_years_end = NULL, last_6_years_start = NULL, min_end_cervical = NULL, last_6_years_end_NHS = NULL, last_6_years_start_NHS = NULL)]
  return(dataset)
}
flu_dt_1 <- process_flu_data(obs_dt_1)
flu_dt_2 <- process_flu_data(obs_dt_2)
flu_dt_3 <- process_flu_data(obs_dt_3)

#unique
flu_dt_1 <- unique(flu_dt_1)
flu_dt_2 <- unique(flu_dt_2)
flu_dt_3 <- unique(flu_dt_3)

#then turn cols to char
col_char <- function(data){
  data[, fluVE_breast_narrow := as.character(fluVE_breast_narrow)]
  data[, fluVE_breast_broad := as.character(fluVE_breast_broad)]
  data[, fluVE_bowel := as.character(fluVE_bowel)]
  data[, fluVE_cervical_narrow := as.character(fluVE_cervical_narrow)]
  data[, fluVE_cervical_broad := as.character(fluVE_cervical_broad)]
  data[, fluVE_PSA := as.character(fluVE_PSA)]
  data[, fluVE_bone := as.character(fluVE_bone)]
  data[, fluVE_NHS := as.character(fluVE_NHS)]
}
flu_dt_1 <- col_char(flu_dt_1)
flu_dt_2 <- col_char(flu_dt_2)
flu_dt_3 <- col_char(flu_dt_3)

#then create one row per patid
flu_dt_1 <- flu_dt_1[, lapply(.SD, collapse_rows), by = patid]
#length(unique(flu_dt_1$patid))
flu_dt_2 <- flu_dt_2[, lapply(.SD, collapse_rows), by = patid]
#length(unique(flu_dt_2$patid))
flu_dt_3 <- flu_dt_3[, lapply(.SD, collapse_rows), by = patid]
#length(unique(flu_dt_3$patid))

#change the format of the data
format_data <- function(data){
  data[, fluVE_breast_narrow := ifelse(fluVE_breast_narrow == TRUE, 1, 0)]
  data[, fluVE_breast_broad := ifelse(fluVE_breast_broad == TRUE, 1, 0)]
  data[, fluVE_bowel  := ifelse(fluVE_bowel  == TRUE, 1, 0)]
  data[, fluVE_cervical_narrow := ifelse(fluVE_cervical_narrow == TRUE, 1, 0)]
  data[, fluVE_cervical_broad := ifelse(fluVE_cervical_broad == TRUE, 1, 0)]
  data[, fluVE_PSA := ifelse(fluVE_PSA == TRUE, 1, 0)]
  data[, fluVE_bone := ifelse(fluVE_bone == TRUE, 1, 0)]
  data[, fluVE_NHS := ifelse(fluVE_NHS == TRUE, 1, 0)]
  data[is.na(data)] <- 0
  return(data)
}
flu_dt_1 <- format_data(flu_dt_1)
flu_dt_2 <- format_data(flu_dt_2)
flu_dt_3 <- format_data(flu_dt_3)

#rbind for flu
flu_dt_all <- rbind(flu_dt_1, flu_dt_2, flu_dt_3)
rm(flu_dt_1, flu_dt_2, flu_dt_3)

#write these for later use
setwd(datafiles)
write_parquet(covid_dt_all, paste0(datafiles, "covid_pt1markers.parquet"))
write_parquet(neg_dt_all, paste0(datafiles, "neg_pt1markers.parquet"))
write_parquet(flu_dt_all, paste0(datafiles, "flu_pt1markers.parquet"))
############################################################

################SHORT LB RESTRICTIVE########################
setwd(datafiles)
obs_dt_1 <- read_parquet("observations_markers_short_lookback_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_short_lookback_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_short_lookback_tmp3.parquet")

#deduplicate each patient for each marker num keeping only the latest event 
process_dataset <- function(dataset) {
  dataset <- dataset[order(patid)]
  dataset[, duplicates := .N > 1, by = c("marker_num", "patid")]
  dataset[duplicates == TRUE, last_entry := lapply(.SD, max), by = c("patid", "marker_num"), .SDcols = "eventdate"]
  dataset[duplicates == TRUE, keep_row := eventdate == last_entry]
  dataset[duplicates == FALSE, keep_row := TRUE]
  dataset <- dataset[keep_row == TRUE]
  dataset[, duplicates := NULL]
  dataset[, last_entry := NULL]
  dataset[, keep_row := NULL]
  return(dataset)
}
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#covid - restrictive definition
process_covid_data <- function(dataset) {
  dataset[, dob := as.Date(dob, format = "%d/%m/%Y")]
  dataset[, age := as.numeric(difftime(as.Date("08/12/2020", format = "%d/%m/%Y"), dob, units = "days")) / 365.25]
  eligible_age_end_breast <- 71
  dataset[, min_end_breast := pmin(eligible_age_end_breast, age)]
  dataset[, last_4_years_end := dob + (min_end_breast*365.25)]
  dataset[, last_4_years_start := last_4_years_end - 1461]
  dataset[, last_4_years_end := as.IDate(last_4_years_end, format = "%d/%m/%Y")]
  dataset[, last_4_years_start := as.IDate(last_4_years_start, format = "%d/%m/%Y")]
  eligible_age_end_cervical <- 64
  dataset[, min_end_cervical := pmin(eligible_age_end_cervical, age)]
  dataset[, last_6_years_end := dob + (min_end_cervical*365.25)]
  dataset[, last_6_years_start := last_6_years_end - 1461]
  dataset[, last_6_years_end := as.IDate(last_6_years_end, format = "%d/%m/%Y")]
  dataset[, last_6_years_start := as.IDate(last_6_years_start, format = "%d/%m/%Y")]
  eligible_age_end <- 74
  dataset[, min_end_bowel := pmin(eligible_age_end, age)]
  dataset[, last_3_years_end := dob + (min_end_bowel*365.25)]
  dataset[, last_3_years_start := last_3_years_end - 1095.75]
  dataset[, last_3_years_end := as.IDate(last_3_years_end, format = "%d/%m/%Y")]
  dataset[, last_3_years_start := as.IDate(last_3_years_start, format = "%d/%m/%Y")]
  dataset[, min_end_NHS := pmin(eligible_age_end, age)]
  dataset[, last_6_years_end_NHS := dob + (min_end_NHS*365.25)]
  dataset[, last_6_years_start_NHS := last_6_years_end_NHS - 2191.5]
  dataset[, last_6_years_end_NHS := as.IDate(last_6_years_end_NHS, format = "%d/%m/%Y")]
  dataset[, last_6_years_start_NHS := as.IDate(last_6_years_start_NHS, format = "%d/%m/%Y")]
  dataset <- dataset[covidVE == TRUE & (marker_num == 2 | marker_num == 16 | marker_num == 3 | marker_num == 4 | marker_num == 17 | marker_num == 8)]
  dataset[marker_num == 2, COVIDVE_breast_narrow_res := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 16, COVIDVE_breast_broad_res := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 3, COVIDVE_bowel_res := eventdate >= last_3_years_start & eventdate <= last_3_years_end, by = .(patid)]
  dataset[marker_num == 4, COVIDVE_cervical_narrow_res := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 17, COVIDVE_cervical_broad_res := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 8, COVIDVE_NHS_res := eventdate >= last_6_years_start_NHS & eventdate <= last_6_years_end_NHS, by = .(patid)]
  dataset <- dataset[COVIDVE_breast_narrow_res == TRUE | COVIDVE_breast_broad_res == TRUE | COVIDVE_bowel_res == TRUE | COVIDVE_cervical_narrow_res == TRUE | COVIDVE_cervical_broad_res == TRUE | COVIDVE_NHS_res == TRUE]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, enterdate = NULL, eventdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL, age = NULL, min_end_breast = NULL, min_end_bowel = NULL, min_end_NHS = NULL, last_4_years_end = NULL, last_4_years_start = NULL, last_3_years_end = NULL, last_3_years_start = NULL, last_6_years_end = NULL, last_6_years_start = NULL, min_end_cervical = NULL, last_6_years_end_NHS = NULL, last_6_years_start_NHS = NULL)]
  return(dataset)
}
covid_dt_1 <- process_covid_data(obs_dt_1)
covid_dt_2 <- process_covid_data(obs_dt_2)
covid_dt_3 <- process_covid_data(obs_dt_3)

#unique columns 
covid_dt_1 <- unique(covid_dt_1)
covid_dt_2 <- unique(covid_dt_2)
covid_dt_3 <- unique(covid_dt_3)

#then turn cols to char
col_char <- function(data){
  data[, COVIDVE_breast_narrow_res := as.character(COVIDVE_breast_narrow_res)]
  data[, COVIDVE_breast_broad_res := as.character(COVIDVE_breast_broad_res)]
  data[, COVIDVE_bowel_res := as.character(COVIDVE_bowel_res)]
  data[, COVIDVE_cervical_narrow_res := as.character(COVIDVE_cervical_narrow_res)]
  data[, COVIDVE_cervical_broad_res := as.character(COVIDVE_cervical_broad_res)]
  data[, COVIDVE_NHS_res := as.character(COVIDVE_NHS_res)]
}
covid_dt_1 <- col_char(covid_dt_1)
covid_dt_2 <- col_char(covid_dt_2)
covid_dt_3 <- col_char(covid_dt_3)

#then create one row per patid
covid_dt_1 <- covid_dt_1[, lapply(.SD, collapse_rows), by = patid]
#length(unique(covid_dt_1$patid))
covid_dt_2 <- covid_dt_2[, lapply(.SD, collapse_rows), by = patid]
#length(unique(covid_dt_2$patid))
covid_dt_3 <- covid_dt_3[, lapply(.SD, collapse_rows), by = patid]
#length(unique(covid_dt_3$patid))

#change the format of the data
format_data <- function(data){
  data[, COVIDVE_breast_narrow_res := ifelse(COVIDVE_breast_narrow_res == TRUE, 1, 0)]
  data[, COVIDVE_breast_broad_res := ifelse(COVIDVE_breast_broad_res == TRUE, 1, 0)]
  data[, COVIDVE_bowel_res  := ifelse(COVIDVE_bowel_res == TRUE, 1, 0)]
  data[, COVIDVE_cervical_narrow_res := ifelse(COVIDVE_cervical_narrow_res == TRUE, 1, 0)]
  data[, COVIDVE_cervical_broad_res := ifelse(COVIDVE_cervical_broad_res == TRUE, 1, 0)]
  data[, COVIDVE_NHS_res := ifelse(COVIDVE_NHS_res == TRUE, 1, 0)]
  data[is.na(data)] <- 0
  return(data)
}
covid_dt_1 <- format_data(covid_dt_1)
covid_dt_2 <- format_data(covid_dt_2)
covid_dt_3 <- format_data(covid_dt_3)

#rbind for covid
covid_dt_all <- rbind(covid_dt_1, covid_dt_2, covid_dt_3)
rm(covid_dt_1, covid_dt_2, covid_dt_3)

#now do for neg control
obs_dt_1 <- read_parquet("observations_markers_short_lookback_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_short_lookback_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_short_lookback_tmp3.parquet")

#remove events that occur after 1 jan 2020
obs_dt_1 <- obs_dt_1[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
obs_dt_2 <- obs_dt_2[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
obs_dt_3 <- obs_dt_3[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]

#then run process data function 
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#neg - restrictive definition
process_neg_data <- function(dataset) {
  dataset[, dob := as.Date(dob, format = "%d/%m/%Y")]
  dataset[, age := as.numeric(difftime(as.Date("01/01/2020", format = "%d/%m/%Y"), dob, units = "days")) / 365.25]
  eligible_age_end_breast <- 71
  dataset[, min_end_breast := pmin(eligible_age_end_breast, age)]
  dataset[, last_4_years_end := dob + (min_end_breast*365.25)]
  dataset[, last_4_years_start := last_4_years_end - 1461]
  dataset[, last_4_years_end := as.IDate(last_4_years_end, format = "%d/%m/%Y")]
  dataset[, last_4_years_start := as.IDate(last_4_years_start, format = "%d/%m/%Y")]
  eligible_age_end_cervical <- 64
  dataset[, min_end_cervical := pmin(eligible_age_end_cervical, age)]
  dataset[, last_6_years_end := dob + (min_end_cervical*365.25)]
  dataset[, last_6_years_start := last_6_years_end - 1461]
  dataset[, last_6_years_end := as.IDate(last_6_years_end, format = "%d/%m/%Y")]
  dataset[, last_6_years_start := as.IDate(last_6_years_start, format = "%d/%m/%Y")]
  eligible_age_end <- 74
  dataset[, min_end_bowel := pmin(eligible_age_end, age)]
  dataset[, last_3_years_end := dob + (min_end_bowel*365.25)]
  dataset[, last_3_years_start := last_3_years_end - 1095.75]
  dataset[, last_3_years_end := as.IDate(last_3_years_end, format = "%d/%m/%Y")]
  dataset[, last_3_years_start := as.IDate(last_3_years_start, format = "%d/%m/%Y")]
  dataset[, min_end_NHS := pmin(eligible_age_end, age)]
  dataset[, last_6_years_end_NHS := dob + (min_end_NHS*365.25)]
  dataset[, last_6_years_start_NHS := last_6_years_end_NHS - 2191.5]
  dataset[, last_6_years_end_NHS := as.IDate(last_6_years_end_NHS, format = "%d/%m/%Y")]
  dataset[, last_6_years_start_NHS := as.IDate(last_6_years_start_NHS, format = "%d/%m/%Y")]
  dataset <- dataset[negcontrolVE == TRUE & (marker_num == 2 | marker_num == 16 | marker_num == 3 | marker_num == 4 | marker_num == 17 | marker_num == 8)]
  dataset[marker_num == 2, negVE_breast_narrow_res := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 16, negVE_breast_broad_res := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 3, negVE_bowel_res := eventdate >= last_3_years_start & eventdate <= last_3_years_end, by = .(patid)]
  dataset[marker_num == 4, negVE_cervical_narrow_res := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 17, negVE_cervical_broad_res := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 8, negVE_NHS_res := eventdate >= last_6_years_start_NHS & eventdate <= last_6_years_end_NHS, by = .(patid)]
  dataset <- dataset[negVE_breast_narrow_res == TRUE | negVE_breast_broad_res == TRUE | negVE_bowel_res == TRUE | negVE_cervical_narrow_res == TRUE | negVE_cervical_broad_res == TRUE | negVE_NHS_res == TRUE]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, enterdate = NULL, eventdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL, age = NULL, min_end_breast = NULL, min_end_bowel = NULL, min_end_NHS = NULL, last_4_years_end = NULL, last_4_years_start = NULL, last_3_years_end = NULL, last_3_years_start = NULL, last_6_years_end = NULL, last_6_years_start = NULL, min_end_cervical = NULL, last_6_years_end_NHS = NULL, last_6_years_start_NHS = NULL)]
  return(dataset)
}
neg_dt_1 <- process_neg_data(obs_dt_1)
neg_dt_2 <- process_neg_data(obs_dt_2)
neg_dt_3 <- process_neg_data(obs_dt_3)

#unique
neg_dt_1 <- unique(neg_dt_1)
neg_dt_2 <- unique(neg_dt_2)
neg_dt_3 <- unique(neg_dt_3)

#then turn cols to char
col_char <- function(data){
  data[, negVE_breast_narrow_res := as.character(negVE_breast_narrow_res)]
  data[, negVE_breast_broad_res := as.character(negVE_breast_broad_res)]
  data[, negVE_bowel_res := as.character(negVE_bowel_res)]
  data[, negVE_cervical_narrow_res := as.character(negVE_cervical_narrow_res)]
  data[, negVE_cervical_broad_res := as.character(negVE_cervical_broad_res)]
  data[, negVE_NHS_res := as.character(negVE_NHS_res)]
}
neg_dt_1 <- col_char(neg_dt_1)
neg_dt_2 <- col_char(neg_dt_2)
neg_dt_3 <- col_char(neg_dt_3)

#then create one row per patid
neg_dt_1 <- neg_dt_1[, lapply(.SD, collapse_rows), by = patid]
#length(unique(neg_dt_1$patid))
neg_dt_2 <- neg_dt_2[, lapply(.SD, collapse_rows), by = patid]
#length(unique(neg_dt_2$patid))
neg_dt_3 <- neg_dt_3[, lapply(.SD, collapse_rows), by = patid]
#length(unique(neg_dt_3$patid))

#change the format of the data
format_data <- function(data){
  data[, negVE_breast_narrow_res := ifelse(negVE_breast_narrow_res == TRUE, 1, 0)]
  data[, negVE_breast_broad_res := ifelse(negVE_breast_broad_res == TRUE, 1, 0)]
  data[, negVE_bowel_res  := ifelse(negVE_bowel_res == TRUE, 1, 0)]
  data[, negVE_cervical_narrow_res := ifelse(negVE_cervical_narrow_res == TRUE, 1, 0)]
  data[, negVE_cervical_broad_res := ifelse(negVE_cervical_broad_res == TRUE, 1, 0)]
  data[, negVE_NHS_res := ifelse(negVE_NHS_res == TRUE, 1, 0)]
  data[is.na(data)] <- 0
  return(data)
}
neg_dt_1 <- format_data(neg_dt_1)
neg_dt_2 <- format_data(neg_dt_2)
neg_dt_3 <- format_data(neg_dt_3)

#rbind for neg
neg_dt_all <- rbind(neg_dt_1, neg_dt_2, neg_dt_3)
rm(neg_dt_1, neg_dt_2, neg_dt_3)

#now do for flu
obs_dt_1 <- read_parquet("observations_markers_short_lookback_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_short_lookback_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_short_lookback_tmp3.parquet")

#restrict to observations before 1 september 2019
obs_dt_1 <- obs_dt_1[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
obs_dt_2 <- obs_dt_2[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
obs_dt_3 <- obs_dt_3[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]

#then run process data function 
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#flu - restrictive 
process_flu_data <- function(dataset) {
  dataset[, dob := as.Date(dob, format = "%d/%m/%Y")]
  dataset[, age := as.numeric(difftime(as.Date("01/09/2019", format = "%d/%m/%Y"), dob, units = "days")) / 365.25]
  eligible_age_end_breast <- 71
  dataset[, min_end_breast := pmin(eligible_age_end_breast, age)]
  dataset[, last_4_years_end := dob + (min_end_breast*365.25)]
  dataset[, last_4_years_start := last_4_years_end - 1461]
  dataset[, last_4_years_end := as.IDate(last_4_years_end, format = "%d/%m/%Y")]
  dataset[, last_4_years_start := as.IDate(last_4_years_start, format = "%d/%m/%Y")]
  eligible_age_end_cervical <- 64
  dataset[, min_end_cervical := pmin(eligible_age_end_cervical, age)]
  dataset[, last_6_years_end := dob + (min_end_cervical*365.25)]
  dataset[, last_6_years_start := last_6_years_end - 1461]
  dataset[, last_6_years_end := as.IDate(last_6_years_end, format = "%d/%m/%Y")]
  dataset[, last_6_years_start := as.IDate(last_6_years_start, format = "%d/%m/%Y")]
  eligible_age_end <- 74
  dataset[, min_end_bowel := pmin(eligible_age_end, age)]
  dataset[, last_3_years_end := dob + (min_end_bowel*365.25)]
  dataset[, last_3_years_start := last_3_years_end - 1095.75]
  dataset[, last_3_years_end := as.IDate(last_3_years_end, format = "%d/%m/%Y")]
  dataset[, last_3_years_start := as.IDate(last_3_years_start, format = "%d/%m/%Y")]
  dataset[, min_end_NHS := pmin(eligible_age_end, age)]
  dataset[, last_6_years_end_NHS := dob + (min_end_NHS*365.25)]
  dataset[, last_6_years_start_NHS := last_6_years_end_NHS - 2191.5]
  dataset[, last_6_years_end_NHS := as.IDate(last_6_years_end_NHS, format = "%d/%m/%Y")]
  dataset[, last_6_years_start_NHS := as.IDate(last_6_years_start_NHS, format = "%d/%m/%Y")]
  dataset <- dataset[influenzaVE == TRUE & (marker_num == 2 | marker_num == 16 | marker_num == 3 | marker_num == 4 | marker_num == 17 | marker_num == 8)]
  dataset[marker_num == 2, fluVE_breast_narrow_res := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 16, fluVE_breast_broad_res := eventdate >= last_4_years_start & eventdate <= last_4_years_end, by = .(patid)]
  dataset[marker_num == 3, fluVE_bowel_res := eventdate >= last_3_years_start & eventdate <= last_3_years_end, by = .(patid)]
  dataset[marker_num == 4, fluVE_cervical_narrow_res := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 17, fluVE_cervical_broad_res := eventdate >= last_6_years_start & eventdate <= last_6_years_end, by = .(patid)]
  dataset[marker_num == 8, fluVE_NHS_res := eventdate >= last_6_years_start_NHS & eventdate <= last_6_years_end_NHS, by = .(patid)]
  dataset <- dataset[fluVE_breast_narrow_res == TRUE | fluVE_breast_broad_res == TRUE | fluVE_bowel_res == TRUE | fluVE_cervical_narrow_res == TRUE | fluVE_cervical_broad_res == TRUE | fluVE_NHS_res == TRUE]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, enterdate = NULL, eventdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL, age = NULL, min_end_breast = NULL, min_end_bowel = NULL, min_end_NHS = NULL, last_4_years_end = NULL, last_4_years_start = NULL, last_3_years_end = NULL, last_3_years_start = NULL, last_6_years_end = NULL, last_6_years_start = NULL, min_end_cervical = NULL, last_6_years_end_NHS = NULL, last_6_years_start_NHS = NULL)]
  return(dataset)
}
flu_dt_1 <- process_flu_data(obs_dt_1)
flu_dt_2 <- process_flu_data(obs_dt_2)
flu_dt_3 <- process_flu_data(obs_dt_3)

#unique
flu_dt_1 <- unique(flu_dt_1)
flu_dt_2 <- unique(flu_dt_2)
flu_dt_3 <- unique(flu_dt_3)

#then turn cols to char
col_char <- function(data){
  data[, fluVE_breast_narrow_res := as.character(fluVE_breast_narrow_res)]
  data[, fluVE_breast_broad_res := as.character(fluVE_breast_broad_res)]
  data[, fluVE_bowel_res := as.character(fluVE_bowel_res)]
  data[, fluVE_cervical_narrow_res := as.character(fluVE_cervical_narrow_res)]
  data[, fluVE_cervical_broad_res := as.character(fluVE_cervical_broad_res)]
  data[, fluVE_NHS_res := as.character(fluVE_NHS_res)]
}
flu_dt_1 <- col_char(flu_dt_1)
flu_dt_2 <- col_char(flu_dt_2)
flu_dt_3 <- col_char(flu_dt_3)

#then create one row per patid
flu_dt_1 <- flu_dt_1[, lapply(.SD, collapse_rows), by = patid]
#length(unique(flu_dt_1$patid))
flu_dt_2 <- flu_dt_2[, lapply(.SD, collapse_rows), by = patid]
#length(unique(flu_dt_2$patid))
flu_dt_3 <- flu_dt_3[, lapply(.SD, collapse_rows), by = patid]
#length(unique(flu_dt_3$patid))

#change the format of the data
format_data <- function(data){
  data[, fluVE_breast_narrow_res := ifelse(fluVE_breast_narrow_res == TRUE, 1, 0)]
  data[, fluVE_breast_broad_res := ifelse(fluVE_breast_broad_res == TRUE, 1, 0)]
  data[, fluVE_bowel_res := ifelse(fluVE_bowel_res == TRUE, 1, 0)]
  data[, fluVE_cervical_narrow_res := ifelse(fluVE_cervical_narrow_res == TRUE, 1, 0)]
  data[, fluVE_cervical_broad_res := ifelse(fluVE_cervical_broad_res == TRUE, 1, 0)]
  data[, fluVE_NHS_res := ifelse(fluVE_NHS_res == TRUE, 1, 0)]
  data[is.na(data)] <- 0
  return(data)
}
flu_dt_1 <- format_data(flu_dt_1)
flu_dt_2 <- format_data(flu_dt_2)
flu_dt_3 <- format_data(flu_dt_3)

#rbind for flu
flu_dt_all <- rbind(flu_dt_1, flu_dt_2, flu_dt_3)
rm(flu_dt_1, flu_dt_2, flu_dt_3)

#write these for later use
setwd(datafiles)
write_parquet(covid_dt_all, paste0(datafiles, "covid_pt1markers_res.parquet"))
write_parquet(neg_dt_all, paste0(datafiles, "neg_pt1markers_res.parquet"))
write_parquet(flu_dt_all, paste0(datafiles, "flu_pt1markers_res.parquet"))
############################################################

################PRIMARY DNA################
obs_dt_1 <- read_parquet("observations_markers_pt2_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_pt2_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_pt2_tmp3.parquet")

#remove events that occur after 8 dec 2020 
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  return(dt)
}
obs_dt_1 <- remove_rows(obs_dt_1)
obs_dt_2 <- remove_rows(obs_dt_2)
obs_dt_3 <- remove_rows(obs_dt_3)

#add in analysis population from patient file
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
setkey(patient, patid)

#add info from patient file into observational files 
obs_dt_1[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
obs_dt_2[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
obs_dt_3[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
rm(patient)

#marker num is 9 for all
obs_dt_1[, marker_num := 9]
obs_dt_2[, marker_num := 9]
obs_dt_3[, marker_num := 9]

#write for later use
setwd(datafiles)
data_tables <- list (obs_dt_1 , obs_dt_2, obs_dt_3)

write_data <- function(data, datafiles) {
  for (i in seq_along(data)) {
    write_parquet(data[[i]], paste0(datafiles, "/observations_markers_primaryDNA_tmp", i, ".parquet"))
  }
}

write_data(data_tables, datafiles)
rm(write_data, data_tables)

#obs_dt_1 <- read_parquet("observations_markers_primaryDNA_tmp1.parquet")
#obs_dt_2 <- read_parquet("observations_markers_primaryDNA_tmp2.parquet")
#obs_dt_3 <- read_parquet("observations_markers_primaryDNA_tmp3.parquet")

#apply the process data function
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#load primary DNA code list
setwd(codelists)

#read in code lists from multiple sheets  
markers <- map(set_names(excel_sheets("Markers_medcodes_part2.xlsx"), c("Aurum_primary_DNA")),
               read_excel, path = "Markers_medcodes_part2.xlsx"
)

#merge code lists 
for (i in markers){
  allmarkers <- rbindlist(markers, use.names = TRUE, idcol = TRUE, fill=TRUE)
}
rm(i)

#then run for COVID
process_covid_data <- function(dataset) {
  dataset[covidVE == TRUE & marker_num == 9, COVIDVE_DNA := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset <- dataset[COVIDVE_DNA == TRUE]
  dataset[allmarkers, `:=`(flag = i.flag), on = "medcodeid"]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, eventdate = NULL, enterdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL)]
  return(dataset)
}
covid_dt_1 <- process_covid_data(obs_dt_1)
covid_dt_2 <- process_covid_data(obs_dt_2)
covid_dt_3 <- process_covid_data(obs_dt_3)

#unique columns 
covid_dt_1 <- unique(covid_dt_1)
covid_dt_2 <- unique(covid_dt_2)
covid_dt_3 <- unique(covid_dt_3)

#bind for covid
covid_dt_all <- rbind(covid_dt_1, covid_dt_2, covid_dt_3)

#need to change the format of the data
covid_dt_all[, COVID_VE_Narrow_DNA := ifelse(COVIDVE_DNA == TRUE & flag == c("Narrow"), 1, 0)]
covid_dt_all[, COVIDVE_DNA := NULL]
covid_dt_all[, flag := NULL]

#now do for negative control
obs_dt_1 <- read_parquet("observations_markers_primaryDNA_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_primaryDNA_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_primaryDNA_tmp3.parquet")

#remove events that occur after 1 jan 2020
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  return(dt)
}
obs_dt_1 <- remove_rows(obs_dt_1)
obs_dt_2 <- remove_rows(obs_dt_2)
obs_dt_3 <- remove_rows(obs_dt_3)

#apply the process data function
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#then run for neg
process_neg_data <- function(dataset) {
  dataset[negcontrolVE == TRUE & marker_num == 9, negVE_DNA := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset <- dataset[negVE_DNA == TRUE]
  dataset[allmarkers, `:=`(flag = i.flag), on = "medcodeid"]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, eventdate = NULL, enterdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL)]
  return(dataset)
}
neg_dt_1 <- process_neg_data(obs_dt_1)
neg_dt_2 <- process_neg_data(obs_dt_2)
neg_dt_3 <- process_neg_data(obs_dt_3)

#unique
neg_dt_1 <- unique(neg_dt_1)
neg_dt_2 <- unique(neg_dt_2)
neg_dt_3 <- unique(neg_dt_3)

#rbind
neg_dt_all <- rbind(neg_dt_1, neg_dt_2, neg_dt_3)

#need to change the format of the data
neg_dt_all[, neg_VE_Narrow_DNA := ifelse(negVE_DNA == TRUE & flag == c("Narrow"), 1, 0)]
neg_dt_all[, negVE_DNA := NULL]
neg_dt_all[, flag := NULL]

#now do for flu
obs_dt_1 <- read_parquet("observations_markers_primaryDNA_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_primaryDNA_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_primaryDNA_tmp3.parquet")

#remove events that occur after 1 sept 2019
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  return(dt)
}
obs_dt_1 <- remove_rows(obs_dt_1)
obs_dt_2 <- remove_rows(obs_dt_2)
obs_dt_3 <- remove_rows(obs_dt_3)

#apply the process data function
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#then run for flu
process_flu_data <- function(dataset) {
  dataset[influenzaVE == TRUE & marker_num == 9, fluVE_DNA := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset <- dataset[fluVE_DNA == TRUE]
  dataset[allmarkers, `:=`(flag = i.flag), on = "medcodeid"]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, eventdate = NULL, enterdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL)]
  return(dataset)
}
flu_dt_1 <- process_flu_data(obs_dt_1)
flu_dt_2 <- process_flu_data(obs_dt_2)
flu_dt_3 <- process_flu_data(obs_dt_3)

#unique
flu_dt_1 <- unique(flu_dt_1)
flu_dt_2 <- unique(flu_dt_2)
flu_dt_3 <- unique(flu_dt_3)

#rbind
flu_dt_all <- rbind(flu_dt_1, flu_dt_2, flu_dt_3)

#need to change the format of the data
flu_dt_all[, flu_VE_Narrow_DNA := ifelse(fluVE_DNA == TRUE & flag == c("Narrow"), 1, 0)]
flu_dt_all[, fluVE_DNA := NULL]
flu_dt_all[, flag := NULL]

#write these for later use
setwd(datafiles)
write_parquet(covid_dt_all, paste0(datafiles, "covid_DNAmarkers.parquet"))
write_parquet(neg_dt_all, paste0(datafiles, "neg_DNAmarkers.parquet"))
write_parquet(flu_dt_all, paste0(datafiles, "flu_DNAmarkers.parquet"))
############################################################

################BLOOD PRESSURE################
obs_dt_1 <- read_parquet("observations_markers_pt3_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_pt3_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_pt3_tmp3.parquet")

#remove events that occur after 8 dec 2020 
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  return(dt)
}
obs_dt_1 <- remove_rows(obs_dt_1)
obs_dt_2 <- remove_rows(obs_dt_2)
obs_dt_3 <- remove_rows(obs_dt_3)

#add in analysis population from patient file
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
setkey(patient, patid)

#add info from patient file into observational files 
obs_dt_1[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
obs_dt_2[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
obs_dt_3[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
rm(patient)

#marker num is 10 for all
obs_dt_1[, marker_num := 10]
obs_dt_2[, marker_num := 10]
obs_dt_3[, marker_num := 10]

#write for later use
setwd(datafiles)
data_tables <- list (obs_dt_1 , obs_dt_2, obs_dt_3)

write_data <- function(data, datafiles) {
  for (i in seq_along(data)) {
    write_parquet(data[[i]], paste0(datafiles, "/observations_markers_BP_tmp", i, ".parquet"))
  }
}

write_data(data_tables, datafiles)
rm(write_data, data_tables) 

#apply the process data function
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#read in code lists from multiple sheets  
setwd(codelists)
markers <- map(set_names(excel_sheets("Markers_medcodes_part3.xlsx"), c("Aurum_BP")),
               read_excel, path = "Markers_medcodes_part3.xlsx"
)

#merge code lists 
for (i in markers){
  allmarkers <- rbindlist(markers, use.names = TRUE, idcol = TRUE, fill=TRUE)
}
rm(i)

#then run for COVID
process_covid_data <- function(dataset) {
  dataset[covidVE == TRUE & marker_num == 10, COVIDVE_BP := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset <- dataset[COVIDVE_BP == TRUE]
  dataset[allmarkers, `:=`(flag = i.flag), on = "medcodeid"]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, eventdate = NULL, enterdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL)]
  return(dataset)
}
covid_dt_1 <- process_covid_data(obs_dt_1)
covid_dt_2 <- process_covid_data(obs_dt_2)
covid_dt_3 <- process_covid_data(obs_dt_3)

#unique columns 
covid_dt_1 <- unique(covid_dt_1)
covid_dt_2 <- unique(covid_dt_2)
covid_dt_3 <- unique(covid_dt_3)

#bind for covid
covid_dt_all <- rbind(covid_dt_1, covid_dt_2, covid_dt_3)

#then need to change the format of the data
covid_dt_all[, COVID_VE_Narrow_BP := ifelse(COVIDVE_BP == TRUE & flag == c("Narrow"), 1, 0)]
covid_dt_all[, COVIDVE_BP := NULL]
covid_dt_all[, flag := NULL]

#now do for negative control
obs_dt_1 <- read_parquet("observations_markers_BP_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_BP_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_BP_tmp3.parquet")

#remove events that occur after 1 jan 2020
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  return(dt)
}
obs_dt_1 <- remove_rows(obs_dt_1)
obs_dt_2 <- remove_rows(obs_dt_2)
obs_dt_3 <- remove_rows(obs_dt_3)

#apply the process data function
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#then run for neg
process_neg_data <- function(dataset) {
  dataset[negcontrolVE == TRUE & marker_num == 10, negVE_BP := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset <- dataset[negVE_BP == TRUE]
  dataset[allmarkers, `:=`(flag = i.flag), on = "medcodeid"]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, eventdate = NULL, enterdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL)]
  return(dataset)
}
neg_dt_1 <- process_neg_data(obs_dt_1)
neg_dt_2 <- process_neg_data(obs_dt_2)
neg_dt_3 <- process_neg_data(obs_dt_3)

#unique
neg_dt_1 <- unique(neg_dt_1)
neg_dt_2 <- unique(neg_dt_2)
neg_dt_3 <- unique(neg_dt_3)

#rbind
neg_dt_all <- rbind(neg_dt_1, neg_dt_2, neg_dt_3)

#then need to change the format of the data
neg_dt_all[, neg_VE_Narrow_BP := ifelse(negVE_BP == TRUE & flag == c("Narrow"), 1, 0)]
neg_dt_all[, negVE_BP := NULL]
neg_dt_all[, flag := NULL]

#now do for flu
obs_dt_1 <- read_parquet("observations_markers_BP_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_BP_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_BP_tmp3.parquet")

#remove events that occur after 1 sept 2019
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  return(dt)
}
obs_dt_1 <- remove_rows(obs_dt_1)
obs_dt_2 <- remove_rows(obs_dt_2)
obs_dt_3 <- remove_rows(obs_dt_3)

#apply the process data function
obs_dt_1 <- process_dataset(obs_dt_1)
obs_dt_2 <- process_dataset(obs_dt_2)
obs_dt_3 <- process_dataset(obs_dt_3)

#then run for flu
process_flu_data <- function(dataset) {
  dataset[influenzaVE == TRUE & marker_num == 10, fluVE_BP := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset <- dataset[fluVE_BP == TRUE]
  dataset[allmarkers, `:=`(flag = i.flag), on = "medcodeid"]
  dataset[, `:=` (obsid = NULL, obsdate = NULL, eventdate = NULL, enterdate = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, medcodeid = NULL, marker_num = NULL, obstypeid = NULL, dob = NULL)]
  return(dataset)
}
flu_dt_1 <- process_flu_data(obs_dt_1)
flu_dt_2 <- process_flu_data(obs_dt_2)
flu_dt_3 <- process_flu_data(obs_dt_3)

#unique
flu_dt_1 <- unique(flu_dt_1)
flu_dt_2 <- unique(flu_dt_2)
flu_dt_3 <- unique(flu_dt_3)

#rbind
flu_dt_all <- rbind(flu_dt_1, flu_dt_2, flu_dt_3)

#then need to change the format of the data
flu_dt_all[, flu_VE_Narrow_BP := ifelse(fluVE_BP == TRUE & flag == c("Narrow"), 1, 0)]
flu_dt_all[, fluVE_BP := NULL]
flu_dt_all[, flag := NULL]

#write these for later use
setwd(datafiles)
write_parquet(covid_dt_all, paste0(datafiles, "covid_BPmarkers.parquet"))
write_parquet(neg_dt_all, paste0(datafiles, "neg_BPmarkers.parquet"))
write_parquet(flu_dt_all, paste0(datafiles, "flu_BPmarkers.parquet"))
############################################################

################PNEUNOCOCCAL VACCINE########################
#read in medcodes
setwd(datafiles)
obs_dt_1 <- read_parquet("observations_markers_pt1_tmp1.parquet")
obs_dt_2 <- read_parquet("observations_markers_pt1_tmp2.parquet")
obs_dt_3 <- read_parquet("observations_markers_pt1_tmp3.parquet")

#restrict to pneunococcal vaccine
obs_dt_1 <- obs_dt_1[marker_num == 5]
obs_dt_2 <- obs_dt_2[marker_num == 5]
obs_dt_3 <- obs_dt_3[marker_num == 5]

#remove unwanted columns
obs_dt_1[, `:=` (medcodeid = NULL, marker_num = NULL, obsid = NULL, obsdate = NULL, enterdate = NULL, obstypeid = NULL, dob = NULL)]
obs_dt_2[, `:=` (medcodeid = NULL, marker_num = NULL, obsid = NULL, obsdate = NULL, enterdate = NULL, obstypeid = NULL, dob = NULL)]
obs_dt_3[, `:=` (medcodeid = NULL, marker_num = NULL, obsid = NULL, obsdate = NULL, enterdate = NULL, obstypeid = NULL, dob = NULL)]

#read in prodcodes 
prod_dt_1 <- read_parquet("drugissue_tmp1.parquet")
prod_dt_2 <- read_parquet("drugissue_tmp2.parquet")
prod_dt_3 <- read_parquet("drugissue_tmp3.parquet")

#pull in prodcode list 
setwd(codelists)
prodcodes <- map(set_names(excel_sheets("Prodcodes.xlsx"), c("Aurum_covid_vacc", "Aurum_flu_vacc", "Aurum_pneu_vacc", "Aurum_immuno")),
                 read_excel, path = "Prodcodes.xlsx"
)

#merge code lists 
for (i in prodcodes){
  allprodcodes <- rbindlist(prodcodes, use.names = TRUE, idcol = TRUE, fill=TRUE)
}
rm(i)

#set key for the patient file as this is the smallest 
setkey(allprodcodes, prodcodeid)

#add info from prodcode lists
prod_dt_1[allprodcodes, `:=`(.id = i..id), on = "prodcodeid"]
prod_dt_2[allprodcodes, `:=`(.id = i..id), on = "prodcodeid"]
prod_dt_3[allprodcodes, `:=`(.id = i..id), on = "prodcodeid"]
rm(allprodcodes, prodcodes)

#identify pneu vaccinations
prod_dt_1 <- prod_dt_1[.id == c("Aurum_pneu_vacc")]
prod_dt_2 <- prod_dt_2[.id == c("Aurum_pneu_vacc")]
prod_dt_3 <- prod_dt_3[.id == c("Aurum_pneu_vacc")]

#remove unwanted columns
prod_dt_1[, `:=` (prodcodeid = NULL, enterdate = NULL, quantity = NULL, quantunitid = NULL, dosageid = NULL, dob = NULL, .id = NULL)]
prod_dt_2[, `:=` (prodcodeid = NULL, enterdate = NULL, quantity = NULL, quantunitid = NULL, dosageid = NULL, dob = NULL, .id = NULL)]
prod_dt_3[, `:=` (prodcodeid = NULL, enterdate = NULL, quantity = NULL, quantunitid = NULL, dosageid = NULL, dob = NULL, .id = NULL)]

#issuedate is enterdate
prod_dt_1[, eventdate := issuedate]
prod_dt_2[, eventdate := issuedate]
prod_dt_3[, eventdate := issuedate]

#remove issuedate
prod_dt_1[, `:=` (issuedate = NULL)]
prod_dt_2[, `:=` (issuedate = NULL)]
prod_dt_3[, `:=` (issuedate = NULL)]

#rbind all the files 
pneu_dt_all <- rbind(obs_dt_1, obs_dt_2, obs_dt_3, prod_dt_1, prod_dt_2, prod_dt_3)#2,179,001

#unique
pneu_dt_all <- unique(pneu_dt_all)

#add in marker_num for all
pneu_dt_all[, marker_num := 5]

#deduplicate keeping only the earliest 
process_dataset <- function(dataset) {
  dataset <- dataset[order(patid)]
  dataset[, duplicates := .N > 1, by = c("marker_num", "patid")]
  dataset[duplicates == TRUE, first_entry := lapply(.SD, min), by = c("patid", "marker_num"), .SDcols = "eventdate"]
  dataset[duplicates == TRUE, keep_row := eventdate == first_entry]
  dataset[duplicates == FALSE, keep_row := TRUE]
  dataset <- dataset[keep_row == TRUE]
  dataset[, c("duplicates", "first_entry", "keep_row") := NULL]
  return(dataset)
}
pneu_dt_all <- process_dataset(pneu_dt_all)

#add in analysis population from patient file
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
setkey(patient, patid)

#add info from patient file into observational files 
pneu_dt_all[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
rm(patient)

#then apply condition for each 
apply_conditions <- function(data) {
  data[covidVE == TRUE & marker_num == 5, COVIDVE_pneu := eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  data[influenzaVE == TRUE & marker_num == 5, fluVE_pneu := eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  data[negcontrolVE == TRUE & marker_num == 5, negVE_pneu := eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  return(data)
}
pneu_dt_all <- apply_conditions(pneu_dt_all)

#all columns are narrow 
pneu_dt_all[, flag := c("Narrow")]

#change the format of the data
pneu_dt_all[, COVID_VE_Narrow_pneu := ifelse(COVIDVE_pneu == TRUE & flag == c("Narrow"), 1, 0)]
pneu_dt_all[, flu_VE_Narrow_pneu := ifelse(fluVE_pneu == TRUE & flag == c("Narrow"), 1, 0)]
pneu_dt_all[, neg_VE_Narrow_pneu := ifelse(negVE_pneu == TRUE & flag == c("Narrow"), 1, 0)]

#remove unwanted columns
pneu_dt_all[, `:=` (eventdate = NULL, marker_num = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, COVIDVE_pneu = NULL, fluVE_pneu = NULL, negVE_pneu = NULL, flag = NULL)]

#write for later use
setwd(datafiles)
write_parquet(pneu_dt_all, paste0(datafiles, "observations_markers_pneu.parquet"))
############################################################

################INFLUENZA VACCINE###########################
setwd(datafiles)
prod_dt_1 <- read_parquet("drugissue_tmp1.parquet")
prod_dt_2 <- read_parquet("drugissue_tmp2.parquet")
prod_dt_3 <- read_parquet("drugissue_tmp3.parquet")

#pull in prodcode list so that we can identify covid-19 and flu vaccinations
setwd(codelists)
prodcodes <- map(set_names(excel_sheets("Prodcodes.xlsx"), c("Aurum_covid_vacc", "Aurum_flu_vacc", "Aurum_pneu_vacc", "Aurum_immuno")),
                 read_excel, path = "Prodcodes.xlsx"
)

#merge code lists 
for (i in prodcodes){
  allprodcodes <- rbindlist(prodcodes, use.names = TRUE, idcol = TRUE, fill=TRUE)
}
rm(i)

#set key for the patient file as this is the smallest 
setkey(allprodcodes, prodcodeid)

#add info from prodcode lists 
prod_dt_1[allprodcodes, `:=`(.id = i..id), on = "prodcodeid"]
prod_dt_2[allprodcodes, `:=`(.id = i..id), on = "prodcodeid"]
prod_dt_3[allprodcodes, `:=`(.id = i..id), on = "prodcodeid"]
rm(allprodcodes, prodcodes)

#identify flu vaccinations
flu_prod_dt_1 <- prod_dt_1[.id == c("Aurum_flu_vacc")]
flu_prod_dt_2 <- prod_dt_2[.id == c("Aurum_flu_vacc")]
flu_prod_dt_3 <- prod_dt_3[.id == c("Aurum_flu_vacc")]

#clean up the file - exclude all vaccinations before 01.09.2018 and after 31.03.2020
clean_file <- function(data){
  data <- data[issuedate >= as.Date("01/09/2018", format = "%d/%m/%Y") & issuedate <= as.Date("31/03/2020", format = "%d/%m/%Y")]
  data[, eventdate := issuedate]
  data[, `:=` (prodcodeid = NULL, issuedate = NULL, enterdate = NULL, quantity = NULL, quantunitid = NULL, dosageid = NULL, duration = NULL, dob = NULL, .id = NULL)]
  colnames(data) <- c("patid",
                      "vacc_date")
  return(data)
}
flu_prod_dt_1 <- clean_file(flu_prod_dt_1)
flu_prod_dt_2 <- clean_file(flu_prod_dt_2)
flu_prod_dt_3 <- clean_file(flu_prod_dt_3)

#check for unique rows 
flu_prod_dt_1 <- unique(flu_prod_dt_1)
flu_prod_dt_2 <- unique(flu_prod_dt_2)
flu_prod_dt_3 <- unique(flu_prod_dt_3)

#rbind
flu_prod_dt <- rbind(flu_prod_dt_1, flu_prod_dt_2, flu_prod_dt_3)
rm(flu_prod_dt_1, flu_prod_dt_2, flu_prod_dt_3)

#for all prodcodes flag is prod
flu_prod_dt[, flag := c("Product")]

#now want to read in the medcodes for flu
setwd(datafiles)
med_dt_1 <- read_parquet("observations_vaccinemed_tmp1.1.parquet")
med_dt_2 <- read_parquet("observations_vaccinemed_tmp2.1.parquet")
med_dt_3 <- read_parquet("observations_vaccinemed_tmp3.1.parquet")

#want to merge this with code lists so that we can identify flu and covid vaccinations
#read in code lists from multiple sheets 
setwd(codelists)
exposuremed <- map(set_names(excel_sheets("Exposure_medcodes.xlsx"), c("Aurum_covid_vacc", "Aurum_flu_vacc")),
                   read_excel, path = "Exposure_medcodes.xlsx"
)

#merge code lists 
for (i in exposuremed){
  allexposuremed <- rbindlist(exposuremed, use.names = TRUE, idcol = TRUE, fill=TRUE)
}
rm(i)

#merge
setkey(allexposuremed, medcodeid)

#add info from patient file into observational files - do one at a time so doesn't crash
med_dt_1[allexposuremed, `:=`(.id = i..id, flag = i.flag), on = "medcodeid"]
med_dt_2[allexposuremed, `:=`(.id = i..id, flag = i.flag), on = "medcodeid"]
med_dt_3[allexposuremed, `:=`(.id = i..id, flag = i.flag), on = "medcodeid"]
rm(allexposuremed, exposuremed)

#identify flu vaccinations
flu_med_dt_1 <- med_dt_1[.id == c("Aurum_flu_vacc")]
flu_med_dt_2 <- med_dt_2[.id == c("Aurum_flu_vacc")]
flu_med_dt_3 <- med_dt_3[.id == c("Aurum_flu_vacc")]

#clean up the file - exclude all vaccinations before 01.09.2018 and after 31.03.2020
clean_file <- function(data){
  data <- data[eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("31/03/2020", format = "%d/%m/%Y")]
  data[, `:=` (obsid = NULL, obsdate = NULL, enterdate = NULL, medcodeid = NULL, value = NULL, numunitid = NULL, numrangehigh = NULL, numrangelow = NULL, obstypeid = NULL, dob = NULL, .id = NULL)]
  colnames(data) <- c("patid",
                      "vacc_date",
                      "flag")
  return(data)
}
flu_med_dt_1 <- clean_file(flu_med_dt_1)
flu_med_dt_2 <- clean_file(flu_med_dt_2)
flu_med_dt_3 <- clean_file(flu_med_dt_3)

#check for unique rows 
flu_med_dt_1 <- unique(flu_med_dt_1)
flu_med_dt_2 <- unique(flu_med_dt_2)
flu_med_dt_3 <- unique(flu_med_dt_3)

#rbind
flu_med_dt <- rbind(flu_med_dt_1, flu_med_dt_2, flu_med_dt_3)
rm(flu_med_dt_1, flu_med_dt_2, flu_med_dt_3)

#merge the med and prodcodes
flu_dt <- rbind(flu_med_dt, flu_prod_dt)

#create a wide dataset
flu_wide <- dcast(flu_dt, patid + vacc_date ~ flag, fun.aggregate = length,
                  value.var = "flag")
flu_wide <- as.data.table(flu_wide)
flu_wide[, .(sum_given = sum(Given), 
             sum_given_lag = sum(Given_lag),
             sum_neutral = sum(Neutral),
             sum_absent = sum(Absent),
             sum_adverse = sum(Adverse),
             sum_product = sum(Product))]

#summing up the same codes per day 
flu_filtered <- flu_wide[, by=c("patid", "vacc_date"), lapply(.SD, sum),
                         .SDcols = c("Given", "Given_lag", "Neutral", "Absent", "Adverse", "Product")]

###########ALGORITHM###########  
#algorithm - if only a med code you could also censor / exclude if there are very few 
#1. all codes count as vaccinated
flu_filtered[, flu_vacc := "vaccinated"]
#2. Just an absent code on one day => remove 
flu_filtered[Absent >0 & Given == 0 & Given_lag == 0 & Neutral == 0 & Adverse == 0 & Product ==0, flu_vacc:= "remove"]
#3. Just a neutral and absent on same day => remove 
flu_filtered[Absent >0 & Given == 0 & Given_lag == 0 & Neutral > 0 & Adverse == 0 & Product ==0, flu_vacc:= "remove"]
#4. Just a given and absent on same day => remove 
flu_filtered[Absent > 0 & Given > 0 & Given_lag == 0 & Neutral == 0 & Adverse == 0 & Product ==0, flu_vacc:= "remove"]
#5. Just a product and absent on same day => remove 
flu_filtered[Absent > 0 & Given == 0 & Given_lag == 0 & Neutral == 0 & Adverse == 0 & Product > 0, flu_vacc:= "remove"]
#6. Just med codes and no product code 
flu_filtered[Product ==0 & (Given > 0 | Given_lag >0 | Neutral > 0 | Adverse > 0), flu_vacc:= "med only"]

#check how many of each
janitor::tabyl(flu_filtered$flu_vacc) %>% adorn_pct_formatting()

#remove those with remove 
tmp <- flu_filtered[!flu_vacc == c("remove")]
rm(flu_filtered)
janitor::tabyl(tmp$flu_vacc) %>% adorn_pct_formatting()

#check number of doses
tmp2 <- tmp
tmp2[, n_dose := seq_len(.N), by = "patid"]
tmp2[, doses := .N, by = c("patid")]
janitor::tabyl(tmp2$doses) %>% adorn_pct_formatting() 

#if n_doses =2 then remove first row if neutral
tmp2[doses == 2 & n_dose == 1 & (Neutral > 0 & Given == 0 & Given_lag == 0 & Adverse ==0 & Product ==0), remove := c("yes")]#1,908,845
tmp3 <- tmp2[is.na(remove)]#1,763,451
tmp3[, n_dose := seq_len(.N), by = "patid"]
tmp3[, doses := .N, by = c("patid")]
tmp3[doses == 1 & (Neutral > 0 & Given == 0 & Given_lag == 0 & Adverse ==0 & Product ==0), .N]#30,955 / 1,763,451 = 1.8% - should we remove these? 

#create a wide dataset per patient
tmp3[, vacc_date := as.IDate(vacc_date, format = "%d/%m/%Y")]
tmp_wide <- dcast(tmp2, patid + doses ~ n_dose,
                  fun.aggregate = NULL,
                  value.var = "vacc_date")
tmp_wide <- as.data.table(tmp_wide)

#rename the columns
colnames(tmp_wide)[3] = "firstvacc"
colnames(tmp_wide)[4] = "secondvacc"
colnames(tmp_wide)[5] = "thirdvacc"
colnames(tmp_wide)[6] = "fourthvacc"
colnames(tmp_wide)[7] = "fifthvacc"
colnames(tmp_wide)[8] = "sixthvacc"
colnames(tmp_wide)[9] = "seventhvacc"
colnames(tmp_wide)[10] = "eighthvacc"
colnames(tmp_wide)[11] = "ninthvacc"
colnames(tmp_wide)[12] = "tenthvacc"
colnames(tmp_wide)[13] = "eleventhvacc"
colnames(tmp_wide)[14] = "twelthvacc"
colnames(tmp_wide)[15] = "thirtieththvacc"
colnames(tmp_wide)[16] = "fourteenthvacc"
colnames(tmp_wide)[17] = "fifteenthvacc"
colnames(tmp_wide)[18] = "sixteenthvacc"
colnames(tmp_wide)[19] = "seventeenthvacc"
colnames(tmp_wide)[20] = "eighteenthvacc"

tmp_wide[, firstvacc := as.IDate(firstvacc, format = "%d/%m/%Y")]
tmp_wide[, secondvacc := as.IDate(secondvacc, format = "%d/%m/%Y")]
tmp_wide[, thirdvacc := as.IDate(thirdvacc, format = "%d/%m/%Y")]
tmp_wide[, fourthvacc := as.IDate(fourthvacc, format = "%d/%m/%Y")]
tmp_wide[, fifthvacc := as.IDate(fifthvacc, format = "%d/%m/%Y")]
tmp_wide[, sixthvacc := as.IDate(sixthvacc, format = "%d/%m/%Y")]
tmp_wide[, seventhvacc := as.IDate(seventhvacc, format = "%d/%m/%Y")]
tmp_wide[, eighthvacc := as.IDate(eighthvacc, format = "%d/%m/%Y")]
tmp_wide[, ninthvacc := as.IDate(ninthvacc, format = "%d/%m/%Y")]
tmp_wide[, tenthvacc := as.IDate(tenthvacc, format = "%d/%m/%Y")]
tmp_wide[, eleventhvacc := as.IDate(eleventhvacc, format = "%d/%m/%Y")]
tmp_wide[, twelthvacc := as.IDate(twelthvacc, format = "%d/%m/%Y")]
tmp_wide[, thirtieththvacc := as.IDate(thirtieththvacc, format = "%d/%m/%Y")]
tmp_wide[, fourteenthvacc := as.IDate(fourteenthvacc, format = "%d/%m/%Y")]
tmp_wide[, fifteenthvacc := as.IDate(fifteenthvacc, format = "%d/%m/%Y")]
tmp_wide[, sixteenthvacc := as.IDate(sixteenthvacc, format = "%d/%m/%Y")]
tmp_wide[, seventeenthvacc := as.IDate(seventeenthvacc, format = "%d/%m/%Y")]
tmp_wide[, eighteenthvacc := as.IDate(eighteenthvacc, format = "%d/%m/%Y")]
tmp_wide[, doses := NULL]

#add in analysis population from patient file
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
setkey(patient, patid)

#add info from patient file into observational files 
tmp_wide[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]

#identify 
apply_conditions <- function(data) {
  covid_date_range_start <- as.Date("01/09/2019", format = "%d/%m/%Y")
  covid_date_range_end <- as.Date("31/03/2020", format = "%d/%m/%Y")
  flu_date_range_start <- as.Date("01/09/2018", format = "%d/%m/%Y")
  flu_date_range_end <- as.Date("31/03/2019", format = "%d/%m/%Y")
  data[, COVIDVE_flu := (covidVE == TRUE & firstvacc >= covid_date_range_start & firstvacc <= covid_date_range_end) |
         (covidVE == TRUE & secondvacc >= covid_date_range_start & secondvacc <= covid_date_range_end) |
         (covidVE == TRUE & thirdvacc >= covid_date_range_start & thirdvacc <= covid_date_range_end) |
         (covidVE == TRUE & fourthvacc >= covid_date_range_start & fourthvacc <= covid_date_range_end) |
         (covidVE == TRUE & fifthvacc >= covid_date_range_start & fifthvacc <= covid_date_range_end) |
         (covidVE == TRUE & sixthvacc >= covid_date_range_start & sixthvacc <= covid_date_range_end) |
         (covidVE == TRUE & seventhvacc >= covid_date_range_start & seventhvacc <= covid_date_range_end) |
         (covidVE == TRUE & eighthvacc >= covid_date_range_start & eighthvacc <= covid_date_range_end) |
         (covidVE == TRUE & ninthvacc >= covid_date_range_start & ninthvacc <= covid_date_range_end) |
         (covidVE == TRUE & tenthvacc >= covid_date_range_start & tenthvacc <= covid_date_range_end) |
         (covidVE == TRUE & eleventhvacc >= covid_date_range_start & eleventhvacc <= covid_date_range_end) |
         (covidVE == TRUE & twelthvacc >= covid_date_range_start & twelthvacc <= covid_date_range_end) |
         (covidVE == TRUE & thirtieththvacc >= covid_date_range_start & thirtieththvacc <= covid_date_range_end) |
         (covidVE == TRUE & fourteenthvacc >= covid_date_range_start & fourteenthvacc <= covid_date_range_end) |
         (covidVE == TRUE & fifteenthvacc >= covid_date_range_start & fifteenthvacc <= covid_date_range_end) |
         (covidVE == TRUE & sixteenthvacc >= covid_date_range_start & sixteenthvacc <= covid_date_range_end) |
         (covidVE == TRUE & seventeenthvacc >= covid_date_range_start & seventeenthvacc <= covid_date_range_end) |
         (covidVE == TRUE & eighteenthvacc >= covid_date_range_start & eighteenthvacc <= covid_date_range_end)]
  data[, fluVE_flu := (influenzaVE == TRUE & firstvacc >= flu_date_range_start & firstvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & secondvacc >= flu_date_range_start & secondvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & thirdvacc >= flu_date_range_start & thirdvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & fourthvacc >= flu_date_range_start & fourthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & fifthvacc >= flu_date_range_start & fifthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & sixthvacc >= flu_date_range_start & sixthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & seventhvacc >= flu_date_range_start & seventhvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & eighthvacc >= flu_date_range_start & eighthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & ninthvacc >= flu_date_range_start & ninthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & tenthvacc >= flu_date_range_start & tenthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & eleventhvacc >= flu_date_range_start & eleventhvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & twelthvacc >= flu_date_range_start & twelthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & thirtieththvacc >= flu_date_range_start & thirtieththvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & fourteenthvacc >= flu_date_range_start & fourteenthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & fifteenthvacc >= flu_date_range_start & fifteenthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & sixteenthvacc >= flu_date_range_start & sixteenthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & seventeenthvacc >= flu_date_range_start & seventeenthvacc <= flu_date_range_end) |
         (influenzaVE == TRUE & eighteenthvacc >= flu_date_range_start & eighteenthvacc <= flu_date_range_end)]
  data[, negVE_flu := (negcontrolVE == TRUE & firstvacc >= flu_date_range_start & firstvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & secondvacc >= flu_date_range_start & secondvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & thirdvacc >= flu_date_range_start & thirdvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & fourthvacc >= flu_date_range_start & fourthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & fifthvacc >= flu_date_range_start & fifthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & sixthvacc >= flu_date_range_start & sixthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & seventhvacc >= flu_date_range_start & seventhvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & eighthvacc >= flu_date_range_start & eighthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & ninthvacc >= flu_date_range_start & ninthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & tenthvacc >= flu_date_range_start & tenthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & eleventhvacc >= flu_date_range_start & eleventhvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & twelthvacc >= flu_date_range_start & twelthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & thirtieththvacc >= flu_date_range_start & thirtieththvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & fourteenthvacc >= flu_date_range_start & fourteenthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & fifteenthvacc >= flu_date_range_start & fifteenthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & sixteenthvacc >= flu_date_range_start & sixteenthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & seventeenthvacc >= flu_date_range_start & seventeenthvacc <= flu_date_range_end) |
         (negcontrolVE == TRUE & eighteenthvacc >= flu_date_range_start & eighteenthvacc <= flu_date_range_end)]
  data <- data[COVIDVE_flu | fluVE_flu | negVE_flu]
  data[, `:=` (firstvacc = NULL, secondvacc = NULL, thirdvacc = NULL, fourthvacc = NULL, fifthvacc = NULL, sixthvacc = NULL, seventhvacc = NULL, eighthvacc = NULL, ninthvacc = NULL, tenthvacc = NULL, eleventhvacc = NULL, twelthvacc = NULL, thirtieththvacc = NULL, fourteenthvacc = NULL, fifteenthvacc = NULL, sixteenthvacc = NULL, seventeenthvacc = NULL, eighteenthvacc = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL)]
  return(data)
}
tmp_wide <- apply_conditions(tmp_wide)

#change the format
tmp_wide[, COVIDVE_flu := ifelse(COVIDVE_flu == TRUE, 1, 0)]
tmp_wide[, fluVE_flu := ifelse(fluVE_flu == TRUE, 1, 0)]
tmp_wide[, negVE_flu  := ifelse(negVE_flu  == TRUE, 1, 0)]

#save this for later
setwd(datafiles)
write_parquet(tmp_wide, paste0(datafiles, "observations_markers_flu.parquet"))
################################################################################

################LOW VALUE PRESCRIPTIONS######################################### 
#read in the data
prod_dt_1 <- read_parquet("observations_markers_lowvaluepres_tmp1.parquet")
prod_dt_2 <- read_parquet("observations_markers_lowvaluepres_tmp2.parquet")
prod_dt_3 <- read_parquet("observations_markers_lowvaluepres_tmp3.parquet")

#remove events that occur after 8 dec 2020 
remove_rows <- function(dt){
  dt <- dt[issuedate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  return(dt)
}
prod_dt_1 <- remove_rows(prod_dt_1)
prod_dt_2 <- remove_rows(prod_dt_2)
prod_dt_3 <- remove_rows(prod_dt_3)
rm(remove_rows)

#add in analysis population from patient file
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
setkey(patient, patid)

#add info from patient file into observational files 
prod_dt_1[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
prod_dt_2[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
prod_dt_3[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
rm(patient)

#remove unwanted columns
prod_dt_1[, `:=` (prodcodeid = NULL, enterdate = NULL, quantity = NULL, quantunitid = NULL, dosageid = NULL, dob = NULL)]
prod_dt_2[, `:=` (prodcodeid = NULL, enterdate = NULL, quantity = NULL, quantunitid = NULL, dosageid = NULL, dob = NULL)]
prod_dt_3[, `:=` (prodcodeid = NULL, enterdate = NULL, quantity = NULL, quantunitid = NULL, dosageid = NULL, dob = NULL)]

#issuedate is enterdate
prod_dt_1[, eventdate := issuedate]
prod_dt_2[, eventdate := issuedate]
prod_dt_3[, eventdate := issuedate]

#remove issuedate
prod_dt_1[, `:=` (issuedate = NULL)]
prod_dt_2[, `:=` (issuedate = NULL)]
prod_dt_3[, `:=` (issuedate = NULL)]

#marker num is 12 for all
prod_dt_1[, marker_num := 12]
prod_dt_2[, marker_num := 12]
prod_dt_3[, marker_num := 12]

#write for later use
setwd(datafiles)
data_tables <- list (prod_dt_1 , prod_dt_2, prod_dt_3)

write_data <- function(data, datafiles) {
  for (i in seq_along(data)) {
    write_parquet(data[[i]], paste0(datafiles, "/observations_markers_lowpres_tmp", i, ".parquet"))
  }
}
write_data(data_tables, datafiles)
rm(write_data, data_tables) 

#prod_dt_1 <- read_parquet("observations_markers_lowpres_tmp1.parquet")
#prod_dt_2 <- read_parquet("observations_markers_lowpres_tmp2.parquet")
#prod_dt_3 <- read_parquet("observations_markers_lowpres_tmp3.parquet")

#deduplicate each patient for each marker num keeping only the latest event 
process_dataset <- function(dataset) {
  dataset <- dataset[order(patid)]
  dataset[, duplicates := .N > 1, by = c("marker_num", "patid")]
  dataset[duplicates == TRUE, last_entry := lapply(.SD, max), by = c("patid", "marker_num"), .SDcols = "eventdate"]
  dataset[duplicates == TRUE, keep_row := eventdate == last_entry]
  dataset[duplicates == FALSE, keep_row := TRUE]
  dataset <- dataset[keep_row == TRUE]
  dataset[, duplicates := NULL]
  dataset[, last_entry := NULL]
  dataset[, keep_row := NULL]
  return(dataset)
}
prod_dt_1 <- process_dataset(prod_dt_1)
prod_dt_2 <- process_dataset(prod_dt_2)
prod_dt_3 <- process_dataset(prod_dt_3)

#then run for COVID
process_covid_data <- function(dataset) {
  dataset[covidVE == TRUE & marker_num == 12, COVIDVE_lowpres := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset <- dataset[COVIDVE_lowpres == TRUE]
  dataset[, `:=` (flag = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, marker_num = NULL, eventdate = NULL)]
  return(dataset)
}
covid_dt_1 <- process_covid_data(prod_dt_1)
covid_dt_2 <- process_covid_data(prod_dt_2)
covid_dt_3 <- process_covid_data(prod_dt_3)

#unique
covid_dt_1 <- unique(covid_dt_1)
covid_dt_2 <- unique(covid_dt_2)
covid_dt_3 <- unique(covid_dt_3)

#bind for covid
covid_dt_all <- rbind(covid_dt_1, covid_dt_2, covid_dt_3)

#narrow for all
covid_dt_all[, flag := c("Narrow")]

#then need to change the format of the data
covid_dt_all[, COVID_VE_Narrow_lowpres := ifelse(COVIDVE_lowpres == TRUE & flag == c("Narrow"), 1, 0)]
covid_dt_all[, COVIDVE_lowpres := NULL]
covid_dt_all[, flag := NULL]

#now need to do for neg control
prod_dt_1 <- read_parquet("observations_markers_lowpres_tmp1.parquet")
prod_dt_2 <- read_parquet("observations_markers_lowpres_tmp2.parquet")
prod_dt_3 <- read_parquet("observations_markers_lowpres_tmp3.parquet")

#remove observations after 1 jan 2020 
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  return(dt)
}
prod_dt_1 <- remove_rows(prod_dt_1)
prod_dt_2 <- remove_rows(prod_dt_2)
prod_dt_3 <- remove_rows(prod_dt_3)

#process data 
prod_dt_1 <- process_dataset(prod_dt_1)
prod_dt_2 <- process_dataset(prod_dt_2)
prod_dt_3 <- process_dataset(prod_dt_3)

#then process for neg control
process_neg_data <- function(dataset) {
  dataset[negcontrolVE == TRUE & marker_num == 12, negVE_lowpres := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset <- dataset[negVE_lowpres == TRUE]
  dataset[, `:=` (flag = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, marker_num = NULL, eventdate = NULL)]
  return(dataset)
}
neg_dt_1 <- process_neg_data(prod_dt_1)
neg_dt_2 <- process_neg_data(prod_dt_2)
neg_dt_3 <- process_neg_data(prod_dt_3)

#unique
neg_dt_1 <- unique(neg_dt_1)
neg_dt_2 <- unique(neg_dt_2)
neg_dt_3 <- unique(neg_dt_3)

#rbind
neg_dt_all <- rbind(neg_dt_1, neg_dt_2, neg_dt_3)

#narrow for all
neg_dt_all[, flag := c("Narrow")]

#then need to change the format of the data
neg_dt_all[, neg_VE_Narrow_lowpres := ifelse(negVE_lowpres == TRUE & flag == c("Narrow"), 1, 0)]
neg_dt_all[, negVE_lowpres := NULL]
neg_dt_all[, flag := NULL]

#then need to do for flu
prod_dt_1 <- read_parquet("observations_markers_lowpres_tmp1.parquet")
prod_dt_2 <- read_parquet("observations_markers_lowpres_tmp2.parquet")
prod_dt_3 <- read_parquet("observations_markers_lowpres_tmp3.parquet")

#remove observations after 1 sept 2019 
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  return(dt)
}
prod_dt_1 <- remove_rows(prod_dt_1)
prod_dt_2 <- remove_rows(prod_dt_2)
prod_dt_3 <- remove_rows(prod_dt_3)

#process data 
prod_dt_1 <- process_dataset(prod_dt_1)
prod_dt_2 <- process_dataset(prod_dt_2)
prod_dt_3 <- process_dataset(prod_dt_3)

#now process for flu
process_flu_data <- function(dataset) {
  dataset[influenzaVE == TRUE & marker_num == 12, fluVE_lowpres := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset <- dataset[fluVE_lowpres == TRUE]
  dataset[, `:=` (flag = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, marker_num = NULL, eventdate = NULL)]
  return(dataset)
}
flu_dt_1 <- process_flu_data(prod_dt_1)
flu_dt_2 <- process_flu_data(prod_dt_2)
flu_dt_3 <- process_flu_data(prod_dt_3)

#unique
flu_dt_1 <- unique(flu_dt_1)
flu_dt_2 <- unique(flu_dt_2)
flu_dt_3 <- unique(flu_dt_3)

#rbind
flu_dt_all <- rbind(flu_dt_1, flu_dt_2, flu_dt_3)

#narrow for all
flu_dt_all[, flag := c("Narrow")]

#then need to change the format of the data
flu_dt_all[, flu_VE_Narrow_lowpres := ifelse(fluVE_lowpres == TRUE & flag == c("Narrow"), 1, 0)]
flu_dt_all[, fluVE_lowpres := NULL]
flu_dt_all[, flag := NULL]

#write these for later use
setwd(datafiles)
write_parquet(covid_dt_all, paste0(datafiles, "covid_lowpresmarkers.parquet"))
write_parquet(neg_dt_all, paste0(datafiles, "neg_lowpresmarkers.parquet"))
write_parquet(flu_dt_all, paste0(datafiles, "flu_lowpresmarkers.parquet"))
############################################################

################ACS CONDITIONS##############################  
hes_dt_1 <- read_parquet("HES_markers_ACS.parquet")

#remove events that occur after 8 dec 2020 
remove_rows <- function(dt){
  dt <- dt[admidate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  return(dt)
}
hes_dt_1 <- remove_rows(hes_dt_1)

#add in analysis population from patient file
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
setkey(patient, patid)

#add info from patient file into observational files 
hes_dt_1[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]

#admidate is enterdate
hes_dt_1[, eventdate := admidate]

#remove admidate
hes_dt_1[, `:=` (admidate = NULL)]

#marker num is 13 for all
hes_dt_1[, marker_num := 13]

#write this for later use 
setwd(datafiles)
write_parquet(hes_dt_1, paste0(datafiles, "HES_markers_ACS.1.parquet"))

#hes_dt_1 <- read_parquet("HES_markers_ACS.1.parquet")

#process the data 
process_dataset <- function(dataset) {
dataset <- dataset[order(patid)]
dataset[, duplicates := .N > 1, by = c("marker_num", "patid")]
dataset[duplicates == TRUE, last_entry := lapply(.SD, max), by = c("patid", "marker_num"), .SDcols = "eventdate"]
dataset[duplicates == TRUE, keep_row := eventdate == last_entry]
dataset[duplicates == FALSE, keep_row := TRUE]
dataset <- dataset[keep_row == TRUE]
dataset[, duplicates := NULL]
dataset[, last_entry := NULL]
dataset[, keep_row := NULL]
return(dataset)
}
hes_dt_1 <- process_dataset(hes_dt_1)

#then run for COVID
process_covid_data <- function(dataset) {
  dataset[covidVE == TRUE & marker_num == 13, COVIDVE_ACS := eventdate >= as.Date("01/09/2014", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset <- dataset[COVIDVE_ACS == TRUE]
  dataset[, `:=` (ICD_PRIMARY = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, eventdate = NULL, marker_num = NULL)]
  return(dataset)
}
covid_dt_1 <- process_covid_data(hes_dt_1)

#unique
covid_dt_1 <- unique(covid_dt_1)

#narrow for all
covid_dt_1[, flag := c("Narrow")]

#change the format of the data
covid_dt_1[, COVID_VE_Narrow_ACS := ifelse(COVIDVE_ACS == TRUE & flag == c("Narrow"), 1, 0)]
covid_dt_1[, COVIDVE_ACS := NULL]
covid_dt_1[, flag := NULL]

#now run for neg control
hes_dt_1 <- read_parquet("HES_markers_ACS.1.parquet")

#exclude those after 1 jan 2020
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  return(dt)
}
hes_dt_1 <- remove_rows(hes_dt_1)

#run the process data 
hes_dt_1 <- process_dataset(hes_dt_1)

#then run for neg control
process_neg_data <- function(dataset) {
  dataset[negcontrolVE == TRUE & marker_num == 13, negVE_ACS := eventdate >= as.Date("01/09/2014", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset <- dataset[negVE_ACS == TRUE]
  dataset[, `:=` (ICD_PRIMARY = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, eventdate = NULL, marker_num = NULL)]
  return(dataset)
}
neg_dt_1 <- process_neg_data(hes_dt_1)

#unique
neg_dt_1 <- unique(neg_dt_1)

#narrow for all
neg_dt_1[, flag := c("Narrow")]

#change the format of the data
neg_dt_1[, neg_VE_Narrow_ACS := ifelse(negVE_ACS == TRUE & flag == c("Narrow"), 1, 0)]
neg_dt_1[, negVE_ACS := NULL]
neg_dt_1[, flag := NULL]

#now run for flu
hes_dt_1 <- read_parquet("HES_markers_ACS.1.parquet")

#exclude those after 1 sept 2019
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  return(dt)
}
hes_dt_1 <- remove_rows(hes_dt_1)

#run the process data 
hes_dt_1 <- process_dataset(hes_dt_1)

#then run for neg control
process_flu_data <- function(dataset) {
  dataset[influenzaVE == TRUE & marker_num == 13, fluVE_ACS := eventdate >= as.Date("01/09/2014", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset <- dataset[fluVE_ACS == TRUE]
  dataset[, `:=` (ICD_PRIMARY = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, eventdate = NULL, marker_num = NULL)]
  return(dataset)
}
flu_dt_1 <- process_flu_data(hes_dt_1)

#unique
flu_dt_1 <- unique(flu_dt_1)

#narrow for all
flu_dt_1[, flag := c("Narrow")]

#change the format of the data
flu_dt_1[, flu_VE_Narrow_ACS := ifelse(fluVE_ACS == TRUE & flag == c("Narrow"), 1, 0)]
flu_dt_1[, fluVE_ACS := NULL]
flu_dt_1[, flag := NULL]

#write these for later use
setwd(datafiles)
write_parquet(covid_dt_1, paste0(datafiles, "covid_ACS.parquet"))
write_parquet(neg_dt_1, paste0(datafiles, "neg_ACS.parquet"))
write_parquet(flu_dt_1, paste0(datafiles, "flu_ACS.parquet"))
############################################################

################LOW VALUE PROCEDURES########################
setwd(datafiles)
proc_dt_1 <- read_parquet("HES_markers_lowvalueproc.parquet")

#remove events that occur after 8 dec 2020 
remove_rows <- function(dt){
  dt <- dt[admidate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  return(dt)
}
proc_dt_1 <- remove_rows(proc_dt_1)

#add in analysis population from patient file
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
setkey(patient, patid)

#add info from patient file into observational files 
proc_dt_1[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
rm(patient)

#admidate is enterdate
proc_dt_1[, eventdate := admidate]

#remove admidate
proc_dt_1[, `:=` (admidate = NULL)]

#marker num is 14 for all
proc_dt_1[, marker_num := 14]

#write this for later use 
setwd(datafiles)
write_parquet(proc_dt_1, paste0(datafiles, "HES_markers_lowvalueproc.1.parquet"))

#proc_dt_1 <- read_parquet("HES_markers_lowvalueproc.1.parquet")

#process the data 
process_dataset <- function(dataset) {
  dataset <- dataset[order(patid)]
  dataset[, duplicates := .N > 1, by = c("marker_num", "patid")]
  dataset[duplicates == TRUE, last_entry := lapply(.SD, max), by = c("patid", "marker_num"), .SDcols = "eventdate"]
  dataset[duplicates == TRUE, keep_row := eventdate == last_entry]
  dataset[duplicates == FALSE, keep_row := TRUE]
  dataset <- dataset[keep_row == TRUE]
  dataset[, duplicates := NULL]
  dataset[, last_entry := NULL]
  dataset[, keep_row := NULL]
  return(dataset)
}
proc_dt_1 <- process_dataset(proc_dt_1)

#then run for COVID
process_covid_data <- function(dataset) {
  dataset[covidVE == TRUE & marker_num == 14, COVIDVE_lowvalueproc := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset <- dataset[COVIDVE_lowvalueproc == TRUE]
  dataset[, `:=` (OPCS = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, eventdate = NULL, marker_num = NULL, evdate = NULL)]
  return(dataset)
}
covid_dt_1 <- process_covid_data(proc_dt_1)

#unique
covid_dt_1 <- unique(covid_dt_1)

#narrow for all
covid_dt_1[, flag := c("Narrow")]

#change the format of the data
covid_dt_1[, COVID_VE_Narrow_lowvalueproc := ifelse(COVIDVE_lowvalueproc == TRUE & flag == c("Narrow"), 1, 0)]
covid_dt_1[, COVIDVE_lowvalueproc := NULL]
covid_dt_1[, flag := NULL]

#now run for neg control
proc_dt_1 <- read_parquet("HES_markers_lowvalueproc.1.parquet")

#exclude those after 1 jan 2020
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  return(dt)
}
proc_dt_1 <- remove_rows(proc_dt_1)

#run the process data 
proc_dt_1 <- process_dataset(proc_dt_1)

#then run for neg control
process_neg_data <- function(dataset) {
  dataset[negcontrolVE == TRUE & marker_num == 14, negVE_lowvalueproc := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset <- dataset[negVE_lowvalueproc == TRUE]
  dataset[, `:=` (OPCS = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, eventdate = NULL, marker_num = NULL, evdate = NULL)]
  return(dataset)
}
neg_dt_1 <- process_neg_data(proc_dt_1)

#unique
neg_dt_1 <- unique(neg_dt_1)

#narrow for all
neg_dt_1[, flag := c("Narrow")]

#change the format of the data
neg_dt_1[, neg_VE_Narrow_lowvalueproc := ifelse(negVE_lowvalueproc == TRUE & flag == c("Narrow"), 1, 0)]
neg_dt_1[, negVE_lowvalueproc := NULL]
neg_dt_1[, flag := NULL]

#now run for flu
proc_dt_1 <- read_parquet("HES_markers_lowvalueproc.1.parquet")

#exclude those after 1 sept 2019
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  return(dt)
}
proc_dt_1 <- remove_rows(proc_dt_1)

#run the process data 
proc_dt_1 <- process_dataset(proc_dt_1)

#then run for neg control
process_flu_data <- function(dataset) {
  dataset[influenzaVE == TRUE & marker_num == 14, fluVE_lowvalueproc := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset <- dataset[fluVE_lowvalueproc == TRUE]
  dataset[, `:=` (OPCS = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL, eventdate = NULL, marker_num = NULL, evdate = NULL)]
  return(dataset)
}
flu_dt_1 <- process_flu_data(proc_dt_1)

#unique
flu_dt_1 <- unique(flu_dt_1)

#narrow for all
flu_dt_1[, flag := c("Narrow")]

#change the format of the data
flu_dt_1[, flu_VE_Narrow_lowvalueproc := ifelse(fluVE_lowvalueproc == TRUE & flag == c("Narrow"), 1, 0)]
flu_dt_1[, fluVE_lowvalueproc := NULL]
flu_dt_1[, flag := NULL]

#write these for later use
setwd(datafiles)
write_parquet(covid_dt_1, paste0(datafiles, "covid_lowvalueproc.parquet"))
write_parquet(neg_dt_1, paste0(datafiles, "neg_lowvalueproc.parquet"))
write_parquet(flu_dt_1, paste0(datafiles, "flu_lowvalueproc.parquet"))
############################################################

################GP VISITS####################
#read in the files
gp_dt_1 <- read_parquet("observations_gp_visits_tmp1.1.parquet")
gp_dt_2 <- read_parquet("observations_gp_visits_tmp2.1.parquet")
gp_dt_3 <- read_parquet("observations_gp_visits_tmp3.1.parquet")

#add in analysis population from patient file
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
setkey(patient, patid)

#add info from patient file into observational files 
gp_dt_1[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
gp_dt_2[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
gp_dt_3[patient, `:=`(covidVE = i.covidVE, influenzaVE = i.influenzaVE, negcontrolVE = i.negcontrolVE), on = "patid"]
rm(patient)

#write this because will need to do separately for diff analyses
data_tables <- list (gp_dt_1, gp_dt_2, gp_dt_3)

write_data <- function(data, datafiles) {
  iteration <- seq(length(data)) + 0.2
  for (i in seq_along(data)) {
    write_parquet(data[[i]], paste0(datafiles, "/observations_gp_visits_tmp", iteration[i], ".parquet"))
  }
}
write_data(data_tables, datafiles)
rm(data_tables)

#then process the data
process_dataset <- function(dataset) {
  dataset <- dataset[order(patid)]
  dataset[, duplicates := .N > 1, by = c("patid")]
  dataset[duplicates == TRUE, last_entry := lapply(.SD, max), by = c("patid"), .SDcols = "eventdate"]
  dataset[duplicates == TRUE, keep_row := eventdate == last_entry]
  dataset[duplicates == FALSE, keep_row := TRUE]
  dataset <- dataset[keep_row == TRUE]
  dataset[, duplicates := NULL]
  dataset[, last_entry := NULL]
  dataset[, keep_row := NULL]
  return(dataset)
}
gp_dt_1 <- process_dataset(gp_dt_1)
gp_dt_2 <- process_dataset(gp_dt_2)
gp_dt_3 <- process_dataset(gp_dt_3)

#then process for covid
process_covid_data <- function(dataset) {
  dataset[covidVE == TRUE, COVIDVE_gpvisits := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("08/12/2020", format = "%d/%m/%Y")]
  dataset <- dataset[COVIDVE_gpvisits == TRUE]
  dataset[, `:=` (staffid = NULL, conssourceid = NULL, consmedcodeid = NULL, cprdconstype = NULL, consdate = NULL, enterdate = NULL, eventdate = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL)]
  return(dataset)
}
covid_dt_1 <- process_covid_data(gp_dt_1)
covid_dt_2 <- process_covid_data(gp_dt_2)
covid_dt_3 <- process_covid_data(gp_dt_3)

#then identify unique 
covid_dt_1 <- unique(covid_dt_1)
covid_dt_2 <- unique(covid_dt_2)
covid_dt_3 <- unique(covid_dt_3)
#length(unique(covid_dt_1$patid))

#then format the data
format_data <- function(data){
  data[, COVIDVE_gpvisits := ifelse(COVIDVE_gpvisits == TRUE, 1, 0)]
  data[is.na(data)] <- 0
  return(data)
}
covid_dt_1 <- format_data(covid_dt_1)
covid_dt_2 <- format_data(covid_dt_2)
covid_dt_3 <- format_data(covid_dt_3)

#then rbind
covid_dt <- rbind(covid_dt_1, covid_dt_2, covid_dt_3)
rm(covid_dt_1, covid_dt_2, covid_dt_3)

#now do for neg - read in the data
gp_dt_1 <- read_parquet("observations_gp_visits_tmp1.2.parquet")
gp_dt_2 <- read_parquet("observations_gp_visits_tmp2.2.parquet")
gp_dt_3 <- read_parquet("observations_gp_visits_tmp3.2.parquet")

#exclude visits after 1 jan 2020
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  return(dt)
}
gp_dt_1 <- remove_rows(gp_dt_1)
gp_dt_2 <- remove_rows(gp_dt_2)
gp_dt_3 <- remove_rows(gp_dt_3)

#then process the data
gp_dt_1 <- process_dataset(gp_dt_1)
gp_dt_2 <- process_dataset(gp_dt_2)
gp_dt_3 <- process_dataset(gp_dt_3)

#then process for neg
process_neg_data <- function(dataset) {
  dataset[negcontrolVE == TRUE, negVE_gpvisits := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/01/2020", format = "%d/%m/%Y")]
  dataset <- dataset[negVE_gpvisits == TRUE]
  dataset[, `:=` (staffid = NULL, conssourceid = NULL, consmedcodeid = NULL, cprdconstype = NULL, consdate = NULL, enterdate = NULL, eventdate = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL)]
  return(dataset)
}
neg_dt_1 <- process_neg_data(gp_dt_1)
neg_dt_2 <- process_neg_data(gp_dt_2)
neg_dt_3 <- process_neg_data(gp_dt_3)

#then unique
neg_dt_1 <- unique(neg_dt_1)
neg_dt_2 <- unique(neg_dt_2)
neg_dt_3 <- unique(neg_dt_3)

#then format the data
format_data <- function(data){
  data[, negVE_gpvisits := ifelse(negVE_gpvisits == TRUE, 1, 0)]
  data[is.na(data)] <- 0
  return(data)
}
neg_dt_1 <- format_data(neg_dt_1)
neg_dt_2 <- format_data(neg_dt_2)
neg_dt_3 <- format_data(neg_dt_3)

#then rbind
neg_dt <- rbind(neg_dt_1, neg_dt_2, neg_dt_3)
rm(neg_dt_1, neg_dt_2, neg_dt_3)

#now do for flu - read in the data
gp_dt_1 <- read_parquet("observations_gp_visits_tmp1.2.parquet")
gp_dt_2 <- read_parquet("observations_gp_visits_tmp2.2.parquet")
gp_dt_3 <- read_parquet("observations_gp_visits_tmp3.2.parquet")

#exclude visits after 1 jan 2020
remove_rows <- function(dt){
  dt <- dt[eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  return(dt)
}
gp_dt_1 <- remove_rows(gp_dt_1)
gp_dt_2 <- remove_rows(gp_dt_2)
gp_dt_3 <- remove_rows(gp_dt_3)

#then process the data
gp_dt_1 <- process_dataset(gp_dt_1)
gp_dt_2 <- process_dataset(gp_dt_2)
gp_dt_3 <- process_dataset(gp_dt_3)

#then process for neg
process_flu_data <- function(dataset) {
  dataset[influenzaVE == TRUE, fluVE_gpvisits := eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
  dataset <- dataset[fluVE_gpvisits == TRUE]
  dataset[, `:=` (staffid = NULL, conssourceid = NULL, consmedcodeid = NULL, cprdconstype = NULL, consdate = NULL, enterdate = NULL, eventdate = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL)]
  return(dataset)
}
flu_dt_1 <- process_flu_data(gp_dt_1)
flu_dt_2 <- process_flu_data(gp_dt_2)
flu_dt_3 <- process_flu_data(gp_dt_3)

#then unique
flu_dt_1 <- unique(flu_dt_1)
flu_dt_2 <- unique(flu_dt_2)
flu_dt_3 <- unique(flu_dt_3)

#then format the data
format_data <- function(data){
  data[, fluVE_gpvisits := ifelse(fluVE_gpvisits == TRUE, 1, 0)]
  data[is.na(data)] <- 0
  return(data)
}
flu_dt_1 <- format_data(flu_dt_1)
flu_dt_2 <- format_data(flu_dt_2)
flu_dt_3 <- format_data(flu_dt_3)

#then rbind
flu_dt <- rbind(flu_dt_1, flu_dt_2, flu_dt_3)
rm(flu_dt_1, flu_dt_2, flu_dt_3)

#write these for merging
setwd(datafiles)
write_parquet(covid_dt, paste0(datafiles, "covid_gp_visits.parquet"))
write_parquet(neg_dt, paste0(datafiles, "neg_gp_visits.parquet"))
write_parquet(flu_dt, paste0(datafiles, "flu_gp_visits.parquet"))
#################GP VISITS - PAPER 2 ADDITIONAL INVESTIGATION###################
setwd(datafiles)
gp_dt_1 <- read_parquet("observations_gp_visits_tmp1.2.parquet")
gp_dt_2 <- read_parquet("observations_gp_visits_tmp2.2.parquet")
gp_dt_3 <- read_parquet("observations_gp_visits_tmp3.2.parquet")

#remove unwanted columns 
gp_dt_1[, `:=` (staffid = NULL, conssourceid = NULL, consmedcodeid = NULL, cprdconstype = NULL, consdate = NULL, enterdate = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL)]
gp_dt_2[, `:=` (staffid = NULL, conssourceid = NULL, consmedcodeid = NULL, cprdconstype = NULL, consdate = NULL, enterdate = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL)]
gp_dt_3[, `:=` (staffid = NULL, conssourceid = NULL, consmedcodeid = NULL, cprdconstype = NULL, consdate = NULL, enterdate = NULL, dob = NULL, covidVE = NULL, influenzaVE = NULL, negcontrolVE = NULL)]

#unique
gp_dt_1 <- unique(gp_dt_1)
gp_dt_2 <- unique(gp_dt_2)
gp_dt_3 <- unique(gp_dt_3)

#create flu specific and merge
flu_dt_1 <- gp_dt_1[eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
flu_dt_2 <- gp_dt_2[eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
flu_dt_3 <- gp_dt_3[eventdate >= as.Date("01/09/2018", format = "%d/%m/%Y") & eventdate <= as.Date("01/09/2019", format = "%d/%m/%Y")]
flu_dt <- rbind(flu_dt_1, flu_dt_2, flu_dt_3)
rm(flu_dt_1, flu_dt_2, flu_dt_3)

#how many visits per patient
flu_dt[, n_visit := seq_len(.N), by = "patid"]
flu_dt[, totalvisit := .N, by = c("patid")]
janitor::tabyl(flu_dt$totalvisit) %>% adorn_pct_formatting()

#make this into one row per patient
flu_wide <- dcast(flu_dt, patid + totalvisit ~ n_visit,
                  fun.aggregate = NULL,
                  value.var = "eventdate")

#remove the dates - not needed
flu_wide <- flu_wide[, -c(3:233)]

#check number of visits per patient 
janitor::tabyl(flu_wide$totalvisit) %>% adorn_pct_formatting()

#write for later use
write_parquet(flu_wide, paste0(datafiles, "flu_gp_visits_paper2.parquet"))
################################################################################

################MERGE WITH PATIENT FILE####################
#now merge with the patient file - will create three populations 
covid_population <- read_parquet("study_population_covid.parquet")
flu_population <- read_parquet("study_population_flu.parquet")
neg_population <- read_parquet("study_population_neg.parquet")

#pull in AA data
aa_dt <- read_parquet("aa_markers.parquet")

#add the right columns to each
covid_population[aa_dt, `:=`(COVIDVE_AAA_narrow = i.COVIDVE_AAA_narrow, COVIDVE_AAA_broad = i.COVIDVE_AAA_broad), on = "patid"]
flu_population[aa_dt, `:=`(fluVE_AAA_narrow = i.fluVE_AAA_narrow, fluVE_AAA_broad = i.fluVE_AAA_broad), on = "patid"]
neg_population[aa_dt, `:=`(negcon_AAA_narrow = i.negcon_AAA_narrow, negcon_AAA_broad = i.negcon_AAA_broad), on = "patid"]
rm(aa_dt)

#pull in short look back 
covid_dt <- read_parquet("covid_pt1markers.parquet")
flu_dt <- read_parquet("flu_pt1markers.parquet")
neg_dt <- read_parquet("neg_pt1markers.parquet")

#add the right columns to each
covid_population[covid_dt, `:=`(COVIDVE_breast_narrow = i.COVIDVE_breast_narrow, COVIDVE_breast_broad = i.COVIDVE_breast_broad, COVIDVE_bowel = i.COVIDVE_bowel, COVIDVE_cervical_narrow = i.COVIDVE_cervical_narrow, COVIDVE_cervical_broad = i.COVIDVE_cervical_broad, COVIDVE_PSA = i.COVIDVE_PSA, COVIDVE_bone = i.COVIDVE_bone, COVIDVE_NHS = i.COVIDVE_NHS), on = "patid"]
flu_population[flu_dt, `:=`(fluVE_breast_narrow = i.fluVE_breast_narrow, fluVE_breast_broad = i.fluVE_breast_broad, fluVE_bowel = i.fluVE_bowel, fluVE_cervical_narrow = i.fluVE_cervical_narrow, fluVE_cervical_broad = i.fluVE_cervical_broad, fluVE_PSA  = i.fluVE_PSA, fluVE_bone = i.fluVE_bone, fluVE_NHS = i.fluVE_NHS), on = "patid"]
neg_population[neg_dt, `:=`(negVE_breast_narrow = i.negVE_breast_narrow, negVE_breast_broad = i.negVE_breast_broad, negVE_bowel = i.negVE_bowel, negVE_cervical_narrow = i.negVE_cervical_narrow, negVE_cervical_broad = i.negVE_cervical_broad, negVE_PSA = i.negVE_PSA, negVE_bone = i.negVE_bone, negVE_NHS = i.negVE_NHS), on = "patid"]
rm(covid_dt, flu_dt, neg_dt)

#pull in all short back - restrictive 
covid_dt <- read_parquet("covid_pt1markers_res.parquet")
flu_dt <- read_parquet("flu_pt1markers_res.parquet")
neg_dt <- read_parquet("neg_pt1markers_res.parquet")

#add the right columns to each
covid_population[covid_dt, `:=`(COVIDVE_breast_narrow_res = i.COVIDVE_breast_narrow_res, COVIDVE_breast_broad_res = i.COVIDVE_breast_broad_res, COVIDVE_bowel_res = i.COVIDVE_bowel_res, COVIDVE_cervical_narrow_res = i.COVIDVE_cervical_narrow_res, COVIDVE_cervical_broad_res = i.COVIDVE_cervical_broad_res, COVIDVE_NHS_res = i.COVIDVE_NHS_res), on = "patid"]
flu_population[flu_dt, `:=`(fluVE_breast_narrow_res = i.fluVE_breast_narrow_res, fluVE_breast_broad_res = i.fluVE_breast_broad_res, fluVE_bowel_res = i.fluVE_bowel_res, fluVE_cervical_narrow_res = i.fluVE_cervical_narrow_res, fluVE_cervical_broad_res = i.fluVE_cervical_broad_res, fluVE_NHS_res = i.fluVE_NHS_res), on = "patid"]
neg_population[neg_dt, `:=`(negVE_breast_narrow_res = i.negVE_breast_narrow_res, negVE_breast_broad_res = i.negVE_breast_broad_res, negVE_bowel_res = i.negVE_bowel_res, negVE_cervical_narrow_res = i.negVE_cervical_narrow_res, negVE_cervical_broad_res = i.negVE_cervical_broad_res, negVE_NHS_res = i.negVE_NHS_res), on = "patid"]
rm(covid_dt, flu_dt, neg_dt)

#pull in primary DNA
covid_dt <- read_parquet("covid_DNAmarkers.parquet")
flu_dt <- read_parquet("flu_DNAmarkers.parquet")
neg_dt <- read_parquet("neg_DNAmarkers.parquet")

#add the right columns to each 
covid_population[covid_dt, `:=`(COVID_VE_Narrow_DNA = i.COVID_VE_Narrow_DNA), on = "patid"]
flu_population[flu_dt, `:=`(flu_VE_Narrow_DNA = i.flu_VE_Narrow_DNA), on = "patid"]
neg_population[neg_dt, `:=`(neg_VE_Narrow_DNA = i.neg_VE_Narrow_DNA), on = "patid"]
rm(covid_dt, flu_dt, neg_dt)

#pull in blood pressure 
covid_dt <- read_parquet("covid_BPmarkers.parquet")
flu_dt <- read_parquet("flu_BPmarkers.parquet")
neg_dt <- read_parquet("neg_BPmarkers.parquet")

#add the right columns to each
covid_population[covid_dt, `:=`(COVID_VE_Narrow_BP = i.COVID_VE_Narrow_BP), on = "patid"]
flu_population[flu_dt, `:=`(flu_VE_Narrow_BP = i.flu_VE_Narrow_BP), on = "patid"]
neg_population[neg_dt, `:=`(neg_VE_Narrow_BP = i.neg_VE_Narrow_BP), on = "patid"]
rm(covid_dt, flu_dt, neg_dt)

#pull in pneunococcal vac
pneu_dt <- read_parquet("observations_markers_pneu.parquet")

#add the right columns to each
covid_population[pneu_dt, `:=`(COVID_VE_Narrow_pneu = i.COVID_VE_Narrow_pneu), on = "patid"]
flu_population[pneu_dt, `:=`(flu_VE_Narrow_pneu = i.flu_VE_Narrow_pneu), on = "patid"]
neg_population[pneu_dt, `:=`(neg_VE_Narrow_pneu = i.neg_VE_Narrow_pneu), on = "patid"]
rm(pneu_dt)

#add in flu 
flu_dt <- read_parquet("observations_markers_flu.parquet")

#add to populations
covid_population[flu_dt, `:=`(COVIDVE_flu = i.COVIDVE_flu), on = "patid"]
flu_population[flu_dt, `:=`(fluVE_flu = i.fluVE_flu), on = "patid"]
neg_population[flu_dt, `:=`(negVE_flu = i.negVE_flu), on = "patid"]
rm(flu_dt)

#pull in low value prescriptions 
covid_dt <- read_parquet("covid_lowpresmarkers.parquet")
flu_dt <- read_parquet("flu_lowpresmarkers.parquet")
neg_dt <- read_parquet("neg_lowpresmarkers.parquet")

#add the right columns to each 
covid_population[covid_dt, `:=`(COVID_VE_Narrow_lowpres = i.COVID_VE_Narrow_lowpres), on = "patid"]
flu_population[flu_dt, `:=`(flu_VE_Narrow_lowpres = i.flu_VE_Narrow_lowpres), on = "patid"]
neg_population[neg_dt, `:=`(neg_VE_Narrow_lowpres = i.neg_VE_Narrow_lowpres), on = "patid"]
rm(covid_dt, flu_dt, neg_dt)

#pull in ACS
covid_dt <- read_parquet("covid_ACS.parquet")
flu_dt <- read_parquet("flu_ACS.parquet")
neg_dt <- read_parquet("neg_ACS.parquet")

#add the right columns to each 
covid_population[covid_dt, `:=`(COVID_VE_Narrow_ACS = i.COVID_VE_Narrow_ACS), on = "patid"]
flu_population[flu_dt, `:=`(flu_VE_Narrow_ACS = i.flu_VE_Narrow_ACS), on = "patid"]
neg_population[neg_dt, `:=`(neg_VE_Narrow_ACS = i.neg_VE_Narrow_ACS), on = "patid"]
rm(covid_dt, flu_dt, neg_dt)

#pull in low value procedures 
covid_dt <- read_parquet("covid_lowvalueproc.parquet")
flu_dt <- read_parquet("flu_lowvalueproc.parquet")
neg_dt <- read_parquet("neg_lowvalueproc.parquet")

#add the right columns to each 
covid_population[covid_dt, `:=`(COVID_VE_Narrow_lowvalueproc = i.COVID_VE_Narrow_lowvalueproc), on = "patid"]
flu_population[flu_dt, `:=`(flu_VE_Narrow_lowvalueproc = i.flu_VE_Narrow_lowvalueproc), on = "patid"]
neg_population[neg_dt, `:=`(neg_VE_Narrow_lowvalueproc= i.neg_VE_Narrow_lowvalueproc), on = "patid"]
rm(covid_dt, flu_dt, neg_dt)

#pull in gp visits 
covid_dt <- read_parquet("covid_gp_visits.parquet")
neg_dt <- read_parquet("neg_gp_visits.parquet")
flu_dt <- read_parquet("flu_gp_visits.parquet")

#add the right columns to each
covid_population[covid_dt, `:=`(COVIDVE_gpvisits = i.COVIDVE_gpvisits), on = "patid"]
flu_population[flu_dt, `:=`(fluVE_gpvisits = i.fluVE_gpvisits), on = "patid"]
neg_population[neg_dt, `:=`(negVE_gpvisits= i.negVE_gpvisits), on = "patid"]
rm(covid_dt, flu_dt, neg_dt)

#need to edit the names of the columns for reading in 
colnames(covid_population) <- gsub("COVID_VE_", "", colnames(covid_population))
colnames(covid_population) <- gsub("COVIDVE_", "", colnames(covid_population))
colnames(flu_population) <- gsub("flu_VE_", "", colnames(flu_population))
colnames(flu_population) <- gsub("fluVE_", "", colnames(flu_population))
colnames(neg_population) <- gsub("neg_VE_", "", colnames(neg_population))
colnames(neg_population) <- gsub("negVE_", "", colnames(neg_population))
colnames(neg_population) <- gsub("negcon_", "", colnames(neg_population))

#update NA to 0 in marker columns 
na_to_zero <- function(dataset) {
  dataset$AAA_narrow[is.na(dataset$AAA_narrow)] <- 0
  dataset$AAA_broad[is.na(dataset$AAA_broad)] <- 0
  dataset$breast_narrow[is.na(dataset$breast_narrow)] <- 0
  dataset$cervical_broad[is.na(dataset$cervical_broad)] <- 0
  dataset$NHS[is.na(dataset$NHS)] <- 0
  dataset$breast_narrow_res[is.na(dataset$breast_narrow_res)] <- 0
  dataset$cervical_broad_res[is.na(dataset$cervical_broad_res)] <- 0
  dataset$NHS_res[is.na(dataset$NHS_res)] <- 0
  dataset$PSA[is.na(dataset$PSA)] <- 0
  dataset$bone[is.na(dataset$bone)] <- 0
  dataset$bowel[is.na(dataset$bowel)] <- 0
  dataset$bowel_res[is.na(dataset$bowel_res)] <- 0
  dataset$gpvisits[is.na(dataset$gpvisits)] <- 0
  dataset$breast_broad[is.na(dataset$breast_broad)] <- 0
  dataset$cervical_narrow[is.na(dataset$cervical_narrow)] <- 0
  dataset$breast_broad_res[is.na(dataset$breast_broad_res)] <- 0
  dataset$cervical_narrow_res[is.na(dataset$cervical_narrow_res)] <- 0
  dataset$Narrow_DNA[is.na(dataset$Narrow_DNA)] <- 0
  dataset$Narrow_BP[is.na(dataset$Narrow_BP)] <- 0
  dataset$Narrow_pneu[is.na(dataset$Narrow_pneu)] <- 0
  dataset$flu[is.na(dataset$flu)] <- 0
  dataset$Narrow_lowpres[is.na(dataset$Narrow_lowpres)] <- 0
  dataset$Narrow_ACS[is.na(dataset$Narrow_ACS)] <- 0
  dataset$Narrow_lowvalueproc[is.na(dataset$Narrow_lowvalueproc)] <- 0
  return(dataset)
}
covid_population <- na_to_zero(covid_population)
flu_population <- na_to_zero(flu_population)
neg_population <- na_to_zero(neg_population)

#before we save restrict to those with HES and ONS linkage 
covid_population <- covid_population[hes_apc_e == 1]#
flu_population <- flu_population[hes_apc_e == 1]#
neg_population <- neg_population[hes_apc_e == 1]#

#save these populations for descriptive tables 
setwd(datafiles)
write_parquet(covid_population, paste0(datafiles, "covid_marker_pop.parquet"))
write_parquet(flu_population, paste0(datafiles, "flu_marker_pop.parquet"))
write_parquet(neg_population, paste0(datafiles, "neg_marker_pop.parquet"))