################################################################################
# Author: Sophie Graham
# Date: 04/07/2023
# Version: R 4.2.2
# File name: healthseekingbehaviour_descriptive_markers.R
# Status: in progress
# CPRD version: May 2022
# Data sets used: CPRD Aurum and HES APC
# R scripts needed: global (not included as has links to secure severs), data 
# prep new (not included as merging code lists with raw data files) and markers 
# Description of file: create a prevalence table of all markers and also 
# stratifies this by gender. Creates dendrograms to show clusters of markers in 
# the data. Creates plots that stratify prevalence of each marker by age 
# category. 
################################################################################

#set rounding and redaction thresholds
rounding_threshold = 1 
redaction_threshold = 5

#create output directory
fs::dir_create(here::here(results))

#Import data 
setwd(datafiles)
covid_population <- read_parquet("covid_marker_pop.parquet")
flu_population <- read_parquet("flu_marker_pop.parquet")
neg_population <- read_parquet("neg_marker_pop.parquet")

#from the covid population we need to exclude those that started their vaccination before 8 December 2020
covid_exclude <- readr::read_tsv("covid_exclude.txt", col_type = cols(.default = col_character())) 
excluded_patids <- as.character(covid_exclude$patid)
covid_population <- covid_population[!patid %in% excluded_patids]

#relabel gender
relabel_gender <- function(dataset) {
  dataset$gender <- recode(dataset$gender,
                           "1" = "male",
                           "2" = "female",
                           "3" = "indeterminate",
                           "4" = "unknown")
  return(dataset)
}
covid_population <- relabel_gender(covid_population)
flu_population <- relabel_gender(flu_population)
neg_population <- relabel_gender(neg_population)

#split the populations into male and female
covid_population_male <- covid_population[gender == c("male")]
covid_population_female <- covid_population[gender == c("female")]
flu_population_male <- flu_population[gender == c("male")]
flu_population_female <- flu_population[gender == c("female")]
neg_population_male <- neg_population[gender == c("male")]
neg_population_female <- neg_population[gender == c("female")]

#Format data 
format_data <- function(data) {
  data <- data %>%
    mutate(
      N = 1,
      allpop = "All",
    )
  return(data)
}
covid_population <- format_data(covid_population)
flu_population <- format_data(flu_population)
neg_population <- format_data(neg_population)
covid_population_male <- format_data(covid_population_male)
covid_population_female <- format_data(covid_population_female)
flu_population_male <- format_data(flu_population_male)
flu_population_female <- format_data(flu_population_female)
neg_population_male <- format_data(neg_population_male)
neg_population_female <- format_data(neg_population_female)

#define variables of interest
define_variables <- function(data) {
  counts <- data %>% 
    select(
      N,
      allpop,
      ## Demographics
      gender,
      ## Markers
      AAA_narrow,
      AAA_broad,
      breast_narrow,
      breast_narrow_res,
      breast_broad,
      breast_broad_res,
      cervical_narrow,
      cervical_narrow_res,
      cervical_broad,
      cervical_broad_res,
      NHS,
      NHS_res,
      PSA,
      bone,
      bowel,
      bowel_res,
      Narrow_DNA,
      Narrow_BP,
      Narrow_pneu,
      flu, 
      Narrow_lowpres,
      Narrow_ACS,
      Narrow_lowvalueproc,
      gpvisits
    )
  return(counts)
}
counts_covid <- define_variables(covid_population)
counts_flu <- define_variables(flu_population)
counts_neg <- define_variables(neg_population)
counts_covid_male <- define_variables(covid_population_male)
counts_covid_female <- define_variables(covid_population_female)
counts_flu_male <- define_variables(flu_population_male)
counts_flu_female <- define_variables(flu_population_female)
counts_neg_male <- define_variables(neg_population_male)
counts_neg_female <- define_variables(neg_population_female)

#clean table names
clean_table_names <- function(input_table) {
  # Relabel variables for plotting
  input_table$Variable[input_table$Variable=="AAA_narrow"] = "AAA screen narrow"
  input_table$Variable[input_table$Variable=="AAA_broad"] = "AAA screen broad"
  input_table$Variable[input_table$Variable=="breast_narrow"] = "Breast screen narrow"
  input_table$Variable[input_table$Variable=="breast_broad"] = "Breast screen broad"
  input_table$Variable[input_table$Variable=="cervical_narrow"] = "Cervical screen narrow"
  input_table$Variable[input_table$Variable=="cervical_broad"] = "Cervical screen broad"
  input_table$Variable[input_table$Variable=="NHS"] = "NHS health checks"
  input_table$Variable[input_table$Variable=="AAA_narrow_res"] = "AAA screen narrow restrictive"
  input_table$Variable[input_table$Variable=="AAA_broad_res"] = "AAA screen broad restrictive"
  input_table$Variable[input_table$Variable=="breast_narrow_res"] = "Breast screen narrow restrictive"
  input_table$Variable[input_table$Variable=="breast_broad_res"] = "Breast screen broad restrictive"
  input_table$Variable[input_table$Variable=="cervical_narrow_res"] = "Cervical screen narrow restrictive"
  input_table$Variable[input_table$Variable=="cervical_broad_res"] = "Cervical screen broad restrictive"
  input_table$Variable[input_table$Variable=="NHS_res"] = "NHS health checks restrictive"
  input_table$Variable[input_table$Variable=="PSA"] = "PSA test"
  input_table$Variable[input_table$Variable=="bowel"] = "Bowel screen"
  input_table$Variable[input_table$Variable=="bowel_res"] = "Bowel screen restrictive"
  input_table$Variable[input_table$Variable=="Narrow_DNA"] = "Primary care DNA"
  input_table$Variable[input_table$Variable=="Narrow_BP"] = "Blood pressure test"
  input_table$Variable[input_table$Variable=="Narrow_pneu"] = "Pneumococcal vaccine"
  input_table$Variable[input_table$Variable=="flu"] = "Flu vaccine"
  input_table$Variable[input_table$Variable=="Narrow_lowpres"] = "Low value prescriptions"
  input_table$Variable[input_table$Variable=="Narrow_ACS"] = "Ambulatory care sensitive"
  input_table$Variable[input_table$Variable=="Narrow_lowvalueproc"] = "Low value procedures"
  input_table$Variable[input_table$Variable=="bone"] = "Bone scan"
  input_table$Variable[input_table$Variable=="gpvisits"] = "GP visits"
  input_table$Group[input_table$Group=="gender"] = "gender"
  input_table$Group[input_table$Variable=="N"] = "N"
  return(input_table)
}
#Generate full and stratified table
levels = c("All")

#Generate table - full and stratified populations
generate_summary_table <- function(counts, levels, gender_var, allpop_var, rounding_threshold, redaction_threshold) {
  data_cohort <- counts
  collated_table <- NULL
  for (i in 1:length(levels)) {
    if (i == 1) { 
      data_subset <- counts
      counts_summary <- data_subset %>% 
        tbl_summary(by = allpop_var)
      counts_summary$inputs$data <- NULL
    } else { 
      data_subset <- subset(counts, !!sym(gender_var) == levels[i]) 
      if (nrow(data_subset) == 0) {
        next 
      }
      counts_summary <- data_subset %>% 
        select(-!!sym(gender_var)) %>% 
        tbl_summary(by = allpop_var)
      counts_summary$inputs$data <- NULL
    }
    table1 <- counts_summary$table_body %>%
      select(group = variable, variable = label, count = stat_1) %>%
      separate(count, c("count","perc"), sep = "([(])") %>%
      mutate(count = gsub(" ", "", count)) %>%
      mutate(count = as.numeric(gsub(",", "", count))) %>%
      filter(!(is.na(count))) %>%
      select(-perc)
    table1$percent <- round(table1$count/nrow(data_cohort)*100, 1)
    colnames(table1) <- c("Group", "Variable", "Count", "Percent")
    ## Clean names
    table1_clean <- clean_table_names(table1)
    ## Calculate rounded total
    rounded_n <- plyr::round_any(nrow(data_subset), rounding_threshold)
    ## Round individual values to rounding threshold
    table1_redacted <- table1_clean %>%
      mutate(Count = plyr::round_any(Count, rounding_threshold))
    table1_redacted$Percent <- round(table1_redacted$Count/rounded_n*100, 1)
    table1_redacted$Non_Count <- rounded_n - table1_redacted$Count
    ## Redact any rows with rounded cell counts or non-counts <= redaction threshold 
    table1_redacted$Summary <- paste0(prettyNum(table1_redacted$Count, big.mark = ","), " (", format(table1_redacted$Percent, nsmall = 1), "%)")
    table1_redacted$Summary <- gsub(" ", "", table1_redacted$Summary, fixed = TRUE) # Remove spaces generated by decimal formatting
    table1_redacted$Summary <- gsub("(", " (", table1_redacted$Summary, fixed = TRUE) # Add first space before (
    table1_redacted$Summary[(table1_redacted$Count > 0 & table1_redacted$Count <= redaction_threshold) | (table1_redacted$Non_Count > 0 & table1_redacted$Non_Count <= redaction_threshold)] <- "[Redacted]"
    table1_redacted$Summary[table1_redacted$Variable == "N"] <- prettyNum(table1_redacted$Count[table1_redacted$Variable == "N"], big.mark = ",")
    table1_redacted <- table1_redacted %>% select(-Non_Count, -Count, -Percent)
    names(table1_redacted)[3] <- levels[i]
    if (i == 1) { 
      collated_table <- table1_redacted 
    } else { 
      collated_table <- collated_table %>% left_join(table1_redacted[,2:3], by = "Variable") 
      collated_table[,i+2][is.na(collated_table[,i+2])] <- "--"
    }
  }
  return(collated_table)
}
summary_table_covid <- generate_summary_table(counts_covid, levels, "gender", "allpop", rounding_threshold, redaction_threshold)
summary_table_flu <- generate_summary_table(counts_flu, levels, "gender", "allpop", rounding_threshold, redaction_threshold)
summary_table_neg <- generate_summary_table(counts_neg, levels, "gender", "allpop", rounding_threshold, redaction_threshold)
summary_table_covid_male <- generate_summary_table(counts_covid_male, levels, "gender", "allpop", rounding_threshold, redaction_threshold)
summary_table_covid_female <- generate_summary_table(counts_covid_female, levels, "gender", "allpop", rounding_threshold, redaction_threshold)
summary_table_flu_male <- generate_summary_table(counts_flu_male, levels, "gender", "allpop", rounding_threshold, redaction_threshold)
summary_table_flu_female <- generate_summary_table(counts_flu_female, levels, "gender", "allpop", rounding_threshold, redaction_threshold)
summary_table_neg_male <- generate_summary_table(counts_neg_male, levels, "gender", "allpop", rounding_threshold, redaction_threshold)
summary_table_neg_female <- generate_summary_table(counts_neg_female, levels, "gender", "allpop", rounding_threshold, redaction_threshold)

#perform a left join
perform_left_join <- function(data, data_male, data_female) {
  result <- data %>%
    left_join(data_male, by = "Variable") %>%
    left_join(data_female, by = "Variable") %>%
    select(Variable, All_Total = All.x, Male = All.y, Female = All) %>%
    rename(All = All_Total)
  
  return(result)
}
covid_results <- perform_left_join(summary_table_covid, summary_table_covid_male, summary_table_covid_female)
flu_results <- perform_left_join(summary_table_flu, summary_table_flu_male, summary_table_flu_female)
neg_results <- perform_left_join(summary_table_neg, summary_table_neg_male, summary_table_neg_female)

#save as csv 
setwd(results)
data_tables <- list (covid_results, flu_results, neg_results)
write_csv <- function(data, results) {
  for (i in seq_along(data)) {
    file_name <- switch(i,
                        "1" = "covid",
                        "2" = "flu",
                        "3" = "neg")
    write.csv(data[[i]], paste0(results, "/table1_markers", file_name, ".csv"))
  }
}
write_csv(data_tables, results)

#Correlation matrix - have not included restrictive 
#prep the data
setwd(datafiles)
covid_population <- read_parquet("covid_marker_pop.parquet")
flu_population <- read_parquet("flu_marker_pop.parquet")
neg_population <- read_parquet("neg_marker_pop.parquet")
flu_population_male <- flu_population[gender == c("1")]
flu_population_female <- flu_population[gender == c("2")]

#counts
define_variables <- function(data) {
  counts <- data %>% 
    select(
      ## Markers
      AAA_narrow,
      AAA_broad,
      breast_narrow,
      breast_broad,
      cervical_narrow,
      cervical_broad,
      NHS, 
      PSA,
      bone,
      bowel,
      Narrow_DNA,
      Narrow_BP,
      Narrow_pneu,
      flu, 
      Narrow_lowpres,
      Narrow_ACS,
      Narrow_lowvalueproc,
      gpvisits
    )
  return(counts)
}
counts_covid <- define_variables(covid_population)
counts_flu <- define_variables(flu_population)
counts_neg <- define_variables(neg_population)
counts_flu_male <- define_variables(flu_population_male)
counts_flu_female <- define_variables(flu_population_female)

#ggcorrplot for manuscript
#remove broad columns in counts
counts_flu[, AAA_broad := NULL]
counts_flu[, AAA_narrow := NULL]
counts_flu[, breast_broad := NULL]
counts_flu[, breast_narrow := NULL]
counts_flu[, cervical_broad := NULL]
counts_flu[, cervical_narrow := NULL]
counts_flu[, PSA := NULL]
counts_flu[, Narrow_lowpres := NULL]
counts_flu_female[, AAA_broad := NULL]
counts_flu_female[, AAA_narrow := NULL]
counts_flu_female[, PSA := NULL]
counts_flu_female[, breast_broad := NULL]
counts_flu_female[, cervical_broad := NULL]
counts_flu_female[, Narrow_lowpres := NULL]
counts_flu_male[, AAA_broad := NULL]
counts_flu_male[, breast_broad := NULL]
counts_flu_male[, breast_narrow := NULL]
counts_flu_male[, cervical_broad := NULL]
counts_flu_male[, cervical_narrow := NULL]
counts_flu_male[, Narrow_lowpres := NULL]

#Dendrogram
phi_matrix_flu <- rcorr(as.matrix(counts_flu))
phi_matrix_flu_male <- rcorr(as.matrix(counts_flu_male))
phi_matrix_flu_female <- rcorr(as.matrix(counts_flu_female))
a <- cor(counts_flu)
b <- cor(counts_flu_male)
c <- cor(counts_flu_female)
anorm <- abs(a)
bnorm <- abs(b)
cnorm <- abs(c)

dendr_flu <- heatmaply_cor(
  a,
  node_type = "scatter",
  point_size_mat = anorm, 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation")
)

dendr_flu_males <- heatmaply_cor(
  b,
  node_type = "scatter",
  point_size_mat = bnorm, 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation")
)

dendr_flu_females <- heatmaply_cor(
  c,
  node_type = "scatter",
  point_size_mat = cnorm, 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation")
)

#print
#render and save at 900 by 600 pixels using R studio viewer

#write phi coefficients for supplement in paper 2
phi_flu <- as.data.frame(a)
phi_flu_male <- as.data.frame(b)
phi_flu_female <- as.data.frame(c)

phi_flu[1, 2:10] <- NA
phi_flu[2, 3:10] <- NA
phi_flu[3, 4:10] <- NA
phi_flu[4, 5:10] <- NA
phi_flu[5, 6:10] <- NA
phi_flu[6, 7:10] <- NA
phi_flu[7, 8:10] <- NA
phi_flu[8, 9:10] <- NA
phi_flu[9, 10] <- NA

phi_flu_male[1, 2:12] <- NA
phi_flu_male[2, 3:12] <- NA
phi_flu_male[3, 4:12] <- NA
phi_flu_male[4, 5:12] <- NA
phi_flu_male[5, 6:12] <- NA
phi_flu_male[6, 7:12] <- NA
phi_flu_male[7, 8:12] <- NA
phi_flu_male[8, 9:12] <- NA
phi_flu_male[9, 10:12] <- NA
phi_flu_male[10, 11:12] <- NA
phi_flu_male[11, 12] <- NA

phi_flu_female[1, 2:12] <- NA
phi_flu_female[2, 3:12] <- NA
phi_flu_female[3, 4:12] <- NA
phi_flu_female[4, 5:12] <- NA
phi_flu_female[5, 6:12] <- NA
phi_flu_female[6, 7:12] <- NA
phi_flu_female[7, 8:12] <- NA
phi_flu_female[8, 9:12] <- NA
phi_flu_female[9, 10:12] <- NA
phi_flu_female[10, 11:12] <- NA
phi_flu_female[11, 12] <- NA

phi_flu <- round(phi_flu, digits = 2)
phi_flu_male <- round(phi_flu_male, digits = 2)
phi_flu_female <- round(phi_flu_female, digits = 2)

#write
setwd(results)
write.csv(phi_flu, "Supplementary_material_paper2_phi_flu.csv")
write.csv(phi_flu_male, "Supplementary_material_paper2_phi_flu_male.csv")
write.csv(phi_flu_female, "Supplementary_material_paper2_phi_flu_female.csv")

######PAPER 2 SUPPLEMENTARY MATERIALS ADDITIONAL ANALYSES#######################
#prevalence by age cat at index

#Import data 
setwd(datafiles)
flu_population <- read_parquet("flu_marker_pop.parquet")

#add in dob
patient <- read_parquet("study_population.parquet")
flu_population[patient, `:=`(dob = i.dob), on = "patid"]

#create age at index variable
calculate_age <- function(dob, reference_date) {
  age <- as.integer(difftime(reference_date, dob, units = "days") / 365.25)
  return(age)
}
reference_date_flu <- as.Date("2019-09-01", format = "%Y-%m-%d")
flu_population[, age := calculate_age(dob, reference_date_flu)]

age_cat <- function(data){
  data[, age_category := cut(
    age,
    breaks = c(65, 70, 75, 80, 85, 90, 95, Inf),
    labels = c("65-69", "70-74", "75-79", "80-84", "85-89", "90-95", "95+"),
    right = FALSE
  )]
  return(data)
}
flu_population <- age_cat(flu_population)

#relabel gender
relabel_gender <- function(dataset) {
  dataset$gender <- recode(dataset$gender,
                           "1" = "male",
                           "2" = "female",
                           "3" = "indeterminate",
                           "4" = "unknown")
  return(dataset)
}
flu_population <- relabel_gender(flu_population)

#split the populations into male and female
flu_population_male <- flu_population[gender == c("male")]
flu_population_female <- flu_population[gender == c("female")]

options(scipen=999)

#calculate prevalence
prev_flu <- bind_rows(
  flu_population_male %>%
    select(AAA_narrow, AAA_broad, PSA, age_category) %>%
    pivot_longer(cols = c('AAA_narrow', 'AAA_broad', 'PSA'),
                 names_to = 'markers', 
                 values_to = 'counts') %>%
    group_by(age_category, markers) %>%
    summarise(count = sum(counts)) %>%
    ungroup() %>%
    group_by(age_category, markers) %>%
    mutate(nrows = sum(flu_population_male$age_category == age_category),
           prev = (count / nrows) * 100),
  
  flu_population_female %>%
    select(breast_narrow, breast_narrow_res, breast_broad, breast_broad_res, cervical_narrow, cervical_narrow_res, cervical_broad, cervical_broad_res, age_category) %>%
    pivot_longer(cols = c('breast_narrow', 'breast_narrow_res', 'breast_broad', 'breast_broad_res', 'cervical_narrow', 'cervical_narrow_res', 'cervical_broad', 'cervical_broad_res'),
                 names_to = 'markers', 
                 values_to = 'counts') %>%
    group_by(age_category, markers) %>%
    summarise(count = sum(counts)) %>%
    ungroup() %>%
    group_by(age_category, markers) %>%
    mutate(nrows = sum(flu_population_female$age_category == age_category),
           prev = (count / nrows) * 100),
  flu_population %>%
    select(bowel, bowel_res, flu, Narrow_pneu, NHS, NHS_res, bone, Narrow_lowvalueproc, Narrow_lowpres, gpvisits, Narrow_DNA, Narrow_ACS, Narrow_BP, age_category) %>%
    pivot_longer(cols = c('bowel', 'bowel_res', 'flu', 'Narrow_pneu', 'NHS', 'NHS_res', 'bone', 'Narrow_lowvalueproc', 'Narrow_lowpres', 'gpvisits', 'Narrow_DNA', 'Narrow_ACS', 'Narrow_BP'),
                 names_to = 'markers',
                 values_to = 'counts') %>%
    group_by(age_category, markers) %>%
    summarise(count = sum(counts)) %>%
    ungroup() %>%
    group_by(age_category, markers) %>%
    mutate(nrows = sum(flu_population$age_category == age_category),
           prev = (count / nrows) * 100)
) %>%
  mutate(
    label = case_when(
      markers %in% c('breast_broad_res', 'cervical_broad_res') ~ 'Broad restrictive',
      markers %in% c('breast_broad', 'cervical_broad', 'AAA_broad') ~ 'Broad standard',
      markers %in% c('breast_narrow_res', 'cervical_narrow_res', 'NHS_res', 'bowel_res') ~ 'Narrow restrictive',
  TRUE ~ 'Narrow standard'))

# Create a grouping variable for "narrow" and "broad" markers 
desired_order <- c(
  "AAA Screening Males", 
  "Breast Cancer Screening Females", 
  "Cervical Cancer Screening Females", 
  "Bowel Cancer Screening", 
  "Influenza Vaccination", 
  "Pneumococcal Vaccination", 
  "NHS Health Checks", 
  "PSA Test Males", 
  "Bone Density Scan", 
  "Low-Value Procedures", 
  "Low-Value Prescriptions", 
  "GP Visits", 
  "DNA Primary Care Visit", 
  "ACS Condition Hospital Visit", 
  "Blood Pressure Test"
)

plot_flu <- prev_flu %>%
  mutate(marker_group = case_when(
    markers %in% c("AAA_narrow", "AAA_broad") ~ "AAA Screening Males",
    markers %in% c("breast_narrow", "breast_narrow_res", "breast_broad", "breast_broad_res") ~ "Breast Cancer Screening Females",
    markers %in% c("cervical_narrow", "cervical_narrow_res", "cervical_broad", "cervical_broad_res") ~ "Cervical Cancer Screening Females",
    markers %in% c("bowel", "bowel_res") ~ "Bowel Cancer Screening",
    markers == "flu" ~ "Influenza Vaccination",
    markers == "Narrow_pneu" ~ "Pneumococcal Vaccination",
    markers %in% c("NHS", "NHS_res") ~ "NHS Health Checks",
    markers == "PSA" ~ "PSA Test Males",
    markers == "bone" ~ "Bone Density Scan",
    markers == "Narrow_lowvalueproc" ~ "Low-Value Procedures",
    markers == "Narrow_lowpres" ~ "Low-Value Prescriptions",
    markers == "gpvisits" ~ "GP Visits",
    markers == "Narrow_DNA" ~ "DNA Primary Care Visit",
    markers == "Narrow_ACS" ~ "ACS Condition Hospital Visit",
    markers == "Narrow_BP" ~ "Blood Pressure Test",
    TRUE ~ markers
  )) %>% 
  mutate(marker_group = factor(marker_group, levels = desired_order)) %>%
  ggplot(., aes(x = age_category, y = prev, fill = label)) +  # Change fill to label
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ marker_group, ncol = 3, scales = "free_y") +
  scale_fill_manual("Legend", values = c("Broad standard" = "#003078", "Broad restrictive" = "#0b0c0c", "Narrow standard" = "#5694ca", "Narrow restrictive" = "#1d70b8")) +
  labs(x = "Age Category",
       y = "Prevalence",
       fill = "Marker") +  
  scale_y_continuous(expand = c(0, 0), breaks = scales::pretty_breaks()) +
  theme_minimal() +
  theme(legend.position = "right",
        text = element_text(size = 10),
        axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
        axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
        panel.grid = element_blank())

#write to pdf
setwd(results)
pdf("bar_chart_by_age_flu_paper_2.pdf", width = 11, height = 7)
print(plot_flu)
dev.off()

#for paper 2 also want to extract raw data for prevalence by age
prev_flu <- as.data.table(prev_flu)
prev_flu <- prev_flu[order(markers)]
marker_names <- data.table(
  old_name = c("AAA_narrow", "AAA_broad", "cervical_narrow", "cervical_narrow_res", "cervical_broad", "cervical_broad_res", "breast_narrow", "breast_narrow_res", "breast_broad", "breast_broad_res", "bowel", "bowel_res", "flu", "Narrow_pneu", "NHS", "NHS_res", "Narrow_BP", "bone", "Narrow_lowvalueproc", "Narrow_lowpres", "gpvisits", "Narrow_DNA", "Narrow_ACS"),
  new_name = c("AAA screening narrow", "AAA screening broad", "Cervical screening narrow standard", "Cervical screening narrow restrictive", "Cervical screening broad standard", "Cervical screening broad restrictive", "Breast screening narrow standard", "Breast screening narrow restrictive", "Breast screening broad standard", "Breast screening broad restrictive", "Bowel cancer screening standard", "Bowel cancer screening restrictive", "Influenza vaccination", "Pneumococcal vaccination", "NHS health checks standard", "NHS health checks restrictive", "Blood pressure test", "Bone density scans", "Low value procedures", "Low value prescriptions", "GP visits", "DNA primary care", "Hospital visit ACS condition")
)
prev_flu[marker_names, markers := i.new_name, on = .(markers = old_name)]
prev_flu[, `:=` (nrows = NULL, label = NULL, marker_group = NULL)]
colnames(prev_flu)[1] = "Age category"
colnames(prev_flu)[2] = "Marker"
colnames(prev_flu)[3] = "Count"
colnames(prev_flu)[4] = "Prevalence"
prev_flu$Prevalence <- round(prev_flu$Prevalence, 1)
prev_flu <- prev_flu[, c("Marker", "Age category", "Count", "Prevalence")]
setwd(results)
write.csv(prev_flu, "Supplementary_material_paper2_prevalence_markers_flu.csv")
################################################################################

#####################PAPER 2 GP ADDITIONAL ANALYSES#############################
#Import data 
setwd(datafiles)
flu_population <- read_parquet("flu_marker_pop.parquet")

#remove marker info not needed
flu_population <- flu_population[, -c(2, 4:33)]

#import number of visits
gp_dt <- read_parquet("flu_gp_visits_paper2.parquet")

#merge total visits into flu pop
flu_population[gp_dt, `:=`(totalvisit = i.totalvisit), on = "patid"]
flu_population[!is.na(totalvisit), .N]
janitor::tabyl(flu_population$gpvisits) %>% adorn_pct_formatting()

#median number of visits
median <- median(flu_population$totalvisit, na.rm = TRUE)
iqr <- quantile(flu_population$totalvisit, na.rm = TRUE)

#histogram GP visits
hist_gp <- ggplot(flu_population, aes(x = totalvisit)) +
  geom_histogram(fill = "#5694ca", color = "white", binwidth = 1, boundary = 20) + 
  labs(x = "Total Visits Per Patient", 
       y = "Frequency") +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), breaks = scales::pretty_breaks()) +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 10),
        axis.line = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
        panel.grid = element_blank(),
        panel.border = element_blank())

# Set the desired x-axis limits
hist_gp <- hist_gp + coord_cartesian(xlim = c(0, 50))

hist_gp

#export
setwd(results)
ggsave("hist_gp_visits_flu_paper_2.pdf", hist_gp, width = 11, height = 7)

#create a variable for over median (7) visits per year
flu_population[, totalvisit := as.numeric(totalvisit)]
flu_population[, overmedianvisits := totalvisit > 7]
flu_population[, overmedianvisits := ifelse(overmedianvisits == TRUE, 1, 0)]

#relabel gender
relabel_gender <- function(dataset) {
  dataset$gender <- recode(dataset$gender,
                           "1" = "male",
                           "2" = "female",
                           "3" = "indeterminate",
                           "4" = "unknown")
  return(dataset)
}
flu_population <- relabel_gender(flu_population)

#split the data
flu_population_male <- flu_population[gender == c("male")]
flu_population_female <- flu_population[gender == c("female")]

#format data
format_data <- function(data) {
  data <- data %>%
    mutate(
      N = 1,
      allpop = "All",
    )
  return(data)
}
flu_population <- format_data(flu_population)
flu_population_male <- format_data(flu_population_male)
flu_population_female <- format_data(flu_population_female)

#Define variables of interest
define_variables <- function(data) {
  counts <- data %>% 
    select(
      N,
      allpop,
      ## Demographics
      gender,
      ## Markers
      overmedianvisits
    )
  return(counts)
}
counts_flu <- define_variables(flu_population)
counts_flu_male <- define_variables(flu_population_male)
counts_flu_female <- define_variables(flu_population_female)

#clean table names
clean_table_names <- function(input_table) {
  # Relabel variables for plotting
  input_table$Variable[input_table$Variable=="overmedianvisits"] = "Over median number of GP visits per year"
  input_table$Group[input_table$Group=="gender"] = "gender"
  input_table$Group[input_table$Variable=="N"] = "N"
  return(input_table)
}

#Generate full and stratified table
levels = c("All")

#Generate table - full and stratified populations
generate_summary_table <- function(counts, levels, gender_var, allpop_var, rounding_threshold, redaction_threshold) {
  data_cohort <- counts
  collated_table <- NULL
  for (i in 1:length(levels)) {
    if (i == 1) { 
      data_subset <- counts
      counts_summary <- data_subset %>% 
        tbl_summary(by = allpop_var)
      counts_summary$inputs$data <- NULL
    } else { 
      data_subset <- subset(counts, !!sym(gender_var) == levels[i]) 
      if (nrow(data_subset) == 0) {
        next 
      }
      counts_summary <- data_subset %>% 
        select(-!!sym(gender_var)) %>% 
        tbl_summary(by = allpop_var)
      counts_summary$inputs$data <- NULL
    }
    table1 <- counts_summary$table_body %>%
      select(group = variable, variable = label, count = stat_1) %>%
      separate(count, c("count","perc"), sep = "([(])") %>%
      mutate(count = gsub(" ", "", count)) %>%
      mutate(count = as.numeric(gsub(",", "", count))) %>%
      filter(!(is.na(count))) %>%
      select(-perc)
    table1$percent <- round(table1$count/nrow(data_cohort)*100, 1)
    colnames(table1) <- c("Group", "Variable", "Count", "Percent")
    ## Clean names
    table1_clean <- clean_table_names(table1)
    ## Calculate rounded total
    rounded_n <- plyr::round_any(nrow(data_subset), rounding_threshold)
    ## Round individual values to rounding threshold
    table1_redacted <- table1_clean %>%
      mutate(Count = plyr::round_any(Count, rounding_threshold))
    table1_redacted$Percent <- round(table1_redacted$Count/rounded_n*100, 1)
    table1_redacted$Non_Count <- rounded_n - table1_redacted$Count
    ## Redact any rows with rounded cell counts or non-counts <= redaction threshold 
    table1_redacted$Summary <- paste0(prettyNum(table1_redacted$Count, big.mark = ","), " (", format(table1_redacted$Percent, nsmall = 1), "%)")
    table1_redacted$Summary <- gsub(" ", "", table1_redacted$Summary, fixed = TRUE) # Remove spaces generated by decimal formatting
    table1_redacted$Summary <- gsub("(", " (", table1_redacted$Summary, fixed = TRUE) # Add first space before (
    table1_redacted$Summary[(table1_redacted$Count > 0 & table1_redacted$Count <= redaction_threshold) | (table1_redacted$Non_Count > 0 & table1_redacted$Non_Count <= redaction_threshold)] <- "[Redacted]"
    table1_redacted$Summary[table1_redacted$Variable == "N"] <- prettyNum(table1_redacted$Count[table1_redacted$Variable == "N"], big.mark = ",")
    table1_redacted <- table1_redacted %>% select(-Non_Count, -Count, -Percent)
    names(table1_redacted)[3] <- levels[i]
    if (i == 1) { 
      collated_table <- table1_redacted 
    } else { 
      collated_table <- collated_table %>% left_join(table1_redacted[,2:3], by = "Variable") 
      collated_table[,i+2][is.na(collated_table[,i+2])] <- "--"
    }
  }
  return(collated_table)
}
summary_table_flu <- generate_summary_table(counts_flu, levels, "gender", "allpop", rounding_threshold, redaction_threshold)
summary_table_flu_male <- generate_summary_table(counts_flu_male, levels, "gender", "allpop", rounding_threshold, redaction_threshold)
summary_table_flu_female <- generate_summary_table(counts_flu_female, levels, "gender", "allpop", rounding_threshold, redaction_threshold)

#perform a left join
perform_left_join <- function(data, data_male, data_female) {
  result <- data %>%
    left_join(data_male, by = "Variable") %>%
    left_join(data_female, by = "Variable") %>%
    select(Variable, All_Total = All.x, Male = All.y, Female = All) %>%
    rename(All = All_Total)
  
  return(result)
}
flu_results <- perform_left_join(summary_table_flu, summary_table_flu_male, summary_table_flu_female)

#save as csv 
setwd(results)
data_tables <- list (flu_results)
write_csv <- function(data, results) {
  for (i in seq_along(data)) {
    file_name <- switch(i,
                        "1" = "flu")
    write.csv(data[[i]], paste0(results, "/table1_markers_gp_visits_paper2_", file_name, ".csv"))
  }
}
write_csv(data_tables, results)

#stratify over median number of visits by age
#add in dob
setwd(datafiles)
patient <- read_parquet("study_population.parquet")
flu_population[patient, `:=`(dob = i.dob), on = "patid"]

#create age at index variable
calculate_age <- function(dob, reference_date) {
  age <- as.integer(difftime(reference_date, dob, units = "days") / 365.25)
  return(age)
}
reference_date_flu <- as.Date("2019-09-01", format = "%Y-%m-%d")
flu_population[, age := calculate_age(dob, reference_date_flu)]

age_cat <- function(data){
  data[, age_category := cut(
    age,
    breaks = c(65, 70, 75, 80, 85, 90, 95, Inf),
    labels = c("65-69", "70-74", "75-79", "80-84", "85-89", "90-95", "95+"),
    right = FALSE
  )]
  return(data)
}
flu_population <- age_cat(flu_population)

#calculate prevalence 
prev_flu_gp <- flu_population %>%
  group_by(age_category) %>%
  summarise(count = sum(overmedianvisits, na.rm = TRUE), 
            nrows = n(),  
            prev = (count / nrows) * 100) %>%
  ungroup()

#bar plot
plot_flu_gp <- prev_flu_gp %>%
  ggplot(., aes(x = age_category, y = prev)) +  
  geom_bar(stat = "identity", position = "dodge", fill = "#5694ca") +
  labs(x = "Age Category",
       y = "Prevalence") +  
  scale_y_continuous(expand = c(0, 0), breaks = scales::pretty_breaks()) +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 10),
        axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
        axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
        panel.grid = element_blank())

#export
setwd(results)
ggsave("median_gp_visits_by_age_paper_2.pdf", plot_flu_gp, width = 11, height = 7)