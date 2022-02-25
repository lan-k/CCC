# 1_ecmo_ich_data_prep.R
# used to prepare data for 1_ecmo_ich_part_x.Rmd
# march 2021
#updated Aug 2021

rm(list=ls())
source('99_functions.R')
library(dplyr)
library(naniar)
options(dplyr.summarise.inform = FALSE) # turn off warning
library(stringr)
library(tidyr)
library(broom)
library(hablar)
# library(R2WinBUGS)
# bugs_location = "c:/Program Files/WinBUGS14" # location of the file WinBUGS14.exe
# library(mstate)
# library(cmprsk) # for competing risks
# library(survminer) # for nice plots of competing risks; also test of influential values
library(flextable) # for nice tables
library(captioner)
table_nums <- captioner(prefix = "Table", levels=2, type = c("n", "c")) # number and character
figure_nums <- captioner(prefix = "Figure", levels=2, type = c("n", "c"))
#source("Martin/ext_mstate.R")

# graphics things:
library(diagram)
library(ggplot2)
library(janitor)
theme_set(theme_bw())
library(ggupset)  # to plot combination of symptoms
g.theme = theme_bw() + theme(panel.grid.minor = element_blank())
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999", "#CC79A7")


fup_stroke = 90  ##days of follow-up for stroke study


##may need to exclude spain  - date of stroke not recorded 

## get the data
WB <- read.csv("data/World_Bank_Income_List.csv", stringsAsFactors = F) %>%
  janitor::clean_names() %>%
  select(country, income_group)

spain_pins <- read.csv("data/Spain_pin.csv", stringsAsFactors = F)

load(file='data/COVID_Data_stroke.RData') # from 0_read_data.R  
start_date = as.Date(data_date) - 90 #as.Date("2021-08-11")
# if(source == 'dummy_data'){
#   load(file='data/COVID_Data_stroke_dummy.RData') # from 0_read_data.R, dummy data
#   # fixes with stroke
#   patients = mutate(patients,
#                     stroke_group = ifelse(is.na(complication_stroke_date)==TRUE, 'None', stroke_group), # blank stroke group for those with no date
#                     stroke_group = ifelse(is.na(complication_stroke_date)==FALSE, sample(c('Hemorrhage','Ischemic','Other/Unknown','None'), size=n(), replace=TRUE), stroke_group) # randomly add stroke group for those with a date
#   )
# }
# if(source == 'real_data'){
# load(file='data/COVID_Data_stroke.RData') # from 0_read_data.R
# }


#date for the eras
date1 <- as.Date("2020-01-01", format="%Y-%m-%d")
date2 <- as.Date("2020-06-30", format="%Y-%m-%d")
date3 <- as.Date("2020-12-31", format="%Y-%m-%d")
date4 <- as.Date(data_date) #as.Date("2021-11-30", format="%Y-%m-%d")



patients <- patients %>% left_join(WB, by="country") 

##remove patients from Spanish database since date of stroke is not recorded
patients <- patients %>% anti_join(spain_pins)
daily <- daily %>% anti_join(spain_pins)


  
# make last possible date based on database freeze
last_possible = as.Date(data_date, format='%Y-%m-%d')
#last_possible = as.Date(paste(data_date$day, data_date$month, '2020', sep='-'), format='%d-%b-%Y')
censor_day = 90 # maximum day for plots and estimates

## data edits
# a) patients
# slim down variables:
patients = select(patients, 'pin', 'site_name','country',income_group,
                  'age', 'sex', 'bmi','ethnicity', 'pregnant',
                  'diagnosis_coronavirus','diagnosis_coronavirus_type',
                  contains("date"),
                  #' 'date_icu',
                  #' #'date_ecmo','date_ecmo_discontinued', 
                  #' 'date_death','date_discharge','last_date','date_first_symptom',
                  #' 'date_mech_vent_discontinued', 'date_mechanical_ventilation',
                  "height","weight",'date_admission','date_hospital_discharge',
                  starts_with('complication_'),
                  starts_with('comorbidity_'),
                  'sofa','apache_ii',
                  contains("ecmo"),
                  # 'ecmo_type',
                  # 'ecmo_drainage_cannula_insertion_site',"ecmo_drainage_cannula_size",
                  # 'ecmo_return_cannula_insertion_site',"ecmo_return_cannula_size",
                  # "ecmo_cardiac_arrest_before", "ecmo_worst_pa_o2_6hr_before",
                  # "ecmo_worst_pa_co2_6hr_before", 'ecmo_continuous_rpt_before',
                  contains("stroke"),
                  #'stroke_during_treatment','stroke_type', 'stroke_group', 
                  'outcome','date_outcome','date_transfer',
                  'treatment_trachoestomy','cause_of_death',
                  #'stroke_during_treatment_side','ecmo_in_elso_registry',
                  'ability_to_self_care', medication_anticoagulation,
                  'complication_cardiac_arrest', 'mech_vent_cardiac_arrest',
                  #'ecmo_cardiac_arrest_before','ecmo_neuromuscolar_blockage_before',
                 # ecmo_vasoactive_drugs_before, 
                  mech_vent_vasoactive_drugs_before, treatment_inotropes_vasopressors) %>%
  mutate(pregnant = factor(case_when(pregnant == 1 ~ "Yes",
                              (pregnant == 0) | (sex == "Male") | (!is.na(age) & (age >= 55 | age <= 12))~ "No",
                              pregnant == 2 ~ "Unknown",
                              TRUE ~ "Missing"
          ))  ,# nicer label
         date_ecmo_discontinued = ifelse(is.na(date_ecmo), NA, date_ecmo_discontinued), # if no start date for ecmo, then remove any end date
         date_ecmo_discontinued = as.Date(date_ecmo_discontinued, origin='1970-01-01'),
         ecmo_type = case_when(ecmo_type == "V-A"~"Venous-Arterial",
                               ecmo_type == "V-V" ~ "Venous-Venous",
                               ecmo_type == "Veno-arterial" ~ "Venous-Arterial",
                               TRUE ~ ecmo_type
                               ),
         comorbidity_obesity = case_when(is.na(comorbidity_obesity) & !is.na(bmi) & bmi >= 30 ~ 1,
                                         is.na(comorbidity_obesity) & !is.na(bmi) & bmi < 30 ~ 0,
                                         TRUE ~ comorbidity_obesity),
         cannula_lumen = case_when(
           !is.na(ecmo_return_cannula_insertion_site) & 
             !is.na(ecmo_drainage_cannula_insertion_site) &
             trimws(ecmo_return_cannula_insertion_site) ==
             trimws(ecmo_drainage_cannula_insertion_site) ~ "Double lumen",
           !is.na(ecmo_return_cannula_insertion_site) & 
             !is.na(ecmo_drainage_cannula_insertion_site) &
             trimws(ecmo_return_cannula_insertion_site) !=
             trimws(ecmo_drainage_cannula_insertion_site) ~ "Single lumen",
           !is.na(ecmo_return_cannula_insertion_site) | 
             !is.na(ecmo_drainage_cannula_insertion_site) ~ NA_character_
         ))

tapply(patients$age, patients$pregnant, summary)
## inclusion
table(patients$country, patients$stroke_group, useNA = "always")
# patients = filter(patients,
#                   !country %in% c('Australia','United Kingdom')) # remove countries that did not collect stroke data

## exclusions
# age, must be adults
starting_n = nrow(patients)
(n_exclude_age_stroke = nrow(filter(patients, age < 18, !is.na(complication_stroke_date)))) # count those excluded in stroke group
patients = filter(patients,
                  !is.na(age),
                  age >= 18) # 
n_exclude_age = starting_n - nrow(patients) 
# no ICU date
previous_n = nrow(patients)
(n_exclude_icu_stroke = nrow(filter(patients, is.na(date_icu), !is.na(complication_stroke_date)))) # count those excluded in stroke
patients = filter(patients,
                  !is.na(date_icu))
n_exclude_icu = previous_n -nrow(patients) 
# missing sex - none
previous_n = nrow(patients)
##no patients with missing sex, 3 with 'Not specified'
patients = filter(patients,
                  sex != 'Missing', # remove those with missing sex as they are missing almost all their data
                  !is.na(sex))
n_exclude_missing = nrow(patients) - previous_n  #none, missing sex fixed

# further data edits

## change 0/1 variables to yes/no

# set_2labels <- function(var) {
#   v <- factor(var, levels = c(0,1), labels=c("No","Yes"))
# }


patients  = mutate(patients,
                   # impossible admission date set to missing:
                   date_admission = ifelse(date_admission >= last_possible, NA, date_admission), 
                   date_admission = as.Date(date_admission, origin='1970-01-01'),
                   # update stroke to add EOT
                   complication_stroke = ifelse(!is.na(stroke_during_treatment) & stroke_during_treatment==1, 1, complication_stroke),
                   stroke_before_ecmo = as.integer(if_else(!is.na(complication_stroke_date) &
                                                  !is.na(date_ecmo),
                                                difftime(complication_stroke_date, date_ecmo) < 0, NA)),
                   stroke_after_ecmo = 1 - stroke_before_ecmo,
                   # labels for factors
                   comorbidity_smoking = factor(comorbidity_smoking, levels=0:2, labels=c('Never Smoked','Current Smoker','Former Smoker')),
                   # ,complication_stroke = factor(complication_stroke, levels=0:1, labels=c('No','Yes')),
                   # complication_seizure = factor(complication_seizure, levels=0:1, labels=c('No','Yes')),
                   # complication_heart_failure = factor(complication_heart_failure, levels=0:1, labels=c('No','Yes')),
                   # comorbidity_chronic_kidney_disease = factor(comorbidity_chronic_kidney_disease, levels=0:1, labels=c('No','Yes')),
                   # comorbidity_malignant_neoplasm = factor(comorbidity_malignant_neoplasm, levels=0:1, labels=c('No','Yes')),
                   # comorbidity_smoking = factor(comorbidity_smoking, levels=0:1, labels=c('No','Yes')),
                   # comorbidity_severe_liver_disease = factor(comorbidity_severe_liver_disease, levels=0:1, labels=c('No','Yes')),
                   # comorbidity_diabetes = factor(comorbidity_diabetes, levels=0:1, labels=c('No','Yes')),
                   # comorbidity_hypertension = factor(comorbidity_hypertension, levels=0:1, labels=c('No','Yes')),
                   # comorbidity_chronic_neurological_disorder = factor(comorbidity_chronic_neurological_disorder, levels=0:1, labels=c('No','Yes')),
                   # comorbidity_chronic_cardiac_disease = factor(comorbidity_chronic_cardiac_disease, levels=0:1, labels=c('No','Yes')),
                   # ecmo_cardiac_arrest_before = factor(ecmo_cardiac_arrest_before, levels=0:1, labels=c('No','Yes')),
                   # ecmo_continuous_rpt_before = factor(ecmo_continuous_rpt_before, levels=0:1, labels=c('No','Yes')),
                    any_mv = as.numeric(!is.na(date_mechanical_ventilation)),
                   # any_mv = factor(any_mv, levels=0:1, labels=c('No','Yes')),
                    any_ecmo = as.numeric(!is.na(ecmo_type)),
                    #complication_stroke = ifelse(is.na(complication_stroke), 0,complication_stroke)
                   # any_ecmo = factor(any_ecmo, levels=0:1, labels=c('No','Yes'))
                     )
##remove 2 from yes/no variables  - does this mean unknown?
patients$complication_stroke[patients$complication_stroke == 2] <- NA
patients$complication_bacterial_pneumonia[patients$complication_bacterial_pneumonia  == 2] <- NA
patients$complication_cardiac_arrest[patients$complication_cardiac_arrest  == 2] <- NA
patients$complication_bacteraemia[patients$complication_bacteraemia  == 2] <- NA

yesnovars <- sapply(patients,function(x) { all(na.omit(x) %in% 0:1) }) #binary variables

yesnolist <- as.vector(names(patients %>% select_if(yesnovars)))


patients[yesnolist] <- as.data.frame(
  lapply(patients[yesnolist], 
         FUN=function(x) {factor(x, levels = c(0,1),
                                 labels=c("No","Yes"))}))



patients <- patients %>%
  mutate(
    days_ecmo = as.numeric(date_ecmo_discontinued - date_ecmo + 0.5),
    days_hosp = as.numeric(date_hospital_discharge - date_admission + 0.5),
    days_icu = as.numeric(date_hospital_discharge - date_icu + 0.5),
    days_vent = as.numeric(date_mech_vent_discontinued - date_mechanical_ventilation + 0.5),
    days_vent_ecmo = as.numeric(date_ecmo - date_mechanical_ventilation + 0.5),
    days_sympt_hosp = as.numeric(date_admission - date_first_symptom + 0.5),
    day = as.numeric(date_ecmo - date_icu))



# b) daily data

daily = select(daily, pin, date_daily, ecmo, platelet_count, "p_h", "aptt", "inr", 
               alt_sgpt, ldh, ferritin,
               serum_creatinine, sodium, potassium, crp, fibrinogen, 
               glucose, neutrophil_count,
               eotd_respiratory_rate, respiratory_rate,
               eotd_ecmo_circuit_change,
               #lactate,
               "ast_sgot", 'bilirubin','blood_urea_nitrogen', 
               "lymphocyte_count", 'eotd_anticoagulants', "eotd_anticoagulants_type", 'd_dimer', 'aptt_aptr',
               "pa_o2", "pa_co2", 'fi_o2','wbc', 'crp',  'troponin_i','il_6','haemoglobin',
               "eotd_haemorrhagic_complication", "eotd_haemorrhagic_complication_source",
               'eotd_peep','pa_o2_fi_o2',
               #
               'tracheostomy_inserted','neuromuscular_blocking_agents', 'eotd_neuromuscolar_blockage',
               inotropic_support, eotd_vasoactive_drugs
) %>%
  mutate(date_daily = ifelse(date_daily >= last_possible, NA, date_daily), # impossible date set to missing
         date_daily = as.Date(date_daily, origin='1970-01-01'),
         il_6 = ifelse(il_6 > 10000, NA_real_, il_6),
         d_dimer=ifelse(d_dimer>70000, NA_real_, d_dimer),
         aptt=ifelse(aptt>200, NA_real_, aptt)) #remove implausible values



# change 'missing' to NA  -none for ethnicity and outcome
for_missing = select(patients, -pin,  -stroke_type) %>% # remove ID; any_ecmo is never missing; stroke_type is a list variable
  mutate(
    ethnicity = ifelse(ethnicity=='Missing', NA, ethnicity),
    outcome = ifelse(outcome %in% c('','unknown'), NA, outcome)) %>% # tidy up missing outcome
  mutate_if(~is.character(.x), ~ifelse(.x=='', NA, .x)) # change '' to missing for character variables


## consistent colours 
stroke_i = 'purple'
stroke = 'pink'
stroke_h = 'dark red'
discharge = 'forestgreen'
death = 'black'
icu = 'lightcyan3'
ecmo_start = 'dodgerblue'
ecmo_end = 'goldenrod1'


#### Essential work for fixing stroke group ####
# a) no stroke date
f = filter(patients, complication_stroke=='Yes') %>%
  select(pin, date_icu, complication_stroke_date, last_date, date_death, date_discharge, 
         date_hospital_discharge, date_ecmo_discontinued, stroke_group) 
n_start = nrow(f) # starting number
(n_missing_stroke_date = sum(is.na(f$complication_stroke_date))) 
# number excluded for missing stroke, 116 on 19/11/2021
# b) stroke before ICU admission
f2 = filter(f, 
            !is.na(complication_stroke_date)) # now use just complete
n_update = nrow(f2)
f2 = filter(f2, complication_stroke_date >= date_icu)
(n_stroke_prior_icu = n_update - nrow(f2)) # number excluded, 16 on 19/11/2021
# c) stroke after hospital discharge or still on ECMO
f3 = mutate(f2, 
            diff = as.numeric(complication_stroke_date - date_hospital_discharge)) %>%
  filter(is.na(diff) | diff<=0 )
(n_stroke_post_discharge = nrow(f2) - nrow(f3))
n_final = nrow(f3)
# find patients knocked out by `c`
s = f2$pin[f2$pin %in% f3$pin==FALSE]
filter(f, pin%in%s)

###patients with implausible dates
##find outliers in ecmo dates
outofrange <- function(x) {
  x < -1 | x > 200
}

timevars <- patients %>% select(pin, contains("day"), date_admission,
                                any_ecmo, 
                                date_mechanical_ventilation, date_ecmo, 
                                date_icu,date_ecmo_discontinued, 
                                date_mech_vent_discontinued,
                                date_hospital_discharge,complication_stroke_date,
                                stroke_before_ecmo, stroke_after_ecmo)
  
timechecks <-  timevars %>%
  select(!c(days_hosp, days_sympt_hosp) )%>%
  filter(if_any(contains("day"), outofrange) ) 
  # filter(across(contains("day"), ~ . < 0)  |
  #          across(contains("day"), ~  . > 200) ) 

#write.csv(timechecks %>% filter(any_ecmo == "Yes", !is.na(date_ecmo_discontinued)), file="Data/Checks/Check dates.csv", na="", row.names=F)



# now exclude stroke patients with missing date or pre-ICU stroke
ex1 = filter(patients, 
             complication_stroke=='Yes',
             is.na(complication_stroke_date)) %>% pull(pin)
ex2 = filter(patients, 
             complication_stroke=='Yes',
             !is.na(complication_stroke_date),
             complication_stroke_date < date_icu) %>% pull(pin)
ex3 <- timechecks  %>%
  pull(pin)
ex4 <- patients %>%filter(!is.na(date_ecmo_discontinued)) %>% pull(pin)

ex = union(ex1, ex2); length(ex)
#ex = union(union(union(ex1, ex2), ex3),ex4); length(ex)
patients = filter(patients, !pin %in% ex)
daily = filter(daily, !pin %in% ex)






### consort flow chart (March 2021)
###Aug 2021 still need to exclude strokes before ecmo date, include only ecmo patients 
# jpeg('figures/stroke/flow.jpg', width=5.5, height=3.7, units='in', res=300, quality=100)
# par(mai=c(0,0,0,0))
# labels = c(paste('All patients\n n = ', format(starting_n, big.mark = ','), sep=''),
#            paste('Excluded\n- Age under 18 or missing, n = ', n_exclude_age,
#                  '\n- No ICU date, n = ', n_exclude_icu, 
#                  '\n- Stroke with no date, n = ', n_missing_stroke_date, 
#                  '\n- Stroke with date prior to ICU admission, n = ', n_stroke_prior_icu, 
#                  sep=''),
#            paste('Analysed\n n =', format(nrow(patients), big.mark = ',')))
# n_labels = length(labels)
# M = matrix(nrow=n_labels, ncol=n_labels)
# M[3,1] = "' '" 
# pos = matrix(data=c(0.15,0.8,
#                     0.65,0.5,
#                     0.15,0.2), ncol=2, byrow=TRUE)
# sizes=c(1.3,3.3,1.3) / 10
# props = c(0.36,0.32,0.36) # narrower for first and last
# plotmat(M, name=labels, pos=pos, box.type = 'rect', box.size=sizes, box.prop = props, curve = 0, arr.pos=0.85)
# shape::Arrows(x0=0.15, x1=0.3, y0=0.5, y1=0.5, arr.width=0.2, arr.length=0.22, arr.type='triangle')
# dev.off()





###fix values for blood gases
daily$pa_o2[daily$pa_o2 < 40] <- NA  #remove outliers -SCM 20/8/2021
patients$ecmo_worst_pa_co2_6hr_before[patients$ecmo_worst_pa_co2_6hr_before <= 10] <- NA
patients$ecmo_highest_fi_o2_6hr_before[patients$ecmo_highest_fi_o2_6hr_before <= 20] <- NA



##fix some outliers
##assume unit mismatch
daily <- daily %>%
  mutate(d_dimer = ifelse(d_dimer > 50000, d_dimer/1000, d_dimer))




# combined
combined = full_join(patients, daily, by='pin') #%>% unique()

##add variables from daily to patients
combined_ecmo <- combined %>% 
  filter(any_ecmo=='Yes') 



##worst PF-ratio in the day before ECMO
worst_PF_day_before <- combined_ecmo %>%
  filter(any_ecmo=='Yes') %>%
  mutate(diff=as.numeric(date_daily - date_ecmo)) %>%
  filter(diff == -1) %>%
  group_by(pin) %>%
  slice_min(diff, n=1) %>%
  summarise(ecmo_worst_pa_o2_fi_o2_day_before = mymin(pa_o2_fi_o2)) %>%
  ungroup() %>%
  select(pin, ecmo_worst_pa_o2_fi_o2_day_before) %>%
  unique()


PF_before <- combined_ecmo %>%
  filter(any_ecmo=='Yes') %>%
  mutate(diff=as.numeric(date_daily - date_ecmo),
         ecmo_highest_fi_o2_6hr_before=ifelse(ecmo_highest_fi_o2_6hr_before<21|
                                                ecmo_highest_fi_o2_6hr_before>100,NA,
                                              ecmo_highest_fi_o2_6hr_before/100),
         ecmo_worst_pa_o2_fi_o2_6hr_before =  ecmo_worst_pa_o2_6hr_before/ecmo_highest_fi_o2_6hr_before) %>%
  filter(diff < 0, 
         !is.na(pa_o2_fi_o2) | !is.na(ecmo_worst_pa_o2_fi_o2_6hr_before)) %>%
  group_by(pin) %>%
  slice_max(diff, n=1) %>%
  mutate(ecmo_pa_o2_fi_o2_before = pa_o2_fi_o2) %>%
  select(pin,  ecmo_pa_o2_fi_o2_before,pa_o2_fi_o2, diff, 
         ecmo_worst_pa_o2_6hr_before, ecmo_highest_fi_o2_6hr_before,
         ecmo_worst_pa_o2_fi_o2_6hr_before) %>% #pa_o2, fi_o2,
  unique() %>%
  filter(diff >= -7 |
           !is.na(ecmo_worst_pa_o2_fi_o2_6hr_before)) %>%  #restrict to 10 days before the start of ECMO
  mutate(ecmo_worst_pa_o2_fi_o2_before = ifelse(!is.na(ecmo_worst_pa_o2_fi_o2_6hr_before),
                                                min(ecmo_worst_pa_o2_fi_o2_6hr_before,
                                                    ecmo_pa_o2_fi_o2_before, na.rm=T),
                                                ecmo_pa_o2_fi_o2_before))  %>%
  select(pin,  ecmo_worst_pa_o2_fi_o2_before)  %>%  
  ungroup()


patients <- patients %>%
  left_join(worst_PF_day_before, by="pin") %>%
  left_join(PF_before, by="pin")
  


##restricted to patients with ecmo
# restrict to ecmo patients only
ecmo_patients <- patients %>%
  filter(any_ecmo=='Yes') %>%
  mutate(
    days_ecmo = as.numeric(date_ecmo_discontinued - date_ecmo + 0.5),
    days_hosp = as.numeric(date_hospital_discharge - date_admission + 0.5),
    days_icu = as.numeric(date_hospital_discharge - date_icu + 0.5),
    days_vent = as.numeric(date_mech_vent_discontinued - date_mechanical_ventilation + 0.5),
    days_vent_ecmo = as.numeric(date_ecmo - date_mechanical_ventilation + 0.5),
    days_sympt_hosp = as.numeric(date_admission - date_first_symptom + 0.5),
    day = as.numeric(date_ecmo - date_icu),
    stroke_group4 = case_when(stroke_during_treatment_type_intraparenchymal_haemorrhage=="Yes" |
                                stroke_during_treatment_type_subarachnoid_haemorrhage=="Yes" |
                                stroke_during_treatment_type_cerebral_haemorrhage == "Yes" ~
                                "Haemorrhagic stroke",
                              stroke_during_treatment_type_ischemic_stroke=="Yes" |
                                stroke_during_treatment_type_hypoxic_ischemic_brain_injury=="Yes" |
                                stroke_during_treatment_type_other=="Yes" |
                                complication_stroke == "Yes" ~ 
                                'Cerebral ischemia/Undetermined',
                              complication_stroke == "No" ~ "No Stroke",
                              TRUE ~ "Unknown/Missing"))


##fix negative values
daylist <- as.vector(names(ecmo_patients %>% select(starts_with("day"))))
a <-replace_with_na_all(ecmo_patients[daylist], condition = ~.x <0)
ecmo_patients[daylist] <- a



daily_ecmo <- daily %>% filter(ecmo==1)

stroke_ecmo <- ecmo_patients %>% 
  select(pin, any_ecmo, ecmo_type, date_ecmo, 
         contains("stroke"), 
         complication_other_value) %>%
  filter(complication_stroke == "Yes"| stroke_during_treatment ==1 )  
  



###CHECK MISSING VALUES
ecmovars <- combined_ecmo %>%
  mutate(haemorrhage=ifelse(date_daily >= date_ecmo , eotd_haemorrhagic_complication, NA),
         rr_vent = ifelse(date_daily == date_mechanical_ventilation, 
                              eotd_respiratory_rate, NA),
         hyperbilirunemia = ifelse(date_daily >= date_ecmo & !is.na(bilirubin),  bilirubin > 1.2, NA)
         ) %>%
  group_by(pin) %>%
  summarise(
            ac_before_ecmo = any(!is.na(eotd_anticoagulants) & eotd_anticoagulants ==1 
                                 & date_daily < date_ecmo),
            ac_during_ecmo = any(!is.na(eotd_anticoagulants) & eotd_anticoagulants ==1 
                                 & date_daily >= date_ecmo & date_daily <= date_ecmo_discontinued),
            ac_after_ecmo = any(!is.na(eotd_anticoagulants) & eotd_anticoagulants ==1 
                                & date_daily > date_ecmo_discontinued),
            anticoagulants = case_when(
                          all(is.na(eotd_anticoagulants)) ~ NA_character_,
                          ac_before_ecmo & !ac_during_ecmo & !ac_after_ecmo ~ "AC before ECMO only",
                          ac_before_ecmo & ac_during_ecmo ~ "AC before & during ECMO",
                          !ac_before_ecmo & ac_during_ecmo ~ "AC during but not before ECMO",
                          !ac_before_ecmo & !ac_during_ecmo & ac_after_ecmo ~ "AC after ECMO only",
                          !ac_before_ecmo & !ac_during_ecmo & !ac_after_ecmo ~ "None"),
         rr_vent = mymax(rr_vent),
         ecmo_circuit_change = factor(mymax(eotd_ecmo_circuit_change), 
                                      levels=0:1,labels=c("No","Yes")),
         complication_haemorrhage = factor(mymax(haemorrhage), 
                                           levels=0:1,labels=c("No","Yes")),
         # complication_hyperbilirunemia = factor(as.integer(any(hyperbilirunemia, na.rm=T)), 
         #                                        levels=0:1, labels=c("No","Yes"),
         complication_hyperbilirunemia = factor(mymax(as.integer(hyperbilirunemia)), 
                                                levels=0:1, labels=c("No","Yes"))
         ) %>%
  ungroup()

##NEED to also create lactic acidosis lactate > 4 mmol/L

##worst values pre_ecmo
pre_ecmo <- combined_ecmo %>%
  filter(date_daily < date_ecmo, any_ecmo=='Yes') %>%
  group_by(pin) %>%
  summarise(pre_ecmo_high_LDH=mymax(ldh),
            pre_ecmo_high_CRP=mymax(crp),
            pre_ecmo_high_d_dimer=mymax(d_dimer),
            pre_ecmo_high_ferritin=mymax(ferritin),
            pre_ecmo_high_il_6=mymax(il_6),
            pre_ecmo_high_fibrinogen=mymax(fibrinogen),
            pre_ecmo_high_WBC=mymax(wbc),
            pre_ecmo_high_lymphocyte_count=mymax(lymphocyte_count)) %>%
  ungroup() %>%
  select(pin, contains("pre_ecmo_high")) %>%
  unique()



#worst values at ecmo initiation
ecmo_24 <- combined_ecmo %>%
  filter(date_daily >= date_ecmo, any_ecmo=='Yes',
         !is.na(date_ecmo)) %>% # must be at least one day after ECMO date
  mutate(diff = as.numeric(date_daily - date_ecmo )) %>%
  filter(diff <= 1 & diff >= 0) %>%
  group_by(pin) %>%
  arrange(pin, diff) %>%
  slice_head(n=1) %>%
  mutate(ecmo_high_LDH=mymax(ldh),
         ecmo_high_CRP=mymax(crp),
         ecmo_high_d_dimer=mymax(d_dimer),
         ecmo_high_ferritin=mymax(ferritin),
         ecmo_high_il_6=mymax(il_6),
         ecmo_high_fibrinogen=mymax(fibrinogen),
         ecmo_high_WBC=mymax(wbc),
         ecmo_high_lymphocyte_count=mymax(lymphocyte_count),
         ecmo_high_aptt=mymax(aptt),
         ecmo_high_inr=mymax(inr),
         ecmo_start_worst_platelet_count=mymin(platelet_count),
         ecmo_start_worst_p_h = mymin(p_h),
         ecmo_start_worst_pa_o2 = mymin(pa_o2),
         ecmo_start_worst_pa_co2 = mymax(pa_co2)) %>%
  ungroup() %>%
  select(pin, matches("ecmo_high_|ecmo_start_worst_"))

#highest values at ecmo discontinuation
ecmo_discont <- combined_ecmo %>%
  filter(date_daily >= date_ecmo_discontinued, any_ecmo=='Yes',
         !is.na(date_ecmo_discontinued)) %>% # must be at least one day after ECMO discontinuation date
  mutate(diff = as.numeric(date_daily - date_ecmo_discontinued )) %>%
  filter(diff <= 1) %>%
  group_by(pin) %>%
  arrange(pin, diff) %>%
  slice_head(n=1) %>%
  mutate(ecmo_discont_LDH=mymax(ldh),
         ecmo_discont_CRP=mymax(crp),
         ecmo_discont_d_dimer=mymax(d_dimer),
         ecmo_discont_ferritin=mymax(ferritin),
         ecmo_discont_il_6=mymax(il_6),
         ecmo_discont_fibrinogen=mymax(fibrinogen),
         ecmo_discont_WBC=mymax(wbc),
         ecmo_discont_lymphocyte_count=mymax(lymphocyte_count),
         ecmo_discont_p_h = mymin(p_h),
         ecmo_discont_pa_o2 = mymin(pa_o2),
         ecmo_discont_pa_co2 = mymax(pa_co2)) %>%
  ungroup() %>%
  select(pin, contains("ecmo_discont_"))


ecmo_delta_o2 <- combined_ecmo %>%
  mutate(diff = as.numeric(date_daily - date_ecmo),
         diff2 = diff - 2) %>%  #want last value to be 48 hours after ECMO started
  filter( any_ecmo=='Yes',
          !is.na(pa_o2),
         !is.na(date_ecmo),
         diff <= 4, diff > 0) %>% # must be at least one day after ECMO start date
  group_by(pin) %>%
  arrange( diff2) %>%
  slice_min(order_by = diff2, n=1, with_ties = F) %>%  #take the day closest to 48 hours after start
  ungroup() %>%
  rename(ecmo_48_pa_o2 = pa_o2) %>%
  select(pin, ecmo_48_pa_o2)


ecmo_delta_co2 <- combined_ecmo %>%
  mutate(diff = as.numeric(date_daily - date_ecmo),
         diff2 = diff - 2) %>%
  filter( any_ecmo=='Yes',
          !is.na(pa_co2),
          !is.na(date_ecmo),
          diff <= 4, diff > 0) %>% # must be at least one day after ECMO start date
  group_by(pin) %>%
  arrange(pin, diff2) %>%
  slice_min(order_by = diff2, n=1, with_ties = F) %>%  #take the day closest to 48 hours after start
  ungroup() %>%
  rename(ecmo_48_pa_co2 = pa_co2) %>%
  select(pin, ecmo_48_pa_co2)


  
  
ecmo_patients <- left_join(ecmo_patients, ecmovars %>% 
                             select(pin, anticoagulants,ecmo_circuit_change,
                                    rr_vent,complication_haemorrhage,
                                    complication_hyperbilirunemia) %>%
                             group_by(pin) %>%
                             slice_head(n=1) %>%
                             ungroup(),by="pin") %>%
  filter(any_ecmo=='Yes') %>%
  left_join(pre_ecmo, by="pin") %>%
  left_join(ecmo_24, by="pin") %>%
  left_join(ecmo_discont, by="pin") %>%
  left_join(ecmo_delta_o2, by="pin") %>%
  left_join(ecmo_delta_co2, by="pin") %>%
  unique()
 


ecmo_patients <- ecmo_patients %>%
  mutate(
    ecmo_IJ_cannula_size_DL = ifelse(cannula_lumen == "Double lumen" &
                                       grepl("internal jugular",ecmo_drainage_cannula_insertion_site),
                                     ecmo_drainage_cannula_size, NA),
    ecmo_IJ_cannula_size_SL = case_when(
      cannula_lumen == "Single lumen" &
        grepl("internal jugular",ecmo_drainage_cannula_insertion_site) ~ 
        ecmo_drainage_cannula_size,
      cannula_lumen == "Single lumen" &
        grepl("internal jugular",ecmo_return_cannula_insertion_site) ~ 
        ecmo_return_cannula_size),
    delta_o2 = ecmo_48_pa_o2 - ecmo_worst_pa_o2_6hr_before ,# ecmo_start_worst_pa_o2,
    delta_co2 = ecmo_48_pa_co2 - ecmo_worst_pa_co2_6hr_before,#ecmo_start_worst_pa_co2,
    delta_o2 = ifelse(delta_o2 > 200 | delta_o2 < -150, NA, delta_o2),
    delta_co2 = ifelse(delta_co2 > 200 | delta_co2 < -150, NA, delta_co2),
    rel_delta_o2 = delta_o2/ecmo_worst_pa_o2_6hr_before ,# ecmo_start_worst_pa_o2,
    rel_delta_co2 = delta_co2/ecmo_worst_pa_co2_6hr_before
    )




###outcomes for competing risks

ecmo_patients <- ecmo_patients %>%
  mutate(
         days_ecmo_death = as.numeric(date_death-date_ecmo),
         days_ecmo_stroke = as.numeric(complication_stroke_date-date_ecmo),
         days_ecmo_disch = as.numeric(date_hospital_discharge-date_ecmo),
         max_days_fup=as.numeric(last_date - date_ecmo),
         max_days_fup=ifelse(max_days_fup > fup_stroke, fup_stroke, max_days_fup),
         max_days_fup=ifelse(max_days_fup <=0, NA, max_days_fup),
         ##remove implausible values
         days_vent_ecmo = ifelse(days_vent_ecmo < 0 | days_vent_ecmo > 100, NA, days_vent_ecmo),
         stroke_ICH = as.numeric(stroke_group3 == "ICH" & days_ecmo_stroke <= fup_stroke),
         stroke = as.numeric(complication_stroke == "Yes" & days_ecmo_stroke <= fup_stroke),
         stroke_death = factor(case_when(
           stroke_group3 == "ICH" ~ "Hemorrhagic stroke",
           stroke_group3 %in% c("Ischaemic", "Undetermined type") ~ "Other stroke",
           !is.na(date_death) & is.na(complication_stroke_date) & days_ecmo_death <= fup_stroke
             & days_ecmo_death <= days_ecmo_disch ~ "Death",
           is.na(date_death) & is.na(complication_stroke_date) & 
             days_ecmo_death <= days_ecmo & days_ecmo_death <= fup_stroke ~ "Death",
           is.na(date_death) & is.na(complication_stroke_date) &  
             !is.na(complication_stroke) & !is.na(date_discharge) &
              days_ecmo_disch <= fup_stroke ~ "Discharged",
           is.na(date_death) & is.na(complication_stroke_date) & 
             !is.na(complication_stroke) ~ "No Stroke",
           stroke_group4 == "No Stroke" ~ "No Stroke",
           TRUE ~ "Missing"),
           levels=c("Hemorrhagic stroke","Other stroke","Death","Discharged",
                    "No Stroke","Missing")),
           #stroke_death=relevel(stroke_death, ref="None"),
         days_fup=case_when(stroke_death == "Death"  & days_ecmo_death <= fup_stroke ~ days_ecmo_death,
                            stroke_death %in% c("Hemorrhagic stroke","Other stroke") & 
                              days_ecmo_stroke <= fup_stroke ~ days_ecmo_stroke,
                            stroke_death == "None" & days_ecmo_disch <= fup_stroke 
                            & !is.na(date_hospital_discharge) ~ days_ecmo_disch,
                            TRUE ~ max_days_fup),
         days_fup=ifelse(days_fup <0, NA, days_fup))



ecmo_patients <- ecmo_patients %>%
  mutate(era=case_when(as.numeric(date_admission - date1) >= 0 &
                       as.numeric(date2 - date_admission) >= 0 ~ "Jan-Jun 2020",
                       as.numeric(date_admission - date2) > 0 &
                       as.numeric(date3 - date_admission) >= 0
                       ~ "Jul-Dec 2020",
                       as.numeric(date_admission - date3) > 0 &
                       as.numeric(date4 - date_admission) >= 0
                       ~ "Jan-Sep 2021"
                         ),
         era=ordered(era, levels=c("Jan-Jun 2020","Jul-Dec 2020","Jan-Sep 2021")))


#check OX_00581-0005

# pa_o2_data = blood_gases(before_var='ecmo_worst_pa_o2_6hr_before',
#                          after_var = 'pa_o2',
#                          etype="Venous-Venous")$to_plot
# 
# pa_o2_data <- pa_o2_data %>% 
#   group_by(pin) %>%
#   mutate(t = ifelse(time=="before", 0,1)) %>%
#   arrange(t) %>%
#   mutate(delta_o2 = diff(result)) %>%
#   ungroup()
# 
# delta_o2 <- pa_o2_data %>%
#   filter(t==1) %>%
#   select(pin, delta_o2)



remove(a, ecmo_24, ecmo_delta_co2, ecmo_delta_o2,
       ecmo_discont, ecmovars, PF_before, pre_ecmo,
       spain_pins, timechecks, timevars, WB, worst_PF_day_before)



### flow chart (Feb 2022)
# jpeg('figures/flow.jpg', width=5.5, height=3.7, units='in', res=300, quality=100)
# 
# 
# f = patients %>%
#   filter(date_ecmo < start_date | is.na(date_ecmo)) %>%
#   select(pin, age, date_icu, complication_stroke, complication_stroke_date, last_date, 
#          date_death, date_discharge, date_ecmo, date_hospital_discharge, 
#          date_ecmo_discontinued, any_ecmo,ecmo_type) 
# n_start = nrow(f) # starting number
# 
# ## did not start ECMO
# f2 = f %>% filter(any_ecmo == "Yes")
# nf2=nrow(f2)
# 
# # V-A ECMO
# f3 <- f2 %>%
#   filter(ecmo_type == "Venous-Venous" )
# 
# # stroke before  ECMO start
# so = filter(f3, complication_stroke == "Yes")
# 
# n_update = nrow(so)
# so2 = filter(so, 
#              complication_stroke_date >= date_ecmo )
# 
# 
# final <- ecmo_patients %>% filter(
#   !is.na(date_ecmo),
#   !is.na(days_fup),
#   date_ecmo < start_date,
#   ecmo_type == "Venous-Venous",
#   is.na(stroke_before_ecmo) | stroke_before_ecmo != "Yes")
# 
# n_final = nrow(final)
# 
# #exclusions
# n_no_ecmo = n_start-nf2
# n_VA =  nf2 - nrow(f3)
# n_missing_stroke_date = sum(is.na(so$complication_stroke_date)) 
# n_stroke_prior_ecmo = n_update - nrow(so2) # number excluded
# n_missing_date= nrow(ecmo_patients %>% filter(
#   !is.na(date_ecmo),
#   date_ecmo < start_date,
#   ecmo_type == "Venous-Venous",
#   is.na(stroke_before_ecmo) | stroke_before_ecmo != "Yes")) - n_final
# 
# ##plot the flow chart
# par(mai=c(0,0,0,0))
# labels = c(paste('All patients\n n = ', format(n_start, big.mark = ','), sep=''),
#            paste(
#              # 'Excluded\n- Age under 18, n = ', n_exclude_age,
#              'Excluded\n- No ECMO support, n = ',format(n_no_ecmo, big.mark = ',') ,
#              '\n- Veno-arterial ECMO, n = ', n_VA,
#              '\n- Missing discharge or ECMO start date, n = ', n_missing_date,
#              '\n- Stroke with date prior ECMO support, n = ', n_stroke_prior_ecmo,
#              sep=''),
#            paste('Analysed\n n =', format(n_final, big.mark = ',')))
# n_labels = length(labels)
# M = matrix(nrow=n_labels, ncol=n_labels)
# M[3,1] = "' '"
# pos = matrix(data=c(0.15,0.8,
#                     0.65,0.5,
#                     0.15,0.2), ncol=2, byrow=TRUE)
# sizes=c(1.3,3.3,1.3) / 10
# props = c(0.36,0.32,0.32) # narrower for first and last
# 
# plotmat(M, name=labels, pos=pos, box.type = 'rect', box.size=sizes, box.prop = props, curve = 0, arr.pos=0.85)
# shape::Arrows(x0=0.15, x1=0.3, y0=0.5, y1=0.5, arr.width=0.2, arr.length=0.22, arr.type='triangle')
# dev.off()
# 



