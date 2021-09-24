# 0_read_data.R
# read the data for: i) the stroke paper, ii) the ecmo papers, iii) obesity paper
# Jan 2021
library(dplyr)
options(dplyr.summarise.inform=FALSE) # suppress annoying warning
library(janitor)
library(stringr)
library(naniar)
library(UpSetR)

# date of this data
rm(list=ls())
data_date <- '2021-09-24'
admin_censor_date = as.Date(data_date, format='%Y-%m-%d') # using date of data update
folder = paste0("R:/data_deidentified/v1.0.0_",data_date)
pfile = paste(folder, '/patients.csv', sep='')
dfile = paste(folder, '/daily.csv', sep='')

### Section 1 - get the latest data ###

## b) read the patient data ##
patients = read.csv(pfile, stringsAsFactors = FALSE) 
#remove columns with all NA
patients <- patients %>% 
  janitor::clean_names() %>%
  unique() %>% # safety net for duplicates
  mutate_if(~is.character(.x), ~ifelse(trimws(.x)=='', NA, .x)) %>%
  janitor::remove_empty(which = "cols") #%>%
  #select(where(~!all(is.na(.x))))
  #select(where(~sum(is.na(.x)) < 0.9*nrow(patients))) #remove columns with > 90% missing

names(patients)



patients <- patients %>%
  select(pin, site_name, country, sex, ethnicity, bmi, age, outcome, diagnosis_coronavirus, diagnosis_coronavirus_type,
         starts_with("comorbidity"),
         contains("stroke"),
         starts_with("med"),
         contains("ecmo"),
         contains("complic"),
         medication_anticoagulation,
         # 'stroke_during_treatment','stroke_during_treatment_side',
         # starts_with("stroke_during_treatment_type_"),
         'transfer_from_other_facility',
         "diagnosis_coronavirus",
         "healthcare_worker", "lab_worker",
         "pregnant",
         'height','weight',
         # ecmo
         'date_outcome',
         "final_outcome",
         # 'date_ecmo','date_ecmo_discontinued','any_ecmo','ecmo_type',
         # 'ecmo_drainage_cannula_insertion_site','ecmo_drainage_cannula_size',
         # 'ecmo_return_cannula_insertion_site',"ecmo_return_cannula_size",
         # "ecmo_cardiac_arrest_before", "ecmo_worst_pa_o2_6hr_before",
         # 'ecmo_highest_fi_o2_6hr_before',
         # "ecmo_worst_pa_co2_6hr_before", 
         # 'ecmo_worst_pa_o2_6hr_before',
         # "ecmo_prone_before",
         # "complication_ards", "complication_viral_pheumonitis",
         # "complication_heart_failure", # to compare missingness
         # "complication_seizure", # to compare missingness
         # "complication_stroke", 
         # "complication_stroke_date", 
         # "neurology_diag_stroke" , # all missing
         # "neurology_diag_stroke_date" ,
         'ecmo_continuous_rpt_before',
         # clinical:
         'apache_ii','sofa','admission_temperature','admission_heart_rate',
         'admission_respiratory_rate','admission_o2_saturation',
         'admis_worst_pa_o2_6hr',
         # other additions (jan 2021)
         'treatment_trachoestomy','cause_of_death',
         'ability_to_self_care','mech_vent_cardiac_arrest', 'complication_cardiac_arrest' ,
         # 'ecmo_in_elso_registry',
         # 'ecmo_neuromuscolar_blockage_before',
         # 'ecmo_vasoactive_drugs_before', 
         'mech_vent_vasoactive_drugs_before', 'treatment_inotropes_vasopressors',
         # dates:
         contains('date')) %>% # slim down variables for now
  mutate(healthlab_worker = pmax(healthcare_worker, lab_worker, na.rm = TRUE) # if either NA then use non-NA information
         # ,date_did = NA, # set up empty dates for below
         # date_dis = NA,
         # date_d = NA,
         # date_a = NA,
         # date_i = NA,
         # date_f = NA,
         # date_o = NA,
         # date_mv = NA,
         # date_mv_end = NA,
         # date_e = NA,
         # date_s = NA
         )

# # had to loop because these dates would not convert
# for(j in 1:nrow(patients)){
#   patients$date_did[j] = ifelse(patients$date_hospital_discharge[j] !='', as.Date(patients$date_hospital_discharge[j]), NA)
#   patients$date_dis[j] = ifelse(patients$date_discharge[j] !='', as.Date(patients$date_discharge[j]), NA)
#   patients$date_d[j] = ifelse(patients$date_death[j] !='', as.Date(patients$date_death[j]), NA)
#   patients$date_a[j] = ifelse(patients$date_admission[j] !='', as.Date(patients$date_admission[j]), NA)
#   patients$date_i[j] = ifelse(patients$date_icu[j] !='', as.Date(patients$date_icu[j]), NA)
#   patients$date_o[j] = ifelse(patients$date_outcome[j] !='', as.Date(patients$date_outcome[j]), NA)
#   patients$date_f[j] = ifelse(patients$date_first_symptom[j] !='', as.Date(patients$date_first_symptom[j]), NA)
#   patients$date_mv[j] = ifelse(patients$date_mechanical_ventilation[j] !='', as.Date(patients$date_mechanical_ventilation[j]), NA)
#   patients$date_mv_end[j] = ifelse(patients$date_mech_vent_discontinued[j] !='', as.Date(patients$date_mech_vent_discontinued[j]), NA)
#   patients$date_e[j] = ifelse(patients$date_ecmo[j] !='', as.Date(patients$date_ecmo[j]), NA)
#   patients$date_ed[j] = ifelse(patients$date_ecmo_discontinued[j] !='', as.Date(patients$date_ecmo_discontinued[j]), NA)
#   patients$date_s[j] = ifelse(patients$complication_stroke_date[j] !='', as.Date(patients$complication_stroke_date[j]), NA)
# }


#convert to date and fix future dates
timelist = as.vector(names(patients %>% select(contains("date"))))
# b <-as.data.frame(lapply(patients[timelist],  FUN=as.Date, format="%Y-%m-%d", origin='1970-01-01'))
# patients[timelist] <- b
# class(patients$date_admission)

fixdate <- function(x, end_date) {
  y = as.Date(x, format="%Y-%m-%d",origin='1970-01-01')
  t = ifelse(y > end_date, NA, y)
  return(as.Date(t, format="%Y-%m-%d",origin='1970-01-01'))
}

b <-as.data.frame(lapply(patients[timelist],  FUN=fixdate, admin_censor_date))      
                                                              

patients[timelist] <- b
class(patients$date_admission)
remove(b)


## combine duplicate comorbidity fields: comorbidity_chronic_hematologic_disease, comorbidity_severe_liver_disease 
patients = patients %>% mutate(
  comorbidity_chronic_hematologic_disease_comb = pmax(comorbidity_chronic_hematologic_disease,
                                                      comorbidity_chronic_hematologic_disease_eot,na.rm=T),
  comorbidity_severe_liver_disease_comb = pmax(comorbidity_severe_liver_disease,
                                               comorbidity_severe_liver_disease_eot,na.rm=T)) %>%
  select(-comorbidity_chronic_hematologic_disease, -comorbidity_severe_liver_disease) %>%
  # rename new variables to old to make life easier:
  rename('comorbidity_chronic_hematologic_disease' = 'comorbidity_chronic_hematologic_disease_comb',
         'comorbidity_severe_liver_disease' = 'comorbidity_severe_liver_disease_comb') 


##find stroke patients from complications_other field that haven't already been counted
patients$complication_stroke[patients$complication_stroke != 1 & 
                               patients$stroke_during_treatment == 1] <- 1


#these are in other complications - put these in subarachnoid stroke?
# strokelist <- c("ICH with midline shift","Intracerebral bleeding",
#                 "intracerebral hemorrhage, subarachnoid hemorrhage,·")
# 
# p<- patients %>% 
#   filter(trimws(complication_other_value) %in% strokelist | 
#            grepl("cerebral", complication_other_value, ignore.case=T),
#          complication_stroke !=1) %>% 
#   select(pin, complication_stroke, stroke_during_treatment, complication_other_value)
# 
# patients$stroke_during_treatment_type_subarachnoid_haemorrhage[
#   trimws(patients$complication_other_value) %in% strokelist] <- 1
# 
# patients <- patients %>%
#   mutate(complication_stroke = case_when (
#     complication_stroke == 1 ~ 1,
#     trimws(complication_other_value) %in% strokelist ~ 1,
#     TRUE ~ complication_stroke
#   ))



patients <- patients %>% mutate(
  stroke_group2 = case_when(
    stroke_during_treatment_type_intraparenchymal_haemorrhage==1 ~ 'Intraparenchymal haemorrhage',
    stroke_during_treatment_type_subarachnoid_haemorrhage==1 ~ 'Subarachnoid haemorrhage',
    stroke_during_treatment_type_ischemic_stroke==1 ~ 'Ischemic stroke',
    stroke_during_treatment_type_hypoxic_ischemic_brain_injury==1 ~ 'Hypoxic ischemic brain injury',
    stroke_during_treatment_type_cerebral_venous_sinus_thrombosis==1 ~ 'Cerebral venous sinus thrombosis',
    stroke_during_treatment_type_other==1 | complication_stroke == 1 ~ 'Undetermined type',
    stroke_during_treatment_type_unknown==1 ~ 'Unknown'
  ),
  stroke_group3 = case_when(
    stroke_group2 %in% c('Intraparenchymal haemorrhage',
                         'Subarachnoid haemorrhage') ~ "ICH",
    stroke_group2 %in% c('Ischemic stroke','Hypoxic ischemic brain injury') ~ "Ischaemic",
    TRUE ~ stroke_group2
  )  #check if there is a subdural haemorrhage in future - would be ICH
)

##15/9/2021 one patient had both SH and hypoxic ischemic brain injury, allocate to ICH
table( patients$stroke_during_treatment_type_hypoxic_ischemic_brain_injury, 
       patients$stroke_during_treatment_type_subarachnoid_haemorrhage, useNA = "always")
table(patients$stroke_group2, patients$stroke_group3, useNA="always")

## make new stroke type variable into list (August 2020) - CRF instructions: "please select up to two (2) options"
yes_responses <- dplyr::select(patients, pin, starts_with('stroke_during_treatment_type_'))  %>%
  mutate(stroke_during_treatment_type_ischemic_stroke = ifelse(stroke_during_treatment_type_ischemic_stroke==1, 'Ischemic stroke', ''), # convert 0,1 numbers to characters
         stroke_during_treatment_type_intraparenchymal_haemorrhage = ifelse(stroke_during_treatment_type_intraparenchymal_haemorrhage==1, 'Intraparenchymal haemorrhage', ''),
         stroke_during_treatment_type_subarachnoid_haemorrhage = ifelse(stroke_during_treatment_type_subarachnoid_haemorrhage==1, 'Subarachnoid haemorrhage', ''),
         stroke_during_treatment_type_hypoxic_ischemic_brain_injury = ifelse(stroke_during_treatment_type_hypoxic_ischemic_brain_injury==1, 'Hypoxic ischemic brain injury', ''),
         stroke_during_treatment_type_cerebral_venous_sinus_thrombosis  = ifelse(stroke_during_treatment_type_cerebral_venous_sinus_thrombosis ==1, 'Cerebral venous sinus thrombosis', ''),
         stroke_during_treatment_type_other  = ifelse(stroke_during_treatment_type_other ==1, 'Other', ''),
         stroke_during_treatment_type_unknown  = ifelse(stroke_during_treatment_type_unknown ==1, 'Unknown', '')
  ) %>%
  as_tibble() %>%
  tidyr::gather(key='type', value='response', -pin) %>%
  filter(!is.na(response),
         !response =='') # remove missing and empty
# now make list from group responses
list_resp = group_by(yes_responses, pin) %>%
  summarise(listed = list(response)) %>%
  ungroup()
names(list_resp)[2] = 'stroke_type'
#patients = select(patients, -starts_with('stroke_during_treatment_type_')) # remove from data
patients = left_join(patients, list_resp, by='pin')

# # more work on dates
# patients = mutate(patients,
#                   date_did = as.Date(date_did, origin='1970-01-01'),
#                   date_dis = as.Date(date_dis, origin='1970-01-01'),
#                   date_a = as.Date(date_a, origin='1970-01-01'),
#                   date_i = as.Date(date_i, origin='1970-01-01'),
#                   date_o = as.Date(date_o, origin='1970-01-01'),
#                   date_d = as.Date(date_d, origin='1970-01-01'),
#                   date_f = as.Date(date_f, origin='1970-01-01'),
#                   date_mv = as.Date(date_mv, origin='1970-01-01'),
#                   date_mv_end = as.Date(date_mv_end, origin='1970-01-01'),
#                   date_e = as.Date(date_e, origin='1970-01-01'),
#                   date_ed = as.Date(date_ed, origin='1970-01-01'),
#                   date_s = as.Date(date_s, origin='1970-01-01')) %>%
#   select(-date_discharge, -date_hospital_discharge, -date_death, -date_outcome, -date_first_symptom, -date_icu, -date_admission, -date_mechanical_ventilation, -date_mech_vent_discontinued, -date_ecmo, -date_ecmo_discontinued, -complication_stroke_date) %>% # use previous date names
#   rename('date_discharge' = 'date_dis', # ICU discharge
#          'date_hospital_discharge' = 'date_did',
#          'date_admission' = 'date_a',
#          'date_icu' = 'date_i',
#          'date_outcome' = 'date_o',
#          'date_death' = 'date_d',
#          'date_first_symptom' = 'date_f',
#          'date_mechanical_ventilation' = 'date_mv',
#          'date_mech_vent_discontinued' = 'date_mv_end',
#          'date_ecmo' = 'date_e',
#          'date_ecmo_discontinued' = 'date_ed',
#          'complication_stroke_date' = 'date_s'
#   )
# quick check
# select(patients, pin, date_admission, date_icu, date_discharge, date_death) %>% sample_n(size=10)


## tidying: i) centre continuous patient variables, ii) missing categories
patients = mutate(patients,
                  bmi5 = (bmi - 27)/5, # standardise BMI
                  age5 = (age - 60)/5, # standardise age
                  # other tidying:
                  sex = ifelse(sex=='', 'Missing', sex),
                  ethnicity  = ifelse(ethnicity =='', 'Missing', ethnicity ),
                  ethnicity  = ifelse(ethnicity =='notavail', 'Missing', ethnicity ))

## c) read the daily data ##
daily = read.csv(dfile, stringsAsFactors = FALSE) %>%
  janitor::clean_names() %>%
  unique() %>% # safety net for duplicates
  mutate_if(~is.character(.x), ~ifelse(trimws(.x)=='', NA, .x)) %>%
  janitor::remove_empty(which = "cols") #%>%
  #select(where(~!all(is.na(.x))))
 # select(where(~sum(is.na(.x)) < 0.9*nrow(daily)))  

names(daily)
#vis_miss(daily, warn_large_data = F)
 
daily <- daily %>%  
  mutate(date_daily = as.Date(date_daily)) %>%
  select(pin, day, date_daily, in_icu, mean_arterial_pressure, compliance_respiratory_system,
         wbc, pa_co2, serum_creatinine, sodium, potassium, crp, fibrinogen, 
         troponin_i, haemoglobin, il_6, d_dimer,
         glucose, neutrophil_count, lymphocyte_count,
         ldh, ferritin,
         'glasgow_coma_score','avpu',
         'prone_positioning',
         contains("respiratory_rate"),
         'ecmo',  #ecmo_type is in patients dataframe
         platelet_count, 'p_h', "aptt", "aptr", "aptt_aptr", "inr", 
         alt_sgpt, "ast_sgot", 'bilirubin','d_dimer',
         'blood_urea_nitrogen', 'hco3', 'pa_co2','pa_o2', 'pa_o2_fi_o2', 'fi_o2',
         lactate,
         'eotd_anticoagulants', "eotd_anticoagulants_type",
         "eotd_haemorrhagic_complication", # 4.52 Haemorrhagic Complication 1
         'eotd_haemorrhagic_complication_source', # 4.53 Source Of Haemorrhagic Complication 1
         'eotd_peep','eotd_tidal_volume','eotd_ventilatory_mode',
         'eotd_fi_o2','eotd_respiratory_rate','eotd_ecmo_circuit_change',
         'eotd_tidal_volume_ideal','eotd_peep','eotd_airway_plateau_pressure',
         'tracheostomy_inserted',
         #
         'neuromuscular_blocking_agents', 'eotd_neuromuscolar_blockage',
         'inotropic_support', 'eotd_vasoactive_drugs',
         contains('infection')) %>% # slim down variables
  filter(day >= 0) %>% # remove clear errors in negative days
  arrange(pin, day)

## last date

# add maximum date for censoring
max_date = group_by(daily, pin) %>%
  summarise(last_date = max(date_daily)) %>% 
  # if no last date then use admin censoring date
  mutate(last_date = ifelse(is.na(last_date) == TRUE, admin_censor_date, 
                            last_date), # fix a few impossible dates
         last_date = as.Date(last_date, origin='1970-01-01')) %>%
  ungroup()
# add date for administrative censoring (people still in hospital)

# add last date back to patients:
patients = left_join(patients, max_date, by='pin') %>% # add last date to patients
  # make a date of transfer if the outcome is transfer, censor patients at this time (later)
  mutate(date_transfer = ifelse(str_detect(string=outcome, pattern='^Transfer'), 
                                date_outcome, NA), # not for Hospitalization
         date_transfer = as.Date(date_transfer, origin='1970-01-01'))




## add dates of prone positioning and ecmo end
# a) dates ecmo 
dates_ecmo = filter(daily, ecmo==1) %>%
  select(pin, date_daily) %>%
  group_by(pin) %>%
  summarise(ecmo_start = min(date_daily),
            ecmo_end = max(date_daily)) %>%
  ungroup() 
# b) prone positioning
dates_prone = filter(daily, prone_positioning==1) %>%
  select(pin, date_daily) %>%
  group_by(pin) %>%
  summarise(prone_start = min(date_daily),
            prone_end = max(date_daily)) %>%
  ungroup() 
# add prone/ecmo dates to patients
patients = left_join(patients, dates_ecmo, by='pin')
patients = left_join(patients, dates_prone, by='pin')

# add stroke groups for those with a stroke, need to use loop because of lists
table(unlist(patients$stroke_type))

patients$stroke_group = ''

for (k in 1:nrow(patients)){
  types = unlist(patients$stroke_type[k])
  if(is.null(types)==FALSE){
    ##combined groups
    i = max(str_detect(pattern="Ischemic stroke|Hypoxic ischemic brain injury",
                       string=types))
    h = max(str_detect(pattern="Intraparenchymal haemorrhage|Subarachnoid haemorrhage", 
                       string=types))
    o = max(str_detect(pattern="Other|Unknown", string=types))
    if(o == 1){
      patients$stroke_group[k] = 'Other/Unknown'}
    if(h == 1){patients$stroke_group[k] = 'Hemorrhage'}
    if(i == 1){patients$stroke_group[k] = 'Ischemic'} # if more than one pick this one
    #    cat(k, i, h, o, patients$stroke_group[k], '\n', sep=' ')
  }
}

table(patients$stroke_group)
table(patients$stroke_group2)
# other fixes
patients = mutate(patients,
                  # change missing to `no stroke`
                  stroke_group = ifelse(stroke_group=='', 'None', stroke_group),
                  # if there is a date then make sure they are not in "None" group
                  stroke_group = ifelse(!is.na(complication_stroke_date) & 
                                          stroke_group=='None', 'Other/Unknown', stroke_group))

### Section 2: save ###

# dates for data
# month = month.name[as.numeric(stringr::str_split(data_date, '_')[[1]][2])]
# day = as.numeric(stringr::str_split(data_date, '_')[[1]][3])
# data_date = list(month=month, day=day)

# save
filename = 'data/COVID_Data_stroke.RData' # stroke
# if(source == 'dummy_data'){filename = 'data/COVID_Data_stroke_dummy.RData'}
save(daily, patients, data_date, file=filename)

# for finding variables:
# see also https://github.com/Samreay/COVID19/blob/master/definitions/oxford.yml
#names(patients)[grep('ecmo', names(patients))] 
#names(daily)[grep('_x', names(daily))] 


# check dates
select(patients, date_admission, date_icu, date_discharge,
       date_hospital_discharge, outcome, last_date) %>%
  filter(!is.na(date_discharge)) %>%
  sample_n(10)

##check dates
max(patients$last_date, na.rm = T)
summary(patients$date_admission)  #max date > today
summary(patients$date_icu)


max(patients$date_discharge, na.rm = T)  > Sys.Date()
max(patients$date_admission, na.rm = T)  > Sys.Date()
max(patients$date_icu, na.rm = T)  > Sys.Date()


date_check <- patients %>% 
  filter(date_discharge > Sys.Date() |
         date_admission > Sys.Date() |
         date_icu > Sys.Date() |
         date_hospital_discharge > Sys.Date() |
         ecmo_start > Sys.Date()) %>%
  select(pin, site_name, country,date_admission,date_icu, ecmo_start,date_discharge, date_hospital_discharge)

date_check_daily <- daily %>%
  filter(date_daily > Sys.Date()) %>%
  select(pin, date_daily)

fn1 = paste0("Data/Checks/patient date checks ",data_date,".csv")
fn2 = paste0("Data/Checks/daily date checks ",data_date,".csv")


check <- patients %>% 
  filter(stroke_during_treatment_type_other=="Yes" 
         | stroke_during_treatment_type_unknown =="Yes") %>% 
  select(pin, stroke_during_treatment_type_other, stroke_during_treatment_type_unknown,
         any_ecmo)

write.csv(check, file="Data/Checks/Missing stroke type.csv", na="", row.names=F)

# write.csv(date_check, file=fn1,na='', row.names = F)
# write.csv(date_check_daily, file=fn2,na='', row.names = F)
