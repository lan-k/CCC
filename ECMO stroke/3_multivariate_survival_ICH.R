###joint modelling



# ---- get_data ----
rm(list=ls())
source('2_univariate_survival.R') # prepares the data and runs exclusions
library(dplyr)
library(splines)
library(zoo)
library(runner)
library(survminer)
library(cmprsk)
library(goftest)
library(visdat)
options(scipen=999) # avoid scientific presentation






# ---- survival_data ----

start_date = as.Date(data_date) - 90 #as.Date("2021-08-11")



covars <- c("age", "sex", "comorbidity_obesity","comorbidity_hypertension",
            "comorbidity_diabetes","comorbidity_chronic_neurological_disorder",
            "ecmo_worst_pa_o2_fi_o2_before",
            # "ecmo_prone_before",
            "income_region", "site_name", "days_vent_ecmo5", "era")

#start on the first day of ecmo up to 90 days days_fup

surv_vars <- c("age", "sex", #"ecmo",
               # "pin", "site_name",
               "ethnic_white",
               "current_smoker",
               "comorbidity_obesity",
               "comorbidity_hypertension",
               "comorbidity_diabetes",
               "comorbidity_chronic_cardiac_disease",
               # "comorbidity_chronic_neurological_disorder",
               "ecmo_vasoactive_drugs_before",
               "eotd_anticoagulants",
               "cannula_lumen",
               "days_vent_ecmo5",
               "income_region", "era")



###need to include anticoagulants as a time varying exposure
##may not be able to include ECMo as a time varying covariate?
data.ids0 <- ecmo_daily %>% # jmdat0  %>% mutate(age_day = age + day) %>%
  filter(!is.na(ecmo), !is.na(eotd_anticoagulants)) %>%
  group_by(pin, ecmo, eotd_anticoagulants) %>%
  # group_by(pin, eotd_anticoagulants) %>%
  arrange(pin, day) %>%
  #slice_max(day,n=1) %>%
  ungroup() %>%
  arrange(pin, day) %>%
  group_by(pin) %>%
  mutate(next_day=lead(day),
         prev_day=lag(day),
         # fup=ifelse(is.na(prev_day), day, days_fup-prev_day),
         tstart = ifelse(is.na(prev_day),0,day),
         tstop = ifelse(!is.na(next_day), next_day,days_fup),
         # tstart=ifelse(is.na(prev_day), 0 , prev_day),
         # tstop = ifelse(!is.na(lead(prev_day)),day,days_fup) ,
         status = factor(ifelse(day != max(day), "Alive", as.character(status))),
         any_stroke = ifelse(day != max(day), 0, any_stroke),
         stroke_ICH = ifelse(day != max(day), 0, stroke_ICH),
         status_ICH = case_when(stroke_death %in% c("Missing","No Stroke",
                                                    "Discharged") ~ 0,
                                stroke_death  == "Hemorrhagic stroke" ~ 1,
                                stroke_death  == "Other stroke" ~ 2,
                                stroke_death == "Death" ~ 3),
         status_ICH = ifelse(day != max(day), 0, status_ICH),
         status_ICH = factor(status_ICH, 0:3, 
                             labels=c("Censor", "Hemorrhagic Stroke", 
                                      "Other stroke","Death"))) %>%
  ungroup() %>%
  mutate(tstop = ifelse(tstop == tstart, tstop + 0.5, tstop), #for same day stop and start
         age_start = age+tstart,
         age_stop=age + tstop,
         status2 = factor(ifelse(as.character(status) == "Discharged", "Alive", as.character(status))),
         status3 = as.numeric(status %in% c("Stroke","Death")),
         age50=age-50) %>%
  select(!c(p_h, pa_o2, pa_co2,platelet_count,d_dimer, il_6, aptt,inr))


data.ids0a <- ecmo_daily %>% # ecmo_daily  %>% mutate(age_day = age + day) %>%
  filter( !is.na(eotd_anticoagulants)) %>%
  group_by(pin,  eotd_anticoagulants) %>%
  arrange(pin, day) %>%
  #slice_max(day,n=1) %>%
  ungroup() %>%
  arrange(pin, day) %>%
  group_by(pin) %>%
  mutate(next_day=lead(day),
         prev_day=lag(day),
         # fup=ifelse(is.na(prev_day), day, days_fup-prev_day),
         tstart = ifelse(is.na(prev_day),0,day),
         tstop = ifelse(!is.na(next_day), next_day,days_fup),
         # tstart=ifelse(is.na(prev_day), 0 , prev_day),
         # tstop = ifelse(!is.na(lead(prev_day)),day,days_fup) , 
         status = factor(ifelse(day != max(day), "Alive", as.character(status))),
         any_stroke = ifelse(day != max(day), 0, any_stroke),
         stroke_ICH = ifelse(day != max(day), 0, stroke_ICH),
         status_ICH = case_when(stroke_death %in% c("Missing","No Stroke",
                                                    "Discharged") ~ 0,
                                stroke_death  == "Hemorrhagic stroke" ~ 1,
                                stroke_death  == "Other stroke" ~ 2,
                                stroke_death == "Death" ~ 3),
         status_ICH = ifelse(day != max(day), 0, status_ICH),
         status_ICH = factor(status_ICH, 0:3, 
                             labels=c("Censor", "Hemorrhagic Stroke", 
                                      "Other stroke","Death"))) %>%
  ungroup() %>%
  mutate(tstop = ifelse(tstop == tstart, tstop + 0.5, tstop), #for same day stop and start
         age_start = age+tstart,
         age_stop=age + tstop,
         status2 = case_when(as.character(status) %in% c("Discharged","Alive") ~ 0,
                             as.character(status) == "Stroke" ~ 1,
                             as.character(status) == "Death" ~ 2),
         status2 = factor(status2, 0:2, labels=c("Censor", "Stroke", "Death")),
         status3 = as.numeric(status %in% c("Stroke","Death"))) %>% 
  select(!c(p_h, pa_o2, pa_co2,platelet_count,d_dimer, il_6, aptt,inr, haemoglobin)) 




##the Cox model  - any stroke or death as the outcome

# ---- multi_surv_ICH ----

sform = as.formula(
  paste0("Surv(tstart, tstop, stroke_ICH) ~ cluster(pin, site_name) + age^2 + sex + days_vent_ecmo5^2 + "
         , paste0(surv_vars, collapse="+"))) #+ ecmo

survFit.p0 <- coxph(sform ,
                    data = data.ids0a, id=pin)
# cox.zph(survFit.p0) #all PH
coxsum <- summary(survFit.p0)
cox_data <- data.frame(coxsum$conf.int[,c(1,3:4)])
cox_data$Variable=rownames(cox_data)
rownames(cox_data) <- c()

npat <- data.ids0a %>% select(surv_vars, ecmo, pin) %>% na.omit() %>% pull(pin) %>% unique()
cox_data <- cox_data %>%
  rename(mean=exp.coef.,
         lower=lower..95,
         upper=upper..95)%>%
  mutate(
         HR=paste0(roundz(mean, digits=2), " (",
                   roundz(lower, digits=2),", ",
                   roundz(upper, digits=2),")")
         ) %>%
  mutate(Variable=case_when(#Variable == "ecmo" ~"During vs post ECMO",
                            Variable == "age" ~"Age",
                            Variable == "sexMale" ~"Male vs Female",
                            Variable == "ethnic_whiteYes" ~"White Ethnicity Yes vs No",
                            Variable == "current_smokerYes" ~"Current Smoker Yes vs No",
                            Variable == "income_regionHigh Income" ~"High vs Middle Income Region",
                            Variable == "comorbidity_obesityYes" ~"Obese Yes vs No",
                            Variable == "comorbidity_diabetesYes" ~"Co-morbid Diabetes Yes vs No",
                            Variable == "comorbidity_hypertensionYes" ~"Co-morbid Hypertension Yes vs No",
                            Variable == "comorbidity_chronic_cardiac_diseaseYes" ~
                              "Co-morbid Cardiac Disease Yes vs No",
                            Variable == "eotd_anticoagulants" ~"Anticoagulant use during ECMO Yes vs No",
                            Variable =="ecmo_vasoactive_drugs_beforeYes" ~
                              "Pre-ECMO vasoactive medicine use Yes vs No",
                            Variable == "cannula_lumenSingle lumen" ~"Single vs double lumen cannula",
                            Variable == "days_vent_ecmo5" ~"Days ventilated pre-ECMO",
                            Variable == "era.L" ~ "Pandemic era Jul-Dec 2020 vs Jan-Jun 2020",
                            Variable == "era.Q" ~ "Pandemic era Jan-Sep 2021 vs Jan-Dec 2020",
                            TRUE ~ Variable))

# ---- multi_CR_ICH ----
# the CR Fine-Gray model

etime1 <- data.ids0a$tstart
etime2 <- data.ids0a$tstop
event <- data.ids0a$status_ICH #(0=censored, 1 = stroke, 2= death)

crdata <- data.ids0a %>% select(pin, status_ICH,tstart,tstop,  site_name, all_of(surv_vars))
crdata <- finegray(Surv(etime1, etime2, event) ~ ., data=crdata,id=pin,
                   etype="Hemorrhagic Stroke")  
#creates the Fine-Gray weights for stroke as the event

sform2 = as.formula(
  paste0("Surv(fgstart, fgstop, fgstatus) ~ cluster(pin, site_name)  + age^2 + sex + days_vent_ecmo5^2 + "
         , paste0(surv_vars, collapse="+"))) #+ ecmo

fitFG <- coxph( sform2, id=pin,weight=fgwt,
              data=crdata)

mean <- fitFG$coefficients
se <- sqrt(diag(fitFG$var))
lower = mean - 1.96 * se
upper = mean + 1.96 * se
fg_data <- data.frame(exp(cbind(mean, lower, upper)))

fg_data$Variable=rownames(fg_data)

rownames(fg_data) <- c()

fg_data <- fg_data %>%
  mutate(
    sHR=paste0(roundz(mean, digits=2), " (",
              roundz(lower, digits=2),", ",
              roundz(upper, digits=2),")")
  ) %>%
  mutate(Variable=case_when(#Variable == "ecmo" ~"During vs post ECMO",
                            Variable == "age" ~"Age",
                            Variable == "sexMale" ~"Male vs Female",
                            Variable == "ethnic_whiteYes" ~"White Ethnicity Yes vs No",
                            Variable == "current_smokerYes" ~"Current Smoker Yes vs No",
                            Variable == "income_regionHigh Income" ~"High vs Middle Income Region",
                            Variable == "comorbidity_obesityYes" ~"Obese Yes vs No",
                            Variable == "comorbidity_diabetesYes" ~"Co-morbid Diabetes Yes vs No",
                            Variable == "comorbidity_hypertensionYes" ~"Co-morbid Hypertension Yes vs No",
                            Variable == "comorbidity_chronic_cardiac_diseaseYes" ~
                              "Co-morbid Cardiac Disease Yes vs No",
                            Variable == "eotd_anticoagulants" ~"Anticoagulant use during ECMO Yes vs No",
                            Variable =="ecmo_vasoactive_drugs_beforeYes" ~
                              "Pre-ECMO vasoactive medicine use Yes vs No",
                            Variable == "cannula_lumenSingle lumen" ~"Single vs double lumen cannula",
                            Variable == "days_vent_ecmo5" ~"Days ventilated pre-ECMO",
                            Variable == "era.L" ~ "Pandemic era Jul-Dec 2020 vs Jan-Jun 2020",
                            Variable == "era.Q" ~ "Pandemic era Jan-Sep 2021 vs Jan-Dec 2020",
                            TRUE ~ Variable))


# ---- multi_cox_CR_ICH ----
##combine the cox and Cr forestplots
xticks <- c(0,1,2,4,6,8,10)
xtlab <- rep(TRUE, length.out = length(xticks))
attr(xticks, "labels") <- xtlab

headerc <- tibble(Variable = c( "Variable","Variable"),
                  HR = c("subHR/HR (95% CI)","subHR/HR (95% CI)"),
                  Model = c("Competing Risks","Cox"))



fg_data <- fg_data %>%
  mutate(HR=paste0(" ",sHR))

multi_hr_data=bind_rows(fg_data %>% select(!sHR) %>% mutate(Model = "Competing Risks"),
                       cox_data  %>% mutate(Model = "Cox") #%>% select(!HR)
) %>%
  filter(!is.na(mean))

multi_hr_data=bind_rows(headerc,multi_hr_data)

multi_forest_ICH <- multi_hr_data %>% 
  group_by(Model) %>%
  forestplot(labeltext = c(Variable, HR), #
             graph.pos=2,
             title = "Haemorrhagic stroke within 90 days of ECMO (multivariate models)",
             fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
             xticks=xticks,
             xlog = F, 
             zero=1,
             line.margin = 0.5, #.2,
             boxsize = .25,
             txt_gp = fpTxtGp(label = gpar(cex = 0.8),
                              ticks = gpar(cex = 0.8),
                              xlab  = gpar(cex = 0.8)),
             col = fpColors(box = c("cyan4","darkorange1"),
                            line = c("cyan4","darkorange1")),
             # col = fpColors(box = c("royalblue","darkred"),
             #                line = c("darkblue","darkred")),
             vertices = TRUE,
             xlab="subHR/HR")

multi_forest_ICH




