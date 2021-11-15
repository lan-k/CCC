###joint modelling
rm(list=ls())

source('1_ecmo_ich_data_prep.R') # prepares the data and runs exclusions
library(dplyr)
library(splines)
library(JMbayes2)
library(zoo)
library(runner)
library(survminer)
library(cmprsk)
library(goftest)
library(visdat)
options(scipen=999) # avoid scientific presentation

(start_date = as.Date(data_date) - 90) #as.Date("2021-08-11")

biomarkers <- c("ecmo_high_d_dimer","ecmo_high_il_6",
                "ecmo_start_worst_platelet_count",
                "ecmo_worst_arterial_p_h_6hr_before","ecmo_high_aptt",
                "ecmo_high_inr",
                "d_dimer","il_6", "platelet_count",
                "pa_o2_fi_o2",
                "haemoglobin",
                "eotd_anticoagulants")

covars <- c("age", "sex", "comorbidity_obesity","comorbidity_hypertension",
            "comorbidity_diabetes",
            "ecmo_worst_pa_o2_fi_o2_day_before",
             # "ecmo_prone_before",
            "income_region", "site_name", "days_vent_ecmo", "era")

#start on the first day of ecmo up to 90 days days_fup

ecmo_patients <- ecmo_patients %>%
  filter(
    #!is.na(date_ecmo_discontinued),
    !is.na(date_ecmo),
    !is.na(days_fup),
    date_ecmo < start_date,
    ecmo_type == "Venous-Venous",
    is.na(stroke_before_ecmo) | stroke_before_ecmo != "Yes") %>%
  mutate(stroke_death=as.character(stroke_death),
         status = case_when(stroke_death %in% c("Missing","No Stroke") ~ "Alive", 
                                  stroke_death %in% c("Hemorrhagic stroke","Other stroke") ~ "Stroke",
                                TRUE ~ stroke_death),
         #status = factor(ifelse(day == days_fup, stroke_death, "Alive")),
         status = relevel(factor(status), ref="Alive"),
         ac_before = ifelse(!is.na(anticoagulants), 
                            as.numeric(anticoagulants == "AC before ECMO only"), NA),
         any_stroke = as.numeric(stroke_group3 
                                 %in% c("ICH","Ischaemic","Undetermined type")),
         income_region = relevel(factor(ifelse(income_group == "High income",
                                               "High Income","Middle Income")),
                                 ref="Middle Income"))


table(ecmo_patients$stroke_group2, ecmo_patients$stroke_group3, useNA="always")
table(ecmo_patients$stroke_death, ecmo_patients$stroke_group3, useNA="always")

ecmo_daily <- ecmo_patients %>% 
  select(pin,days_fup,income_region, stroke_death,status, ac_before, era,
         any_stroke,delta_o2, delta_co2, days_vent_ecmo, ecmo_worst_pa_o2_fi_o2_day_before)  %>% 
  left_join(combined_ecmo %>% select(!days_vent_ecmo) %>%
              filter(any_ecmo=='Yes', !is.na(date_daily)), by="pin") %>%
  mutate(last_date_bio = date_ecmo+days_fup-1,
         diff = as.numeric(last_date_bio - date_daily),
         day = as.numeric(date_daily - date_ecmo),
         ecmo=case_when(date_daily == date_ecmo~ 1,
                        is.na(ecmo) ~ as.numeric(date_daily >= date_ecmo & 
                                                 date_daily <= date_ecmo_discontinued),
                        TRUE ~ ecmo),
         days_vent_ecmo = ifelse(days_vent_ecmo <0 | days_vent_ecmo > 100,NA, days_vent_ecmo)) %>%
  filter(day >= 0, day <= days_fup)  
  #select(pin, date_daily, date_ecmo, days_fup) %>%



##test joint model
hist(ecmo_daily$p_h)  #don't need to log-transform
summary(ecmo_daily$p_h)
summary(ecmo_daily$days_vent_ecmo)
##remove missing data

jmdat0 <- ecmo_daily %>%
  mutate(age_day = age + day) %>%
  group_by(pin) %>%
  # mutate(p_h = rollapply(p_h,3,median,align='center',fill=NA),
  #        pa_co2 = rollapply(pa_co2,3,median,align='center',fill=NA),
  #        pa_o2 = rollapply(pa_o2,3,median,align='center',fill=NA),
  #        pa_o2_fi_o2 = rollapply(pa_o2_fi_o2,3,median,align='center',fill=NA),
  #        p_h = ifelse(day==0 & is.na(p_h), ecmo_worst_arterial_p_h_6hr_before, p_h),
  #        pa_o2 = ifelse(day==0 & is.na(pa_o2), ecmo_worst_pa_o2_6hr_before, pa_o2),
  #        pa_co2 = ifelse(day==0 & is.na(pa_co2), ecmo_worst_pa_co2_6hr_before, pa_co2)
  #        ) %>%
  # filter(!all(is.na(p_h))
  #        ,!all(is.na(ecmo))
  #        ,!all(is.na(pa_co2))
  #        ,!all(is.na(pa_o2))
  #        ) %>%
  ungroup() %>%
  select(pin, date_daily, date_ecmo, cannula_lumen,
         ecmo_vasoactive_drugs_before, eotd_anticoagulants, #,ac_before,  survival not converging with ac_before
         p_h, pa_o2, pa_co2, pa_o2_fi_o2,haemoglobin,
         # ecmo_worst_pa_o2_6hr_before, ecmo_worst_pa_co2_6hr_before,
         # ecmo_worst_arterial_p_h_6hr_before,
         # ecmo_highest_fi_o2_6hr_before,
         delta_o2, delta_co2,ecmo_worst_pa_o2_fi_o2_day_before,
         platelet_count,d_dimer, il_6, aptt,inr,
         day, ecmo,days_fup, any_stroke, stroke_death,
          age_day,all_of(covars),
          status) 




jmdat0 %>% vis_dat()

###need to include anticoagulants as a time varying exposure
data.ids0 <- jmdat0 %>% # ecmo_daily  %>% mutate(age_day = age + day) %>%
  filter(!is.na(ecmo), !is.na(eotd_anticoagulants)) %>%
  group_by(pin, ecmo, eotd_anticoagulants) %>%
  arrange(pin, day) %>%
  #slice_max(day,n=1) %>%
  ungroup() %>%
  arrange(pin, day) %>%
  group_by(pin) %>%
  mutate(prev_day=lag(day),
         fup=ifelse(is.na(prev_day), day, days_fup-prev_day),
         tstart=ifelse(is.na(prev_day), 0 , prev_day),
         tstop = ifelse(!is.na(lead(prev_day)),day,days_fup) , 
         status = factor(ifelse(day != max(day), "Alive", as.character(status))),
         any_stroke = ifelse(day != max(day), 0, any_stroke)) %>%
  ungroup() %>%
  mutate(tstop = ifelse(tstop == tstart, tstop + 0.5, tstop), #for same day stop and start
         age_start = age+tstart,
         age_stop=age + tstop,
         status2 = factor(ifelse(as.character(status) == "Discharged", "Alive", as.character(status))),
         status3 = as.numeric(status %in% c("Stroke","Death"))) %>% 
  select(!c(p_h, pa_o2, pa_co2,platelet_count,d_dimer, il_6, aptt,inr))  



#remove rows with missing data
jmdat <- jmdat0  %>% 
  select(!c(d_dimer, il_6, aptt,inr,
            delta_o2, delta_co2, ecmo_worst_pa_o2_fi_o2_day_before)) %>% #too much missing data  pa_o2_fi_o2,
  na.omit()

ids <- jmdat$pin %>% unique()
data.ids <- data.ids0[data.ids0$pin %in% ids,]


length(unique(jmdat$pin))  #338
length(unique(data.ids$pin))  

table(data.ids$any_stroke, useNA = "always")


# the linear mixed model
table(data.ids$status, data.ids$status2) #41 strokes reduced to 30 after removing missing data


lmeFit.p1 <- lme(p_h ~ day + ecmo:day  ,
                 data = jmdat, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p2 <- lme(log(pa_o2) ~ day + ecmo:day  ,
                 data = jmdat, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p3 <- lme(log(pa_co2) ~ day + ecmo:day  ,
                 data = jmdat, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p4 <- lme(log(platelet_count) ~ day + ecmo:day  ,
                 data = jmdat, random = ~ day | pin)

lmeFit.p5 <- lme(log(pa_o2_fi_o2) ~ day + ecmo:day  ,
                 data = jmdat, random = ~ day | pin)


lmeFit.p6 <- lme(haemoglobin ~ day + ecmo:day  ,
                 data = jmdat, random = ~ day | pin)

# lmeFit.p5 <- mixed_model(eotd_anticoagulants ~ day + ecmo:day  ,
#                           data = jmdat,
#                           family = binomial(), random = ~ 1 | pin)


##the Cox model  - any stroke or death as the outcome

##test linearity of continuous variables
data.ids <- data.ids %>% mutate(age_std = (age - 50)/12, age2=age+age^2)
tfit <- survFit.p1 <- coxph(Surv(tstart, tstop, any_stroke) ~ age2 #log(age) 
                            +  days_vent_ecmo + days_vent_ecmo^2
                            #+ ecmo_worst_pa_o2_fi_o2_day_before
                            #+ delta_o2 + delta_co2
                            +  cluster(pin), 
                            data = data.ids, id=pin)

ggcoxfunctional(tfit, data=data.ids)  #linearity assumption


##formula Cox model
sform = as.formula(Surv(tstart, tstop, any_stroke) ~ age + age^2 + sex + ecmo  #Surv(fup, any_stroke) , status2
                      + days_vent_ecmo + days_vent_ecmo^2  + cannula_lumen
                      + income_region + era  + eotd_anticoagulants #+ ac_before 
                      + ecmo_vasoactive_drugs_before
                      + comorbidity_obesity + comorbidity_hypertension
                      + comorbidity_diabetes #+ ecmo_prone_before 
                      + cluster(pin, site_name)
                   )



survFit.p1 <- coxph(sform , 
                    data = data.ids, id=pin)
summary(survFit.p1)


cox.zph(survFit.p1) #all PH 



# the joint model

# fForms <- list(
#   "log(serBilir)" = ~ slope(log(serBilir)) + slope(log(serBilir)):sex,
#   "prothrombin" = ~ area(prothrombin)
# )
jointFit1.p1 <- jm(survFit.p1, list(lmeFit.p2,lmeFit.p3, lmeFit.p4,  lmeFit.p5),   # lmeFit.p1, lmeFit.p6,
                   time_var = "day",
                   # functional_forms = fForms,
                   
                   functional_forms = ~ 
                         #value(p_h)  #+ slope(p_h, eps=1, direction = "back")
                        + value(log(pa_o2)) #+ slope(log(pa_o2), eps=1, direction = "back") 
                        + value(log(pa_co2))  # + slope(log(pa_co2), eps=1, direction = "back")
                       # + lag(log(platelet_count)) + slope(log(platelet_count), eps=1, direction = "back")
                       # + vexpit(value(eotd_anticoagulants))
                      
                   ,n_iter = 1000L, n_burnin = 100L) #,n_iter = 10000L, n_burnin = 1000L)

summary(jointFit1.p1)


# Hazard Ratios
(stab <- summary(jointFit1.p1)$Survival)
round(exp(stab[c(1,3,4)]),digits=3)


##shrink coefficients
jointFit2 <- update(jointFit1.p1, priors = list("penalty_alphas" = "ridge")) #horseshoe
coefs <- cbind("un-penalized" = unlist(coef(jointFit1.p1)), 
      "penalized" = unlist(coef(jointFit2)))

exp(coefs[,c(1,2)]) #almost no difference



####competing risks version

##cumulative incidence curves
crdata <- data.ids %>%
  group_by(pin) %>%
  arrange(pin, day) %>%
  slice_tail(n=1) %>% 
  filter(status2 != "Alive")


fit <- cuminc(ftime = crdata$days_fup, fstatus = crdata$status2)

ggcompetingrisks(fit, palette = "Dark2",
                 legend = "top",
                 ggtheme = theme_bw(),
                 xlab = "Days since ECMO",
                 # multiple_panels = F,
                 conf.int = T)
# the CR model

data.CR <- crLong(data.ids, statusVar = "status2",
                  censLevel = "Alive", nameStrata = "CR")

##version for CR model
sform2 = as.formula(Surv(tstart, 
                         tstop, any_stroke) ~ (age + age^2 + sex + ecmo  #Surv(fup, any_stroke) , status2
                        + days_vent_ecmo + days_vent_ecmo^2  + cannula_lumen
                         + income_region + era  + eotd_anticoagulants #+ ac_before + cannula_lumen
                         + ecmo_vasoactive_drugs_before
                         + comorbidity_obesity + comorbidity_hypertension
                      + comorbidity_diabetes)*strata(CR) #+ ecmo_prone_before 
                    + cluster(pin, site_name)
)

survFit.p2 <- coxph(sform2 ,id=pin, data = data.CR)
summary(survFit.p2)
                   

# the CR joint model
jointFit1.p2 <- jm(survFit.p2, list(lmeFit.p2,lmeFit.p3, lmeFit.p4,  lmeFit.p5), 
                   time_var = "day",
                   functional_forms = ~ lag(y) + slope(y, eps=1, direction = "back")
                   ,n_iter = 1000L, n_burnin = 100L) #,n_iter = 10000L, n_burnin = 1000L)

summary(jointFit1.p2)
(stab2 <- summary(jointFit1.p2)$Survival)
round(exp(stab2[c(1,3,4)]),digits=3)

##shrink CR coefficients
jointFit.CR2 <- update(jointFit1.p2, priors = list("penalty_alphas" = "ridge")) #horseshoe


coefs.CR <- cbind("un-penalized" = unlist(coef(jointFit1.p2)), 
               "penalized" = unlist(coef(jointFit.CR2)))

exp(coefs.CR[,c(1,2)]) #almost no difference
exp(unlist(confint(jointFit1.p2)))

