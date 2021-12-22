###joint modelling

# ---- get_data ----
rm(list=ls())
source('3_multivariate_survival_ICH.R') # prepares the data and runs exclusions
library(JMbayes2)

biomarkers <- c(
                # "ecmo_high_d_dimer","ecmo_high_il_6",
                # "ecmo_start_worst_platelet_count",
                # "ecmo_worst_arterial_p_h_6hr_before","ecmo_high_aptt",
                # "ecmo_high_inr",
                "d_dimer","il_6", "platelet_count",
                "pa_o2_fi_o2", "eotd_peep",
                "pa_co2",
                "pa_o2", "p_h",
                "haemoglobin")

biomarkers2 <- c("platelet_count"
                 ,"pa_co2"
                ,"pa_o2"
                # , "p_h"
                ,"haemoglobin"
                )


surv_vars2 <- c("age", "sex", #"ecmo",
                # "pin", "site_name",
                "ethnic_white",
                "current_smoker",
                "comorbidity_obesity",
                "comorbidity_hypertension",
                "comorbidity_diabetes",
                # "comorbidity_chronic_cardiac_disease",
                # "comorbidity_chronic_neurological_disorder",
                "ecmo_vasoactive_drugs_before",
                "eotd_anticoagulants"
                ,"days_vent_ecmo5"
                # ,"cannula_lumen",  #CI for cannula very wide, so leave out
                # "sofa",
                ,"income_region"
                , "era"
                )



surv_vars3 <- c("age", "sex", 
               # "pin", "site_name",
               "comorbidity_obesity",
               "comorbidity_hypertension",
               "comorbidity_diabetes",
               # "comorbidity_chronic_cardiac_disease",
               # "comorbidity_chronic_neurological_disorder",
               "ecmo_vasoactive_drugs_before",
               "eotd_anticoagulants",
               "days_vent_ecmo5",
               # "cannula_lumen",
               # "income_region", 
               "era")
##reduced set of covariates
testvars <- c("age","days_vent_ecmo5","sex", #"ecmo",
              "comorbidity_obesity", "ecmo_vasoactive_drugs_before",
              "eotd_anticoagulants")




# ---- JM_data ----


jmdat0 <- ecmo_daily %>%
  mutate(age_day = age + day) %>%  
  group_by(pin) %>%
  mutate(
         # p_h = rollapply(p_h,3,median,align='center',fill=NA),
         # pa_co2 = rollapply(pa_co2,3,median,align='center',fill=NA),
         # pa_o2 = rollapply(pa_o2,3,median,align='center',fill=NA),
         # haemoglobin = rollapply(haemoglobin,3,median,align='center',fill=NA),
         p_h = ifelse(day==0 & is.na(p_h), ecmo_worst_arterial_p_h_6hr_before, p_h),
         pa_o2 = ifelse(day==0 & is.na(pa_o2), ecmo_worst_pa_o2_6hr_before, pa_o2),
         pa_co2 = ifelse(day==0 & is.na(pa_co2), ecmo_worst_pa_co2_6hr_before, pa_co2)
         ) %>%
  # filter(!all(is.na(p_h))
  #        ,!all(is.na(ecmo))
  #        ,!all(is.na(pa_co2))
#        ,!all(is.na(pa_o2))
#        ) %>%
ungroup() %>%
  select(pin,  day, ecmo,sofa,
         any_stroke, stroke_death, status,stroke_ICH,status_ICH,
          all_of(surv_vars),
         all_of(biomarkers))
  # select(pin, date_daily, date_ecmo, cannula_lumen,
  #        ecmo_vasoactive_drugs_before, eotd_anticoagulants, #,ac_before,  survival not converging with ac_before
  #        p_h, p_h_cent,pa_o2, pa_co2, pa_o2_fi_o2, eotd_peep, haemoglobin,abshb,
  #        # ecmo_worst_pa_o2_6hr_before, ecmo_worst_pa_co2_6hr_before,
  #        # ecmo_worst_arterial_p_h_6hr_before,
  #        # ecmo_highest_fi_o2_6hr_before,
  #        delta_o2, delta_co2,ecmo_worst_pa_o2_fi_o2_before,
  #        platelet_count,d_dimer, il_6, aptt,logappt60, inr,
  #        day, ecmo,days_fup, any_stroke, stroke_death,
  #        age_day,all_of(covars),
  #        status)




# jmdat0 %>% vis_dat()





##remove missing data
#remove rows with missing data
##choose either testvars or surv_vars2
jmdat <- jmdat0  %>% 
  # filter(status_ICH != "Other stroke") %>%
  select(pin,  all_of(biomarkers2), all_of(testvars), #surv_vars2
         any_stroke,  status, #ecmo,
          day)  %>% 
  # select(!c(d_dimer, il_6, aptt,logappt60, inr,eotd_peep,pa_o2_fi_o2,
  #           delta_o2, delta_co2, ecmo_worst_pa_o2_fi_o2_before)) %>% #too much missing data  pa_o2_fi_o2,
  na.omit()

ids <- jmdat$pin %>% unique()
data.ids <- data.ids0a[data.ids0a$pin %in% ids,]
check <- data.ids %>% select(pin, day, prev_day, next_day, tstart, tstop,  days_fup)

length(unique(jmdat$pin))  #389

length(unique(data.ids$pin))  

table(data.ids$any_stroke, useNA = "always")

table(data.ids$status, data.ids$status_ICH) #29 ICH and 6 other stroke after removing missing data


# the joint model

# fForms <- list(
#   "log(serBilir)" = ~ slope(log(serBilir)) + slope(log(serBilir)):sex,
#   "prothrombin" = ~ area(prothrombin)
# )


##longitudinal models
# the linear mixed model
# lmeFit.p1 <- lme(p_h ~ day + ecmo:day  ,
#                  data = jmdat, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p2 <- lme(log2(pa_o2) ~ day #+ ecmo:day  
                 ,data = jmdat, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p3 <- lme(log2(pa_co2) ~ day #+ ecmo:day  
                 ,data = jmdat, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p4 <- lme(log2(platelet_count) ~ day #+ ecmo:day  
                 ,data = jmdat, random = ~ day | pin)

# lmeFit.p5 <- lme(abshb ~ day + ecmo:day  ,  #haemoglobin
#                  data = jmdat, random = ~ day | pin)

# lmeFit.p6 <- lme(log(pa_o2_fi_o2) ~ day + ecmo:day  ,
#                  data = jmdat, random = ~ day | pin)




# lmeFit.p5 <- mixed_model(eotd_anticoagulants ~ day + ecmo:day  ,
#                           data = jmdat,
#                           family = binomial(), random = ~ 1 | pin)


# ---- jm_model_CR ----
####competing risks version

n_other = nrow(data.ids %>% filter(status_ICH == "Other stroke"))
n_ICH = nrow(data.ids %>% filter(status_ICH == "Hemorrhagic Stroke"))
n_death = nrow(data.ids %>% filter(status_ICH == "Death"))

data.CR_ICH <- data.ids %>% select(all_of(testvars), tstop, tstart, status_ICH, pin, site_name) %>%
  rename(status=status_ICH)  #surv_vars2
data.CR_ICH <- crLong(data.CR_ICH, 
                  statusVar = "status",
                  censLevel = "Censor", nameStrata = "CR")

##version for CR model


sform.CR_ICH = as.formula(
  paste0("Surv(tstart, tstop, status2) ~ cluster(pin, site_name) + ( age^2 + days_vent_ecmo5^2 + "
         , paste0(testvars, collapse="+"),")*strata(CR)"))  #surv_vars2 ecmo +

survFit.CR_ICH <- coxph(sform.CR_ICH ,id=pin, data = data.CR_ICH)
summary(survFit.CR_ICH)


                   

# the CR joint model

CR_forms <- list(
  "log2(pa_o2)" = ~ value(log2(pa_o2)):CR
  ,"log2(pa_co2)" = ~ value(log2(pa_co2)):CR
  ,"log2(platelet_count)" = ~ value(log2(platelet_count)):CR
  # ,"haemoglobin" = ~ value(haemoglobin):CR
)
jointFit1.CR_ICH <- jm(survFit.CR_ICH, list(lmeFit.p2,lmeFit.p3, lmeFit.p4),   #,  lmeFit.p5
                   time_var = "day",
                   functional_forms = CR_forms
                   ,n_iter = 25000L, n_burnin = 5000L, n_thin = 5L) #,n_iter = 10000L, n_burnin = 1000L)

summary(jointFit1.CR_ICH)
(stab.CR_ICH <- summary(jointFit1.CR_ICH)$Survival)
(jtab.CR_ICH <- round(exp(stab.CR_ICH[c(1,3,4)]),digits=3))

save(jointFit1.CR_ICH, file="Data/JM_CR_red_ICH.Rdata") #uses testvars

##shrink CR coefficients
# jointFit.CR2 <- update(jointFit1.CR, priors = list("penalty_alphas" = "ridge")) #horseshoe
# 
# 
# coefs.CR <- cbind("un-penalized" = unlist(coef(jointFit1.CR)), 
#                "penalized" = unlist(coef(jointFit.CR2)))
# 
# exp(coefs.CR[,c(1,2)]) #almost no difference
# exp(unlist(confint(jointFit1.CR)))


