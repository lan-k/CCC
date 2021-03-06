###joint modelling

# ---- get_data ----
rm(list=ls())
source('3_multivariate_survival.R') # prepares the data and runs exclusions
library(JMbayes2)
library(dplyr)

biomarkers <- c(
                # "ecmo_high_d_dimer","ecmo_high_il_6",
                # "ecmo_start_worst_platelet_count",
                # "ecmo_worst_arterial_p_h_6hr_before","ecmo_high_aptt",
                # "ecmo_high_inr",
                "d_dimer","il_6", "platelet_count",
                "pa_o2_fi_o2", "eotd_peep","p_h",
                "pa_co2", "pa_o2",
                # "eotd_pa_co2", "eotd_pa_o2",
                "haemoglobin")

biomarkers2 <- c("platelet_count"
                 ,"pa_co2"
                ,"pa_o2"
                # , "eotd_peep"
                # ,"haemoglobin"
                )


surv_vars2 <- c("age", "sex", "ecmo",
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
               "income_region",
               "era")
##reduced set of covariates
testvars <- c("age","days_vent_ecmo5","sex","ecmo",
              "comorbidity_obesity", "ecmo_vasoactive_drugs_before",
              "eotd_anticoagulants")

jm_models <- function(
                      svars,lvar, lfunc,fForm=NULL,
                      n_iter = 10000L, n_burnin = 1000L,
                      survdf=data.ids0) {
  
  ##function for JM with single longitudinal marker
  
  sdata <- survdf %>% select(site_name, pin, tstart, tstop, any_stroke,
                             all_of(svars))
  
  jdata <- ecmo_daily  %>% 
    select(day, ecmo, all_of(lvar),
           site_name, pin,  all_of(svars)) %>%  
    na.omit()
  
  ids <- jdata$pin %>% unique()
  sdata <- sdata[sdata$pin %in% ids,]
  
  sf <- as.formula(
    paste0("Surv(tstart, tstop, any_stroke) ~ cluster(pin, site_name) + age^2 + days_vent_ecmo5^2 + "
           , paste0(svars, collapse="+")))
  
  sFit <- coxph(sf, data = sdata, id=pin)
  
  lf <- as.formula(
    paste0(lfunc,"~ day + ecmo:day"))
  
  lFit <- lme(lf,data = jdata, random = ~ day | pin)
  
  jFit <- jm(sFit, lFit,time_var = "day",
             functional_forms = fForm)
  
  stab <- summary(jFit)$Survival
  jtab <- round(exp(stab[c(1,3,4)]),digits=3)
  colnames(jtab) <- c("mean","lower","upper")
  jtab$N <- jFit$model_data$n
  jtab$nevent<- sFit$nevent
  
  jtab$Variable = rownames(jtab)
  rownames(jtab) <- c()
  
  # return(jtab)
  return(jFit)
                      
}


# ---- JM_data ----
ecmo_daily <- ecmo_daily %>%
  mutate(logappt60 = log2(aptt) - log2(60),
         abshb=sqrt(abs(haemoglobin - 11)),
         eotd_peep = ifelse(eotd_peep == 0.0, 0.01, eotd_peep))



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
         any_stroke = ifelse(day != max(day), 0, any_stroke)) %>%
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
  select(pin,  day, ecmo,
         any_stroke, stroke_death, status,status2, sofa,
          all_of(surv_vars),
         all_of(biomarkers)) #, abshb, logappt60
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





# ---- plots ----
# plot(lowess(jmdat0$p_h_cent^2 ~ jmdat0$day))
plot(lowess(log(jmdat0$pa_o2) ~ jmdat0$day))
plot(lowess(log(jmdat0$pa_co2) ~ jmdat0$day))
plot(lowess(log(jmdat0$platelet_count) ~ jmdat0$day))
plot(lowess(jmdat0$haemoglobin ~ jmdat0$day))

# ---- univ_jm_model ----


# (jm_ph <- jm_models( svars = testvars,"p_h", "p_h", #surv_vars2
#                       n_iter = 1000L, n_burnin = 100L) ) #wide CI for pH
# exp(summary(jm_ph)$Survival[c(1,3,4)])

# (jm_ph2 <- jm_models(svars = surv_vars2,"p_h_cent", "p_h_cent^2",
#                     n_iter = 1000L, n_burnin = 100L) ) #doesn't converge

# (jm_ph_cent <- jm_models(svars = surv_vars2,"p_h_cent", "abs(p_h_cent)",
#                      n_iter = 1000L, n_burnin = 100L) ) 


(jm_o2 <- jm_models(svars = testvars,"pa_o2", "log(pa_o2)",
                   n_iter = 1000L, n_burnin = 100L))
exp(summary(jm_o2)$Survival[c(1,3,4)])

(jm_o2_sl <- jm_models(svars = testvars,"pa_o2", "log2(pa_o2)",
                    fForm =   ~ value(log2(pa_o2))
                     + slope(log2(pa_o2), eps=1, direction = "back" )
                    ,n_iter = 1000L, n_burnin = 100L))
exp(summary(jFit)$Survival[c(1,3,4)])

(jm_co2 <- jm_models(svars = testvars,"pa_co2", "log2(pa_co2)",
                    n_iter = 1000L, n_burnin = 100L))
exp(summary(jm_co2)$Survival[c(1,3,4)])

# (jm_co2_sl <- jm_models(svars = testvars,"pa_co2", "log2(pa_co2)",
#                      fForm =   ~ value(log2(pa_co2)) + 
#                          slope(log2(pa_co2), eps=1, direction = "back")
#                      ,n_iter = 1000L, n_burnin = 100L))
# exp(summary(jm_co2_sl)$Survival[c(1,3,4)])


(jm_pc <- jm_models(svars = testvars,"platelet_count", "log2(platelet_count)",
                    # fForm = ~ value(platelet_count)  + slope(platelet_count, eps=1, direction = "back"),
                    n_iter = 1000L, n_burnin = 100L))
exp(summary(jm_pc)$Survival[c(1,3,4)])



(jm_hb <- jm_models(svars = testvars,"haemoglobin", "haemoglobin",
                   n_iter = 1000L, n_burnin = 100L))
exp(summary(jm_hb)$Survival[c(1,3,4)])


(jm_hbabs <- jm_models(svars = testvars,"abshb", "abshb",
                    n_iter = 1000L, n_burnin = 100L))
exp(summary(jm_hbabs)$Survival[c(1,3,4)])


# (jm_aptt <- jm_models( svars = surv_vars3,"appt", "log2(appt)",
#                        n_iter = 1000L, n_burnin = 100L,
#                        survdf=data.ids0a)) #most patients are on ecmo, high income, country etc
# 
# (jm_aptt2 <- jm_models( svars = testvars,"logappt60", "logappt60^2",
#                      n_iter = 1000L, n_burnin = 100L,
#                      survdf=data.ids0a)) #most patients are on ecmo, high income, country etc


# ---- multi_jm_model ----

##remove missing data
#remove rows with missing data
##choose either testvars or surv_vars2
jmdat <- jmdat0  %>% 
  select(pin,  all_of(biomarkers2), all_of(testvars ), # surv_vars3,abshb
         any_stroke,  status2,
         ecmo, day)  %>%
  # select(!c(d_dimer, il_6, aptt,logappt60, inr,eotd_peep,pa_o2_fi_o2,
  #           delta_o2, delta_co2, ecmo_worst_pa_o2_fi_o2_before)) %>% #too much missing data  pa_o2_fi_o2,
  na.omit()

ids <- jmdat$pin %>% unique()
data.ids <- data.ids0[data.ids0$pin %in% ids,]
check <- data.ids %>% select(pin, day, prev_day, next_day, tstart, tstop,  days_fup)



length(unique(jmdat$pin))  #393

length(unique(data.ids$pin))  

table(data.ids$any_stroke, useNA = "always")


table(data.ids$status, data.ids$status2) #43 strokes reduced to 35 after removing missing data

check2 <- data.ids %>% 
  group_by(pin) %>%
  mutate(any_post_ecmo = min(ecmo)) %>%
  filter(day==max(day))


table(check2$any_post_ecmo, check2$status2)

# the joint model

# fForms <- list(
#   "log(serBilir)" = ~ slope(log(serBilir)) + slope(log(serBilir)):sex,
#   "prothrombin" = ~ area(prothrombin)
# )


##longitudinal models
# the linear mixed model
# lmeFit.p1 <- lme(p_h ~ day + ecmo:day  ,
#                  data = jmdat, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p2 <- lme(log2(pa_o2) ~ day + ecmo:day  ,
                 data = jmdat, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p3 <- lme(log2(pa_co2) ~ day + ecmo:day  ,
                 data = jmdat, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p4 <- lme(log2(platelet_count) ~ day + ecmo:day  ,
                 data = jmdat, random = ~ day | pin)



# lmeFit.p6 <- lme(log(pa_o2_fi_o2) ~ day + ecmo:day  ,
#                  data = jmdat, random = ~ day | pin)




# lmeFit.p5 <- mixed_model(eotd_anticoagulants ~ day + ecmo:day  ,
#                           data = jmdat,
#                           family = binomial(), random = ~ 1 | pin)

ftable(sex ~ comorbidity_obesity + ecmo_vasoactive_drugs_before, data=ecmo_patients)

##survival model using data from longitudinal model
sform = as.formula(
  paste0("Surv(tstart, tstop, any_stroke) ~ cluster(pin, site_name) + ecmo + age^2 + days_vent_ecmo5^2 + " #
         , paste0(testvars, collapse="+")))  #  surv_vars2 surv_vars3
survFit.p1 <- coxph(sform , 
                    data = data.ids, id=pin)
summary(survFit.p1)

cox.zph(survFit.p1) #all PH 


jointFit1.p1 <- jm(survFit.p1, list( lmeFit.p2, lmeFit.p3, lmeFit.p4),   # lmeFit.p1,,  lmeFit.p5, lmeFit.p6
                   time_var = "day",
                   # functional_forms = fForms,
                   
                   functional_forms = ~ 
                         #value(p_h)  #+ slope(p_h, eps=1, direction = "back")
                        # + value(log2(pa_o2)) #+ slope(log(pa_o2), eps=1, direction = "back") 
                        + value(log2(pa_co2))  # + slope(log(pa_co2), eps=1, direction = "back")
                       # + lag(log(platelet_count)) + slope(log(platelet_count), eps=1, direction = "back")
                       # + vexpit(value(eotd_anticoagulants))
                      
                   ,n_iter = 10000L, n_burnin = 1000L, n_thin = 1L) #,n_iter = 10000L, n_burnin = 1000L)

summary(jointFit1.p1)


# Hazard Ratios
stab <- summary(jointFit1.p1)$Survival
(jmtab <- round(exp(stab[c(1,3,4)]),digits=3))

## check convergence

ggtraceplot(jointFit1.p1, "alphas")
ggdensityplot(jointFit1.p1, "alphas")

##shrink coefficients
jointFitr <- update(jointFit1.p1, priors = list("penalty_alphas" = "ridge")) #horseshoe
coefs <- cbind("un-penalized" = unlist(coef(jointFit1.p1)), 
               "penalized" = unlist(coef(jointFitr)))

exp(coefs[,c(1,2)]) #almost no difference

# save(jointFit1.p1,jointFitr,  file="Data/JM_reduced.Rdata")  #uses testvars as covariates
# save(jointFit1.p1,jointFitr,  file="Data/JM.Rdata")  #uses surv_vars2 as covariates, very wide CIs
# save(jointFit1.p1,jointFitr,  file="Data/JM_nosmoke_ethnic.Rdata")  #uses surv_vars3 as covariates



# ---- jm_model_CR ----
####competing risks version

sform.CR = as.formula(
  paste0("Surv(tstart, tstop, status2) ~ cluster(pin, site_name) + (ecmo +  " #age^2 + days_vent_ecmo5^2 +
         , paste0(testvars, collapse="+"),")*strata(CR)"))  #status2, surv_vars2




##survival only
data.CR0 <- data.ids0 %>%
  select(all_of(testvars), tstop, tstart, status2, pin, site_name) %>% #status2
  na.omit() %>% 
  rename(status=status2)  #surv_vars2


table(data.CR0$status)

npat0 <- length(data.CR0  %>% pull(pin) %>% unique())
data.CR0 <- crLong(data.CR0, 
                  statusVar = "status",
                  censLevel = "Alive", nameStrata = "CR") #c("Alive","Discharged")


survFit.CR0 <- coxph(sform.CR ,id=pin, data = data.CR0)
summary(survFit.CR0)



## joint model
data.CR <- data.ids %>% select(all_of(testvars), tstop, tstart, status2, pin, site_name)  %>%
  rename(status=status2)  #surv_vars2
data.CR <- crLong(data.CR, 
                  statusVar = "status",
                  censLevel = "Alive", nameStrata = "CR") #c("Alive","Discharged")

##version for CR model
# sform2 = as.formula(Surv(tstart, 
#                          tstop, status2) ~ (age + age^2 + sex + ecmo  #Surv(fup, any_stroke) , status2
#                                                + days_vent_ecmo + days_vent_ecmo^2  + cannula_lumen
#                                                + income_region + era  + eotd_anticoagulants #+ ac_before + cannula_lumen
#                                                + ecmo_vasoactive_drugs_before
#                                                + comorbidity_obesity + comorbidity_hypertension
#                                                + comorbidity_diabetes)*strata(CR) #+ ecmo_prone_before 
#                     + cluster(pin, site_name)
# )



survFit.CR <- coxph(sform.CR ,id=pin, data = data.CR)
summary(survFit.CR)


                   

# the CR joint model

CR_forms <- list(
  "log2(pa_o2)" = ~ value(log2(pa_o2)):CR
  ,"log2(pa_co2)" = ~ value(log2(pa_co2)):CR
  ,"log2(platelet_count)" = ~ value(log2(platelet_count)):CR
  # ,"haemoglobin" = ~ value(haemoglobin):CR
)
jointFit1.CR <- jm(survFit.CR, list(lmeFit.p2,lmeFit.p3, lmeFit.p4),   #,  lmeFit.p5
                   time_var = "day",
                   functional_forms = CR_forms
                   # ,n_iter = 25000L, n_burnin = 5000L, n_thin = 5L
                    # ,n_iter = 50000L, n_burnin = 10000L, n_thin = 10L, n_chains=4
                   ,n_iter = 75000L, n_burnin = 15000L, n_thin = 15L, n_chains=4, seed=2022
                   #,n_iter = 100000L, n_burnin = 20000L, n_thin = 20L, n_chains=4
                   ) #,n_iter = 10000L, n_burnin = 1000L)

summary(jointFit1.CR)
(stab.CR <- summary(jointFit1.CR)$Survival)
(jtab.CR <- round(exp(stab.CR[c(1,3,4)]),digits=3))

ggtraceplot(jointFit1.CR, "alphas",grid = TRUE, size=0.5)
ggdensityplot(jointFit1.CR, "alphas",grid = TRUE, size=0.5)
ggdensityplot(jointFit1.CR, "betas",grid = TRUE, size=0.5)
ggdensityplot(jointFit1.CR, "all",grid = TRUE, size=0.5)

# save(jointFit1.CR, file="Data/JM_CR_red.Rdata") #uses testvars and 25k iterations
# save(jointFit1.CR, file="Data/JM_CR_red_50kit.Rdata")
# save(jointFit1.CR, file="Data/JM_CR_red_75kit.Rdata")
# save(jointFit1.CR, file="Data/JM_CR_red_100kit.Rdata")


##different functional forms

CR_forms2 <- list(
  "log2(pa_o2)" = ~ (area(log2(pa_o2)) ):CR  #+slope(log(pa_o2), eps=1, direction = "back")
  ,"log2(pa_co2)" = ~ (area(log2(pa_co2))):CR #+slope(log(pa_co2), eps=1, direction = "back")
  ,"log2(platelet_count)" = ~ value(log2(platelet_count)):CR
  # ,"haemoglobin" = ~ value(haemoglobin):CR
) #area for gases and vlaue for platelets has the best fit


# jointFit1.CR2 <- jm(survFit.CR, list(lmeFit.p2,lmeFit.p3, lmeFit.p4),   #,  lmeFit.p5
#                       time_var = "day",
#                       functional_forms = CR_forms2
#                       ,n_iter = 100000L, n_burnin = 20000L, n_thin = 20L, n_chains=4
#                       )


jointFit1.CR2 <- update(jointFit1.CR, functional_forms = CR_forms2
                        ,n_iter = 100000L, n_burnin = 20000L, n_thin = 20L, n_chains=4)

summary(jointFit1.CR2)

(stab.CR2 <- summary(jointFit1.CR2)$Survival)
(jtab.CR2 <- round(exp(stab.CR2[c(1,3,4)]),digits=3))
# save(jointFit1.CR2, file="Data/JM_CR_red_slope_75kit.Rdata")
# save(jointFit1.CR2, file="Data/JM_CR_red_pc_area_75kit.Rdata")
# save(jointFit1.CR2, file="Data/JM_CR_red_area_75kit.Rdata")
# save(jointFit1.CR2, file="Data/JM_CR_red_gas_area_75kit.Rdata")  # did not converge

save(npat0, survFit.CR0, survFit.CR, jointFit1.CR2, 
     file="Data/JM_CR_red_gas_area_100kit.Rdata") #This has the best fit of unpenalisd models


#marginal DIC               WAIC                                        LPML
#431872456287391318006   27906778285344291444806626440826264288448066   -1331734283880501547608462

##try different functional form of days_vent_ecmo
sform2.CR = as.formula(
  paste0("Surv(tstart, tstop, status) ~ cluster(pin, site_name) + (ecmo + I(days_vent_ecmo5^2) + "
         , paste0(testvars, collapse="+"),")*strata(CR)"))  #surv_vars2

survFit2.CR <- coxph(sform2.CR ,id=pin, data = data.CR)
summary(survFit2.CR)

        


library(lmtest)
lrtest(survFit2.CR,survFit.CR)  #marginal difference

jointFit2.CR <- jm(survFit2.CR, list(lmeFit.p2,lmeFit.p3, lmeFit.p4),   #,  lmeFit.p5
                   time_var = "day",
                   functional_forms = CR_forms
                   ,n_iter = 75000L, n_burnin = 15000L, n_thin = 15L, n_chains=4
) #,n_iter = 10000L, n_burnin = 1000L)
summary(jointFit2.CR)
(stab2.CR <- summary(jointFit2.CR)$Survival)
(jtab2.CR <- round(exp(stab2.CR[c(1,3,4)]),digits=3))




##shrink CR coefficients
## value models
# load( file="Data/JM_CR_red_75kit.Rdata")
load(file="Data/JM_CR_red_100kit.Rdata")
jointFit1.CRh <- update(jointFit1.CR, priors = list("penalty_alphas" = "horseshoe")) 
jointFit1.CRr <- update(jointFit1.CR, priors = list("penalty_alphas" = "ridge")) 

save(jointFit1.CRr, jointFit1.CRh, 
     file="Data/JM_CR_red_100kit_pen.Rdata")

## gas area models
load( file="Data/JM_CR_red_gas_area_100kit.Rdata")
jointFit.CRr <- update(jointFit1.CR2, priors = list("penalty_alphas" = "ridge")) #horseshoe


coefs.CRr <- cbind("un-penalized" = unlist(coef(jointFit1.CR2)),
               "penalized" = unlist(coef(jointFit.CRr)))

exp(coefs.CRr[,c(1,2)]) #almost no difference
# exp(unlist(confint(jointFit.CRr)))

jointFit.CRh <- update(jointFit1.CR2, priors = list("penalty_alphas" = "horseshoe")) #best fit
summary(jointFit.CRh)

##horseshoe gase area model fit stats
# marginal DIC        WAIC             LPML
# 150709.9           220791806352    -212327080

# save(jointFit1.CR,jointFit.CR2, jointFit.CRh, file="Data/JM_CR_red_75kit.Rdata")
save( jointFit.CRh, jointFit.CRr,
     file="Data/JM_CR_red_gas_area_100kit_pen.Rdata")


# ---- jm_model_CR-4 ----
##includes peep

jmdat2 <- jmdat0  %>% 
  select(pin,  all_of(biomarkers2),eotd_peep, all_of(testvars), 
         any_stroke,  status,
         ecmo, day)  %>%
  # select(!c(d_dimer, il_6, aptt,logappt60, inr,eotd_peep,pa_o2_fi_o2,
  #           delta_o2, delta_co2, ecmo_worst_pa_o2_fi_o2_before)) %>% #too much missing data  pa_o2_fi_o2,
  na.omit()

ids <- jmdat2$pin %>% unique()
data.ids2 <- data.ids0[data.ids0$pin %in% ids,]


length(unique(jmdat2$pin))  #370

length(unique(data.ids2$pin))  

table(data.ids2$any_stroke, useNA = "always")

table(data.ids2$status, data.ids2$status2) #43 strokes reduced to 34 after removing missing data



lmeFit.p2.2 <- lme(log2(pa_o2) ~ day + ecmo:day  ,
                   data = jmdat2, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p3.2 <- lme(log2(pa_co2) ~ day + ecmo:day  ,
                   data = jmdat2, random = ~ day | pin)  #, na.action=na.exclude

lmeFit.p4.2 <- lme(log2(platelet_count) ~ day + ecmo:day  ,
                   data = jmdat2, random = ~ day | pin)

lmeFit.p5.2 <- lme(log2(eotd_peep) ~ day + ecmo:day  ,  #peep
                 data = jmdat2, random = ~ day | pin)


sform.CR2 = as.formula(
  paste0("Surv(tstart, tstop, status) ~ cluster(pin, site_name) + (ecmo +  " #age^2 + days_vent_ecmo5^2 +
         , paste0(testvars, collapse="+"),")*strata(CR)"))  #surv_vars2



data.CR2 <- data.ids2 %>% select(all_of(testvars), tstop, tstart, status2, pin, site_name) %>%
  rename(status=status2)  #surv_vars2
data.CR2 <- crLong(data.CR2, 
                  statusVar = "status",
                  censLevel = "Alive", nameStrata = "CR")

survFit.CR2 <- coxph(sform.CR2 ,id=pin, data = data.CR2)
summary(survFit.CR2)




# the CR joint model

CR_forms.2 <- list(
  "log2(pa_o2)" = ~ (area(log2(pa_o2)) ):CR  #+slope(log(pa_o2), eps=1, direction = "back")
  ,"log2(pa_co2)" = ~ (area(log2(pa_co2))):CR #+slope(log(pa_co2), eps=1, direction = "back")
  ,"log2(platelet_count)" = ~ value(log2(platelet_count)):CR
  ,"log2(eotd_peep)" = ~ (area(log2(eotd_peep))):CR
)
jointFit.CR2 <- jm(survFit.CR2, list(lmeFit.p2.2,lmeFit.p3.2, lmeFit.p4.2,  lmeFit.p5.2),   #
                   time_var = "day",
                   functional_forms = CR_forms.2
                   # ,n_iter = 25000L, n_burnin = 5000L, n_thin = 5L
                   # ,n_iter = 50000L, n_burnin = 10000L, n_thin = 10L, n_chains=4
                   ,n_iter = 75000L, n_burnin = 15000L, n_thin = 15L, n_chains=4, seed=2022
                   #,n_iter = 100000L, n_burnin = 20000L, n_thin = 20L, n_chains=4
) #,n_iter = 10000L, n_burnin = 1000L)

summary(jointFit.CR2)
(stab.CR2 <- summary(jointFit.CR2)$Survival)
(jtab.CR2 <- round(exp(stab.CR2[c(1,3,4)]),digits=3))

ggtraceplot(jointFit.CR2, "alphas",grid = TRUE, size=0.5)
ggdensityplot(jointFit.CR2, "alphas",grid = TRUE, size=0.5)
ggdensityplot(jointFit.CR2, "betas",grid = TRUE, size=0.5)
ggdensityplot(jointFit.CR2, "all",grid = TRUE, size=0.5)


##shrinkage
jointFit.CRh2 <- update(jointFit.CR2, priors = list("penalty_alphas" = "horseshoe")) 
summary(jointFit.CRh2)

save(jointFit.CR2, jointFit.CRh2, 
     file="Data/JM_CR_red_gas_peep_100kit_pen.Rdata")
