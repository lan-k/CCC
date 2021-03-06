###univariate Cox and CR models

# ---- get_data ----
rm(list=ls())
source('1_ecmo_ich_data_prep.R') # prepares the data and runs exclusions
library(dplyr)
library(splines)
library(JMbayes2)
library(zoo)
library(runner)
library(survminer)
library(survival)
library(cmprsk)
library(goftest)
library(visdat)
library(ggsci)
library(cowplot)
library(forestplot)
options(scipen=999) # avoid scientific presentation


start_date = as.Date(data_date) - 90 #as.Date("2021-08-11")
#start on the first day of ecmo up to 90 days days_fup

covars <- c("age", "sex", "comorbidity_obesity","comorbidity_hypertension",
            "comorbidity_diabetes", "comorbidity_chronic_neurological_disorder",
            "ecmo_worst_pa_o2_fi_o2_before",
            # "ecmo_prone_before",
            "income_region", "site_name", "days_vent_ecmo", "era")

# surv_vars <- c("days_fup", "any_stroke" ,"age", "sex", 
#                # "pin", "ecmo","site_name",
#                "comorbidity_obesity",
#                "comorbidity_hypertension",
#                "comorbidity_diabetes",
#                "comorbidity_chronic_cardiac_disease",
#                "comorbidity_chronic_neurological_disorder",
#                "ecmo_vasoactive_drugs_before",
#                # "anticoagulants",
#                "days_vent_ecmo5",
#                "cannula_lumen",
#                "income_region",   "era")
# 




univ <- function(var,var2=var, outcome= "any_stroke", df) {
  
  form = as.formula(paste0("Surv(tstart, tstop,", outcome," ) ~ cluster(site_name) + ",var2))
  fit <- coxph( form, id=pin,
                data=df)
  
  est <- fit$coefficients
  se <- sqrt(diag(fit$var))
  lower = est - 1.96 * se
  upper = est + 1.96 * se
  N = length(unique(df %>% filter(!is.na(across(all_of(var)))) %>% pull(pin)))
  hr <- cbind(N,exp(cbind(est, lower, upper)))
  
  return(round(hr, digits=3))
  
}


univ_compete <- function(var,var2=var,stat=status2, stype="Stroke", df) { #stat=status2
  
  ##descriptives
  
  # if (!is.numeric(ecmo_patients %>% pull(var ))) {
  #   obs <- ecmo_patients %>%
  #     filter(!is.na(var)) %>%
  #     group_by(across(all_of(var))) %>%
  #     mutate(N=n()) %>%
  #     group_by(across(all_of(c(var,"any_stroke","N")))) %>%
  #     summarise(n=n()) %>%
  #     ungroup() %>%
  #     filter(any_stroke == 1) %>%
  #     mutate(obs=paste0(n,"/",N)) %>%
  #     select(var, obs)
  #   
  # }
  # 
  
  ##CR model
  
  etime1 <- df$tstart
  etime2 <- df$tstop
  event <- df %>% pull({{stat}})
  # event <- df$status2 
  #(0=censored, 1 = stroke, 2= death) OR
  #(0=censored, 1 = ICH stroke, 2= other stroke, 3= death) OR
  
  df <- df %>% select(pin, {{stat}},tstart,tstop,  site_name, all_of(var))
  crdata <- finegray(Surv(etime1, etime2, event) ~ ., data=df,id=pin,
                     etype=stype)  
  #creates the Fine-Gray weights for stroke as the event
  
  
  form = as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ cluster(site_name) + ",var2))
  fit <- coxph( form, id=pin,weight=fgwt,
                data=crdata)
  
  est <- fit$coefficients
  se <- sqrt(diag(fit$var))
  lower = est - 1.96 * se
  upper = est + 1.96 * se
  N = length(unique(df %>% filter(!is.na(across(all_of(var)))) %>% pull(pin)))
  shr <- cbind(N,exp(cbind(est, lower, upper)))
  
  # maxrow = 2*nrow(shr)/3  ##remove sHRs for Censoring
  # shr <- shr[1:maxrow,]
  shr <- round(shr, digits=3) #subdistribution HRs for Stroke and Death
  
  return(shr)
  
}



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
         status2=ifelse(status == "Discharged", "Alive", status),
         status = relevel(factor(status), ref="Alive"),
         comp_stroke_death = as.numeric(status2 %in% c("Stroke", "Death")), 
         ac_before = ifelse(!is.na(anticoagulants), 
                            as.numeric(anticoagulants == "AC before ECMO only"), NA),
         any_stroke = as.numeric(stroke_group3 
                                 %in% c("ICH","Ischaemic","Undetermined type")),
         income_region = relevel(factor(ifelse(income_group == "High income",
                                               "High Income","Middle Income")),
                                 ref="Middle Income"), 
         ethnic_white = ifelse(ethnicity == "white", "Yes", "No"),
         group="ECMO patients",
         sex=relevel(factor(sex), ref="Female"),
         age50=age-50,
         # delta_o2 = ifelse(delta_o2 > 100 | delta_o2 < -100, NA, delta_o2),
         # delta_co2 = ifelse(delta_co2 > 100 | delta_co2 < -100, NA, delta_co2),
         del_o2_plus50 = ifelse(!is.na(rel_delta_o2),
                                delta_o2 < -25 | delta_o2 > 25, NA),
         del_co2_neg50 = ifelse(!is.na(delta_co2),
                                delta_co2 < -25 | delta_co2 > 25, NA),
         rel_delta_o2t = ifelse(rel_delta_o2 < -1 | rel_delta_o2 > 1, NA, rel_delta_o2),
         rel_delta_co2t = ifelse(rel_delta_co2 < -1 | rel_delta_co2 > 1, NA, rel_delta_co2),
         rel_o2_plus = ifelse(!is.na(rel_delta_o2t),rel_delta_o2t >= 0.5, NA),
         rel_co2_neg = ifelse(!is.na(rel_delta_co2t), rel_delta_co2t < -0.5, NA),
         pc_under200 = ecmo_start_worst_platelet_count < 200,
         # ecmo_start_worst_platelet_count = ecmo_start_worst_platelet_count/100,
         current_smoker = relevel(factor(case_when(comorbidity_smoking== "Current Smoker" ~ "Yes",
                                                   !is.na(comorbidity_smoking) ~ "No")),ref="No"),
         ever_smoker = relevel(factor(case_when(comorbidity_smoking %in% c("Current Smoker",
                                                                           "Former Smoker") ~ "Yes",
                                                !is.na(comorbidity_smoking) ~ "No")),ref="No"),
         pregnant = ifelse(as.character(pregnant) %in% c("Missing","Unknown"), NA_character_,
                           as.character(pregnant)),
         ac_before_during = case_when(anticoagulants %in% c("None","AC after ECMO only")  ~ "No",
                                      !is.na(anticoagulants)   ~ "Yes"),
         era1_0 = if_else(era =="Jan-Jun 2020",0 , 1),
         era2_1 = if_else(era !="Jan-Sep 2021", 0, 1),
         days_vent_ecmo = ifelse(days_vent_ecmo <0 | days_vent_ecmo > 100,NA, days_vent_ecmo),
         days_vent_ecmo5 = days_vent_ecmo - 5,
         tstart = 0, tstop=days_fup,
         tstop = ifelse(tstop == tstart, tstop + 0.5, tstop),
         status2 = case_when(as.character(status) %in% c("Discharged","Alive") ~ 0,
                             as.character(status) == "Stroke" ~ 1,
                             as.character(status) == "Death" ~ 2),
         status2 = factor(status2, 0:2, labels=c("Censor", "Stroke", "Death")),
         status_ICH = case_when(stroke_death %in% c("Missing","No Stroke") ~ 0,
                                stroke_death  == "Discharged" ~ 1,
                                stroke_death  == "Hemorrhagic stroke" ~ 2,
                                stroke_death  == "Other stroke" ~ 3,
                                stroke_death == "Death" ~ 4),
         status_ICH = factor(status_ICH, 0:3, 
                             labels=c("Censor", "Hemorrhagic Stroke", 
                                      "Other stroke","Death")))


ecmo_daily <- ecmo_patients %>% 
  select(pin,days_fup,income_region, stroke_death,
         status, status2,status_ICH,stroke_ICH,
         ac_before, era,
         any_stroke,comp_stroke_death,
         delta_o2, delta_co2,
         rel_delta_o2, rel_delta_co2,
         del_o2_plus50,del_co2_neg50,
         rel_o2_plus,rel_co2_neg,
         days_vent_ecmo, days_vent_ecmo5, 
         ecmo_worst_pa_o2_fi_o2_before, ecmo_start_worst_platelet_count, age50,
         current_smoker, ethnic_white, pregnant,era1_0, era2_1 )  %>% 
  left_join(combined_ecmo %>% select(!c(days_vent_ecmo, pregnant)) %>%
              filter(any_ecmo=='Yes', !is.na(date_daily)), by="pin") %>%
  mutate(last_date_bio = date_ecmo+days_fup-1,
         diff = as.numeric(last_date_bio - date_daily),
         day = as.numeric(date_daily - date_ecmo),
         pa_o2 = ifelse(pa_o2 > 900,pa_o2/10, pa_o2 ), #fix implausible values
         pa_co2 = ifelse(pa_co2 <10,NA, pa_co2 ),
         p_h_cent = p_h-7.4,
         ecmo=case_when(date_daily == date_ecmo~ 1,
                        is.na(ecmo) ~ as.numeric(date_daily >= date_ecmo & 
                                                   date_daily <= date_ecmo_discontinued),
                        TRUE ~ ecmo),
         current_smoker = case_when(comorbidity_smoking == "Current Smoker" ~ "Yes",
                                    !is.na(comorbidity_smoking) ~ "No")) %>%
  filter(day >= 0, day <= days_fup)  


data.ids0e <- ecmo_daily %>% # ecmo_daily  %>% mutate(age_day = age + day) %>%
  filter(!is.na(ecmo)) %>%
  group_by(pin, ecmo) %>%
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





data.ids0 <- ecmo_daily %>% # ecmo_daily  %>% mutate(age_day = age + day) %>%
  filter(!is.na(ecmo), !is.na(eotd_anticoagulants)) %>%
  group_by(pin, ecmo, eotd_anticoagulants) %>%
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




# ---- ci_curve0 ----

##cumulative incidence curves

fitci <- cuminc(ftime = ecmo_patients$days_fup, fstatus = ecmo_patients$status2,
                cencode = "Alive")

ci <- ggcompetingrisks(fitci, palette = "Dark2",
                       legend = "top",
                       ggtheme = theme_bw(),
                       xlab = "Days since ECMO",
                       # multiple_panels = F,
                       conf.int = T)


# ---- ci_curve0 ----
##mstate version

fit0 <- survfit(Surv(days_fup, status, type="mstate") ~ 1, data=ecmo_patients)
# ggcompetingrisks(fit1,
#                  ggtheme = theme_cowplot())  + 
#   scale_fill_jco(labels=c("Alive", "Death","Discharged","Stroke"))

# ---- ci_curve1 ----
fit1 <- survfit(Surv(days_fup, status, type="mstate") ~ group, 
                data=ecmo_patients)

g1 <- ggcompetingrisks(fit1,
                       ggtheme = theme_cowplot(),
                       group=group,
                       xlab = "Days since ECMO initiation",
                       ylab="Cumulative Incidence",
                       title="")  + 
  scale_fill_jco(labels=c("Alive", "Death","Discharged","Stroke"), name="")

# tiff(file="paper/Images/Fig1.tif", res=300, compression="lzw",
     # units="mm", width=170, height=120)
g1
# dev.off()

# ---- ci_curve2 ----
fit2 <- survfit(Surv(days_fup, status2, type="mstate") ~ group, 
                data=ecmo_patients)

g2 <- ggcompetingrisks(fit2,
                       ggtheme = theme_cowplot(),
                       group=group,
                       xlab = "Days since ECMO initiation",
                       title="Cumulative Incidence")  + 
  scale_fill_jco(labels=c("Alive", "Death","Stroke"))


g2


# ---- test_linearity ----
##test linearity of continuous variables
tdat <- data.ids0 %>% 
  select(pin, tstart, tstop, any_stroke,
         # ,ecmo_worst_pa_o2_fi_o2_before,sofa,
         delta_o2, delta_co2,
         age, age50
         # , days_vent_ecmo, days_vent_ecmo5
         ) %>%
  # mutate(days_vent_ecmo6 = days_vent_ecmo - 6,
  #        days_vent_ecmo4 = days_vent_ecmo - 4) %>%
  na.omit()
tfit <-  coxph(Surv(tstart, tstop, any_stroke) ~ 
               #   I(age+age^2 )
               # + I(age50 + age50^2)
               # log(age) + log(age)^2
               # + I(days_vent_ecmo + sqrt(days_vent_ecmo))
               # +  log2(days_vent_ecmo)
               # +  I(days_vent_ecmo5 + days_vent_ecmo5^2)
               # + days_vent_ecmo + days_vent_ecmo^2
               # + ecmo_worst_pa_o2_fi_o2_before  
               # + log2(ecmo_worst_pa_o2_fi_o2_before)
               + pspline(delta_o2) + pspline(delta_co2)
               # + delta_o2 + delta_co2
               # + I(sofa + sofa^2)
               +  cluster(pin, site_name), 
               data = tdat, id=pin,
               na.action = na.omit)
# cox.zph(tfit)
# ggcoxfunctional(tfit, data=tdat)  #linearity assumption
# summary(tfit)


# ---- uni_surv_stroke ----
xticks <- c(0,1,2,4,6,8,10)
xtlab <- rep(TRUE, length.out = length(xticks))
attr(xticks, "labels") <- xtlab

hr_ecmo <- univ("ecmo", df= data.ids0e)
hr_age <- univ(var="age",#var2="age + age^2",
               df= ecmo_patients)
hr_sex <- univ("sex",df= ecmo_patients)
hr_vent_ecmo <- univ(var="days_vent_ecmo5",
                      # var2="days_vent_ecmo5 + I(days_vent_ecmo5^2)",
                     df= ecmo_patients)
hr_ethnic <- univ("ethnic_white",df= ecmo_patients)
hr_smoke <- univ("current_smoker",df= ecmo_patients)
hr_cannula_lumen <- univ("cannula_lumen",df= ecmo_patients)
# (hr_ac_before <- univ("ac_before",df= data.ids0)) #all strokes had no AC before
hr_anticoagulants <- univ("eotd_anticoagulants",df= data.ids0)
hr_obesity <- univ( var= "comorbidity_obesity",df= ecmo_patients)
hr_vasoactive <- univ("ecmo_vasoactive_drugs_before",df= ecmo_patients)
hr_diabetes <- univ( var= "comorbidity_diabetes",df= ecmo_patients)
hr_cardiac <- univ("comorbidity_chronic_cardiac_disease",df= ecmo_patients)
# (hr_neuro <- univ("comorbidity_chronic_neurological_disorder",df= ecmo_patients))
# (hr_pregnant <- univ("pregnant",df= data.ids0))
hr_hypertension <- univ("comorbidity_hypertension",df= ecmo_patients)
hr_pf <- univ(var="ecmo_worst_pa_o2_fi_o2_before",
              var2="log2(1/ecmo_worst_pa_o2_fi_o2_before)",df= ecmo_patients)
hr_sofa <- univ(var="sofa",#var2="sofa + I(sofa^2)",
                df= ecmo_patients)
hr_del_o2_plus50 <- univ(var="del_o2_plus50",df= ecmo_patients)
hr_del_co2_neg50 <- univ(var="del_co2_neg50",df= ecmo_patients)
hr_rel_o2_plus <- univ(var="rel_o2_plus",df= ecmo_patients)
hr_rel_co2_neg <- univ(var="rel_co2_neg",df= ecmo_patients)


hr_pc_under200 <- univ(var="pc_under200",df= ecmo_patients)
hr_pc <- univ(var="ecmo_start_worst_platelet_count",
              var2="log2(1/ecmo_start_worst_platelet_count)", df= ecmo_patients)
hr_region <- univ("income_region",df= ecmo_patients)
hr_era <- univ("era",df= ecmo_patients)
# (hr_era1_0 <- univ("era1_0",df= ecmo_patients))
# (hr_era2_1 <- univ("era2_1",df= ecmo_patients))

# fisher.test(ecmo_patients$pregnant, ecmo_patients$any_stroke)

### forestplot

hr_data <- data.frame(rbind(hr_ecmo,hr_age, hr_sex,
                            hr_vent_ecmo,
                            hr_ethnic,
                            hr_smoke, hr_obesity, 
                            hr_diabetes, hr_cardiac,
                            hr_hypertension, #hr_neuro,
                            
                            hr_vasoactive,
                            hr_anticoagulants,
                            hr_cannula_lumen, 
                            hr_pf,
                            hr_sofa,
                            # hr_del_o2_plus50,
                            # hr_del_co2_neg50,
                            hr_rel_o2_plus,
                            hr_rel_co2_neg,
                            hr_pc,
                            # hr_pc_under200,
                            hr_region, hr_era
)) %>%
  rename(mean=est) %>% #, lower=V2, upper = V3
  mutate(HR=paste0(roundz(mean, digits=2), " (",
                   roundz(lower, digits=2),", ",
                   roundz(upper, digits=2),")"),
         N=as.character(N))
hr_data$Variable = rownames(hr_data)
rownames(hr_data) <- c()




##nice labels

hr_data <- hr_data %>%
  mutate(Variable=case_when(Variable == "ecmo" ~"During vs post ECMO",
                            Variable == "age" ~"Age",
                            Variable == "sexMale" ~"Male vs Female",
                            Variable == "ethnic_whiteYes" ~"White Ethnicity Yes vs No",
                            Variable == "current_smokerYes" ~"Current Smoker Yes vs No",
                            Variable == "income_regionHigh Income" ~"High vs Middle Income Region",
                            Variable == "comorbidity_obesityYes" ~"*Obese Yes vs No",
                            Variable == "comorbidity_diabetesYes" ~"Co-morbid Diabetes Yes vs No",
                            Variable == "comorbidity_hypertensionYes" ~"Co-morbid Hypertension Yes vs No",
                            Variable == "comorbidity_chronic_cardiac_diseaseYes" ~
                              "Co-morbid Cardiac Disease Yes vs No",
                            Variable == "eotd_anticoagulants" ~"Anticoagulant use during ECMO Yes vs No",
                            Variable =="ecmo_vasoactive_drugs_beforeYes" ~
                              "*Pre-ECMO vasoactive medicine use Yes vs No",
                            Variable == "cannula_lumenSingle lumen" ~"Single vs double lumen cannula",
                            Variable == "days_vent_ecmo5" ~"Days ventilated pre-ECMO",
                            Variable == "log2(1/ecmo_worst_pa_o2_fi_o2_before)" ~"Halving of P/F ratio at Baseline",
                            Variable == "sofa" ~ "*SOFA",
                            Variable == "del_o2_plus50TRUE" ~"Delta O2 < -25 or > 25 vs -25 to 25",
                            Variable == "del_co2_neg50TRUE" ~"Delta CO2 < -25 or > 25 vs -25 to 25",
                            Variable == "rel_o2_plusTRUE" ~"Relative Delta O2 \u2265 50% vs < 50%",
                            Variable == "rel_co2_negTRUE" ~"Relative Delta CO2 < -50% vs \u2265 -50%",
                            Variable == "log2(1/ecmo_start_worst_platelet_count)" ~ 
                              "*Halving of Platelet Count at Baseline" ,
                            Variable == "era.L" ~ "Pandemic era Jul-Dec 2020 vs Jan-Jun 2020",
                            Variable == "era.Q" ~ "Pandemic era Jan-Sep 2021 vs Jan-Dec 2020",
                            TRUE ~ Variable))


header <- tibble(Variable = c("", "Variable"),
                 N = c("","N"),
                 HR = c("", "HR (95% CI)"))




hr_data <- bind_rows(header,hr_data)



hr_forest <- hr_data %>% 
  forestplot(labeltext = c(Variable, N, HR), 
             title = "Stroke within 90 days of ECMO (Univariable survival models)",
             graph.pos=3,
             xticks=xticks,
             boxsize = .25,
             xlog = F, 
             zero=1,
             txt_gp = fpTxtGp(label = gpar(cex = 0.8),
                              ticks = gpar(cex = 0.8),
                              xlab  = gpar(cex = 0.8)),

             col = fpColors(box = "royalblue",
                            line = "darkblue"),
             vertices = TRUE,
             xlab="Hazard Ratio")

# hr_forest


# ---- subHR ----
shr_ecmo <- univ_compete("ecmo", df= data.ids0e)
shr_age <- univ_compete(var="age", #var2="age + age^2",
                        df= ecmo_patients)
shr_sex <- univ_compete("sex",df= ecmo_patients)
shr_ethnic <- univ_compete("ethnic_white",df= ecmo_patients)
shr_smoke <- univ_compete("current_smoker",df= ecmo_patients)
shr_vent_ecmo <- univ_compete(var="days_vent_ecmo5",
                              # var2="days_vent_ecmo5 + days_vent_ecmo5^2",
                              df= ecmo_patients)
shr_cannula_lumen <- univ_compete("cannula_lumen",
                                  df= ecmo_patients)
# (shr_ac_before <- univ_compete("ac_before",df= data.ids0)) #all strokes had no AC before
shr_anticoagulants <- univ_compete("eotd_anticoagulants",
                                   df= data.ids0)
shr_obesity <- univ_compete( var= "comorbidity_obesity",df= ecmo_patients)
shr_vasoactive <- univ_compete("ecmo_vasoactive_drugs_before",
                               df= ecmo_patients)
shr_diabetes <- univ_compete( var= "comorbidity_diabetes",df= ecmo_patients)
shr_cardiac <- univ_compete("comorbidity_chronic_cardiac_disease",
                            df= ecmo_patients)
shr_hypertension <- univ_compete("comorbidity_hypertension",
                                 df= ecmo_patients)
# (shr_neuro <- univ_compete("comorbidity_chronic_neurological_disorder",df= data.ids0))
shr_pf <- univ_compete(var= "ecmo_worst_pa_o2_fi_o2_before", 
                       var2="log2(1/ecmo_worst_pa_o2_fi_o2_before)",df= ecmo_patients)
shr_sofa <- univ_compete(var= "sofa",#var2="sofa + sofa^2" ,
                         df= ecmo_patients)
shr_del_o2_plus50 <- univ_compete(var="del_o2_plus50",df= ecmo_patients)
shr_del_co2_neg50 <- univ_compete(var="del_co2_neg50",df= ecmo_patients)
shr_rel_o2_plus <- univ_compete(var="rel_o2_plus",df= ecmo_patients)
shr_rel_co2_neg <- univ_compete(var="rel_co2_neg",df= ecmo_patients)
shr_pc_under200 <- univ_compete(var="pc_under200",df= ecmo_patients)
shr_pc <- univ_compete(var="ecmo_start_worst_platelet_count",
                       var2="log2(1/ecmo_start_worst_platelet_count)", df= ecmo_patients)

shr_region <- univ_compete("income_region",
                           df= ecmo_patients)
shr_era <- univ_compete("era",df= ecmo_patients)



# fisher.test(ecmo_patients$comorbidity_chronic_neurological_disorder, ecmo_patients$status2)

### forestplot

shr_data <- data.frame(rbind(shr_ecmo,shr_age, shr_sex,shr_vent_ecmo,
                             shr_ethnic,
                             shr_smoke, shr_obesity, 
                             shr_diabetes, shr_cardiac,
                             shr_hypertension,
                             
                             shr_vasoactive,
                             shr_anticoagulants,
                             shr_cannula_lumen, 
                             shr_pf,
                             shr_sofa,
                             # shr_del_o2_plus50,
                             # shr_del_co2_neg50,
                             shr_rel_o2_plus,
                             shr_rel_co2_neg,
                             # shr_pc_under200,
                             shr_pc,
                             shr_region, shr_era
)) %>%
  rename(mean=est) %>% #, lower=V2, upper = V3
  mutate(sHR=paste0(roundz(mean, digits=2), " (",
                    roundz(lower, digits=2),", ",
                    roundz(upper, digits=2),")"),
         N=as.character(N))
shr_data$Variable = rownames(shr_data)
rownames(shr_data) <- c()

##nice labels
shr_data <- shr_data %>%
  mutate(Variable=case_when(Variable == "ecmo" ~"During vs post ECMO",
                            Variable == "age" ~"Age",
                            Variable == "sexMale" ~"Male vs Female",
                            Variable == "ethnic_whiteYes" ~"White Ethnicity Yes vs No",
                            Variable == "current_smokerYes" ~"Current Smoker Yes vs No",
                            Variable == "income_regionHigh Income" ~"High vs Middle Income Region",
                            Variable == "comorbidity_obesityYes" ~"*Obese Yes vs No",
                            Variable == "comorbidity_diabetesYes" ~"Co-morbid Diabetes Yes vs No",
                            Variable == "comorbidity_hypertensionYes" ~"Co-morbid Hypertension Yes vs No",
                            Variable == "comorbidity_chronic_cardiac_diseaseYes" ~
                              "Co-morbid Cardiac Disease Yes vs No",
                            Variable == "eotd_anticoagulants" ~"Anticoagulant use during ECMO Yes vs No",
                            Variable =="ecmo_vasoactive_drugs_beforeYes" ~
                              "*Pre-ECMO vasoactive medicine use Yes vs No",
                            Variable == "cannula_lumenSingle lumen" ~"Single vs double lumen cannula",
                            Variable == "days_vent_ecmo5" ~"Days ventilated pre-ECMO",
                            Variable == "log2(1/ecmo_worst_pa_o2_fi_o2_before)" ~"Halving of P/F ratio at Baseline",
                            Variable == "sofa" ~ "*SOFA",
                            Variable == "del_o2_plus50TRUE" ~"Delta O2 < -25 or > 25 vs -25 to 25",
                            Variable == "del_co2_neg50TRUE" ~"Delta CO2 < -25 or > 25 vs -25 to 25",
                            Variable == "rel_o2_plusTRUE" ~"Relative Delta O2 \u2265 50% vs < 50%",
                            Variable == "rel_co2_negTRUE" ~"Relative Delta CO2 < -50% vs \u2265 -50%",
                            Variable == "log2(1/ecmo_start_worst_platelet_count)" ~ 
                              "*Halving of Platelet Count at Baseline" ,
                            Variable == "era.L" ~ "Pandemic era Jul-Dec 2020 vs Jan-Jun 2020",
                            Variable == "era.Q" ~ "Pandemic era Jan-Sep 2021 vs Jan-Dec 2020",
                            TRUE ~ Variable))

headerc <- tibble(Variable = c("", "Variable"),
                  N=c("","N"),
                  sHR = c("", "subHR (95% CI)"))




shr_data <- bind_rows(headerc,shr_data)

shr_forest <- shr_data %>% 
  forestplot(labeltext = c(Variable, N, sHR),
             title = "Stroke within 90 days of ECMO (Univariable survival models)  \nDeath as competing risk",
             xticks=xticks,
             boxsize = .25,
             xlog = F, 
             zero=1,
             graph.pos=3,
             txt_gp = fpTxtGp(label = gpar(cex = 0.8),
                              ticks = gpar(cex = 0.8),
                              xlab  = gpar(cex = 0.8)),
             col = fpColors(box = "royalblue",
                            line = "darkblue"),
             vertices = TRUE,
             xlab="Subdistribution Hazard Ratio")

# shr_forest

# ---- cox_CR ----
##combine the cox and Cr forestplots

headerc <- tibble(Variable = c( "Variable","Variable"),
                  N=c("N","N"),
                  HR = c("Stroke subHR/HR (95% CI)","Stroke subHR/HR (95% CI)"),
                  Model = c("Competing Risks","Cox"))



shr_data2 <- shr_data %>%
  mutate(HR=paste0(" ",sHR))

surv_hr_data=bind_rows(shr_data2 %>% select(!sHR) %>% mutate(Model = "Competing Risks"),
                       hr_data  %>% mutate(Model = "Cox") #%>% select(!HR)
) %>%
  filter(!is.na(mean))

surv_hr_data=bind_rows(headerc,surv_hr_data)

all_forest <- surv_hr_data %>% 
  group_by(Model) %>%
  forestplot(labeltext = c(Variable, N, HR), #
             graph.pos=3,
             title = "Stroke within 90 days of ECMO (univariable survival models)",
             fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
             xticks=xticks,
             xlog = F, 
             zero=1,
             line.margin = 0.5, #.2,
             boxsize = .25,
             txt_gp = fpTxtGp(label = gpar(cex = 0.7),
                              ticks = gpar(cex = 0.7),
                              xlab  = gpar(cex = 0.7)),
             col = fpColors(box = c("cyan4","darkorange1"),
                            line = c("cyan4","darkorange1")),
             vertices = TRUE,
             xlab="subHR/HR")
# tiff(file="paper/Images/Fig2.tif", res=300, compression="lzw",
#      units="mm", width=170, height=170)

all_forest

# dev.off()
# ---- cox_CR_spline ----
library(ggplot2)
library(ggpubr)

testdata <- ecmo_patients %>%
  mutate(delta_o2 = ifelse(delta_o2 > 60 | delta_o2 < -60, NA, delta_o2),
        delta_co2 = ifelse(delta_co2 > 75 | delta_co2 < -75, NA, delta_co2),
        rel_delta_o2 = 100* ifelse(rel_delta_o2> 1, NA, rel_delta_o2),
        rel_delta_co2 = 100 * ifelse(rel_delta_co2> 1 | rel_delta_co2 < -0.75, 
                                     NA, rel_delta_co2),
        ecmo_worst_pa_co2_6hr_before = ifelse(ecmo_worst_pa_co2_6hr_before > 120,
                                              NA, ecmo_worst_pa_co2_6hr_before))

hr_spline <- function(var, xtitle = NULL, degf = 3) {
  
  form <- as.formula(
    paste0("Surv(tstart, tstop, any_stroke) ~ cluster(site_name) + pspline("
           , paste0(var),", df=",as.character(degf), ")"))
  
  
  
  c <-  coxph(form, 
              data = testdata, id=pin,
              na.action = na.omit, x=T)
  
  g = termplot(c, term = 1, se = T,plot = F)
  
  # npat <- data.ids0 %>% select(var, pin, site_name) %>% na.omit() %>% pull(pin) %>% unique()
  xtitle = paste0(xtitle, " (N=",as.character(c$n),")")
  
  
  newdf = data.frame(g) 
  colnames(newdf) = c("var", "y", "se")
  newdf <- newdf%>%
    mutate(hr=exp(y), 
           lci=exp(y-1.96*se),
           uci=exp(y+1.96*se))
  
  
  p <- newdf %>% 
    ggplot(aes(x=var, y=hr)) +
    geom_line(col='darkorange1') +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
    geom_hline(yintercept = 1, linetype="dashed", colour="darkgrey") +
    xlab(xtitle) + ylab("HR") +
    theme_bw()
  return(p)
  
}

plot_o2_pre <- hr_spline("ecmo_worst_pa_o2_6hr_before", "Pre-ECMO O2 (mmHg)",degf = 2) #
plot_co2_pre <- hr_spline("ecmo_worst_pa_co2_6hr_before", "Pre-ECMO CO2 (mmHg)",degf = 2) #

plot_del_o2 <- hr_spline("delta_o2", "Delta O2 (mmHg)",degf = 2) #
plot_del_co2 <- hr_spline("delta_co2", "Delta CO2 (mmHg)",degf = 2) #

plot_rel_del_o2 <- hr_spline("rel_delta_o2", "Relative Delta O2 (%)",degf = 2) #
plot_rel_del_co2 <- hr_spline("rel_delta_co2", "Relative Delta CO2 (%)",degf = 2) #


plot_ph <- hr_spline("ecmo_start_worst_p_h", "pH",degf = 2) #
plot_aptt <- hr_spline("ecmo_high_aptt", "aPTT",degf = 2)
plot_pc <- hr_spline("ecmo_start_worst_platelet_count", "Platelet Count",degf = 2)
# plot_d_dimer <- hr_spline("pre_ecmo_high_d_dimer", "D-Dimer",degf = 5)
plot_age <- hr_spline("age", "Age",degf = 2)
plot_vent <- hr_spline("days_vent_ecmo", "Ventilated days",degf = 2)
plot_sofa <- hr_spline("sofa", "SOFA",degf = 3)


ggarrange(          
  plot_rel_del_o2 + ggtitle("(a)") + theme(plot.title = element_text(hjust = 0.5)), 
  plot_rel_del_co2+ ggtitle("(b)") + theme(plot.title = element_text(hjust = 0.5)), 
  ncol=2, nrow=1)


ggarrange(plot_del_o2 + ggtitle("(a)") + theme(plot.title = element_text(hjust = 0.5)),
          plot_del_co2+ ggtitle("(b)") + theme(plot.title = element_text(hjust = 0.5)),
          plot_rel_del_o2 + ggtitle("(c)") + theme(plot.title = element_text(hjust = 0.5)),
          plot_rel_del_co2+ ggtitle("(d)") + theme(plot.title = element_text(hjust = 0.5)),
          ncol=2, nrow=2)





# 
# plot(x = newdf$delta_co2, y = newdf$hr, type = 'l', col = "red",xlab = "Delta CO2",
#      ylab = "Hazard ratio for Delta CO2 with the Cox model using psplines")
# lines(x = newdf$delta_co2, y = newdf$lci, col = "blue")
# lines(x = newdf$delta_co2, y = newdf$uci, col = "blue")
# abline(h = 1, col = "gray")



AIC(coxph(Surv(tstart, tstop, any_stroke) ~ cluster(site_name) + pspline(sofa,df=3), 
            data = testdata, id=pin,
            na.action = na.omit))


summary(coxph(Surv(tstart, tstop, any_stroke) ~ cluster(site_name) + pspline(delta_o2,df=2)  
              +  ecmo_worst_pa_o2_6hr_before
              ,data = testdata, id=pin))

cox_co2 <- coxph(Surv(tstart, tstop, any_stroke) ~  pspline(delta_co2,df=2) + 
                   pspline(ecmo_worst_pa_co2_6hr_before, df=2) + cluster(site_name),
                 data = testdata, id=pin, x=T)
summary(cox_co2)
termplot(cox_co2, term = 1, se = T,plot = T)


