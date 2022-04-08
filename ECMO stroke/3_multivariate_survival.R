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
library(lmtest)
options(scipen=999) # avoid scientific presentation






# ---- survival_data ----

start_date = as.Date(data_date) - 90 #as.Date("2021-08-11")



covars <- c("age", "sex", "comorbidity_obesity","comorbidity_hypertension",
            "comorbidity_diabetes","comorbidity_chronic_neurological_disorder",
            "ecmo_worst_pa_o2_fi_o2_before",
            # "ecmo_prone_before",
            "income_region", "site_name", "days_vent_ecmo5", "era")

#start on the first day of ecmo up to 90 days days_fup

surv_vars <- c("ecmo","age", "sex", 
               "days_vent_ecmo5",
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
               
               "income_region", "era")


gas_vars <- c("age", 
               "sofa",
               "comorbidity_obesity",
               "ecmo_vasoactive_drugs_before",
              "delta_o2",
              "delta_co2")



###need to include anticoagulants as a time varying exposure
data.ids0 <- ecmo_daily %>% # jmdat0  %>% mutate(age_day = age + day) %>%
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





##the Cox model  - any stroke or death as the outcome

# ---- multi_surv ----

sform = as.formula(
  paste0("Surv(tstart, tstop, any_stroke) ~ cluster(pin, site_name) + ecmo + " 
         , paste0(surv_vars, collapse="+")))

#days_vent_ecmo5^2 + log2(1/ecmo_start_worst_platelet_count) +

survFit.p0 <- coxph(sform ,
                    data = data.ids0, id=pin)
# cox.zph(survFit.p0) #all PH
coxsum <- summary(survFit.p0)
cox_data <- data.frame(coxsum$conf.int[,c(1,3:4)])
cox_data$Variable=rownames(cox_data)
rownames(cox_data) <- c()

npat <- data.ids0 %>% select(surv_vars,  ecmo, pin) %>% na.omit() %>% pull(pin) %>% unique()
cox_data <- cox_data %>%
  rename(mean=exp.coef.,
         lower=lower..95,
         upper=upper..95)%>%
  mutate(
         HR=paste0(roundz(mean, digits=2), " (",
                   roundz(lower, digits=2),", ",
                   roundz(upper, digits=2),")")
         ) %>%
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
                            Variable == "era.L" ~ "Pandemic era Jul-Dec 2020 vs Jan-Jun 2020",
                            Variable == "era.Q" ~ "Pandemic era Jan-Sep 2021 vs Jan-Dec 2020",
                            TRUE ~ Variable))

# ---- multi_CR ----
# the CR Fine-Gray model

etime1 <- data.ids0$tstart
etime2 <- data.ids0$tstop
event <- data.ids0$status2 #(0=censored, 1 = stroke, 2=discharged, 3= death) #was status2

crdata <- data.ids0 %>% select(pin, status2,tstart,tstop,  site_name, all_of(surv_vars))
crdata <- finegray(Surv(etime1, etime2, event) ~ ., data=crdata,id=pin,
                   etype="Stroke")  
#creates the Fine-Gray weights for stroke as the event

sform2 = as.formula(
  paste0("Surv(fgstart, fgstop, fgstatus) ~ cluster(pin, site_name) + ecmo + "
         , paste0(surv_vars, collapse="+")))


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
                            Variable == "era.L" ~ "Pandemic era Jul-Dec 2020 vs Jan-Jun 2020",
                            Variable == "era.Q" ~ "Pandemic era Jan-Sep 2021 vs Jan-Dec 2020",
                            TRUE ~ Variable))


# ---- multi_cox_CR ----
##combine the cox and Cr forestplots
xticks <- c(0,1,2,4,6,8, 10)
xtlab <- rep(TRUE, length.out = length(xticks))
attr(xticks, "labels") <- xtlab

headerc <- tibble(Variable = c( "Variable","Variable"),
                  HR = c("Stroke subHR/HR (95% CI)","Stroke subHR/HR (95% CI)"),
                  Model = c("Competing Risks","Cox"))



fg_data <- fg_data %>%
  mutate(HR=paste0(" ",sHR))

multi_hr_data=bind_rows(fg_data %>% select(!sHR) %>% mutate(Model = "Competing Risks"),
                       cox_data  %>% mutate(Model = "Cox") #%>% select(!HR)
) %>%
  filter(!is.na(mean))

multi_hr_data=bind_rows(headerc,multi_hr_data)

title = paste0("Stroke within 90 days of ECMO (multivariable survival models) \nN=",
               as.character(length(npat)))

multi_forest <- multi_hr_data %>% 
  group_by(Model) %>%
  forestplot(labeltext = c(Variable, HR), #
             graph.pos=2,
             title = title,
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
             # col = fpColors(box = c("royalblue","darkred"),
             #                line = c("darkblue","darkred")),
             vertices = TRUE,
             xlab="subHR/HR")

tiff(file="paper/Images/Fig3.tif", res=300, compression="lzw",
     units="mm", width=170, height=140)

multi_forest

dev.off()

# ---- cox_gas ----
##survival model with the gos and other variables with > 50% missing data
# sformgas = as.formula(
#   paste0("Surv(tstart, tstop, any_stroke) ~ cluster(site_name) + " #days_vent_ecmo5^2 + log2(ecmo_worst_pa_o2_fi_o2_before) +
#          , paste0(gas_vars, collapse="+")))
# 
# survFit.pgas <- coxph(sformgas ,
#                     data = ecmo_patients, id=pin)
# # cox.zph(survFit.p0) #all PH
# coxsumgas <- summary(survFit.pgas)
# cox_datagas <- data.frame(coxsumgas$conf.int[,c(1,3:4)])
# cox_datagas$Variable=rownames(cox_datagas)
# rownames(cox_datagas) <- c()
# 
# npatgas <- data.ids0 %>% select(gas_vars, ecmo, pin) %>% na.omit() %>% pull(pin) %>% unique()
# cox_datagas <- cox_datagas %>%
#   rename(mean=exp.coef.,
#          lower=lower..95,
#          upper=upper..95)%>%
#   mutate(
#     HR=paste0(roundz(mean, digits=2), " (",
#               roundz(lower, digits=2),", ",
#               roundz(upper, digits=2),")")
#   ) %>%
#   mutate(Variable=case_when(
#     Variable == "age" ~"Age",
#     Variable == "comorbidity_obesityYes" ~"Obese Yes vs No",
#     Variable =="ecmo_vasoactive_drugs_beforeYes" ~
#       "Pre-ECMO vasoactive medicine use Yes vs No",
#     Variable == "sofa" ~"SOFA",
#     Variable == "rel_o2_plus50TRUE" ~"Relative Delta O2 > 50% vs \u2264 50%",
#     Variable == "rel_co2_neg50TRUE" ~"Relative Delta CO2 < -50% vs \u2265 -50%",
#     TRUE ~ Variable))
# 
# 
# # ---- gas_CR ----
# # the CR Fine-Gray model for gases
# 
# etime1 <- ecmo_patients$tstart
# etime2 <- ecmo_patients$tstop
# event <- ecmo_patients$status #(0=censored, 1 = stroke, 2= death)
# 
# crdata <- ecmo_patients %>% select(pin, status2,tstart,tstop,  site_name, 
#                                ecmo_worst_pa_o2_fi_o2_before, all_of(gas_vars))
# crdata <- finegray(Surv(etime1, etime2, event) ~ ., data=crdata,id=pin,
#                    etype="Stroke")  
# #creates the Fine-Gray weights for stroke as the event
# 
# sformFGgas = as.formula(
#   paste0("Surv(fgstart, fgstop, fgstatus) ~ cluster(pin, site_name) + " #log2(ecmo_worst_pa_o2_fi_o2_before) + 
#          , paste0(gas_vars, collapse="+")))
# 
# 
# fitFGgas <- coxph( sformFGgas, id=pin,weight=fgwt,
#                 data=crdata)
# 
# mean <- fitFGgas$coefficients
# se <- sqrt(diag(fitFGgas$var))
# lower = mean - 1.96 * se
# upper = mean + 1.96 * se
# fg_datagas <- data.frame(exp(cbind(mean, lower, upper)))
# 
# fg_datagas$Variable=rownames(fg_datagas)
# 
# rownames(fg_datagas) <- c()
# 
# fg_datagas <- fg_datagas %>%
#   mutate(
#     sHR=paste0(roundz(mean, digits=2), " (",
#                roundz(lower, digits=2),", ",
#                roundz(upper, digits=2),")")
#   ) %>%
#   mutate(Variable=case_when(
#                             Variable == "age" ~"Age",
#                             Variable == "comorbidity_obesityYes" ~"Obese Yes vs No",
#                             Variable =="ecmo_vasoactive_drugs_beforeYes" ~
#                               "Pre-ECMO vasoactive medicine use Yes vs No",
#                             Variable == "sofa" ~"SOFA",
#                             Variable == "rel_o2_plus50TRUE" ~"Relative Delta O2 > 50% vs \u2264 50%",
#                             Variable == "rel_co2_neg50TRUE" ~"Relative Delta CO2 < -50% vs \u2265 -50%",
#                             TRUE ~ Variable))
# 
# 
# # ---- multi_cox_CR_gas ----
# ##combine the cox and Cr forestplots
# xticks <- c(0.1,1, 2,10, 20, 100)
# xtlab <- rep(TRUE, length.out = length(xticks))
# attr(xticks, "labels") <- xtlab
# 
# headerc <- tibble(Variable = c( "Variable","Variable"),
#                   HR = c("Stroke subHRHR (95% CI)","Stroke subHR/HR (95% CI)"),
#                   Model = c("Competing Risks","Cox"))
# 
# 
# 
# fg_datagas <- fg_datagas %>%
#   mutate(HR=paste0(" ",sHR))
# 
# multi_hr_datagas=bind_rows(fg_datagas %>% select(!sHR) %>% mutate(Model = "Competing Risks"),
#                         cox_datagas  %>% mutate(Model = "Cox") #%>% select(!HR)
# ) %>%
#   filter(!is.na(mean))
# 
# multi_hr_datagas=bind_rows(headerc,multi_hr_datagas)
# 
# title = paste0("Stroke within 90 days of ECMO (multivariable survival models) \nN=",
#                as.character(length(npatgas)), "\n",
#                as.character(fitFGgas$nevent), " strokes")
# 
# multi_forestgas <- multi_hr_datagas %>% 
#   group_by(Model) %>%
#   forestplot(labeltext = c(Variable, HR), #
#              graph.pos=2,
#              title = title,
#              fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
#              xticks=xticks,
#              xlog = T, 
#              zero=1,
#              line.margin = 0.5, #.2,
#              boxsize = .25,
#              txt_gp = fpTxtGp(label = gpar(cex = 0.8),
#                               ticks = gpar(cex = 0.8),
#                               xlab  = gpar(cex = 0.8)),
#              col = fpColors(box = c("cyan4","darkorange1"),
#                             line = c("cyan4","darkorange1")),
#              # col = fpColors(box = c("royalblue","darkred"),
#              #                line = c("darkblue","darkred")),
#              vertices = TRUE,
#              xlab="subHR/HR")
# 
# # jpeg('figures/gas_forest_plot.jpg', width=8, height=5, units='in', res=300, quality=100)
# # multi_forestgas
# # dev.off()



# ---- glmnet_surv ----

##lasso regression
library(glmnet)

# df <- ecmo_patients %>%
#   select(all_of(surv_vars)) %>%
#   na.omit() %>%
#   mutate(#ecmo=factor(ecmo), 
#     days_fup=ifelse(days_fup == 0, 0.1, days_fup),
#     ac_during = factor(case_when(anticoagulants %in% c("AC before & during ECMO",
#                                                        "AC during but not before ECMO") ~ "Yes",
#                                  !is.na(anticoagulants) ~ "No")),
#     sex=factor(sex),
#     # pin=factor(pin),
#     comorbidity_obesity=factor(comorbidity_obesity),
#     # comorbidity_hypertension = factor(comorbidity_hypertension),
#     comorbidity_diabetes = factor(comorbidity_diabetes),
#     # comorbidity_chronic_cardiac_disease = factor(comorbidity_chronic_cardiac_disease),
#     ecmo_vasoactive_drugs_before = factor(ecmo_vasoactive_drugs_before),
#     # cannula_lumen = factor(cannula_lumen),
#     income_region = factor(income_region),
#     # site_name = factor(site_name),
#     era=factor(era)) %>%
#   select(!anticoagulants)
# 
# 
# x <- model.matrix(~. , df %>% select(!c(days_fup, any_stroke)))
# 
# y <- Surv(df$days_fup,  df$any_stroke)
# 
# #multivariate cox model
# coxph_fit <- coxph(y ~ x, ties = "breslow")
# summary(coxph_fit)
# 
# 
# 
# fitl <- glmnet(x, y, family = "cox", lambda=0) #not converging
# cv.fit <- cv.glmnet(x, y, family = "cox", nfolds = 5)
# cv.fit$lambda.min; cv.fit$lambda.1se  #the same
# plot(cv.fit)
# 
# 
# glmnet_fit <- glmnet(x, y, family = "cox", lambda = 0) #ridge
# 
# plot(coef(glmnet_fit), coef(coxph_fit))
# abline(0, 1)
# 
# coef(fitl, s = "lambda.min")
# coef(coxph_fit)  #they are the same
# 




