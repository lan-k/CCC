###univariate Cox and CR models
##version with discharge and death as competing risks

# ---- get_data_disch ----
rm(list=ls())
source('2_univariate_survival.R') # prepares the data and runs exclusions
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


# ---- disch_funcs ----

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


univ_compete_disch <- function(var,var2=var,stat=status, stype="Stroke", df) { #stat=status2
  
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









# ---- subHR_disch ----

xticks <- c(0,1,2,4,6,8,10)
xtlab <- rep(TRUE, length.out = length(xticks))
attr(xticks, "labels") <- xtlab


shr_ecmo <- univ_compete_disch("ecmo", df= data.ids0e)
shr_age <- univ_compete_disch(var="age", #var2="age + age^2",
                        df= ecmo_patients)
shr_sex <- univ_compete_disch("sex",df= ecmo_patients)
shr_ethnic <- univ_compete_disch("ethnic_white",df= ecmo_patients)
shr_smoke <- univ_compete_disch("current_smoker",df= ecmo_patients)
shr_vent_ecmo <- univ_compete_disch(var="days_vent_ecmo5",
                              # var2="days_vent_ecmo5 + days_vent_ecmo5^2",
                              df= ecmo_patients)
shr_cannula_lumen <- univ_compete_disch("cannula_lumen",
                                  df= ecmo_patients)
# (shr_ac_before <- univ_compete_disch("ac_before",df= data.ids0)) #all strokes had no AC before
shr_anticoagulants <- univ_compete_disch("eotd_anticoagulants",
                                   df= data.ids0)
shr_obesity <- univ_compete_disch( var= "comorbidity_obesity",df= ecmo_patients)
shr_vasoactive <- univ_compete_disch("ecmo_vasoactive_drugs_before",
                               df= ecmo_patients)
shr_diabetes <- univ_compete_disch( var= "comorbidity_diabetes",df= ecmo_patients)
shr_cardiac <- univ_compete_disch("comorbidity_chronic_cardiac_disease",
                            df= ecmo_patients)
shr_hypertension <- univ_compete_disch("comorbidity_hypertension",
                                 df= ecmo_patients)
# (shr_neuro <- univ_compete_disch("comorbidity_chronic_neurological_disorder",df= data.ids0))
shr_pf <- univ_compete_disch(var= "ecmo_worst_pa_o2_fi_o2_before", 
                       var2="log2(1/ecmo_worst_pa_o2_fi_o2_before)",df= ecmo_patients)
shr_sofa <- univ_compete_disch(var= "sofa",#var2="sofa + sofa^2" ,
                         df= ecmo_patients)
shr_del_o2_plus50 <- univ_compete_disch(var="del_o2_plus50",df= ecmo_patients)
shr_del_co2_neg50 <- univ_compete_disch(var="del_co2_neg50",df= ecmo_patients)
shr_rel_o2_plus <- univ_compete_disch(var="rel_o2_plus",df= ecmo_patients)
shr_rel_co2_neg <- univ_compete_disch(var="rel_co2_neg",df= ecmo_patients)
shr_pc_under200 <- univ_compete_disch(var="pc_under200",df= ecmo_patients)
shr_pc <- univ_compete_disch(var="ecmo_start_worst_platelet_count",
                       var2="log2(1/ecmo_start_worst_platelet_count)", df= ecmo_patients)

shr_region <- univ_compete_disch("income_region",
                           df= ecmo_patients)
shr_era <- univ_compete_disch("era",df= ecmo_patients)



# fisher.test(ecmo_patients$comorbidity_chronic_neurological_disorder, ecmo_patients$status2)

### forestplot

shr_data_disch <- data.frame(rbind(shr_ecmo,shr_age, shr_sex,shr_vent_ecmo,
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
shr_data_disch$Variable = rownames(shr_data_disch)
rownames(shr_data_disch) <- c()

##nice labels
shr_data_disch <- shr_data_disch %>%
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
                            Variable == "sofa" ~"*SOFA",
                            Variable == "del_o2_plus50TRUE" ~"Delta O2 < -25 or > 25 vs -25 to 25",
                            Variable == "del_co2_neg50TRUE" ~"Delta CO2 < -25 or > 25 vs -25 to 25",
                            Variable == "rel_o2_plusTRUE" ~"Relative Delta O2 \u2265 50% vs < 50%",
                            Variable == "rel_co2_negTRUE" ~"Relative Delta CO2 < -50% vs \u2265 -50%",
                            Variable == "log2(1/ecmo_start_worst_platelet_count)" ~ "*Halving of Platelet Count at Baseline" ,
                            Variable == "era.L" ~ "Pandemic era Jul-Dec 2020 vs Jan-Jun 2020",
                            Variable == "era.Q" ~ "Pandemic era Jan-Sep 2021 vs Jan-Dec 2020",
                            TRUE ~ Variable))

headerc <- tibble(Variable = c("", "Variable"),
                  N=c("","N"),
                  sHR = c("", "subHR (95% CI)"))




shr_data_disch <- bind_rows(headerc,shr_data_disch)

shr_forest_disch <- shr_data_disch %>% 
  forestplot(labeltext = c(Variable, N, sHR),
             title = "Stroke within 90 days of ECMO (Univariable survival models)  \nDeath and discharge as competing risks",
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

# shr_forest_disch

# ---- cox_CR_disch ----
##combine the cox and Cr forestplots

headerc <- tibble(Variable = c( "Variable","Variable"),
                  N=c("N","N"),
                  HR = c("Stroke subHR/HR (95% CI)","Stroke subHR/HR (95% CI)"),
                  Model = c("Competing risk of death and discharge","Cox"))



shr_data2_disch <- shr_data_disch %>%
  mutate(HR=paste0(" ",sHR))

surv_hr_data_disch=bind_rows(shr_data2_disch %>% select(!sHR) %>% 
                               mutate(Model = "Competing risk of death and discharge"),
                       hr_data  %>% mutate(Model = "Cox") #%>% select(!HR)
) %>%
  filter(!is.na(mean))


surv_hr_data_disch=bind_rows(headerc,surv_hr_data_disch)

all_forest_disch <- surv_hr_data_disch %>% 
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

# tiff(file="paper/Images/eFig1.tif", res=300, compression="lzw",
#      units="mm", width=170, height=170)
all_forest_disch
# dev.off()

