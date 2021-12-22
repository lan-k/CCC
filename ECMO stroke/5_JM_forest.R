#forest plots of the joim model HRs
#---- jm_hr ----

rm(list=ls())
library(forestplot)
library(dplyr)
library(stringr)
source("99_functions.R")


load(file="Data/JM_reduced.Rdata")

jfit_red <- jointFit1.p1
jfit_red_ridge <- jointFitr

# 
# summary(jfit_red)
# summary(jfit_red_ridge)

load(file="Data/JM_nosmoke_ethnic.Rdata") #excludes current smoker and white ethnicity
# summary( jointFit1.p1)

## CR
load(file="Data/JM_CR_red.Rdata") 
load( file="Data/JM_CR_red_ICH.Rdata")

jm_hr <- function(jfit, jfitr = NULL) {
  
  stab <- summary(jfit)$Survival
  jmtab <- data.frame(round(exp(stab[c(1,3,4)]),digits=3))
  colnames(jmtab) = c("mean", "lower","upper")
  jmtab$Variable = rownames(jmtab)
  jmtab$Penalty = "Unpenalized"
  
  if (!is.null(jfitr)) {
    mean <-  exp(unlist(coef(jfitr)))
    statsr <- jfitr$statistics
    lower <- exp(unlist(list(statsr$CI_low$gammas, statsr$CI_low$alphas)))
    upper <- exp(unlist(list(statsr$CI_upp$gammas, statsr$CI_upp$alphas)))
    
    jmtabr <- data.frame(cbind(mean, lower, upper))
    jmtabr$Penalty = "Penalized"
    jmtabr$Variable = jmtab$Variable
    
    hr <- bind_rows(jmtab,jmtabr)
    
  } else {
    hr <- jmtab
  }
  
  rownames(hr) <- c()
  
  hr <- hr %>%
    mutate(Variable=case_when(Variable == "ecmo" ~"During vs post ECMO",
                              Variable == "age" ~"Age",
                              Variable == "sexMale" ~"Male vs Female",
                              Variable == "income_regionHigh Income" ~"High vs Middle Income Region",
                              Variable == "comorbidity_obesityYes" ~"Obese Yes vs No",
                              Variable == "comorbidity_diabetesYes" ~"Co-morbid Diabetes Yes vs No",
                              Variable == "comorbidity_hypertensionYes" ~"Co-morbid Hypertension Yes vs No",
                              Variable == "eotd_anticoagulants" ~"Anticoagulant use during ECMO Yes vs No",
                              Variable =="ecmo_vasoactive_drugs_beforeYes" ~
                                "Pre-ECMO vasoactive medicine use Yes vs No",
                              Variable == "days_vent_ecmo5" ~"Days ventilated pre-ECMO",
                              Variable == "era.L" ~ "Pandemic era Jul-Dec 2020 vs Jan-Jun 2020",
                              Variable == "era.Q" ~ "Pandemic era Jan-Sep 2021 vs Jan-Dec 2020",
                              Variable == "value(log2(pa_o2))" ~ "PaO2 (log2-transformed)",
                              Variable == "value(log2(pa_co2))" ~ "PaCO2 (log2-transformed)",
                              Variable == "value(log2(platelet_count))" ~ "Platelet count (log2-transformed)",
                              TRUE ~ Variable))
  return(hr)
  
}



##reduced model
jmhr_reduced <- jm_hr(jfit_red, jfit_red_ridge)
jmhr <- jm_hr(jointFit1.p1, jointFitr)


#---- jm_forest_red_un ----
##unpenalized
xticks <- c(.01, 0.02, 0.1, 0.2, 1, 2, 10)
xtlab <- rep(TRUE, length.out = length(xticks))
attr(xticks, "labels") <- xtlab

##reduced model
headercu <- tibble(Variable = c( "Variable"),
                  HR = c("HR (95% CI)"))


jm_hr_data_redu=jmhr_reduced %>%
  filter(!is.na(mean), Penalty == "Unpenalized") %>%
  mutate(HR=paste0(roundz(mean, digits=2), " (",
                   roundz(lower, digits=2),", ",
                   roundz(upper, digits=2),")"))


jm_hr_data_redu=bind_rows(headercu,jm_hr_data_redu) 

jm_forest_redu <- jm_hr_data_redu %>% 
  forestplot(labeltext = c(Variable, HR), #
             graph.pos=2,
             clip=c(0.01, 10),
             xticks = xticks, #c(.01, 0.02, 0.1, 0.2, 1, 2, 10),
             title = "Stroke within 90 days of ECMO (joint model 1)",
             xlog = T, 
             zero=1,
             line.margin = 0.5, #.2,
             boxsize = .25,
             txt_gp = fpTxtGp(label = gpar(cex = 0.8),
                              ticks = gpar(cex = 0.8),
                              xlab  = gpar(cex = 0.8)),
             col = fpColors(box = c("royalblue"),
                            line = c("darkblue")),
             # col = fpColors(box = c("royalblue","darkred"),
             #                line = c("darkblue","darkred")),
             vertices = TRUE,
             xlab="HR")

jm_forest_redu



#---- jm_forest_un ----

##unpenalized
##Full model
headercu <- tibble(Variable = c( "Variable"),
                  HR = c("HR (95% CI)"))


jm_hr_data=jmhr %>%
  filter(!is.na(mean), Penalty == "Unpenalized") %>%
  mutate(HR=paste0(roundz(mean, digits=2), " (",
                   roundz(lower, digits=2),", ",
                   roundz(upper, digits=2),")"))

jm_hr_data=bind_rows(headercu,jm_hr_data) 

jm_forestu<- jm_hr_data %>% 
  forestplot(labeltext = c(Variable, HR), #
             graph.pos=2,
             title = "Stroke within 90 days of ECMO (joint model 2)",
             xlog = T, 
             zero=1,
             clip=c(0.01, 10),
             xticks = xticks,
             line.margin = 0.5, #.2,
             boxsize = .25,
             txt_gp = fpTxtGp(label = gpar(cex = 0.8),
                              ticks = gpar(cex = 0.8),
                              xlab  = gpar(cex = 0.8)),
             col = fpColors(box = c("royalblue"),
                            line = c("darkblue")),
             vertices = TRUE,
             xlab="HR")

jm_forestu





#---- jm_forest_red ----
#both types of penalty
xticks <- c(.01, 0.02, 0.1, 0.2, 1, 2, 10)
xtlab <- rep(TRUE, length.out = length(xticks))
attr(xticks, "labels") <- xtlab

##reduced model
headerc <- tibble(Variable = c( "Variable","Variable"),
                  HR = c("HR (95% CI), Unpenalized/penalized","HR (95% CI), Unpenalized/penalized"),
                  Penalty = c("Unpenalized","Penalized"))




jm_hr_data_red=jmhr_reduced %>%
  filter(!is.na(mean)) %>%
  mutate(HR=paste0(roundz(mean, digits=2), " (",
                   roundz(lower, digits=2),", ",
                   roundz(upper, digits=2),")"),
         HR=ifelse(Penalty == "Penalized", paste0(" ", HR) ,HR))


jm_hr_data_red=bind_rows(headerc,jm_hr_data_red) 

jm_forest_red <- jm_hr_data_red %>% 
  group_by(Penalty) %>%
  forestplot(labeltext = c(Variable, HR), #
             graph.pos=2,
             clip=c(0.01, 10),
             xticks = xticks, #c(.01, 0.02, 0.1, 0.2, 1, 2, 10),
             title = "Stroke within 90 days of ECMO (joint model 1)",
             fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
             xlog = T, 
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
             xlab="HR")

jm_forest_red



#---- jm_forest_red ----
##Full model
headerc <- tibble(Variable = c( "Variable","Variable"),
                  HR = c("HR (95% CI), Unpenalized/penalized","HR (95% CI), Unpenalized/penalized"),
                  Penalty = c("Unpenalized","Penalized"))


jm_hr_data=jmhr %>%
  filter(!is.na(mean)) %>%
  mutate(HR=paste0(roundz(mean, digits=2), " (",
                   roundz(lower, digits=2),", ",
                   roundz(upper, digits=2),")"),
         HR=ifelse(Penalty == "Penalized", paste0(" ", HR) ,HR))

jm_hr_data=bind_rows(headerc,jm_hr_data) 

jm_forest<- jm_hr_data %>% 
  group_by(Penalty) %>%
  forestplot(labeltext = c(Variable, HR), #
             graph.pos=2,
             title = "Stroke within 90 days of ECMO (joint model 2)",
             fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
             xlog = T, 
             zero=1,
             clip=c(0.01, 10),
             xticks = xticks,
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

jm_forest


# ---- jm_forest_CR ----
##reduced competing risks model
## Does NOT use a Fine-Gray model
## different strata for different outcomes


jm_redu_CR <- jm_hr(jointFit1.CR)


xticks <- c(0.01,.01, 0.02, 0.1, 0.2, 1, 2, 10, 20)
xtlab <- rep(TRUE, length.out = length(xticks))
attr(xticks, "labels") <- xtlab



headercu <- tibble(Variable = c( "Variable", "Variable"),
                   HR = c("HR (95% CI) Stroke, Death","HR (95% CI) Stroke, Death"),
                   Outcome=c("Stroke","Death"))


jm_hr_data_redu_CR=jm_redu_CR %>%
  filter(!is.na(mean), Penalty == "Unpenalized") %>%
  mutate(HR=paste0(roundz(mean, digits=3), " (",
                   roundz(lower, digits=3),", ",
                   roundz(upper, digits=3),")"),
         Outcome = ifelse(grepl("Stroke", Variable), "Stroke", "Death"),
         Variable = ifelse(grepl(":", Variable),str_extract(Variable, "^.*(?=(:))"),
                           Variable),
         Variable=case_when(Variable == "ecmo" ~"During vs post ECMO",
                            Variable == "age" ~"Age",
                            Variable == "sexMale" ~"Male vs Female",
                            Variable == "income_regionHigh Income" ~"High vs Middle Income Region",
                            Variable == "comorbidity_obesityYes" ~"Obese Yes vs No",
                            Variable == "comorbidity_diabetesYes" ~"Co-morbid Diabetes Yes vs No",
                            Variable == "comorbidity_hypertensionYes" ~"Co-morbid Hypertension Yes vs No",
                            Variable == "eotd_anticoagulants" ~"Anticoagulant use during ECMO Yes vs No",
                            Variable =="ecmo_vasoactive_drugs_beforeYes" ~
                              "Pre-ECMO vasoactive medicine use Yes vs No",
                            Variable == "days_vent_ecmo5" ~"Days ventilated pre-ECMO",
                            Variable == "era.L" ~ "Pandemic era Jul-Dec 2020 vs Jan-Jun 2020",
                            Variable == "era.Q" ~ "Pandemic era Jan-Sep 2021 vs Jan-Dec 2020",
                            Variable == "value(log2(pa_o2))" ~ "PaO2 (log2-transformed)",
                            Variable == "value(log2(pa_co2))" ~ "PaCO2 (log2-transformed)",
                            Variable == "value(log2(platelet_count))" ~ "Platelet count (log2-transformed)",
                            TRUE ~ Variable)) %>%
  select(!Penalty) #%>%
  #arrange(desc(Outcome), Variable)


jm_hr_data_redu_CR=bind_rows(headercu,jm_hr_data_redu_CR) 

jm_forest_redu_CR <- jm_hr_data_redu_CR %>% 
  group_by(Outcome) %>%
  forestplot(labeltext = c(Variable, HR), #
             graph.pos=2,
             clip=c(0.01, 20),
             xticks = xticks, #c(.01, 0.02, 0.1, 0.2, 1, 2, 10),
             title = "Stroke within 90 days of ECMO (death as competing risk)",
             fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
             xlog = T, 
             zero=1,
             # line.margin = 0.5, #.2,
             boxsize = 0.1, #.25,
             txt_gp = fpTxtGp(label = gpar(cex = 0.8),
                              ticks = gpar(cex = 0.8),
                              xlab  = gpar(cex = 0.8)),
             # col = fpColors(box = c("royalblue"),
             #                line = c("darkblue")),
             col = fpColors(box = c("cyan4","darkorange1"),
                            line = c("cyan4","darkorange1")),
             vertices = TRUE,
             xlab="HR")

jm_forest_redu_CR



# ---- jm_forest_CR_ICH ----
##reduced competing risks model
## Does NOT use a Fine-Gray model
## different strata for different outcomes


jm_redu_CR_ICH <- jm_hr(jointFit1.CR_ICH)


xticks <- c(0.01,.01, 0.02, 0.1, 0.2, 1, 2, 10, 20)
xtlab <- rep(TRUE, length.out = length(xticks))
attr(xticks, "labels") <- xtlab



headercu <- tibble(Variable = c( "Variable", "Variable"),
                   HR = c("HR (95% CI) Hemorrhagic Stroke, Death",
                          "HR (95% CI) Hemorrhagic Stroke, Death"),
                   Outcome=c("Hemorrhagic Stroke","Death"))


jm_hr_data_redu_CR_ICH=jm_redu_CR_ICH %>%
  filter(!is.na(mean), Penalty == "Unpenalized",
         !grepl("Other stroke", Variable)) %>%  #exclude other strokes
  mutate(HR=paste0(roundz(mean, digits=3), " (",
                   roundz(lower, digits=3),", ",
                   roundz(upper, digits=3),")"),
         Outcome = ifelse(grepl("Stroke", Variable), "Hemorrhagic Stroke", "Death"),
         Variable = ifelse(grepl(":", Variable),str_extract(Variable, "^.*(?=(:))"),
                           Variable),
         Variable=case_when(Variable == "ecmo" ~"During vs post ECMO",
                            Variable == "age" ~"Age",
                            Variable == "sexMale" ~"Male vs Female",
                            Variable == "income_regionHigh Income" ~"High vs Middle Income Region",
                            Variable == "comorbidity_obesityYes" ~"Obese Yes vs No",
                            Variable == "comorbidity_diabetesYes" ~"Co-morbid Diabetes Yes vs No",
                            Variable == "comorbidity_hypertensionYes" ~"Co-morbid Hypertension Yes vs No",
                            Variable == "eotd_anticoagulants" ~"Anticoagulant use during ECMO Yes vs No",
                            Variable =="ecmo_vasoactive_drugs_beforeYes" ~
                              "Pre-ECMO vasoactive medicine use Yes vs No",
                            Variable == "days_vent_ecmo5" ~"Days ventilated pre-ECMO",
                            Variable == "era.L" ~ "Pandemic era Jul-Dec 2020 vs Jan-Jun 2020",
                            Variable == "era.Q" ~ "Pandemic era Jan-Sep 2021 vs Jan-Dec 2020",
                            Variable == "value(log2(pa_o2))" ~ "PaO2 (log2-transformed)",
                            Variable == "value(log2(pa_co2))" ~ "PaCO2 (log2-transformed)",
                            Variable == "value(log2(platelet_count))" ~ "Platelet count (log2-transformed)",
                            TRUE ~ Variable)) %>%
  select(!Penalty) #%>%
#arrange(desc(Outcome), Variable)


jm_hr_data_redu_CR_ICH=bind_rows(headercu,jm_hr_data_redu_CR_ICH) 

jm_forest_redu_CR_ICH <- jm_hr_data_redu_CR_ICH %>% 
  group_by(Outcome) %>%
  forestplot(labeltext = c(Variable, HR), #
             graph.pos=2,
             clip=c(0.01, 20),
             xticks = xticks, #c(.01, 0.02, 0.1, 0.2, 1, 2, 10),
             title = "Hemorrhagic stroke within 90 days of ECMO (death as competing risk)",
             fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
             xlog = T, 
             zero=1,
             # line.margin = 0.5, #.2,
             boxsize = 0.1, #.25,
             txt_gp = fpTxtGp(label = gpar(cex = 0.8),
                              ticks = gpar(cex = 0.8),
                              xlab  = gpar(cex = 0.8)),
             # col = fpColors(box = c("royalblue"),
             #                line = c("darkblue")),
             col = fpColors(box = c("cyan4","darkorange1"),
                            line = c("cyan4","darkorange1")),
             vertices = TRUE,
             xlab="HR")

jm_forest_redu_CR_ICH









