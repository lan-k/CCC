###descriptives

library(tableone)
library(dplyr)
library(visdat)
library(summarytools)
library(pander)
library(naniar)
library(arsenal)

source('99_functions.R')

filename = 'Data/COVID_Data_stroke.RData' # stroke
load(file=filename)

complic <- patients %>% filter(any_ecmo == 1) %>% 
  select(contains("complic")) 
names(complic)
summary(complic)



source('1_ecmo_ich_data_prep.R') # prepares the data and runs exclusions 
cat_vars = c('sex','ethnicity','country','ecmo_type','comorbidity_chronic_kidney_disease',
             'comorbidity_severe_liver_disease','comorbidity_diabetes','comorbidity_hypertension')
#restrict to just ECMO
table(patients$any_ecmo, useNA = "always")

ecmo_patients <- patients %>% filter(any_ecmo=="Yes")
combined_ecmo <- combined %>% filter(any_ecmo=='Yes')

ecmo_patients <- ecmo_patients %>%
  replace_na(all_of(cat_vars),'')

##missing data
for_missing = select(ecmo_patients, -pin,  -stroke_type) %>% # remove ID; any_ecmo is never missing; stroke_type is a list variable
  filter(any_ecmo == "Yes") %>%
mutate(
  ethnicity = ifelse(ethnicity=='Missing', NA, ethnicity),
  outcome = ifelse(outcome %in% c('','unknown'), NA, outcome)) %>% # tidy up missing outcome
  mutate_if(~is.character(.x), ~ifelse(.x=='', NA, .x)) %>%
  select(-any_ecmo)


ncols = ncol(for_missing) 
nhalf = round(ncols/2)
half1 <- for_missing %>% select(c(1:nhalf))
half2 <- for_missing %>% select(!c(1:nhalf))
half1  %>%
  vis_dat()
half2 %>% 
  vis_dat()



missing_gender = sum(is.na(for_missing$sex))


#missing data for ecmo patients only
for_missing %>% 
  filter(any_ecmo == "Yes") %>% 
  select(-any_ecmo, site_name) %>%
  vis_dat()


###### Table of percentages of missing


# added 18 March 2021
missing_perc = mutate_all(for_missing, is.na) %>% # convert all to binary missing yes or no
  pivot_longer(cols=everything()) %>% # long format
  group_by(name, value) %>% # count missing
  tally() %>%
  group_by(name) %>%
  mutate(p = round(prop.table(n)*100),
         value = ifelse(value==TRUE, 'Missing', 'Not_missing')) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = value, values_from=p) %>%
  arrange(Missing) %>%
  select(-'Not_missing') #simplify
ftab = flextable(missing_perc) %>%
  theme_box() %>%
  autofit()
ftab

##ECMO patients only
missing_perc_ecmo <- for_missing %>%
  filter(any_ecmo == "Yes") %>%
  select(-any_ecmo) %>%
  mutate_all(is.na) %>% # convert all to binary missing yes or no
  pivot_longer(cols=everything()) %>% # long format
  group_by(name, value) %>% # count missing
  tally() %>%
  group_by(name) %>%
  mutate(p = round(prop.table(n)*100),
         value = ifelse(value==TRUE, 'Missing', 'Not_missing')) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = value, values_from=p) %>%
  arrange(Missing) %>%
  select(-'Not_missing') #simplify
ftab_ecmo = flextable(missing_perc_ecmo) %>%
  theme_box() %>%
  autofit()
ftab_ecmo


tab_any = with(patients, freq(any_ecmo, round.digits=1, cumul = FALSE, report.nas	=FALSE))
knitr::kable(tab_any)
cont_vars = c('sofa')
tab = select(ecmo_patients, pin, all_of(cont_vars)) %>%
  tidyr::gather(key='Variable', value='Value', -pin) %>%
  mutate(Variable=toupper(Variable)) %>%
  group_by(Variable) %>%
  summarise(Missing = sum(is.na(Value)),
            Median = roundz(median(Value, na.rm=TRUE),0),
            Q1 = roundz(quantile(Value, probs=0.25, na.rm=TRUE),0),
            Q3 = roundz(quantile(Value, probs=0.75, na.rm=TRUE),0),
            IQR = paste(Q1 , ' to ', Q3, sep=''),
            Min = roundz(min(Value, na.rm=TRUE),0),
            Max = roundz(max(Value, na.rm=TRUE),0)) %>%
  select(Variable, Missing, Median, IQR, Min, Max) %>%
  mutate(Variable = ifelse(Variable=='age', 'Age', Variable))
ftab = flextable(tab) %>% 
  theme_box() %>%
  autofit()
ftab



###plots
to_plot <- patients %>% 
  filter(any_ecmo=="Yes") %>%
  select( pin,  height, weight, any_ecmo) %>%
  gather(key='var', value='res', -pin,-any_ecmo) %>%
  filter(!is.na(res))


vplot = ggplot(data=to_plot, aes(x=any_ecmo, y=res)) +
  geom_boxplot()+
  geom_violin(alpha=0.5)+
  facet_wrap(~var)+
  g.theme+
  xlab('')+
  ylab('Height or weight')+
  #scale_y_continuous(breaks=NULL)+
  coord_flip() +
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
 scale_fill_manual( values=cbPalette)
vplot



to_plot = filter(combined_ecmo, 
                 !is.na(date_icu),
                 !is.na(date_daily)) %>%
  mutate(day = as.numeric(date_daily - date_icu)) %>%
  filter(day <= 30)


tplot = ggplot(data=to_plot, aes(x=factor(day), y=platelet_count))+ #, col=factor(any_ecmo
  geom_boxplot()+
  xlab('Day in ICU')+
  ylab('Platelet count')+
  scale_color_manual( values=c('darkorange1'))+ #'Any ECMO',,'cyan4'
  g.theme
tplot


to_stat = group_by(to_plot, day) %>% #, any_ecmo
  summarise(mean = mean(blood_urea_nitrogen, na.rm = TRUE),
            se = sem(blood_urea_nitrogen)) %>%
  ungroup() %>%
  mutate(lower = mean - (z*se),
         upper = mean + (z*se)
         #day = ifelse(any_ecmo=='Yes', day+0.2, day)
  ) # jitter x slightly to avoid overlap
tplot = ggplot(data=to_stat, aes(x=day, y=mean, ymin=lower, ymax=upper))+ #, col=factor(any_ecmo
  geom_point()+
  geom_errorbar(width=0)+
  xlab('Day in ICU')+
  ylab('Blood urea nitrogen')+
  scale_color_manual( values=c('darkorange1'))+ #'Any ECMO',,'cyan4'
  g.theme
tplot

##blood gases

fi_o2_data = blood_gases(
  before_var='ecmo_highest_fi_o2_6hr_before',
  after_var = 'fi_o2',
  etype="Venous-Venous")

peep_data = blood_gases(
  before_var='ecmo_highest_peep_6hr_before',
  after_var = 'eotd_peep',
  etype="Venous-Venous",
  logt=F) # from 99_functions.R

lplot =  ggplot(peep_data$to_plot, aes(x=xaxis, y=result, group=pin))+
  geom_line(col=grey(0.4))+
  geom_line(data=peep_data$av_line, aes(x=xaxis, y=result, group=pin), col='dark red')+ # add average
  scale_x_discrete(expand = c(0.05,0.05)) + # reduce white space
  theme_bw()+
  xlab('')+
  ylab('PEEP') #+
#scale_y_log10()
lplot



to_stat = combined_ecmo %>% 
  filter(ecmo_type == "Venous-Venous" 
         #,is.na(stroke_before_ECMO) | stroke_before_ECMO != "Yes"
         ) %>%
  group_by(pin) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  select(pin, anticoagulants, complication_stroke, stroke_before_ECMO) %>%
  mutate(complication_stroke = ifelse(is.na(complication_stroke), "Missing",
                                            complication_stroke))
          


tab_anticoag = with(to_stat, ctable(anticoagulants, complication_stroke,  chisq=FALSE, 
                                    prop='c', round.digits=0, dnn=c(' ','.')))
tab_anticoag

labels(to_stat) <- c(anticoagulants = "Anticoagulants",
       complication_stroke = "Stroke")
tab_ac <- tableby(complication_stroke ~ notest(anticoagulants), data=to_stat,
                cat.stats=c("Nmiss","countpct"),
                na.action=na.tableby(T),
                 test=FALSE, total=T)


tab_ac2<- freqlist(~includeNA(complication_stroke,"Missing") + anticoagulants , test=FALSE, data = to_stat)
summary(tab_ac2)


summary(tableby(~complication_other_value, data = ecmo_patients %>% filter(complication_stroke == "No")))



neuro_vars <-c("complication_stroke_date", "date_ecmo",
               "stroke_before_ecmo", "stroke_after_ecmo",
               "anticoagulants",
               "stroke_during_treatment_type_intraparenchymal_haemorrhage",
               "stroke_during_treatment_type_ischemic_stroke" ,
               "stroke_during_treatment_type_subarachnoid_haemorrhage",
               "stroke_during_treatment_type_hypoxic_ischemic_brain_injury"
)

factorvars <- c("stroke_before_ecmo", "stroke_after_ecmo",
                "anticoagulants",
                "stroke_during_treatment_type_intraparenchymal_haemorrhage",
                "stroke_during_treatment_type_ischemic_stroke" ,
                "stroke_during_treatment_type_subarachnoid_haemorrhage",
                "stroke_during_treatment_type_hypoxic_ischemic_brain_injury")

tabdata <- ecmo_patients %>%
  select(all_of(neuro_vars))

