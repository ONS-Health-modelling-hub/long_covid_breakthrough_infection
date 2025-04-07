library(haven)
library(sqldf)
library(tidyverse)
library(MatchIt)
library(stats)
library(PropCIs)
library(sandwich)
library(lmtest)
library(ggplot2)
library(splines)

out_dir = "filepath"
dataset_date = "20211122"
cutoff_date = "2021-11-20"
time_cis_offset = 7

########################################### DATA PREP ###########################################

### read in data
input_file <- paste0("filepath/data_participant_clean.dta")

vars_of_interest <- c(
  "participant_id",
  "hh_id_fake",
  "cis20_samp",
  "visit_id",
  "visit_date",
  "visit_status",
  "long_covid_have_symptoms",
  "reduce_activities_long_covid",
  "long_covid_fever",
  "long_covid_weakness_tiredness",
  "long_covid_diarrhoea",
  "long_covid_loss_of_smell",
  "long_covid_shortness_of_breath",
  "long_covid_vertigo_dizziness",
  "long_covid_trouble_sleeping",
  "long_covid_headache",
  "long_covid_nausea_vomiting",
  "long_covid_loss_of_appetite",
  "long_covid_sore_throat",
  "long_covid_chest_pain",
  "long_covid_worry_anxiety",
  "long_covid_memory_loss_confusion",
  "long_covid_muscle_ache",
  "long_covid_abdominal_pain",
  "long_covid_loss_of_taste",
  "long_covid_cough",
  "long_covid_palpitations",
  "long_covid_low_mood_not_enjoying",
  "long_covid_difficult_concentrate",
  "age_at_visit",
  "sex",
  "ethnicityg",
  "country",
  "gor9d",
  "imd_samp",
  "health_conditions",
  "health_conditions_impact",
  "smoke_ever_regularly",
  "smoke_now_cigarettes",
  "smoke_now_cigar",
  "smoke_now_pipe",
  "smoke_now_vape",
  "smoke_now_nothing",
  "work_sector",
  "work_status",
  "work_status_v1",
  "work_status_v2",
  "work_healthcare",
  "work_healthcare_v1",
  "work_socialcare",
  "work_care_nursing_home",
  "work_direct_contact_patients_etc",
  "result_mk",
  "result_combined",
  #"sympt_now_any",
  "sympt_now_fever",
  "sympt_now_muscle_ache_myalgia",
  "sympt_now_fatigue_weakness",
  "sympt_now_sore_throat",
  "sympt_now_cough",
  "sympt_now_shortness_of_breath",
  "sympt_now_headache",
  "sympt_now_nausea_vomiting",
  "sympt_now_abdominal_pain",
  "sympt_now_diarrhoea",
  "sympt_now_loss_of_taste",
  "sympt_now_loss_of_smell",
  "sympt_now_date",
  "covid_test_swab_result",
  "covid_test_swab_pos_first_date",
  "covid_test_swab_neg_last_date",
  "covid_test_blood",
  "covid_test_blood_result",
  "covid_test_blood_pos_first_date",
  "covid_test_blood_neg_last_date",
  "covid_think_havehad",
  "covid_date",
  "covid_nhs_contact",
  "covid_admitted",
  "sympt_covid_any",
  "sympt_covid_fever",
  "sympt_covid_muscle_ache_myalgia",
  "sympt_covid_fatigue_weakness",
  "sympt_covid_sore_throat",
  "sympt_covid_cough",
  "sympt_covid_shortness_of_breath",
  "sympt_covid_headache",
  "sympt_covid_nausea_vomiting",
  "sympt_covid_diarrhoea",
  "sympt_covid_loss_of_taste",
  "sympt_covid_loss_of_smell",
  "covid_vaccine_havehad",
  "covid_vaccine_date",
  "covid_vaccine_n_doses"
)

dat_all <- haven::read_dta(input_file, col_select=all_of(vars_of_interest))
dat_all <- zap_labels(dat_all)

### drop visits after cut-off date
dat_all$visit_date <- as.numeric(as.Date(dat_all$visit_date))
dat_all <- dat_all[dat_all$visit_date <= cutoff_date,]

### keep all visits belonging to participants who ever responded to LC question
part_select <- dat_all$participant_id[!is.na(dat_all$long_covid_have_symptoms) &
                                        dat_all$long_covid_have_symptoms %in% c(0, 1)]
part_select <- part_select[!duplicated(part_select)]
dat_all <- dat_all[dat_all$participant_id %in% part_select,]

### keep all visits belonging to participants who ever tested positive by PCR in CIS or outside of study
part_select <- dat_all$participant_id[dat_all$result_mk==1 | !is.na(dat_all$covid_test_swab_pos_first_date)]
part_select <- part_select[!duplicated(part_select)]
dat_all <- dat_all[dat_all$participant_id %in% part_select,]


### sort dataset by participant ID, visit date, and visit ID
dat_all <- dat_all[with(dat_all, order(participant_id, visit_date, -xtfrm(visit_id))),]


### vaccination date from vaccination data
dat_vacc <- haven::read_dta(paste0("filepath/data_participant_vaccination.dta"),
                            col_select = c("participant_id",
                                           "covid_vaccine_date1",
                                           "covid_vaccine_date2",
                                           "covid_vaccine_type1",
                                           "covid_vaccine_type2"))

dat_full <- merge(x = dat_all,
                  y = dat_vacc,
                  by.x = "participant_id",
                  by.y = "participant_id",
                  all.x = TRUE,
                  all.y = FALSE)

dat_full$covid_vaccine_date1 <- as.numeric(dat_full$covid_vaccine_date1)
dat_full$covid_vaccine_date2 <- as.numeric(dat_full$covid_vaccine_date2)

rm(dat_vacc); gc()
rm(dat_all); gc()

#################################### JOIN COVID STATUS VARIABLES ####################################

### Create first positive study PCR test date
dat_study_swab_pos1_date <- sqldf("
  select
    participant_id,
    min(visit_date) as study_swab_pos1_date
  from dat_full
  where result_mk=1
  group by participant_id
")

### merge first positive study PCR test date
dat_full <- merge(x = dat_full,
                 y = dat_study_swab_pos1_date,
                 by.x = "participant_id",
                 by.y = "participant_id",
                 all.x = TRUE,
                 all.y = FALSE)

rm(dat_study_swab_pos1_date); gc()

### Create first positive non-study PCR test date
dat_non_study_swab_pos1_date <- sqldf("
  select
    participant_id,
    min(covid_test_swab_pos_first_date) as non_study_swab_pos1_date
  from dat_full
  group by participant_id
")

### merge first positive study PCR test date
dat_full <- merge(x = dat_full,
                  y = dat_non_study_swab_pos1_date,
                  by.x = "participant_id",
                  by.y = "participant_id",
                  all.x = TRUE,
                  all.y = FALSE)

rm(dat_non_study_swab_pos1_date); gc()

### Create first positive study blood test date
dat_study_blood_pos1_date <- sqldf("
  select
    participant_id,
    min(visit_date) as study_blood_pos1_date
  from dat_full
  where result_combined=1
  group by participant_id
  ")

### merge first positive study blood test date
dat_full <- merge(x = dat_full,
                  y = dat_study_blood_pos1_date,
                  by.x = "participant_id",
                  by.y = "participant_id",
                  all.x = TRUE,
                  all.y = FALSE)

rm(dat_study_blood_pos1_date); gc()

### Create first positive non-study blood test date
dat_non_study_blood_pos1_date <- sqldf("
  select
    participant_id,
    min(covid_test_blood_pos_first_date) as non_study_blood_pos1_date
  from dat_full
  group by participant_id
  ")

### merge first positive study blood test date
dat_full <- merge(x = dat_full,
                  y = dat_non_study_blood_pos1_date,
                  by.x = "participant_id",
                  by.y = "participant_id",
                  all.x = TRUE,
                  all.y = FALSE)

rm(dat_non_study_blood_pos1_date); gc()

### Change anyone who tested positive for anti-bodies
### after getting vaccinated to negative

### Study
dat_full$result_combined <- ifelse(with(dat_full,
                                        !is.na(study_blood_pos1_date) &
                                          !is.na(covid_vaccine_date1) &
                                          study_blood_pos1_date >= covid_vaccine_date1),
                                   0, dat_full$result_combined)

### Non-study
dat_full$covid_test_blood_pos_first_date <- as.numeric(dat_full$covid_test_blood_pos_first_date)
dat_full$covid_test_blood_result <- ifelse(with(dat_full,
                                                !is.na(covid_test_blood_pos_first_date) &
                                                  !is.na(covid_vaccine_date1) &
                                                  covid_test_blood_pos_first_date >= covid_vaccine_date1),
                                           0, dat_full$covid_test_blood_result)

### Change study blood first positive date to NA if
### the positive test took place after vaccination
dat_full$study_blood_pos1_date <- ifelse(with(dat_full,
                                                  !is.na(study_blood_pos1_date) &
                                                    !is.na(covid_vaccine_date1) &
                                                    study_blood_pos1_date >= covid_vaccine_date1),
                                             NA, dat_full$study_blood_pos1_date)

### Change non-study blood first positive date to NA if
### the positive test took place after vaccination
dat_full$non_study_blood_pos1_date <- ifelse(with(dat_full,
                                                        !is.na(non_study_blood_pos1_date) &
                                                          !is.na(covid_vaccine_date1) &
                                                    non_study_blood_pos1_date >= covid_vaccine_date1),
                                                   NA, dat_full$non_study_blood_pos1_date)

#create earliest swab and blood dates from study and non-study tests
dat_full <- dat_full %>% mutate(first_pos_swab_date = pmin(study_swab_pos1_date, non_study_swab_pos1_date, na.rm = T),
                                first_pos_blood_date = pmin(study_blood_pos1_date, non_study_blood_pos1_date, na.rm = T))

### coerce necessary variables to date format
dat_full$first_pos_swab_date <- as.Date(dat_full$first_pos_swab_date, origin = "1970-01-01")
dat_full$first_pos_blood_date <- as.Date(dat_full$first_pos_blood_date, origin = "1970-01-01")
dat_full$covid_date <- as.Date(dat_full$covid_date, origin = "1970-01-01")
dat_full$visit_date <- as.Date(dat_full$visit_date, origin = "1970-01-01")
dat_full$covid_vaccine_date1 <- as.Date(dat_full$covid_vaccine_date1, origin = "1970-01-01")
dat_full$covid_vaccine_date2 <- as.Date(dat_full$covid_vaccine_date2, origin = "1970-01-01")

### COVID-19 arrived in the UK on 24 Jan 2020 - set any dates before this to NA
dat_full$first_pos_swab_date[dat_full$first_pos_swab_date < as.Date("2020-01-24")] <- NA
dat_full$first_pos_blood_date[dat_full$first_pos_blood_date < as.Date("2020-01-24")] <- NA
dat_full$covid_date[dat_full$covid_date < as.Date("2020-01-24")] <- NA

### Calculate time between first positive swab and first positive blood test
dat_full <- dat_full %>% mutate(time_since_pos_blood = as.numeric(first_pos_swab_date - first_pos_blood_date),
                              time_since_suspected_covid = as.numeric(first_pos_swab_date - covid_date))

#non_study_swab_dates <- dat_all %>% group_by(participant_id) %>% 
#  summarise(time_since_non_study_previous_positive_swab = max(time_since_non_study_previous_positive_swab, na.rm = T)) %>% 
#  mutate(time_since_non_study_previous_positive_swab = ifelse(time_since_non_study_previous_positive_swab=="-Inf", NA, time_since_non_study_previous_positive_swab))

###remove people who had a positive blood test >120 days before first positive swab test

part_select <- dat_full$participant_id[is.na(dat_full$time_since_pos_blood) | dat_full$time_since_pos_blood<120]
part_select <- part_select[!duplicated(part_select)]
dat_full <- dat_full[dat_full$participant_id %in% part_select,]

#think_covid_dates <- dat_all %>% group_by(participant_id) %>% 
#  summarise(time_since_suspected_covid = max(time_since_suspected_covid, na.rm = T)) %>% 
#  mutate(time_since_suspected_covid = ifelse(time_since_suspected_covid =="-Inf", NA, time_since_suspected_covid))

#part_select <- think_covid_dates$participant_id[is.na(think_covid_dates$time_since_suspected_covid) | think_covid_dates$time_since_suspected_covid<120]
#part_select <- part_select[!duplicated(part_select)]
#dat_all <- dat_all[dat_all$participant_id %in% part_select,]

dat_full <- dat_full %>% mutate(duration_after_infection = as.numeric(visit_date - first_pos_swab_date),
                              visit_28_after_infection = ifelse(duration_after_infection >= 28, 1, 0),
                              visit_84_after_infection = ifelse(duration_after_infection >= 84, 1, 0))

### aggregate to person level
dat_lc_status <- sqldf("
  select
    participant_id,
    sex,
    max(long_covid_have_symptoms) as lc_ever,
    min(first_pos_swab_date) as infection_date,
    covid_vaccine_date1 as covid_vaccine_date1,
    covid_vaccine_date2 as covid_vaccine_date2,
    max(covid_admitted) as covid_admitted_ever,
    imd_samp,
    ethnicityg,
    country
  from dat_full
  group by participant_id
")

dat_lc_date <- sqldf("
  select
    participant_id,
    min(visit_date) as lc_first_date
  from dat_full
  where long_covid_have_symptoms=1
  group by participant_id
")

dat_lc_q <- sqldf("
  select
    participant_id,
    max(visit_date) as lc_last_question_date
  from dat_full
  where long_covid_have_symptoms is not null
  group by participant_id
")

dat_lc_last_yes <- sqldf("
  select
    participant_id,
    max(visit_date) as lc_last_yes
  from dat_full
  where long_covid_have_symptoms = 1
  group by participant_id
")


dat_i_to_r_84 <- sqldf("
  select
    participant_id,
    visit_84_after_infection,
    min(visit_date) as first_response_84_date,
    age_at_visit as age_lc_84,
    long_covid_have_symptoms as lc_response_84,
    reduce_activities_long_covid as lc_activity_84,
    health_conditions as health_conditions_84
  from dat_full
  where visit_84_after_infection = 1 
  AND long_covid_have_symptoms is not null
  group by participant_id
")

dat_i_to_r_28 <- sqldf("
  select
    participant_id,
    visit_28_after_infection,
    min(visit_date) as first_response_28_date,
    age_at_visit as age_lc_28,
    long_covid_have_symptoms as lc_response_28,
    reduce_activities_long_covid as lc_activity_28,
    health_conditions as health_conditions_28
  from dat_full
  where visit_28_after_infection = 1 
  AND long_covid_have_symptoms is not null
  group by participant_id
")


dat_lc_status_date <- merge(x = dat_lc_status,
                         y = dat_lc_date,
                         by.x = "participant_id",
                         by.y = "participant_id",
                         all.x = TRUE,
                         all.y = FALSE)

dat_lc_status_date_q <- merge(x = dat_lc_status_date,
                            y = dat_lc_q,
                            by.x = "participant_id",
                            by.y = "participant_id",
                            all.x = TRUE,
                            all.y = FALSE)

dat_lc_status_date_q_yes <- merge(x = dat_lc_status_date_q,
                              y = dat_lc_last_yes,
                              by.x = "participant_id",
                              by.y = "participant_id",
                              all.x = TRUE,
                              all.y = FALSE)

dat_infected_lc <- merge(x = dat_lc_status_date_q_yes,
                       y = dat_i_to_r_84,
                       by.x = "participant_id",
                       by.y = "participant_id",
                       all.x = TRUE,
                       all.y = FALSE)

dat_infected_lc <- merge(x = dat_infected_lc,
                         y = dat_i_to_r_28,
                         by.x = "participant_id",
                         by.y = "participant_id",
                         all.x = TRUE,
                         all.y = FALSE)

dat_infected_lc <- dat_infected_lc %>% mutate(lc_onset_after_first_infection_date = ifelse(lc_first_date >= infection_date, 1, 0),
                                              diff_lc_first_date_infection_date = as.numeric(as.Date(lc_first_date, origin = "1970-01-01") - as.Date(infection_date, origin="1970-01-01")),
                                              diff_infection_lc_q_last_answered = as.numeric(as.Date(lc_last_question_date, origin = "1970-01-01") - as.Date(infection_date, origin="1970-01-01")),
                                              diff_infection_lc_q_last_yes = as.numeric(as.Date(lc_last_yes, origin = "1970-01-01") - as.Date(infection_date, origin="1970-01-01")))

### define IMD quintile
dat_infected_lc$imd_quintile_eng <- ifelse(dat_infected_lc$imd_samp>(0*32844/5) &
                                             dat_infected_lc$imd_samp<=(1*32844/5), 1,
                                   ifelse(dat_infected_lc$imd_samp>(1*32844/5) &
                                            dat_infected_lc$imd_samp<=(2*32844/5), 2,
                                          ifelse(dat_infected_lc$imd_samp>(2*32844/5) &
                                                   dat_infected_lc$imd_samp<=(3*32844/5), 3,
                                                 ifelse(dat_infected_lc$imd_samp>(3*32844/5) &
                                                          dat_infected_lc$imd_samp<=(4*32844/5), 4,
                                                        ifelse(dat_infected_lc$imd_samp>(4*32844/5) &
                                                                 dat_infected_lc$imd_samp<=(5*32844/5), 5, -999)))))

dat_infected_lc$imd_quintile_wal <- ifelse(dat_infected_lc$imd_samp>(0*1909/5) &
                                             dat_infected_lc$imd_samp<=(1*1909/5), 1,
                                   ifelse(dat_infected_lc$imd_samp>(1*1909/5) &
                                            dat_infected_lc$imd_samp<=(2*1909/5), 2,
                                          ifelse(dat_infected_lc$imd_samp>(2*1909/5) &
                                                   dat_infected_lc$imd_samp<=(3*1909/5), 3,
                                                 ifelse(dat_infected_lc$imd_samp>(3*1909/5) &
                                                          dat_infected_lc$imd_samp<=(4*1909/5), 4,
                                                        ifelse(dat_infected_lc$imd_samp>(4*1909/5) &
                                                                 dat_infected_lc$imd_samp<=(5*1909/5), 5, -999)))))

dat_infected_lc$imd_quintile_sco <- ifelse(dat_infected_lc$imd_samp>(0*6976/5) &
                                             dat_infected_lc$imd_samp<=(1*6976/5), 1,
                                   ifelse(dat_infected_lc$imd_samp>(1*6976/5) &
                                            dat_infected_lc$imd_samp<=(2*6976/5), 2,
                                          ifelse(dat_infected_lc$imd_samp>(2*6976/5) &
                                                   dat_infected_lc$imd_samp<=(3*6976/5), 3,
                                                 ifelse(dat_infected_lc$imd_samp>(3*6976/5) &
                                                          dat_infected_lc$imd_samp<=(4*6976/5), 4,
                                                        ifelse(dat_infected_lc$imd_samp>(4*6976/5) &
                                                                 dat_infected_lc$imd_samp<=(5*6976/5), 5, -999)))))

dat_infected_lc$imd_quintile_ni <- ifelse(dat_infected_lc$imd_samp>(0*890/5) &
                                            dat_infected_lc$imd_samp<=(1*890/5), 1,
                                  ifelse(dat_infected_lc$imd_samp>(1*890/5) &
                                           dat_infected_lc$imd_samp<=(2*890/5), 2,
                                         ifelse(dat_infected_lc$imd_samp>(2*890/5) &
                                                  dat_infected_lc$imd_samp<=(3*890/5), 3,
                                                ifelse(dat_infected_lc$imd_samp>(3*890/5) &
                                                         dat_infected_lc$imd_samp<=(4*890/5), 4,
                                                       ifelse(dat_infected_lc$imd_samp>(4*890/5) &
                                                                dat_infected_lc$imd_samp<=(5*890/5), 5, -999)))))

dat_infected_lc$imd_quintile <- ifelse(dat_infected_lc$country==0, dat_infected_lc$imd_quintile_eng,
                               ifelse(dat_infected_lc$country==1, dat_infected_lc$imd_quintile_wal,
                                      ifelse(dat_infected_lc$country==3, dat_infected_lc$imd_quintile_sco,
                                             ifelse(dat_infected_lc$country==2, dat_infected_lc$imd_quintile_ni, -999))))
### define white/non-white variable
dat_infected_lc$white <- ifelse(dat_infected_lc$ethnicityg==1, 1, 0)

##define activity limitation
dat_infected_lc$lc_activity_28_pooled <- ifelse(dat_infected_lc$lc_activity_28==1 | is.na(dat_infected_lc$lc_activity_28), 0, 1)
dat_infected_lc$lc_activity_84_pooled <- ifelse(dat_infected_lc$lc_activity_84==1 | is.na(dat_infected_lc$lc_activity_84), 0, 1)

#define month/year of infection
dat_infected_lc$infection_month <- format(as.Date(dat_infected_lc$infection_date), "%Y-%m")

##derive vaccinated before infection flags and remove people never vaccinated
dat_infected_lc <- dat_infected_lc %>% mutate(infection_date = as.Date(infection_date, origin = "1970-01-01"),
                                   first_response_84_date = as.Date(first_response_84_date, origin = "1970-01-01" ),
                                   first_response_28_date = as.Date(first_response_28_date, origin = "1970-01-01" )) %>%
                            filter(!is.na(covid_vaccine_date1)) %>%
                            filter(is.na(lc_onset_after_first_infection_date) | lc_onset_after_first_infection_date == 1)
                            #filter(age_at_infection >= 18 & age_at_infection < 70)
dat_infected_lc <- dat_infected_lc %>%  mutate(diff_betweeen_vacc2_infection = as.Date(infection_date, origin = "1970-01-01") - covid_vaccine_date2,
                                   diff_betweeen_vacc1_infection = as.Date(infection_date, origin = "1970-01-01") - covid_vaccine_date1,
                                   vacc2_same_day_first_inf = ifelse(covid_vaccine_date2 == infection_date, 1, 0),
                              duration_between_inf_response_84 = as.numeric(first_response_84_date - infection_date),
                              duration_between_inf_response_28 = as.numeric(first_response_28_date - infection_date),
                              first_vacc_between_infection_lc_response_84 = ifelse(covid_vaccine_date1 >= infection_date & covid_vaccine_date1 <= first_response_84_date, 1, 0),
                              first_vacc_between_infection_lc_response_28 = ifelse(covid_vaccine_date1 >= infection_date & covid_vaccine_date1 <= first_response_28_date, 1, 0),
                              second_vacc_between_infection_lc_response_84 = ifelse(covid_vaccine_date2 >= infection_date & covid_vaccine_date2 <= first_response_84_date, 1, 0),
                              second_vacc_between_infection_lc_response_28 = ifelse(covid_vaccine_date2 >= infection_date & covid_vaccine_date2 <= first_response_28_date, 1, 0),
                              vacc1_1_days_before_infection = ifelse(diff_betweeen_vacc1_infection >=1, 1, 0),
                              vacc1_14_days_before_infection = ifelse(diff_betweeen_vacc1_infection >=14, 1, 0),
                              vacc1_21_days_before_infection = ifelse(diff_betweeen_vacc1_infection >=21, 1, 0),
                              vacc1_28_days_before_infection = ifelse(diff_betweeen_vacc1_infection >=28, 1, 0),
                              vacc1_1_to_13_days_before_infection = ifelse(vacc1_1_days_before_infection == 1 & vacc1_14_days_before_infection == 0, 1, 0),
                              vacc1_1_to_20_days_before_infection = ifelse(vacc1_1_days_before_infection == 1 & vacc1_21_days_before_infection == 0, 1, 0),
                              vacc1_1_to_27_days_before_infection = ifelse(vacc1_1_days_before_infection == 1 & vacc1_28_days_before_infection == 0, 1, 0),
                              vacc2_1_days_before_infection = ifelse(diff_betweeen_vacc2_infection >=1, 1, 0),
                              vacc2_14_days_before_infection = ifelse(diff_betweeen_vacc2_infection >=14, 1, 0),
                              vacc2_21_days_before_infection = ifelse(diff_betweeen_vacc2_infection >=21, 1, 0),
                              vacc2_28_days_before_infection = ifelse(diff_betweeen_vacc2_infection >=28, 1, 0),
                              vacc2_1_to_20_days_before_infection = ifelse(vacc2_1_days_before_infection == 1 & vacc2_21_days_before_infection == 0, 1, 0),
                              vacc2_1_to_13_days_before_infection = ifelse(vacc2_1_days_before_infection == 1 & vacc2_14_days_before_infection == 0, 1, 0),
                              vacc2_1_to_27_days_before_infection = ifelse(vacc2_1_days_before_infection == 1 & vacc2_28_days_before_infection == 0, 1, 0),
                              lc_fu_time_28 = ifelse(diff_infection_lc_q_last_answered>=28, 1, 0),
                              lc_fu_time_84 = ifelse(diff_infection_lc_q_last_answered>=84, 1, 0),
                              lc_yes_at_28_or_more = ifelse(diff_infection_lc_q_last_yes>=28, 1, 0),
                              lc_yes_at_84_or_more = ifelse(diff_infection_lc_q_last_yes>=84, 1, 0))
                              
###single dose###

dat_matching_84_3weeks <- dat_infected_lc %>% filter(age_lc_84 >= 18 & age_lc_84 < 70) %>%
  filter(visit_84_after_infection == 1) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
  filter(first_vacc_between_infection_lc_response_84 == 0) %>%                                 
  filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = ifelse(vacc1_21_days_before_infection == 1, 1, 0),
         duration_between_inf_response_84_grouped_7 = cut(duration_between_inf_response_84, seq(84, 238, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_7 = ifelse(is.na(duration_between_inf_response_84_grouped_7), ">=238", duration_between_inf_response_84_grouped_7),
         duration_between_inf_response_84_grouped_14 = cut(duration_between_inf_response_84, seq(84, 238, 14), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_14 = ifelse(is.na(duration_between_inf_response_84_grouped_14), ">=238", duration_between_inf_response_84_grouped_14),
         health_conditions_84 = ifelse(is.na(health_conditions_84), 0, health_conditions_84),
         covid_admitted_ever = ifelse(is.na(covid_admitted_ever), 0, covid_admitted_ever),
         age_group_10_year = case_when(age_lc_84 >=18 & age_lc_84 <= 29 ~ "18 to 29",
                                       age_lc_84 >=30 & age_lc_84 < 40 ~ "30 to 39",
                                       age_lc_84 >=40 & age_lc_84 < 50 ~ "40 to 49",
                                       age_lc_84 >=50 & age_lc_84 < 70 ~ "50 to 69")) %>%
  select(participant_id, case, age_lc_84, age_group_10_year, covid_admitted_ever, health_conditions_84, imd_quintile, ethnicityg, white, sex, infection_date, infection_month, first_response_84_date, duration_between_inf_response_84, duration_between_inf_response_84_grouped_7, duration_between_inf_response_84_grouped_14, covid_vaccine_date1, covid_vaccine_date2, lc_response_84, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_28_days_before_infection, vacc2_28_days_before_infection, first_vacc_between_infection_lc_response_84, vacc1_1_to_27_days_before_infection, vacc2_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_84)


ggplot(dat_matching_84_3weeks, aes(x = age_lc_84, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Age (years)") +
  scale_fill_discrete("Case")

ggplot(dat_matching_84_3weeks, aes(x = duration_between_inf_response_84, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Days between infection & response") +
  scale_fill_discrete("Case")


##matching
matched_84_3weeks_age_sex <- matchit(case ~ age_group_10_year + sex + duration_between_inf_response_84_grouped_7, 
                                     data = dat_matching_84_3weeks,
                                     method = "nearest",
                                     exact = ~ age_group_10_year + sex + duration_between_inf_response_84_grouped_7) 
summary(matched_84_3weeks_age_sex)
matched_84_data_3weeks_age_sex <- match.data(matched_84_3weeks_age_sex)


dat_matching_28_3weeks <- dat_infected_lc %>% filter(age_lc_28 >= 18 & age_lc_28 < 70) %>%
  filter(visit_28_after_infection == 1) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
  filter(first_vacc_between_infection_lc_response_28 == 0) %>%
  filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = ifelse(vacc1_21_days_before_infection == 1, 1, 0),
         duration_between_inf_response_28_grouped_7 = cut(duration_between_inf_response_28, seq(28, 238, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_28_grouped_7 = ifelse(is.na(duration_between_inf_response_28_grouped_7), ">=238", duration_between_inf_response_28_grouped_7),
         health_conditions_28 = ifelse(is.na(health_conditions_28), 0, health_conditions_28),
         covid_admitted_ever = ifelse(is.na(covid_admitted_ever), 0, covid_admitted_ever),
         age_group_10_year = case_when(age_lc_28 >=18 & age_lc_28 <= 29 ~ "18 to 29",
                                       age_lc_28 >=30 & age_lc_28 < 40 ~ "30 to 39",
                                       age_lc_28 >=40 & age_lc_28 < 50 ~ "40 to 49",
                                       age_lc_28 >=50 & age_lc_28 < 70 ~ "50 to 69")) %>%
  select(participant_id, case, age_lc_28, age_group_10_year, covid_admitted_ever, health_conditions_28, imd_quintile, ethnicityg, white,  sex, infection_date, first_response_28_date, duration_between_inf_response_28, duration_between_inf_response_28_grouped_7, covid_vaccine_date1, covid_vaccine_date2, lc_response_28, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_28_days_before_infection, vacc2_28_days_before_infection, first_vacc_between_infection_lc_response_28, vacc1_1_to_27_days_before_infection, vacc2_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_28)


ggplot(dat_matching_28_3weeks, aes(x = age_lc_28, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Age (years)") +
  scale_fill_discrete("Case")

ggplot(dat_matching_28_3weeks, aes(x = duration_between_inf_response_28, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Days between infection & response") +
  scale_fill_discrete("Case")

##matching

matched_28_3weeks_age_sex <- matchit(case ~ age_group_10_year + sex + duration_between_inf_response_28_grouped_7, 
                                     data = dat_matching_28_3weeks,
                                     method = "nearest",
                                     exact = ~ age_group_10_year + sex + duration_between_inf_response_28_grouped_7) 
summary(matched_28_3weeks_age_sex)
matched_28_data_3weeks_age_sex <- match.data(matched_28_3weeks_age_sex)

##calculate proportions
proportions_28 <- matched_28_data_3weeks_age_sex %>% group_by(case, lc_response_28) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
case_28_counts <- proportions_28 %>% filter(case == 1) %>% select(n)
control_28_counts <- proportions_28 %>% filter(case == 0) %>% select(n)

case_28_ci_lc_yes <- scoreci(x = case_28_counts[2,2], n = sum(case_28_counts[,2]), conf.level = 0.95)
case_28_ci_lc_no <- scoreci(x = case_28_counts[1,2], n = sum(case_28_counts[,2]), conf.level = 0.95)

case_28_confints_lc_yes <- data.frame(case_28_ci_lc_yes$conf.int)
colnames(case_28_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
case_28_confints_lc_yes <- case_28_confints_lc_yes %>% mutate(case = 1, lc_response_28 = 1)

case_28_confints_lc_no <- data.frame(case_28_ci_lc_no$conf.int)
colnames(case_28_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
case_28_confints_lc_no <- case_28_confints_lc_no %>% mutate(case = 1, lc_response_28 = 0)

control_28_ci_lc_yes <- scoreci(x = control_28_counts[2,2], n = sum(control_28_counts[,2]), conf.level = 0.95)
control_28_ci_lc_no <- scoreci(x = control_28_counts[1,2], n = sum(control_28_counts[,2]), conf.level = 0.95)

control_28_confints_lc_yes <- data.frame(control_28_ci_lc_yes$conf.int)
colnames(control_28_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
control_28_confints_lc_yes <- control_28_confints_lc_yes %>% mutate(case = 0, lc_response_28 = 1)

control_28_confints_lc_no <- data.frame(control_28_ci_lc_no$conf.int)
colnames(control_28_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
control_28_confints_lc_no <- control_28_confints_lc_no %>% mutate(case = 0, lc_response_28 = 0)

case_28_confints <- rbind(case_28_confints_lc_no, case_28_confints_lc_yes, control_28_confints_lc_yes, control_28_confints_lc_no)

results_28_single_dose <- left_join(proportions_28, case_28_confints, by = c("case", "lc_response_28"))
write.csv(results_28_single_dose, paste0(out_dir, "\\single_dose_proportions_28.csv"), row.names = F)

##unadjusted 
mod_28_unadj_1_dose <- glm(lc_response_28 ~ case, data = matched_28_data_3weeks_age_sex, family = binomial)
vcov_28_unadj_1_dose <- vcovCL(mod_28_unadj_1_dose, cluster = ~subclass)
output_28_unadj_1_dose <- coeftest(mod_28_unadj_1_dose, vcov = vcov_28_unadj_1_dose)

unadjust_28_results_1_dose <- data.frame(model = "unadjusted",
                                  OR = exp(output_28_unadj_1_dose[2,1]),
                                  lower_ci = exp(output_28_unadj_1_dose[2,1] - 1.96*output_28_unadj_1_dose[2,2]),
                                  upper_ci = exp(output_28_unadj_1_dose[2,1] + 1.96*output_28_unadj_1_dose[2,2]),
                                  p = output_28_unadj_1_dose[2,4],
                                  row.names = NULL,
                                  stringsAsFactors = FALSE)  

##adjusted
mod_28_adj_1_dose <- glm(lc_response_28 ~ case 
                  + ns(age_lc_28, df=2, Boundary.knots=quantile(age_lc_28, c(.10, .90)))
                  + sex 
                  + health_conditions_28 
                  + covid_admitted_ever 
                  + ns(duration_between_inf_response_28, df=2, Boundary.knots=quantile(duration_between_inf_response_28, c(.10, .90))) 
                  + imd_quintile 
                  + white, 
                  data = matched_28_data_3weeks_age_sex, family = binomial)
vcov_28_adj_1_dose <- vcovCL(mod_28_adj_1_dose, cluster = ~subclass)
output_28_adj_1_dose <- coeftest(mod_28_adj_1_dose, vcov = vcov_28_adj_1_dose)

adjust_28_results_1_dose <- data.frame(model = "adjusted",
                                  OR = exp(output_28_adj_1_dose[2,1]),
                                  lower_ci = exp(output_28_adj_1_dose[2,1] - 1.96*output_28_adj_1_dose[2,2]),
                                  upper_ci = exp(output_28_adj_1_dose[2,1] + 1.96*output_28_adj_1_dose[2,2]),
                                  p = output_28_adj_1_dose[2,4],
                                  row.names = NULL,
                                  stringsAsFactors = FALSE)  

write.csv(rbind(unadjust_28_results_1_dose, adjust_28_results_1_dose), paste0(out_dir, "\\single_dose_OR_28.csv"), row.names = F) 

##calculate proportions
proportions_84 <- matched_84_data_3weeks_age_sex %>% group_by(case, lc_response_84) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
case_84_counts <- proportions_84 %>% filter(case == 1) %>% select(n)
control_84_counts <- proportions_84 %>% filter(case == 0) %>% select(n)

case_84_ci_lc_yes <- scoreci(x = case_84_counts[2,2], n = sum(case_84_counts[,2]), conf.level = 0.95)
case_84_ci_lc_no <- scoreci(x = case_84_counts[1,2], n = sum(case_84_counts[,2]), conf.level = 0.95)


case_84_confints_lc_yes <- data.frame(case_84_ci_lc_yes$conf.int)
colnames(case_84_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
case_84_confints_lc_yes <- case_84_confints_lc_yes %>% mutate(case = 1, lc_response_84 = 1)

case_84_confints_lc_no <- data.frame(case_84_ci_lc_no$conf.int)
colnames(case_84_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
case_84_confints_lc_no <- case_84_confints_lc_no %>% mutate(case = 1, lc_response_84 = 0)


control_84_ci_lc_yes <- scoreci(x = control_84_counts[2,2], n = sum(control_84_counts[,2]), conf.level = 0.95)
control_84_ci_lc_no <- scoreci(x = control_84_counts[1,2], n = sum(control_84_counts[,2]), conf.level = 0.95)


control_84_confints_lc_yes <- data.frame(control_84_ci_lc_yes$conf.int)
colnames(control_84_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
control_84_confints_lc_yes <- control_84_confints_lc_yes %>% mutate(case = 0, lc_response_84 = 1)

control_84_confints_lc_no <- data.frame(control_84_ci_lc_no$conf.int)
colnames(control_84_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
control_84_confints_lc_no <- control_84_confints_lc_no %>% mutate(case = 0, lc_response_84 = 0)

case_84_confints <- rbind(case_84_confints_lc_no, case_84_confints_lc_yes, control_84_confints_lc_yes, control_84_confints_lc_no)

results_84_single_dose <- left_join(proportions_84, case_84_confints, by = c("case", "lc_response_84"))
write.csv(results_84_single_dose, paste0(out_dir, "\\single_dose_proportions_84.csv"), row.names = F)

##unadjusted
mod_84_unadj_1_dose <- glm(lc_response_84 ~ case, data = matched_84_data_3weeks_age_sex, family = binomial)
vcov_84_unadj_1_dose <- vcovCL(mod_84_unadj_1_dose, cluster = ~subclass)
output_84_unadj_1_dose <- coeftest(mod_84_unadj_1_dose, vcov = vcov_84_unadj_1_dose)

unadjust_84_results_1_dose <- data.frame(model = "unadjusted",
                                  OR = exp(output_84_unadj_1_dose[2,1]),
                                  lower_ci = exp(output_84_unadj_1_dose[2,1] - 1.96*output_84_unadj_1_dose[2,2]),
                                  upper_ci = exp(output_84_unadj_1_dose[2,1] + 1.96*output_84_unadj_1_dose[2,2]),
                                  p = output_84_unadj_1_dose[2,4],
                                  row.names = NULL,
                                  stringsAsFactors = FALSE)  

##adjusted
mod_84_adj_1_dose <- glm(lc_response_84 ~ case 
                  + ns(age_lc_84, df=2, Boundary.knots=quantile(age_lc_84, c(.10, .90)))
                  + sex 
                  + health_conditions_84 
                  + covid_admitted_ever 
                  + ns(duration_between_inf_response_84, df=2, Boundary.knots=quantile(duration_between_inf_response_84, c(.10, .90)))
                  + imd_quintile 
                  + white, 
                  data = matched_84_data_3weeks_age_sex, family = binomial)
vcov_84_adj_1_dose <- vcovCL(mod_84_adj_1_dose, cluster = ~subclass)
output_84_adj_1_dose <- coeftest(mod_84_adj_1_dose, vcov = vcov_84_adj_1_dose)

adjust_84_results_1_dose <- data.frame(model = "adjusted",
                                OR = exp(output_84_adj_1_dose[2,1]),
                                lower_ci = exp(output_84_adj_1_dose[2,1] - 1.96*output_84_adj_1_dose[2,2]),
                                upper_ci = exp(output_84_adj_1_dose[2,1] + 1.96*output_84_adj_1_dose[2,2]),
                                p = output_84_adj_1_dose[2,4],
                                row.names = NULL,
                                stringsAsFactors = FALSE)  

write.csv(rbind(unadjust_84_results_1_dose, adjust_84_results_1_dose), paste0(out_dir, "\\single_dose_OR_84.csv"), row.names = F) 


###subgroup analysis excluding cases who received their second dose between infection and response
dat_matching_84_3weeks_subgroup <- dat_infected_lc %>% filter(visit_84_after_infection == 1) %>%
  filter(age_lc_84 >= 18 & age_lc_84 < 70) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
  filter(first_vacc_between_infection_lc_response_84 == 0) %>%
  filter(is.na(second_vacc_between_infection_lc_response_84) | second_vacc_between_infection_lc_response_84 == 0) %>%
  filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = ifelse(vacc1_21_days_before_infection == 1, 1, 0),
         duration_between_inf_response_84_grouped_7 = cut(duration_between_inf_response_84, seq(84, 238, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_7 = ifelse(is.na(duration_between_inf_response_84_grouped_7), ">=238", duration_between_inf_response_84_grouped_7),
         health_conditions_84 = ifelse(is.na(health_conditions_84), 0, health_conditions_84),
         covid_admitted_ever = ifelse(is.na(covid_admitted_ever), 0, covid_admitted_ever)) %>%
  select(participant_id, case, age_lc_84, imd_quintile, ethnicityg, white, health_conditions_84, covid_admitted_ever, sex, infection_date, first_response_84_date, duration_between_inf_response_84, duration_between_inf_response_84_grouped_7, covid_vaccine_date1, covid_vaccine_date2, lc_response_84, vacc1_1_days_before_infection, vacc1_21_days_before_infection, first_vacc_between_infection_lc_response_84, second_vacc_between_infection_lc_response_84)


dat_matching_28_3weeks_subgroup <- dat_infected_lc %>% filter(visit_28_after_infection == 1) %>%
  filter(age_lc_28 >= 18 & age_lc_28 < 70) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
  filter(first_vacc_between_infection_lc_response_28 == 0) %>%
  filter(is.na(second_vacc_between_infection_lc_response_28) | second_vacc_between_infection_lc_response_28 == 0) %>%
  filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = ifelse(vacc1_21_days_before_infection == 1, 1, 0),
         duration_between_inf_response_28_grouped_7 = cut(duration_between_inf_response_28, seq(28, 238, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_28_grouped_7 = ifelse(is.na(duration_between_inf_response_28_grouped_7), ">=238", duration_between_inf_response_28_grouped_7),
         health_conditions_28 = ifelse(is.na(health_conditions_28), 0, health_conditions_28),
         covid_admitted_ever = ifelse(is.na(covid_admitted_ever), 0, covid_admitted_ever),
         age_group_10_year = case_when(age_lc_28 >=18 & age_lc_28 <= 29 ~ "18 to 29",
                                       age_lc_28 >=30 & age_lc_28 < 40 ~ "30 to 39",
                                       age_lc_28 >=40 & age_lc_28 < 50 ~ "40 to 49",
                                       age_lc_28 >=50 & age_lc_28 < 70 ~ "50 to 69")) %>%
  select(participant_id, case, age_lc_28, age_group_10_year, imd_quintile, ethnicityg, white, health_conditions_28, covid_admitted_ever, sex, infection_date, first_response_28_date, duration_between_inf_response_28, duration_between_inf_response_28_grouped_7, covid_vaccine_date1, covid_vaccine_date2, lc_response_28, vacc1_1_days_before_infection, vacc1_21_days_before_infection, first_vacc_between_infection_lc_response_28, vacc1_1_to_20_days_before_infection, second_vacc_between_infection_lc_response_28)

##matching
matched_28_3weeks_age_sex_subgroup <- matchit(case ~ age_group_10_year + sex + duration_between_inf_response_28_grouped_7, 
                                     data = dat_matching_28_3weeks_subgroup,
                                     method = "nearest",
                                     exact = ~ age_group_10_year + sex + duration_between_inf_response_28_grouped_7) 
summary(matched_28_3weeks_age_sex_subgroup)
matched_28_data_3weeks_age_sex_subgroup <- match.data(matched_28_3weeks_age_sex_subgroup)


##calculate proportions
proportions_28 <- matched_28_data_3weeks_age_sex_subgroup %>% group_by(case, lc_response_28) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
case_28_counts <- proportions_28 %>% filter(case == 1) %>% select(n)
control_28_counts <- proportions_28 %>% filter(case == 0) %>% select(n)

case_28_ci_lc_yes <- scoreci(x = case_28_counts[2,2], n = sum(case_28_counts[,2]), conf.level = 0.95)
case_28_ci_lc_no <- scoreci(x = case_28_counts[1,2], n = sum(case_28_counts[,2]), conf.level = 0.95)

case_28_confints_lc_yes <- data.frame(case_28_ci_lc_yes$conf.int)
colnames(case_28_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
case_28_confints_lc_yes <- case_28_confints_lc_yes %>% mutate(case = 1, lc_response_28 = 1)

case_28_confints_lc_no <- data.frame(case_28_ci_lc_no$conf.int)
colnames(case_28_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
case_28_confints_lc_no <- case_28_confints_lc_no %>% mutate(case = 1, lc_response_28 = 0)

control_28_ci_lc_yes <- scoreci(x = control_28_counts[2,2], n = sum(control_28_counts[,2]), conf.level = 0.95)
control_28_ci_lc_no <- scoreci(x = control_28_counts[1,2], n = sum(control_28_counts[,2]), conf.level = 0.95)

control_28_confints_lc_yes <- data.frame(control_28_ci_lc_yes$conf.int)
colnames(control_28_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
control_28_confints_lc_yes <- control_28_confints_lc_yes %>% mutate(case = 0, lc_response_28 = 1)

control_28_confints_lc_no <- data.frame(control_28_ci_lc_no$conf.int)
colnames(control_28_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
control_28_confints_lc_no <- control_28_confints_lc_no %>% mutate(case = 0, lc_response_28 = 0)

case_28_confints <- rbind(case_28_confints_lc_no, case_28_confints_lc_yes, control_28_confints_lc_yes, control_28_confints_lc_no)

results_28_single_dose_subgroup <- left_join(proportions_28, case_28_confints, by = c("case", "lc_response_28"))
write.csv(results_28_single_dose_subgroup, paste0(out_dir, "\\single_dose_proportions_28_subgroup.csv"), row.names = F)

##undadjusted 
mod_28_unadj_1_dose_sub <- glm(lc_response_28 ~ case, data = matched_28_data_3weeks_age_sex_subgroup, family = binomial)
vcov_28_unadj_1_dose_sub <- vcovCL(mod_28_unadj_1_dose_sub, cluster = ~subclass)
output_28_unadj_1_dose_sub <- coeftest(mod_28_unadj_1_dose_sub, vcov = vcov_28_unadj_1_dose_sub)

unadjust_28_results_1_dose_sub <- data.frame(model = "unadjusted",
                                         OR = exp(output_28_unadj_1_dose_sub[2,1]),
                                         lower_ci = exp(output_28_unadj_1_dose_sub[2,1] - 1.96*output_28_unadj_1_dose_sub[2,2]),
                                         upper_ci = exp(output_28_unadj_1_dose_sub[2,1] + 1.96*output_28_unadj_1_dose_sub[2,2]),
                                         p = output_28_unadj_1_dose_sub[2,4],
                                         row.names = NULL,
                                         stringsAsFactors = FALSE) 

##adjusted
mod_28_adj_1_dose_sub <- glm(lc_response_28 ~ case 
                             + ns(age_lc_28, df=2, Boundary.knots=quantile(age_lc_28, c(.10, .90)))
                             + sex 
                             + health_conditions_28 
                             + covid_admitted_ever 
                             + ns(duration_between_inf_response_28, df=2, Boundary.knots=quantile(duration_between_inf_response_28, c(.10, .90))) 
                             + imd_quintile 
                             + white, 
                             data = matched_28_data_3weeks_age_sex_subgroup, family = binomial)
vcov_28_adj_1_dose_sub <- vcovCL(mod_28_adj_1_dose_sub, cluster = ~subclass)
output_28_adj_1_dose_sub <- coeftest(mod_28_adj_1_dose_sub, vcov = vcov_28_adj_1_dose_sub)

adjust_28_results_1_dose_sub <- data.frame(model = "adjusted",
                                             OR = exp(output_28_adj_1_dose_sub[2,1]),
                                             lower_ci = exp(output_28_adj_1_dose_sub[2,1] - 1.96*output_28_adj_1_dose_sub[2,2]),
                                             upper_ci = exp(output_28_adj_1_dose_sub[2,1] + 1.96*output_28_adj_1_dose_sub[2,2]),
                                             p = output_28_adj_1_dose_sub[2,4],
                                             row.names = NULL,
                                             stringsAsFactors = FALSE)

write.csv(rbind(unadjust_28_results_1_dose_sub, adjust_28_results_1_dose_sub), paste0(out_dir, "\\single_dose_OR_28_subgroup.csv"), row.names = F) 


###fully vaccinated###
dat_matching_84_2weeks_2nd_dose <- dat_infected_lc %>% filter(visit_84_after_infection == 1) %>%
  filter(age_lc_84 >= 18 & age_lc_84 < 70) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(first_vacc_between_infection_lc_response_84 == 0) %>%
  filter(second_vacc_between_infection_lc_response_84 == 0 | is.na(second_vacc_between_infection_lc_response_84)) %>%                                 
  filter(vacc2_1_to_13_days_before_infection == 0 | is.na(vacc2_1_to_13_days_before_infection)) %>%
  #filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = case_when(vacc2_14_days_before_infection == 1 ~ 1, 
                          vacc2_14_days_before_infection == 0 | is.na(vacc2_14_days_before_infection) ~ 0),
         duration_between_inf_response_84_grouped_7 = cut(duration_between_inf_response_84, seq(84, 238, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_7 = ifelse(is.na(duration_between_inf_response_84_grouped_7), ">=238", duration_between_inf_response_84_grouped_7),
         duration_between_inf_response_84_grouped_14 = cut(duration_between_inf_response_84, seq(84, 238, 14), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_14 = ifelse(is.na(duration_between_inf_response_84_grouped_14), ">=238", duration_between_inf_response_84_grouped_14),
         health_conditions_84 = ifelse(is.na(health_conditions_84), 0, health_conditions_84),
         covid_admitted_ever = ifelse(is.na(covid_admitted_ever), 0, covid_admitted_ever),
         age_group_10_year = case_when(age_lc_84 >=18 & age_lc_84 <= 29 ~ "18 to 29",
                                       age_lc_84 >=30 & age_lc_84 < 40 ~ "30 to 39",
                                       age_lc_84 >=40 & age_lc_84 < 50 ~ "40 to 49",
                                       age_lc_84 >=50 & age_lc_84 < 70 ~ "50 to 69"),
         age_40_plus = ifelse(age_lc_84>=40, 1, 0)) %>%
  filter(case == 1 | (case == 0 & vacc1_1_days_before_infection == 0)) %>%
  select(participant_id, case, age_lc_84, age_group_10_year, age_40_plus, imd_quintile, ethnicityg, white, health_conditions_84, covid_admitted_ever, sex, infection_date, first_response_84_date, duration_between_inf_response_84, duration_between_inf_response_84_grouped_7, duration_between_inf_response_84_grouped_14, covid_vaccine_date1, covid_vaccine_date2, lc_response_84, lc_activity_84_pooled, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_21_days_before_infection, vacc2_21_days_before_infection, first_vacc_between_infection_lc_response_84, vacc1_1_to_20_days_before_infection, vacc2_1_to_20_days_before_infection, second_vacc_between_infection_lc_response_84)


dat_matching_28_2weeks_2nd_dose <- dat_infected_lc %>% filter(visit_28_after_infection == 1) %>%
  filter(age_lc_28 >= 18 & age_lc_28 < 70) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(first_vacc_between_infection_lc_response_28 == 0) %>%
  filter(second_vacc_between_infection_lc_response_28 == 0 | is.na(second_vacc_between_infection_lc_response_28)) %>%                                 
  filter(vacc2_1_to_13_days_before_infection == 0 | is.na(vacc2_1_to_13_days_before_infection)) %>%
  #filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = case_when(vacc2_14_days_before_infection == 1 ~ 1, 
                          vacc2_14_days_before_infection == 0 | is.na(vacc2_14_days_before_infection) ~ 0),
         duration_between_inf_response_28_grouped_7 = cut(duration_between_inf_response_28, seq(28, 238, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_28_grouped_7 = ifelse(is.na(duration_between_inf_response_28_grouped_7), ">=238", duration_between_inf_response_28_grouped_7),
         duration_between_inf_response_28_grouped_14 = cut(duration_between_inf_response_28, seq(84, 238, 14), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_28_grouped_14 = ifelse(is.na(duration_between_inf_response_28_grouped_14), ">=238", duration_between_inf_response_28_grouped_14),
         health_conditions_28 = ifelse(is.na(health_conditions_28), 0, health_conditions_28),
         covid_admitted_ever = ifelse(is.na(covid_admitted_ever), 0, covid_admitted_ever),
         age_group_10_year = case_when(age_lc_28 >=18 & age_lc_28 <= 29 ~ "18 to 29",
                                       age_lc_28 >=30 & age_lc_28 < 40 ~ "30 to 39",
                                       age_lc_28 >=40 & age_lc_28 < 50 ~ "40 to 49",
                                       age_lc_28 >=50 & age_lc_28 < 70 ~ "50 to 69"),
         age_40_plus = ifelse(age_lc_28>=40, 1, 0)) %>%
  filter(case == 1 | (case == 0 & vacc1_1_days_before_infection == 0)) %>%
  select(participant_id, case, age_lc_28, age_group_10_year, age_40_plus, imd_quintile, ethnicityg, white, health_conditions_28, covid_admitted_ever, sex, infection_date, first_response_28_date, duration_between_inf_response_28, duration_between_inf_response_28_grouped_7, duration_between_inf_response_28_grouped_14, covid_vaccine_date1, covid_vaccine_date2, lc_response_28, lc_activity_28_pooled, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_21_days_before_infection, vacc2_21_days_before_infection, first_vacc_between_infection_lc_response_28, vacc1_1_to_20_days_before_infection, vacc2_1_to_20_days_before_infection, second_vacc_between_infection_lc_response_28)


##matching
matched_28_2weeks_age_sex_2nd_dose <- matchit(case ~ ns(age_lc_28, df=2, Boundary.knots=quantile(age_lc_28, c(.10, .90)))
                                              + sex 
                                              + health_conditions_28 
                                              + covid_admitted_ever 
                                              duration_between_inf_response_28_grouped_7 
                                              + imd_quintile 
                                              + white, 
                                              data = dat_matching_28_2weeks_2nd_dose,
                                              method = "nearest",
                                              exact = ~ duration_between_inf_response_28_grouped_7) 
summary(matched_28_2weeks_age_sex_2nd_dose)
matched_28_data_2weeks_age_sex_2nd_dose <- match.data(matched_28_2weeks_age_sex_2nd_dose)

matched_84_2weeks_age_sex_2nd_dose <- matchit(case ~ ns(age_lc_84, df=2, Boundary.knots=quantile(age_lc_84, c(.10, .90)))
                                              + sex 
                                              + health_conditions_84 
                                              + covid_admitted_ever 
                                              + duration_between_inf_response_84_grouped_7 
                                              + imd_quintile 
                                              + white, 
                                              data = dat_matching_84_2weeks_2nd_dose,
                                              method = "nearest",
                                              exact = ~ duration_between_inf_response_84_grouped_7) 
summary(matched_84_2weeks_age_sex_2nd_dose)
matched_84_data_2weeks_age_sex_2nd_dose <- match.data(matched_84_2weeks_age_sex_2nd_dose)




##calculate proportions
proportions_28 <- matched_28_data_2weeks_age_sex_2nd_dose %>% group_by(case, lc_response_28) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
case_28_counts <- proportions_28 %>% filter(case == 1) %>% select(n)
control_28_counts <- proportions_28 %>% filter(case == 0) %>% select(n)

case_28_ci_lc_yes <- scoreci(x = case_28_counts[2,2], n = sum(case_28_counts[,2]), conf.level = 0.95)
case_28_ci_lc_no <- scoreci(x = case_28_counts[1,2], n = sum(case_28_counts[,2]), conf.level = 0.95)

case_28_confints_lc_yes <- data.frame(case_28_ci_lc_yes$conf.int)
colnames(case_28_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
case_28_confints_lc_yes <- case_28_confints_lc_yes %>% mutate(case = 1, lc_response_28 = 1)

case_28_confints_lc_no <- data.frame(case_28_ci_lc_no$conf.int)
colnames(case_28_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
case_28_confints_lc_no <- case_28_confints_lc_no %>% mutate(case = 1, lc_response_28 = 0)

control_28_ci_lc_yes <- scoreci(x = control_28_counts[2,2], n = sum(control_28_counts[,2]), conf.level = 0.95)
control_28_ci_lc_no <- scoreci(x = control_28_counts[1,2], n = sum(control_28_counts[,2]), conf.level = 0.95)

control_28_confints_lc_yes <- data.frame(control_28_ci_lc_yes$conf.int)
colnames(control_28_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
control_28_confints_lc_yes <- control_28_confints_lc_yes %>% mutate(case = 0, lc_response_28 = 1)

control_28_confints_lc_no <- data.frame(control_28_ci_lc_no$conf.int)
colnames(control_28_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
control_28_confints_lc_no <- control_28_confints_lc_no %>% mutate(case = 0, lc_response_28 = 0)

case_28_confints <- rbind(case_28_confints_lc_no, case_28_confints_lc_yes, control_28_confints_lc_yes, control_28_confints_lc_no)

results_28_fully_vacc <- left_join(proportions_28, case_28_confints, by = c("case", "lc_response_28"))
write.csv(results_28_fully_vacc, paste0(out_dir, "\\two_dose_proportions_28.csv"), row.names = F)


proportions_84 <- matched_84_data_2weeks_age_sex_2nd_dose %>% group_by(case, lc_response_84) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
case_84_counts <- proportions_84 %>% filter(case == 1) %>% select(n)
control_84_counts <- proportions_84 %>% filter(case == 0) %>% select(n)

case_84_ci_lc_yes <- scoreci(x = case_84_counts[2,2], n = sum(case_84_counts[,2]), conf.level = 0.95)
case_84_ci_lc_no <- scoreci(x = case_84_counts[1,2], n = sum(case_84_counts[,2]), conf.level = 0.95)


case_84_confints_lc_yes <- data.frame(case_84_ci_lc_yes$conf.int)
colnames(case_84_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
case_84_confints_lc_yes <- case_84_confints_lc_yes %>% mutate(case = 1, lc_response_84 = 1)

case_84_confints_lc_no <- data.frame(case_84_ci_lc_no$conf.int)
colnames(case_84_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
case_84_confints_lc_no <- case_84_confints_lc_no %>% mutate(case = 1, lc_response_84 = 0)


control_84_ci_lc_yes <- scoreci(x = control_84_counts[2,2], n = sum(control_84_counts[,2]), conf.level = 0.95)
control_84_ci_lc_no <- scoreci(x = control_84_counts[1,2], n = sum(control_84_counts[,2]), conf.level = 0.95)


control_84_confints_lc_yes <- data.frame(control_84_ci_lc_yes$conf.int)
colnames(control_84_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
control_84_confints_lc_yes <- control_84_confints_lc_yes %>% mutate(case = 0, lc_response_84 = 1)

control_84_confints_lc_no <- data.frame(control_84_ci_lc_no$conf.int)
colnames(control_84_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
control_84_confints_lc_no <- control_84_confints_lc_no %>% mutate(case = 0, lc_response_84 = 0)

case_84_confints <- rbind(case_84_confints_lc_no, case_84_confints_lc_yes, control_84_confints_lc_yes, control_84_confints_lc_no)

results_84_fully_vacc <- left_join(proportions_84, case_84_confints, by = c("case", "lc_response_84"))
write.csv(results_84_fully_vacc, paste0(out_dir, "\\two_dose_proportions_84.csv"), row.names = F)


##undadjusted 
mod_28_unadj_2_dose <- glm(lc_response_28 ~ case, data = matched_28_data_2weeks_age_sex_2nd_dose, family = binomial)
vcov_28_unadj_2_dose <- vcovCL(mod_28_unadj_2_dose, cluster = ~subclass)
output_28_unadj_2_dose <- coeftest(mod_28_unadj_2_dose, vcov = vcov_28_unadj_2_dose)

unadjust_28_results_2_dose <- data.frame(model = "unadjusted",
                                         OR = exp(output_28_unadj_2_dose[2,1]),
                                         lower_ci = exp(output_28_unadj_2_dose[2,1] - 1.96*output_28_unadj_2_dose[2,2]),
                                         upper_ci = exp(output_28_unadj_2_dose[2,1] + 1.96*output_28_unadj_2_dose[2,2]),
                                         p = output_28_unadj_2_dose[2,4],
                                         row.names = NULL,
                                         stringsAsFactors = FALSE)

##adjusted
mod_28_adj_2_dose <- glm(lc_response_28 ~ case 
                  + ns(age_lc_28, df=2, Boundary.knots=quantile(age_lc_28, c(.10, .90)))
                  + sex 
                  + health_conditions_28 
                  + covid_admitted_ever 
                  + ns(duration_between_inf_response_28, df=2, Boundary.knots=quantile(duration_between_inf_response_28, c(.10, .90))) 
                  + imd_quintile 
                  + white,
                  data = matched_28_data_2weeks_age_sex_2nd_dose, family = binomial)
vcov_28_adj_2_dose <- vcovCL(mod_28_adj_2_dose, cluster = ~subclass)
output_28_adj_2_dose <- coeftest(mod_28_adj_2_dose, vcov = vcov_28_adj_2_dose)

adjust_28_results_2_dose <- data.frame(model = "adjusted",
                                           OR = exp(output_28_adj_2_dose[2,1]),
                                           lower_ci = exp(output_28_adj_2_dose[2,1] - 1.96*output_28_adj_2_dose[2,2]),
                                           upper_ci = exp(output_28_adj_2_dose[2,1] + 1.96*output_28_adj_2_dose[2,2]),
                                           p = output_28_adj_2_dose[2,4],
                                           row.names = NULL,
                                           stringsAsFactors = FALSE)

write.csv(rbind(unadjust_28_results_2_dose, adjust_28_results_2_dose), paste0(out_dir, "\\two_dose_OR_28.csv"), row.names = F) 



##undadjusted 
mod_84_unadj_2_dose <- glm(lc_response_84 ~ case, data = matched_84_data_2weeks_age_sex_2nd_dose, family = binomial)
vcov_84_unadj_2_dose <- vcovCL(mod_84_unadj_2_dose, cluster = ~subclass)
output_84_unadj_2_dose <- coeftest(mod_84_unadj_2_dose, vcov = vcov_84_unadj_2_dose)

unadjust_84_results_2_dose <- data.frame(model = "unadjusted",
                                         OR = exp(output_84_unadj_2_dose[2,1]),
                                         lower_ci = exp(output_84_unadj_2_dose[2,1] - 1.96*output_84_unadj_2_dose[2,2]),
                                         upper_ci = exp(output_84_unadj_2_dose[2,1] + 1.96*output_84_unadj_2_dose[2,2]),
                                         p = output_84_unadj_2_dose[2,4],
                                         row.names = NULL,
                                         stringsAsFactors = FALSE)

##adjusted
mod_84_adj_2_dose <- glm(lc_response_84 ~ case 
                  + ns(age_lc_84, df=2, Boundary.knots=quantile(age_lc_84, c(.10, .90)))
                  + sex 
                  + health_conditions_84 
                  + covid_admitted_ever 
                  + ns(duration_between_inf_response_84, df=2, Boundary.knots=quantile(duration_between_inf_response_84, c(.10, .90))) 
                  + imd_quintile 
                  + white, 
                  data = matched_84_data_2weeks_age_sex_2nd_dose, family = binomial)
vcov_84_adj_2_dose <- vcovCL(mod_84_adj_2_dose, cluster = ~subclass)
output_84_adj_2_dose <- coeftest(mod_84_adj_2_dose, vcov = vcov_84_adj_2_dose)

adjust_84_results_2_dose <- data.frame(model = "adjusted",
                                       OR = exp(output_84_adj_2_dose[2,1]),
                                       lower_ci = exp(output_84_adj_2_dose[2,1] - 1.96*output_84_adj_2_dose[2,2]),
                                       upper_ci = exp(output_84_adj_2_dose[2,1] + 1.96*output_84_adj_2_dose[2,2]),
                                       p = output_84_adj_2_dose[2,4],
                                       row.names = NULL,
                                       stringsAsFactors = FALSE)

write.csv(rbind(unadjust_84_results_2_dose, adjust_84_results_2_dose), paste0(out_dir, "\\two_dose_OR_84.csv"), row.names = F) 


##het effects

mod_28_adj_2_dose <- glm(lc_response_28 ~ case
                         #+ age_40_plus:case
                         + ns(age_lc_28, df=2, Boundary.knots=quantile(age_lc_28, c(.10, .90)))*case
                         + sex 
                         + health_conditions_28 
                         + covid_admitted_ever 
                         + ns(duration_between_inf_response_28, df=2, Boundary.knots=quantile(duration_between_inf_response_28, c(.10, .90))) 
                         + imd_quintile 
                         + white,
                         data = matched_28_data_2weeks_age_sex_2nd_dose, family = binomial)
vcov_28_adj_2_dose <- vcovCL(mod_28_adj_2_dose, cluster = ~subclass)
output_28_adj_2_dose <- coeftest(mod_28_adj_2_dose, vcov = vcov_28_adj_2_dose)

#OR for <55 is exp(case coefficient) and 95% CI can be calculated using SE for case
#OR for >=55 is exp(case coefficient + case:age interaction term)
#Use vcov matrix to calculate SE for >= 55 sqrt(var(case coefficient) + var(case:age) + 2*cov(case coefficient, case:age))
#Var is diagonals [1,1], [2,2], [3,3], [4,4]
#covar for case and case:age is [2,4]

















ggplot(dat_matching_84, aes(x = age_at_infection, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Age (years)") +
  scale_fill_discrete("Case")

ggplot(dat_matching_84, aes(x = duration_between_inf_response_84, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Days between infection & response") +
  scale_fill_discrete("Case")

##remove 3 people who received 2nd jab same day as first infection
##filter to people with 28 days fu time from first infection date 
dat_matching_28 <- dat_lc_vacc %>% filter(visit_28_after_infection == 1) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(first_vacc_between_infection_lc_response_28 == 0 & second_vacc_between_infection_lc_response_28 == 0) %>%
  filter(vacc2_1_to_27_days_before_infection == 0 & vacc1_1_to_27_days_before_infection == 0) %>%
  mutate(case = ifelse(vacc2_28_days_before_infection == 1, 1, 0),
         duration_between_inf_response_28_grouped = cut(duration_between_inf_response_28, seq(28, 119, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_28_grouped = ifelse(is.na(duration_between_inf_response_28_grouped), ">=119", duration_between_inf_response_28_grouped),
         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
  select(participant_id, case, age_at_infection, age_55_plus, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_28_date, duration_between_inf_response_28, duration_between_inf_response_28_grouped, covid_vaccine_date1, covid_vaccine_date2, lc_response_28, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_28_days_before_infection, vacc2_28_days_before_infection, first_vacc_between_infection_lc_response_28, vacc1_1_to_27_days_before_infection, vacc2_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_28)

ggplot(dat_matching_28, aes(x = age_at_infection, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Age (years)") +
  scale_fill_discrete("Case")

ggplot(dat_matching_28, aes(x = duration_between_inf_response_28, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Days between infection & response") +
  scale_fill_discrete("Case")

##matching
  matched_28 <- matchit(case ~ age_group + sex + duration_between_inf_response_28_grouped + health_conditions_at_infection, 
                    data = dat_matching_28,
                    method = "nearest",
                    exact = ~ age_group + sex + duration_between_inf_response_28_grouped + health_conditions_at_infection) 
  summary(matched_28)
  matched_28_data <- match.data(matched_28)
  
 props_28 <- table(matched_28_data$lc_response_28, matched_28_data$case)
 
 prop.test(x = props_28[2,], n = c(colSums(props_28)))
 
 
 proportions_28 <- matched_28_data %>% group_by(case, lc_response_28) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
 case_28_counts <- proportions_28 %>% filter(case == 1) %>% select(n)
 control_28_counts <- proportions_28 %>% filter(case == 0) %>% select(n)
 
 case_28_ci_lc_yes <- scoreci(x = case_28_counts[2,2], n = sum(case_28_counts[,2]), conf.level = 0.95)
 case_28_ci_lc_no <- scoreci(x = case_28_counts[1,2], n = sum(case_28_counts[,2]), conf.level = 0.95)
 
 case_28_confints_lc_yes <- data.frame(case_28_ci_lc_yes$conf.int)
 colnames(case_28_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
 case_28_confints_lc_yes <- case_28_confints_lc_yes %>% mutate(case = 1, lc_response_28 = 1)

 case_28_confints_lc_no <- data.frame(case_28_ci_lc_no$conf.int)
 colnames(case_28_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
 case_28_confints_lc_no <- case_28_confints_lc_no %>% mutate(case = 1, lc_response_28 = 0)
 
 control_28_ci_lc_yes <- scoreci(x = control_28_counts[2,2], n = sum(control_28_counts[,2]), conf.level = 0.95)
 control_28_ci_lc_no <- scoreci(x = control_28_counts[1,2], n = sum(control_28_counts[,2]), conf.level = 0.95)
 
 control_28_confints_lc_yes <- data.frame(control_28_ci_lc_yes$conf.int)
 colnames(control_28_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
 control_28_confints_lc_yes <- control_28_confints_lc_yes %>% mutate(case = 0, lc_response_28 = 1)
 
 control_28_confints_lc_no <- data.frame(control_28_ci_lc_no$conf.int)
 colnames(control_28_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
 control_28_confints_lc_no <- control_28_confints_lc_no %>% mutate(case = 0, lc_response_28 = 0)
 
 case_28_confints <- rbind(case_28_confints_lc_no, case_28_confints_lc_yes, control_28_confints_lc_yes, control_28_confints_lc_no)
 
 results_28 <- left_join(proportions_28, case_28_confints, by = c("case", "lc_response_28"))

 ##undadjusted 
 mod <- glm(lc_response_28 ~ case, data = matched_28_data, family = binomial)
 vcov <- vcovCL(mod, cluster = ~subclass)
 output_28 <- coeftest(mod, vcov = vcov)
 
 ##adjusted
 mod <- glm(lc_response_28 ~ case + age_at_infection + sex + health_conditions_at_infection + duration_between_inf_response_28 + imd_quintile + white, data = matched_28_data, family = binomial)
 vcov <- vcovCL(mod, cluster = ~subclass)
 output_28 <- coeftest(mod, vcov = vcov)

 ##het effects
 ##unadjusted
 mod <- glm(lc_response_28 ~ case + age_55_plus*case, data = matched_28_data, family = binomial)
 vcov <- vcovCL(mod, cluster = ~subclass)
 output_28 <- coeftest(mod, vcov = vcov)
 
 ##adjusted
 mod <- glm(lc_response_28 ~ case + age_55_plus*case + age_at_infection + sex + health_conditions_at_infection + duration_between_inf_response_28, data = matched_28_data, family = binomial)
 vcov <- vcovCL(mod, cluster = ~subclass)
 output_28 <- coeftest(mod, vcov = vcov)
 
 #OR for <55 is exp(case coefficient) and 95% CI can be calculated using SE for case
 #OR for >=55 is exp(case coefficient + case:age interaction term)
 #Use vcov matrix to calculate SE for >= 55 sqrt(var(case coefficient) + var(case:age) + 2*cov(case coefficient, case:age))
 #Var is diagonals [1,1], [2,2], [3,3], [4,4]
 #covar for case and case:age is [2,4]
 
 matched_84 <- matchit(case ~ age_group + sex + duration_between_inf_response_84_grouped_7 + health_conditions_at_infection, 
                        data = dat_matching_84,
                        method = "nearest",
                        exact = ~ age_group + sex + duration_between_inf_response_84_grouped_7 + health_conditions_at_infection) 
  summary(matched_84)
  matched_84_data <- match.data(matched_84)
  
  props_84 <- table(matched_84_data$lc_response_84, matched_84_data$case)
  
  prop.test(x = props_84[2,], n = c(colSums(props_84)))
  
  proportions_84 <- matched_84_data %>% group_by(case, lc_response_84) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
  case_84_counts <- proportions_84 %>% filter(case == 1) %>% select(n)
  control_84_counts <- proportions_84 %>% filter(case == 0) %>% select(n)
  
  case_84_ci_lc_yes <- scoreci(x = case_84_counts[2,2], n = sum(case_84_counts[,2]), conf.level = 0.95)
  case_84_ci_lc_no <- scoreci(x = case_84_counts[1,2], n = sum(case_84_counts[,2]), conf.level = 0.95)
  
  
  case_84_confints_lc_yes <- data.frame(case_84_ci_lc_yes$conf.int)
  colnames(case_84_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
  case_84_confints_lc_yes <- case_84_confints_lc_yes %>% mutate(case = 1, lc_response_84 = 1)
  
  case_84_confints_lc_no <- data.frame(case_84_ci_lc_no$conf.int)
  colnames(case_84_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
  case_84_confints_lc_no <- case_84_confints_lc_no %>% mutate(case = 1, lc_response_84 = 0)
  
  
  control_84_ci_lc_yes <- scoreci(x = control_84_counts[2,2], n = sum(control_84_counts[,2]), conf.level = 0.95)
  control_84_ci_lc_no <- scoreci(x = control_84_counts[1,2], n = sum(control_84_counts[,2]), conf.level = 0.95)
  
  
  control_84_confints_lc_yes <- data.frame(control_84_ci_lc_yes$conf.int)
  colnames(control_84_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
  control_84_confints_lc_yes <- control_84_confints_lc_yes %>% mutate(case = 0, lc_response_84 = 1)
  
  control_84_confints_lc_no <- data.frame(control_84_ci_lc_no$conf.int)
  colnames(control_84_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
  control_84_confints_lc_no <- control_84_confints_lc_no %>% mutate(case = 0, lc_response_84 = 0)
  
  case_84_confints <- rbind(case_84_confints_lc_no, case_84_confints_lc_yes, control_84_confints_lc_yes, control_84_confints_lc_no)
  
  results_84 <- left_join(proportions_84, case_84_confints, by = c("case", "lc_response_84"))
  
  ##unadjusted
  mod <- glm(lc_response_84 ~ case, data = matched_84_data, family = binomial)
  vcov <- vcovCL(mod, cluster = ~subclass)
  output_84 <- coeftest(mod, vcov = vcov)
  
  ##adjusted
  mod <- glm(lc_response_84 ~ case + age_at_infection + sex + health_conditions_at_infection + duration_between_inf_response_84 + imd_quintile + white, data = matched_84_data, family = binomial)
  vcov <- vcovCL(mod, cluster = ~subclass)
  output_84 <- coeftest(mod, vcov = vcov)
  
  ###subgroup analysis excluding cases who received their second dose between infection and response
  dat_matching_84 <- dat_lc_vacc %>% filter(visit_84_after_infection == 1) %>%
    filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
    filter(first_vacc_between_infection_lc_response_84 == 0) %>%
    filter(is.na(second_vacc_between_infection_lc_response_84) | second_vacc_between_infection_lc_response_84 == 0) %>%
    filter(vacc_1_to_27_days_before_infection == 0) %>%
    mutate(case = ifelse(vacc_28_days_before_infection == 1, 1, 0),
           duration_between_inf_response_84_grouped = cut(duration_between_inf_response_84, seq(84, 196, 7), include.lowest = FALSE, right = FALSE),
           duration_between_inf_response_84_grouped = ifelse(is.na(duration_between_inf_response_84_grouped), ">=196", duration_between_inf_response_84_grouped),
           health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
    select(participant_id, case, age_at_infection, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_84_date, duration_between_inf_response_84, duration_between_inf_response_84_grouped, covid_vaccine_date1, covid_vaccine_date2, lc_response_84, vacc_1_days_before_infection, vacc_28_days_before_infection, first_vacc_between_infection_lc_response_84, vacc_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_84)
  
  
  ##remove 3 people who received 2nd jab same day as first infection
  ##filter to people with 28 days fu time from first infection date 
  dat_matching_28 <- dat_lc_vacc %>% filter(visit_28_after_infection == 1) %>%
    filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
    filter(first_vacc_between_infection_lc_response_28 == 0) %>%
    filter(is.na(second_vacc_between_infection_lc_response_28) | second_vacc_between_infection_lc_response_28 == 0) %>%
    filter(vacc_1_to_27_days_before_infection == 0) %>%
    mutate(case = ifelse(vacc_28_days_before_infection == 1, 1, 0),
           duration_between_inf_response_28_grouped = cut(duration_between_inf_response_28, seq(28, 119, 7), include.lowest = FALSE, right = FALSE),
           duration_between_inf_response_28_grouped = ifelse(is.na(duration_between_inf_response_28_grouped), ">=119", duration_between_inf_response_28_grouped),
           health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
    select(participant_id, case, age_at_infection, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_28_date, duration_between_inf_response_28, duration_between_inf_response_28_grouped, covid_vaccine_date1, covid_vaccine_date2, lc_response_28, vacc_1_days_before_infection, vacc_28_days_before_infection, first_vacc_between_infection_lc_response_28, vacc_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_28)
  
  ##matching
  matched_28 <- matchit(case ~ age_group + sex + duration_between_inf_response_28_grouped + health_conditions_at_infection, 
                        data = dat_matching_28,
                        method = "nearest",
                        exact = ~ age_group + sex + duration_between_inf_response_28_grouped + health_conditions_at_infection) 
  summary(matched_28)
  matched_28_data <- match.data(matched_28)
  
  props_28 <- table(matched_28_data$lc_response_28, matched_28_data$case)
  
  prop.test(x = props_28[2,], n = c(colSums(props_28)))
  
  
  proportions_28 <- matched_28_data %>% group_by(case, lc_response_28) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
  case_28_counts <- proportions_28 %>% filter(case == 1) %>% select(n)
  control_28_counts <- proportions_28 %>% filter(case == 0) %>% select(n)
  
  case_28_ci_lc_yes <- scoreci(x = case_28_counts[2,2], n = sum(case_28_counts[,2]), conf.level = 0.95)
  case_28_ci_lc_no <- scoreci(x = case_28_counts[1,2], n = sum(case_28_counts[,2]), conf.level = 0.95)
  
  case_28_confints_lc_yes <- data.frame(case_28_ci_lc_yes$conf.int)
  colnames(case_28_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
  case_28_confints_lc_yes <- case_28_confints_lc_yes %>% mutate(case = 1, lc_response_28 = 1)
  
  case_28_confints_lc_no <- data.frame(case_28_ci_lc_no$conf.int)
  colnames(case_28_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
  case_28_confints_lc_no <- case_28_confints_lc_no %>% mutate(case = 1, lc_response_28 = 0)
  
  control_28_ci_lc_yes <- scoreci(x = control_28_counts[2,2], n = sum(control_28_counts[,2]), conf.level = 0.95)
  control_28_ci_lc_no <- scoreci(x = control_28_counts[1,2], n = sum(control_28_counts[,2]), conf.level = 0.95)
  
  control_28_confints_lc_yes <- data.frame(control_28_ci_lc_yes$conf.int)
  colnames(control_28_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
  control_28_confints_lc_yes <- control_28_confints_lc_yes %>% mutate(case = 0, lc_response_28 = 1)
  
  control_28_confints_lc_no <- data.frame(control_28_ci_lc_no$conf.int)
  colnames(control_28_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
  control_28_confints_lc_no <- control_28_confints_lc_no %>% mutate(case = 0, lc_response_28 = 0)
  
  case_28_confints <- rbind(case_28_confints_lc_no, case_28_confints_lc_yes, control_28_confints_lc_yes, control_28_confints_lc_no)
  
  results_28 <- left_join(proportions_28, case_28_confints, by = c("case", "lc_response_28"))
  
  
  ##unadjusted
  mod <- glm(lc_response_28 ~ case, data = matched_28_data, family = binomial)
  vcov <- vcovCL(mod, cluster = ~subclass)
  output_28 <- coeftest(mod, vcov = vcov)
  
  ##adjusted
  mod <- glm(lc_response_28 ~ case + age_at_infection + sex + health_conditions_at_infection + duration_between_inf_response_28 + imd_quintile + white, data = matched_28_data, family = binomial)
  vcov <- vcovCL(mod, cluster = ~subclass)
  output_28 <- coeftest(mod, vcov = vcov)
  
  matched_84 <- matchit(case ~ age_group + sex + duration_between_inf_response_84_grouped + health_conditions_at_infection, 
                        data = dat_matching_84,
                        method = "nearest",
                        exact = ~ age_group + sex + duration_between_inf_response_84_grouped + health_conditions_at_infection) 
  summary(matched_84)
  matched_84_data <- match.data(matched_84)
  
  props_84 <- table(matched_84_data$lc_response_84, matched_84_data$case)
  
  prop.test(x = props_84[2,], n = c(colSums(props_84)))
  
  proportions_84 <- matched_84_data %>% group_by(case, lc_response_84) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
  case_84_counts <- proportions_84 %>% filter(case == 1) %>% select(n)
  control_84_counts <- proportions_84 %>% filter(case == 0) %>% select(n)
  
  case_84_ci_lc_yes <- scoreci(x = case_84_counts[2,2], n = sum(case_84_counts[,2]), conf.level = 0.95)
  case_84_ci_lc_no <- scoreci(x = case_84_counts[1,2], n = sum(case_84_counts[,2]), conf.level = 0.95)
  
  
  case_84_confints_lc_yes <- data.frame(case_84_ci_lc_yes$conf.int)
  colnames(case_84_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
  case_84_confints_lc_yes <- case_84_confints_lc_yes %>% mutate(case = 1, lc_response_84 = 1)
  
  case_84_confints_lc_no <- data.frame(case_84_ci_lc_no$conf.int)
  colnames(case_84_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
  case_84_confints_lc_no <- case_84_confints_lc_no %>% mutate(case = 1, lc_response_84 = 0)
  
  
  control_84_ci_lc_yes <- scoreci(x = control_84_counts[2,2], n = sum(control_84_counts[,2]), conf.level = 0.95)
  control_84_ci_lc_no <- scoreci(x = control_84_counts[1,2], n = sum(control_84_counts[,2]), conf.level = 0.95)
  
  
  control_84_confints_lc_yes <- data.frame(control_84_ci_lc_yes$conf.int)
  colnames(control_84_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
  control_84_confints_lc_yes <- control_84_confints_lc_yes %>% mutate(case = 0, lc_response_84 = 1)
  
  control_84_confints_lc_no <- data.frame(control_84_ci_lc_no$conf.int)
  colnames(control_84_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
  control_84_confints_lc_no <- control_84_confints_lc_no %>% mutate(case = 0, lc_response_84 = 0)
  
  
  
  case_84_confints <- rbind(case_84_confints_lc_no, case_84_confints_lc_yes, control_84_confints_lc_yes, control_84_confints_lc_no)
  
  results_84 <- left_join(proportions_84, case_84_confints, by = c("case", "lc_response_84"))
  
  mod <- glm(lc_response_84 ~ case, data = matched_84_data, family = binomial)
  vcov <- vcovCL(mod, cluster = ~subclass)
  output_84 <- coeftest(mod, vcov = vcov)