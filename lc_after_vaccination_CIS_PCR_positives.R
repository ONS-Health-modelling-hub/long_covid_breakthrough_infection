library(haven)
library(sqldf)
library(tidyverse)
library(MatchIt)
library(stats)
library(PropCIs)
library(sandwich)
library(lmtest)
library(ggplot2)

out_dir = "filepath"
dataset_date = "20211004"
cutoff_date = "2021-10-01"
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
  "sympt_now_any",
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

### keep all visits belonging to participants who ever tested positive by PCR in CIS
part_select <- dat_all$participant_id[dat_all$result_mk==1]
part_select <- part_select[!duplicated(part_select)]
dat_all <- dat_all[dat_all$participant_id %in% part_select,]

### sort dataset by participant ID, visit date, and visit ID
dat_all <- dat_all[with(dat_all, order(participant_id, visit_date, -xtfrm(visit_id))),]

#################################### JOIN COVID STATUS VARIABLES ####################################

### Create first positive study PCR test date
dat_study_swab_pos1_date <- sqldf("
  select
    participant_id,
    min(age_at_visit) as age_at_infection,
    min(visit_date) as study_swab_pos1_date
  from dat_all
  where result_mk=1
  group by participant_id
")

### merge first positive study PCR test date
dat_all <- merge(x = dat_all,
                 y = dat_study_swab_pos1_date,
                 by.x = "participant_id",
                 by.y = "participant_id",
                 all.x = TRUE,
                 all.y = FALSE)


### coerce necessary variables to date format
dat_all$covid_test_swab_pos_first_date <- as.Date(dat_all$covid_test_swab_pos_first_date)
dat_all$covid_test_blood_pos_first_date <- as.Date(dat_all$covid_test_blood_pos_first_date,
                                                   origin = "1970-01-01")
dat_all$covid_date <- as.Date(dat_all$covid_date)
dat_all$visit_date <- as.Date(dat_all$visit_date, origin = "1970-01-01")
dat_all$study_swab_pos1_date <- as.Date(dat_all$study_swab_pos1_date, origin = "1970-01-01")

### COVID-19 arrived in the UK on 24 Jan 2020 - set any dates before this to NA
dat_all$covid_test_swab_pos_first_date[dat_all$covid_test_swab_pos_first_date < as.Date("2020-01-24")] <- NA
dat_all$covid_test_blood_pos_first_date[dat_all$covid_test_blood_pos_first_date < as.Date("2020-01-24")] <- NA
dat_all$covid_date[dat_all$covid_date < as.Date("2020-01-24")] <- NA

dat_all <- dat_all %>% mutate(time_since_non_study_previous_positive_swab = as.numeric(study_swab_pos1_date - covid_test_swab_pos_first_date),
                              time_since_suspected_covid = as.numeric(study_swab_pos1_date - covid_date))

non_study_swab_dates <- dat_all %>% group_by(participant_id) %>% 
  summarise(time_since_non_study_previous_positive_swab = max(time_since_non_study_previous_positive_swab, na.rm = T)) %>% 
  mutate(time_since_non_study_previous_positive_swab = ifelse(time_since_non_study_previous_positive_swab=="-Inf", NA, time_since_non_study_previous_positive_swab))

###remove people who had a confirmed or suspected infection >120 days before first positive test in CIS

part_select <- non_study_swab_dates$participant_id[is.na(non_study_swab_dates$time_since_non_study_previous_positive_swab) | non_study_swab_dates$time_since_non_study_previous_positive_swab<120]
#part_select <- part_select[!duplicated(part_select)]
dat_all <- dat_all[dat_all$participant_id %in% part_select,]

think_covid_dates <- dat_all %>% group_by(participant_id) %>% 
  summarise(time_since_suspected_covid = max(time_since_suspected_covid, na.rm = T)) %>% 
  mutate(time_since_suspected_covid = ifelse(time_since_suspected_covid =="-Inf", NA, time_since_suspected_covid))

part_select <- think_covid_dates$participant_id[is.na(think_covid_dates$time_since_suspected_covid) | think_covid_dates$time_since_suspected_covid<120]
#part_select <- part_select[!duplicated(part_select)]
dat_all <- dat_all[dat_all$participant_id %in% part_select,]

dat_all <- dat_all %>% mutate(duration_after_infection = as.numeric(visit_date - study_swab_pos1_date),
                              visit_28_after_infection = ifelse(duration_after_infection >= 28, 1, 0),
                              visit_84_after_infection = ifelse(duration_after_infection >= 84, 1, 0))

### aggregate to person level
dat_lc_status <- sqldf("
  select
    participant_id,
    sex,
    max(long_covid_have_symptoms) as lc_ever
  from dat_all
  group by participant_id
")

dat_lc_date <- sqldf("
  select
    participant_id,
    min(visit_date) as lc_first_date
  from dat_all
  where long_covid_have_symptoms=1
  group by participant_id
")

dat_lc_q <- sqldf("
  select
    participant_id,
    max(visit_date) as lc_last_question_date
  from dat_all
  where long_covid_have_symptoms is not null
  group by participant_id
")

dat_lc_last_yes <- sqldf("
  select
    participant_id,
    max(visit_date) as lc_last_yes
  from dat_all
  where long_covid_have_symptoms = 1
  group by participant_id
")


dat_i_to_r_84 <- sqldf("
  select
    participant_id,
    visit_84_after_infection,
    min(visit_date) as first_response_84_date,
    long_covid_have_symptoms as lc_response_84
  from dat_all
  where visit_84_after_infection = 1 
  AND long_covid_have_symptoms is not null
  group by participant_id
")

dat_i_to_r_28 <- sqldf("
  select
    participant_id,
    visit_28_after_infection,
    min(visit_date) as first_response_28_date,
    long_covid_have_symptoms as lc_response_28
  from dat_all
  where visit_28_after_infection = 1 
  AND long_covid_have_symptoms is not null
  group by participant_id
")

dat_health_conditions <- sqldf("
  select
    participant_id,
    health_conditions as health_conditions_at_infection,
    imd_samp,
    ethnicityg,
    country
  from dat_all
  where study_swab_pos1_date = visit_date
  group by participant_id
")

dat_health_conditions_before <- sqldf("
  select
    participant_id,
    max(health_conditions) as health_conditions_before_infection
  from dat_all
  where study_swab_pos1_date >= visit_date
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

dat_infected_lc <- merge(x = dat_study_swab_pos1_date,
                         y = dat_lc_status_date_q_yes,
                         by.x = "participant_id",
                         by.y = "participant_id",
                         all.x = FALSE,
                         all.y = TRUE)

dat_infected_lc <- merge(x = dat_infected_lc,
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

dat_infected_lc <- merge(x = dat_infected_lc,
                         y = dat_health_conditions,
                         by.x = "participant_id",
                         by.y = "participant_id",
                         all.x = TRUE,
                         all.y = FALSE)

dat_infected_lc <- merge(x = dat_infected_lc,
                         y = dat_health_conditions_before,
                         by.x = "participant_id",
                         by.y = "participant_id",
                         all.x = TRUE,
                         all.y = FALSE)

dat_infected_lc <- dat_infected_lc %>% mutate(lc_onset_after_first_infection_date = ifelse(lc_first_date >= study_swab_pos1_date, 1, 0),
                                              diff_lc_first_date_infection_date = as.numeric(as.Date(lc_first_date, origin = "1970-01-01") - as.Date(study_swab_pos1_date, origin="1970-01-01")),
                                              diff_infection_lc_q_last_answered = as.numeric(as.Date(lc_last_question_date, origin = "1970-01-01") - as.Date(study_swab_pos1_date, origin="1970-01-01")),
                                              diff_infection_lc_q_last_yes = as.numeric(as.Date(lc_last_yes, origin = "1970-01-01") - as.Date(study_swab_pos1_date, origin="1970-01-01")))

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

### vaccination date from vaccination data
dat_vacc <- haven::read_dta(paste0("filepath/data_participant_vaccination.dta"),
                            col_select = c("participant_id",
                                           "covid_vaccine_date1",
                                           "covid_vaccine_date2",
                                           "covid_vaccine_type1",
                                           "covid_vaccine_type2"))

dat_positive_vacc <- merge(x = dat_infected_lc,
                           y = dat_vacc,
                           by.x = "participant_id",
                           by.y = "participant_id",
                           all.x = TRUE,
                           all.y = FALSE)

###filter out anyone who has never been vaccinated and 14 people who report having LC before a positive PCR test in the CIS
##restrict analysis to people aged 18 to 70
dat_positive_vacc <- dat_positive_vacc %>%
                            mutate(study_swab_pos1_date = as.Date(study_swab_pos1_date, origin = "1970-01-01"),
                                   first_response_84_date = as.Date(first_response_84_date, origin = "1970-01-01" ),
                                   first_response_28_date = as.Date(first_response_28_date, origin = "1970-01-01" )) %>%
                            filter(!is.na(covid_vaccine_date1)) %>%
                            filter(is.na(lc_onset_after_first_infection_date) | lc_onset_after_first_infection_date == 1) %>%
                            filter(age_at_infection >= 18 & age_at_infection < 70)
dat_lc_vacc <- dat_positive_vacc %>%
                            mutate(diff_betweeen_vacc2_infection = as.Date(study_swab_pos1_date, origin = "1970-01-01") - covid_vaccine_date2,
                                   diff_betweeen_vacc1_infection = as.Date(study_swab_pos1_date, origin = "1970-01-01") - covid_vaccine_date1,
                                   vacc2_same_day_first_inf = ifelse(covid_vaccine_date2 == study_swab_pos1_date, 1, 0),
                              duration_between_inf_response_84 = as.numeric(first_response_84_date - study_swab_pos1_date),
                              duration_between_inf_response_28 = as.numeric(first_response_28_date - study_swab_pos1_date),
                              first_vacc_between_infection_lc_response_84 = ifelse(covid_vaccine_date1 >= study_swab_pos1_date & covid_vaccine_date1 <= first_response_84_date, 1, 0),
                              first_vacc_between_infection_lc_response_28 = ifelse(covid_vaccine_date1 >= study_swab_pos1_date & covid_vaccine_date1 <= first_response_28_date, 1, 0),
                              second_vacc_between_infection_lc_response_84 = ifelse(covid_vaccine_date2 >= study_swab_pos1_date & covid_vaccine_date2 <= first_response_84_date, 1, 0),
                              second_vacc_between_infection_lc_response_28 = ifelse(covid_vaccine_date2 >= study_swab_pos1_date & covid_vaccine_date2 <= first_response_28_date, 1, 0),
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
                              lc_yes_at_84_or_more = ifelse(diff_infection_lc_q_last_yes>=84, 1, 0),
                              age_group = case_when(age_at_infection >=18 & age_at_infection <= 24 ~ "18 to 24",
                                                    age_at_infection >=25 & age_at_infection <= 34 ~ "25 to 34",
                                                    age_at_infection >=35 & age_at_infection <= 49 ~ "35 to 49",
                                                    age_at_infection >=50 & age_at_infection <= 69 ~ "50 to 69",
                                                    age_at_infection >=70 ~ "Over 70"),
                              age_group_5yr = case_when(age_at_infection >=18 & age_at_infection <= 24 ~ "18 to 24",
                                                           age_at_infection >=25 & age_at_infection < 30 ~ "25 to 29",
                                                           age_at_infection >=30 & age_at_infection < 35 ~ "30 to 34",
                                                           age_at_infection >=35 & age_at_infection < 40 ~ "35 to 39",
                                                           age_at_infection >=40 & age_at_infection < 45 ~ "40 to 44",
                                                           age_at_infection >=45 & age_at_infection < 50 ~ "45 to 49",
                                                           age_at_infection >=50 & age_at_infection < 55 ~ "50 to 54",
                                                           age_at_infection >=55 & age_at_infection < 60 ~ "55 to 59",
                                                           age_at_infection >=60 & age_at_infection < 65 ~ "60 to 64",
                                                           age_at_infection >=65 & age_at_infection <= 70 ~ "65 to 69"),
                              age_group_10_year = case_when(age_at_infection >=18 & age_at_infection <= 29 ~ "18 to 29",
                                                        age_at_infection >=30 & age_at_infection < 40 ~ "30 to 39",
                                                        age_at_infection >=40 & age_at_infection < 50 ~ "40 to 49",
                                                        age_at_infection >=50 & age_at_infection <= 70 ~ "50 to 69"),
                              age_55_plus = ifelse(age_at_infection >= 55, 1, 0))


###single dose###
#dat_matching_84_4weeks <- dat_lc_vacc %>% filter(visit_84_after_infection == 1) %>%
#  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
#  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
#  filter(first_vacc_between_infection_lc_response_84 == 0) %>%                                 
#  filter(vacc1_1_to_27_days_before_infection == 0) %>%
#  mutate(case = ifelse(vacc1_28_days_before_infection == 1, 1, 0),
#         duration_between_inf_response_84_grouped_7 = cut(duration_between_inf_response_84, seq(84, 196, 7), include.lowest = FALSE, right = FALSE),
#         duration_between_inf_response_84_grouped_7 = ifelse(is.na(duration_between_inf_response_84_grouped_7), ">=196", duration_between_inf_response_84_grouped_7),
#         duration_between_inf_response_84_grouped_14 = cut(duration_between_inf_response_84, seq(84, 196, 14), include.lowest = FALSE, right = FALSE),
#         duration_between_inf_response_84_grouped_14 = ifelse(is.na(duration_between_inf_response_84_grouped_14), ">=196", duration_between_inf_response_84_grouped_14),
#         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
#  select(participant_id, case, age_at_infection, age_55_plus, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_84_date, duration_between_inf_response_84, duration_between_inf_response_84_grouped_7, duration_between_inf_response_84_grouped_14, covid_vaccine_date1, covid_vaccine_date2, lc_response_84, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_28_days_before_infection, vacc2_28_days_before_infection, first_vacc_between_infection_lc_response_84, vacc1_1_to_27_days_before_infection, vacc2_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_84)

dat_matching_84_3weeks <- dat_lc_vacc %>% filter(visit_84_after_infection == 1) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
  filter(first_vacc_between_infection_lc_response_84 == 0) %>%                                 
  filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = ifelse(vacc1_21_days_before_infection == 1, 1, 0),
         duration_between_inf_response_84_grouped_7 = cut(duration_between_inf_response_84, seq(84, 196, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_7 = ifelse(is.na(duration_between_inf_response_84_grouped_7), ">=196", duration_between_inf_response_84_grouped_7),
         duration_between_inf_response_84_grouped_14 = cut(duration_between_inf_response_84, seq(84, 196, 14), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_14 = ifelse(is.na(duration_between_inf_response_84_grouped_14), ">=196", duration_between_inf_response_84_grouped_14),
         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
  select(participant_id, case, age_at_infection, age_group_10_year, age_55_plus, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_84_date, duration_between_inf_response_84, duration_between_inf_response_84_grouped_7, duration_between_inf_response_84_grouped_14, covid_vaccine_date1, covid_vaccine_date2, lc_response_84, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_28_days_before_infection, vacc2_28_days_before_infection, first_vacc_between_infection_lc_response_84, vacc1_1_to_27_days_before_infection, vacc2_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_84)

#dat_matching_84_2weeks <- dat_lc_vacc %>% filter(visit_84_after_infection == 1) %>%
#  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
#  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
#  filter(first_vacc_between_infection_lc_response_84 == 0) %>%                                 
#  filter(vacc1_1_to_13_days_before_infection == 0) %>%
#  mutate(case = ifelse(vacc1_14_days_before_infection == 1, 1, 0),
#         duration_between_inf_response_84_grouped_7 = cut(duration_between_inf_response_84, seq(84, 196, 7), include.lowest = FALSE, right = FALSE),
#         duration_between_inf_response_84_grouped_7 = ifelse(is.na(duration_between_inf_response_84_grouped_7), ">=196", duration_between_inf_response_84_grouped_7),
#         duration_between_inf_response_84_grouped_14 = cut(duration_between_inf_response_84, seq(84, 196, 14), include.lowest = FALSE, right = FALSE),
#         duration_between_inf_response_84_grouped_14 = ifelse(is.na(duration_between_inf_response_84_grouped_14), ">=196", duration_between_inf_response_84_grouped_14),
#         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
#  select(participant_id, case, age_at_infection, age_55_plus, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_84_date, duration_between_inf_response_84, duration_between_inf_response_84_grouped_7, duration_between_inf_response_84_grouped_14, covid_vaccine_date1, covid_vaccine_date2, lc_response_84, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_28_days_before_infection, vacc2_28_days_before_infection, first_vacc_between_infection_lc_response_84, vacc1_1_to_27_days_before_infection, vacc2_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_84)

ggplot(dat_matching_84_3weeks, aes(x = age_at_infection, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Age (years)") +
  scale_fill_discrete("Case")

ggplot(dat_matching_84_3weeks, aes(x = duration_between_inf_response_84, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Days between infection & response") +
  scale_fill_discrete("Case")


##matching
matched_84_3weeks_age_sex <- matchit(case ~ age_group_10_year + sex + duration_between_inf_response_84, 
                                     data = dat_matching_84_3weeks,
                                     method = "nearest",
                                     exact = ~ age_group_10_year + sex + duration_between_inf_response_84) 
summary(matched_84_3weeks_age_sex)
matched_84_data_3weeks_age_sex <- match.data(matched_84_3weeks_age_sex)

#unmatched <- anti_join(dat_matching_84_3weeks, matched_84_data_3weeks_age_sex, by = "participant_id")
#unmatched_cases <- unmatched %>% filter(case == 1) %>% mutate(matched = 0)
#matched_84_data_3weeks_age_sex_cases <- matched_84_data_3weeks_age_sex %>% filter(case == 1) %>% mutate(matched = 1) %>% select(-distance, -weights, - subclass)
#matched_vs_unmatched <- rbind(matched_84_data_3weeks_age_sex_cases, unmatched_cases)

#ggplot(matched_vs_unmatched, aes(x = age_at_infection, fill = factor(matched))) +
#  geom_density(alpha = 0.4) +
#  xlab("Age (years)") +
#  scale_fill_discrete("Matched") +
#  scale_x_continuous(breaks = seq(20, 100, 5))

#dat_matching_28_4weeks <- dat_lc_vacc %>% filter(visit_28_after_infection == 1) %>%
#  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
#  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
#  filter(first_vacc_between_infection_lc_response_28 == 0) %>%
#  filter(vacc1_1_to_27_days_before_infection == 0) %>%
#  mutate(case = ifelse(vacc1_28_days_before_infection == 1, 1, 0),
#         duration_between_inf_response_28_grouped = cut(duration_between_inf_response_28, seq(28, 119, 7), include.lowest = FALSE, right = FALSE),
#         duration_between_inf_response_28_grouped = ifelse(is.na(duration_between_inf_response_28_grouped), ">=119", duration_between_inf_response_28_grouped),
#         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
#  select(participant_id, case, age_at_infection, age_55_plus, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_28_date, duration_between_inf_response_28, duration_between_inf_response_28_grouped, covid_vaccine_date1, covid_vaccine_date2, lc_response_28, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_28_days_before_infection, vacc2_28_days_before_infection, first_vacc_between_infection_lc_response_28, vacc1_1_to_27_days_before_infection, vacc2_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_28)

dat_matching_28_3weeks <- dat_lc_vacc %>% filter(visit_28_after_infection == 1) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
  filter(first_vacc_between_infection_lc_response_28 == 0) %>%
  filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = ifelse(vacc1_21_days_before_infection == 1, 1, 0),
         duration_between_inf_response_28_grouped = cut(duration_between_inf_response_28, seq(28, 119, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_28_grouped = ifelse(is.na(duration_between_inf_response_28_grouped), ">=119", duration_between_inf_response_28_grouped),
         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
  select(participant_id, case, age_at_infection, age_group_10_year, age_55_plus, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_28_date, duration_between_inf_response_28, duration_between_inf_response_28_grouped, covid_vaccine_date1, covid_vaccine_date2, lc_response_28, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_28_days_before_infection, vacc2_28_days_before_infection, first_vacc_between_infection_lc_response_28, vacc1_1_to_27_days_before_infection, vacc2_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_28)

#dat_matching_28_2weeks <- dat_lc_vacc %>% filter(visit_28_after_infection == 1) %>%
#  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
#  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
#  filter(first_vacc_between_infection_lc_response_28 == 0) %>%
#  filter(vacc1_1_to_13_days_before_infection == 0) %>%
#  mutate(case = ifelse(vacc1_14_days_before_infection == 1, 1, 0),
#         duration_between_inf_response_28_grouped = cut(duration_between_inf_response_28, seq(28, 119, 7), include.lowest = FALSE, right = FALSE),
#         duration_between_inf_response_28_grouped = ifelse(is.na(duration_between_inf_response_28_grouped), ">=119", duration_between_inf_response_28_grouped),
#         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
#  select(participant_id, case, age_at_infection, age_55_plus, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_28_date, duration_between_inf_response_28, duration_between_inf_response_28_grouped, covid_vaccine_date1, covid_vaccine_date2, lc_response_28, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_28_days_before_infection, vacc2_28_days_before_infection, first_vacc_between_infection_lc_response_28, vacc1_1_to_27_days_before_infection, vacc2_1_to_27_days_before_infection, second_vacc_between_infection_lc_response_28)

ggplot(dat_matching_28_3weeks, aes(x = age_at_infection, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Age (years)") +
  scale_fill_discrete("Case")

ggplot(dat_matching_28_3weeks, aes(x = duration_between_inf_response_28, fill = factor(case))) +
  geom_density(alpha = 0.4) +
  xlab("Days between infection & response") +
  scale_fill_discrete("Case")

##matching
matched_28_3weeks_age_sex <- matchit(case ~ age_group_10_year + sex + duration_between_inf_response_28, 
                                     data = dat_matching_28_3weeks,
                                     method = "nearest",
                                     exact = ~ age_group_10_year + sex + duration_between_inf_response_28) 
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

results_28 <- left_join(proportions_28, case_28_confints, by = c("case", "lc_response_28"))

##undadjusted 
mod_28_unadj <- glm(lc_response_28 ~ case, data = matched_28_data_3weeks_age_sex, family = binomial)
vcov_28_unadj <- vcovCL(mod_28_unadj, cluster = ~subclass)
output_28_unadj <- coeftest(mod_28_unadj, vcov = vcov_28_unadj)

##adjusted
mod_28_adj <- glm(lc_response_28 ~ case + age_at_infection + sex + health_conditions_at_infection + duration_between_inf_response_28 + imd_quintile + white, data = matched_28_data_3weeks_age_sex, family = binomial)
vcov_28_adj <- vcovCL(mod_28_adj, cluster = ~subclass)
output_28_adj <- coeftest(mod_28_adj, vcov = vcov_28_adj)

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

results_84 <- left_join(proportions_84, case_84_confints, by = c("case", "lc_response_84"))

##unadjusted
mod_84_unadj <- glm(lc_response_84 ~ case, data = matched_84_data_3weeks_age_sex, family = binomial)
vcov_84_unadj <- vcovCL(mod_84_unadj, cluster = ~subclass)
output_84_unadj <- coeftest(mod_84_unadj, vcov = vcov_84_unadj)

##adjusted
mod_84_adj <- glm(lc_response_84 ~ case + age_at_infection + sex + health_conditions_at_infection + duration_between_inf_response_84 + imd_quintile + white, data = matched_84_data_3weeks_age_sex, family = binomial)
vcov_84_adj <- vcovCL(mod_84_adj, cluster = ~subclass)
output_84_adj <- coeftest(mod_84_adj, vcov = vcov_84_adj)

###subgroup analysis excluding cases who received their second dose between infection and response
dat_matching_84_3weeks_subgroup <- dat_lc_vacc %>% filter(visit_84_after_infection == 1) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
  filter(first_vacc_between_infection_lc_response_84 == 0) %>%
  filter(is.na(second_vacc_between_infection_lc_response_84) | second_vacc_between_infection_lc_response_84 == 0) %>%
  filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = ifelse(vacc1_21_days_before_infection == 1, 1, 0),
         duration_between_inf_response_84_grouped = cut(duration_between_inf_response_84, seq(84, 196, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped = ifelse(is.na(duration_between_inf_response_84_grouped), ">=196", duration_between_inf_response_84_grouped),
         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
  select(participant_id, case, age_at_infection, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_84_date, duration_between_inf_response_84, duration_between_inf_response_84_grouped, covid_vaccine_date1, covid_vaccine_date2, lc_response_84, vacc1_1_days_before_infection, vacc1_21_days_before_infection, first_vacc_between_infection_lc_response_84, second_vacc_between_infection_lc_response_84)


dat_matching_28_3weeks_subgroup <- dat_lc_vacc %>% filter(visit_28_after_infection == 1) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(is.na(vacc2_1_days_before_infection) | vacc2_1_days_before_infection == 0) %>%
  filter(first_vacc_between_infection_lc_response_28 == 0) %>%
  filter(is.na(second_vacc_between_infection_lc_response_28) | second_vacc_between_infection_lc_response_28 == 0) %>%
  filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = ifelse(vacc1_21_days_before_infection == 1, 1, 0),
         duration_between_inf_response_28_grouped = cut(duration_between_inf_response_28, seq(28, 119, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_28_grouped = ifelse(is.na(duration_between_inf_response_28_grouped), ">=119", duration_between_inf_response_28_grouped),
         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
  select(participant_id, case, age_at_infection, age_group_10_year, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_28_date, duration_between_inf_response_28, duration_between_inf_response_28_grouped, covid_vaccine_date1, covid_vaccine_date2, lc_response_28, vacc1_1_days_before_infection, vacc1_21_days_before_infection, first_vacc_between_infection_lc_response_28, vacc1_1_to_20_days_before_infection, second_vacc_between_infection_lc_response_28)

##matching
matched_28_3weeks_age_sex_subgroup <- matchit(case ~ age_group_10_year + sex + duration_between_inf_response_28, 
                                     data = dat_matching_28_3weeks_subgroup,
                                     method = "nearest",
                                     exact = ~ age_group_10_year + sex + duration_between_inf_response_28) 
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

results_28 <- left_join(proportions_28, case_28_confints, by = c("case", "lc_response_28"))

##undadjusted 
mod_28_unadj <- glm(lc_response_28 ~ case, data = matched_28_data_3weeks_age_sex_subgroup, family = binomial)
vcov_28_unadj <- vcovCL(mod_28_unadj, cluster = ~subclass)
output_28_unadj <- coeftest(mod_28_unadj, vcov = vcov_28_unadj)

##adjusted
mod_28_adj <- glm(lc_response_28 ~ case + age_at_infection + sex + health_conditions_at_infection + duration_between_inf_response_28 + imd_quintile + white, data = matched_28_data_3weeks_age_sex_subgroup, family = binomial)
vcov_28_adj <- vcovCL(mod_28_adj, cluster = ~subclass)
output_28_adj <- coeftest(mod_28_adj, vcov = vcov_28_adj)


###fully vaccinated###
dat_matching_84_3weeks_2nd_dose <- dat_lc_vacc %>% filter(visit_84_after_infection == 1) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(first_vacc_between_infection_lc_response_84 == 0) %>%
  filter(second_vacc_between_infection_lc_response_84 == 0 | is.na(second_vacc_between_infection_lc_response_84)) %>%                                 
  filter(vacc2_1_to_20_days_before_infection == 0 | is.na(vacc2_1_to_20_days_before_infection)) %>%
  #filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = case_when(vacc2_21_days_before_infection == 1 ~ 1, 
                          vacc2_21_days_before_infection == 0 | is.na(vacc2_21_days_before_infection) ~ 0),
         duration_between_inf_response_84_grouped_7 = cut(duration_between_inf_response_84, seq(84, 196, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_7 = ifelse(is.na(duration_between_inf_response_84_grouped_7), ">=196", duration_between_inf_response_84_grouped_7),
         duration_between_inf_response_84_grouped_14 = cut(duration_between_inf_response_84, seq(84, 196, 14), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_14 = ifelse(is.na(duration_between_inf_response_84_grouped_14), ">=196", duration_between_inf_response_84_grouped_14),
         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
  filter(case == 1 | (case == 0 & vacc1_1_days_before_infection == 0)) %>%
  select(participant_id, case, age_at_infection, age_group_10_year, age_55_plus, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_84_date, duration_between_inf_response_84, duration_between_inf_response_84_grouped_7, duration_between_inf_response_84_grouped_14, covid_vaccine_date1, covid_vaccine_date2, lc_response_84, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_21_days_before_infection, vacc2_21_days_before_infection, first_vacc_between_infection_lc_response_84, vacc1_1_to_20_days_before_infection, vacc2_1_to_20_days_before_infection, second_vacc_between_infection_lc_response_84)

dat_matching_84_2weeks_2nd_dose <- dat_lc_vacc %>% filter(visit_84_after_infection == 1) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(first_vacc_between_infection_lc_response_84 == 0) %>%
  filter(second_vacc_between_infection_lc_response_84 == 0 | is.na(second_vacc_between_infection_lc_response_84)) %>%                                 
  filter(vacc2_1_to_13_days_before_infection == 0 | is.na(vacc2_1_to_13_days_before_infection)) %>%
  #filter(vacc1_1_to_13_days_before_infection == 0) %>%
  mutate(case = case_when(vacc2_14_days_before_infection == 1 ~ 1, 
                          vacc2_14_days_before_infection == 0 | is.na(vacc2_14_days_before_infection) ~ 0),
         duration_between_inf_response_84_grouped_7 = cut(duration_between_inf_response_84, seq(84, 196, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_7 = ifelse(is.na(duration_between_inf_response_84_grouped_7), ">=196", duration_between_inf_response_84_grouped_7),
         duration_between_inf_response_84_grouped_14 = cut(duration_between_inf_response_84, seq(84, 196, 14), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_14 = ifelse(is.na(duration_between_inf_response_84_grouped_14), ">=196", duration_between_inf_response_84_grouped_14),
         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
  filter(case == 1 | (case == 0 & vacc1_1_days_before_infection == 0)) %>%
  select(participant_id, case, age_at_infection, age_group_10_year, age_55_plus, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_84_date, duration_between_inf_response_84, duration_between_inf_response_84_grouped_7, duration_between_inf_response_84_grouped_14, covid_vaccine_date1, covid_vaccine_date2, lc_response_84, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_21_days_before_infection, vacc2_21_days_before_infection, first_vacc_between_infection_lc_response_84, vacc1_1_to_20_days_before_infection, vacc2_1_to_20_days_before_infection, second_vacc_between_infection_lc_response_84)


dat_matching_28_3weeks_2nd_dose <- dat_lc_vacc %>% filter(visit_28_after_infection == 1) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(first_vacc_between_infection_lc_response_28 == 0) %>%
  filter(second_vacc_between_infection_lc_response_28 == 0 | is.na(second_vacc_between_infection_lc_response_28)) %>%                                 
  filter(vacc2_1_to_20_days_before_infection == 0 | is.na(vacc2_1_to_20_days_before_infection)) %>%
  #filter(vacc1_1_to_20_days_before_infection == 0) %>%
  mutate(case = case_when(vacc2_21_days_before_infection == 1 ~ 1, 
                          vacc2_21_days_before_infection == 0 | is.na(vacc2_21_days_before_infection) ~ 0),
         #duration_between_inf_response_84_grouped_7 = cut(duration_between_inf_response_84, seq(84, 196, 7), include.lowest = FALSE, right = FALSE),
         #duration_between_inf_response_84_grouped_7 = ifelse(is.na(duration_between_inf_response_84_grouped_7), ">=196", duration_between_inf_response_84_grouped_7),
         #duration_between_inf_response_84_grouped_14 = cut(duration_between_inf_response_84, seq(84, 196, 14), include.lowest = FALSE, right = FALSE),
         #duration_between_inf_response_84_grouped_14 = ifelse(is.na(duration_between_inf_response_84_grouped_14), ">=196", duration_between_inf_response_84_grouped_14),
         health_conditions_at_infection = ifelse(is.na(health_conditions_at_infection), 0, health_conditions_at_infection)) %>%
  filter(case == 1 | (case == 0 & vacc1_1_days_before_infection == 0)) %>%
  select(participant_id, case, age_at_infection, age_group_10_year, age_55_plus, imd_quintile, ethnicityg, white, health_conditions_before_infection, health_conditions_at_infection, age_group, sex, study_swab_pos1_date, first_response_28_date, duration_between_inf_response_28, covid_vaccine_date1, covid_vaccine_date2, lc_response_28, vacc1_1_days_before_infection, vacc2_1_days_before_infection, vacc1_21_days_before_infection, vacc2_21_days_before_infection, first_vacc_between_infection_lc_response_28, vacc1_1_to_20_days_before_infection, vacc2_1_to_20_days_before_infection, second_vacc_between_infection_lc_response_28)


##matching
matched_28_3weeks_age_sex_2nd_dose <- matchit(case ~ age_group_10_year + sex + duration_between_inf_response_28, 
                                              data = dat_matching_28_3weeks_2nd_dose,
                                              method = "nearest",
                                              exact = ~ age_group_10_year + sex + duration_between_inf_response_28) 
summary(matched_28_3weeks_age_sex_2nd_dose)
matched_28_data_3weeks_age_sex_2nd_dose <- match.data(matched_28_3weeks_age_sex_2nd_dose)

matched_84_3weeks_age_sex_2nd_dose <- matchit(case ~ age_group_10_year + sex + duration_between_inf_response_84, 
                                              data = dat_matching_84_3weeks_2nd_dose,
                                              method = "nearest",
                                              exact = ~ age_group_10_year + sex + duration_between_inf_response_84) 
summary(matched_84_3weeks_age_sex_2nd_dose)
matched_84_data_3weeks_age_sex_2nd_dose <- match.data(matched_84_3weeks_age_sex_2nd_dose)


##calculate proportions
proportions_28 <- matched_28_data_3weeks_age_sex_2nd_dose %>% group_by(case, lc_response_28) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
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


proportions_84 <- matched_84_data_3weeks_age_sex_2nd_dose %>% group_by(case, lc_response_84) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
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

##undadjusted 
mod_28_unadj <- glm(lc_response_28 ~ case, data = matched_28_data_3weeks_age_sex_2nd_dose, family = binomial)
vcov_28_unadj <- vcovCL(mod_28_unadj, cluster = ~subclass)
output_28_unadj <- coeftest(mod_28_unadj, vcov = vcov_28_unadj)

##adjusted
mod_28_adj <- glm(lc_response_28 ~ case + age_at_infection + sex + health_conditions_at_infection + duration_between_inf_response_28 + imd_quintile + white, data = matched_28_data_3weeks_age_sex_2nd_dose, family = binomial)
vcov_28_adj <- vcovCL(mod_28_adj, cluster = ~subclass)
output_28_adj <- coeftest(mod_28_adj, vcov = vcov_28_adj)

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