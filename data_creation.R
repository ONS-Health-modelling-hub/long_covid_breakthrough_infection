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
dataset_date = "20211222"
cutoff_date = "2021-11-30"

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
dat_all <- dat_all[dat_all$visit_date <= cutoff_date,]
dat_all$visit_date <- as.numeric(as.Date(dat_all$visit_date))

### keep all visits belonging to participants who ever responded to LC question
part_select <- dat_all$participant_id[!is.na(dat_all$long_covid_have_symptoms) &
                                        dat_all$long_covid_have_symptoms %in% c(0, 1)]
part_select <- part_select[!duplicated(part_select)]
dat_all <- dat_all[dat_all$participant_id %in% part_select,]

### keep all visits belonging to participants who ever tested positive by swab in CIS or outside of study
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
dat_full$covid_vaccine_type1 <- as.numeric(dat_full$covid_vaccine_type1)
dat_full$covid_vaccine_type2 <- as.numeric(dat_full$covid_vaccine_type2)

rm(dat_vacc); gc()
rm(dat_all); gc()

#################################### JOIN COVID STATUS VARIABLES ####################################

### Create first suspected COVID date
dat_study_suspected1_date <- sqldf("
  select
    participant_id,
    min(visit_date) as study_suspected1_date
  from dat_full
  where covid_think_havehad=1
  group by participant_id
")

### merge first suspected COVID date
dat_full <- merge(x = dat_full,
                 y = dat_study_suspected1_date,
                 by.x = "participant_id",
                 by.y = "participant_id",
                 all.x = TRUE,
                 all.y = FALSE)

rm(dat_study_suspected1_date); gc()

### Create first long COVID date
dat_study_lc1_date <- sqldf("
  select
    participant_id,
    min(visit_date) as study_lc1_date
  from dat_full
  where long_covid_have_symptoms=1
  group by participant_id
")

### merge first long COVID date
dat_full <- merge(x = dat_full,
                 y = dat_study_lc1_date,
                 by.x = "participant_id",
                 by.y = "participant_id",
                 all.x = TRUE,
                 all.y = FALSE)

rm(dat_study_lc1_date); gc()

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

### Create first positive non-study swab test date
dat_non_study_swab_pos1_date <- sqldf("
  select
    participant_id,
    min(covid_test_swab_pos_first_date) as non_study_swab_pos1_date
  from dat_full
  group by participant_id
")

### merge first positive non-study swab test date
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

### create earliest swab and blood dates from study and non-study tests
dat_full <- dat_full %>% mutate(first_pos_swab_date = pmin(study_swab_pos1_date, non_study_swab_pos1_date, na.rm = T),
                                first_pos_blood_date = pmin(study_blood_pos1_date, non_study_blood_pos1_date, na.rm = T))

### coerce necessary variables to date format
dat_full$study_suspected1_date <- as.Date(dat_full$study_suspected1_date, origin = "1970-01-01")
dat_full$study_lc1_date <- as.Date(dat_full$study_lc1_date, origin = "1970-01-01")
dat_full$first_pos_swab_date <- as.Date(dat_full$first_pos_swab_date, origin = "1970-01-01")
dat_full$first_pos_blood_date <- as.Date(dat_full$first_pos_blood_date, origin = "1970-01-01")
dat_full$covid_date <- as.Date(dat_full$covid_date, origin = "1970-01-01")
dat_full$visit_date <- as.Date(dat_full$visit_date, origin = "1970-01-01")
dat_full$covid_vaccine_date1 <- as.Date(dat_full$covid_vaccine_date1, origin = "1970-01-01")
dat_full$covid_vaccine_date2 <- as.Date(dat_full$covid_vaccine_date2, origin = "1970-01-01")

### COVID-19 arrived in the UK on 24 Jan 2020 - set any dates before this to NA
dat_full$study_suspected1_date[dat_full$study_suspected1_date < as.Date("2020-01-24")] <- NA
dat_full$study_lc1_date[dat_full$study_lc1_date < as.Date("2020-01-24")] <- NA
dat_full$first_pos_swab_date[dat_full$first_pos_swab_date < as.Date("2020-01-24")] <- NA
dat_full$first_pos_blood_date[dat_full$first_pos_blood_date < as.Date("2020-01-24")] <- NA
dat_full$covid_date[dat_full$covid_date < as.Date("2020-01-24")] <- NA

### remove visits where first_pos_swab_date is unknown
dat_full <- dat_full[!is.na(dat_full$first_pos_swab_date),]

### Calculate time between first positive swab and first positive blood test, first suspected COVID, and first LC symptoms
dat_full <- dat_full %>% mutate(time_since_pos_blood = as.numeric(first_pos_swab_date - first_pos_blood_date),
                                time_since_suspected_covid = as.numeric(first_pos_swab_date - study_suspected1_date),
                                time_since_lc = as.numeric(first_pos_swab_date - study_lc1_date))

###remove people who had a positive blood test >14 days before first positive swab test
part_select_blood <- dat_full$participant_id[is.na(dat_full$time_since_pos_blood) | dat_full$time_since_pos_blood<14]
part_select_blood <- part_select_blood[!duplicated(part_select_blood)]
dat_full <- dat_full[dat_full$participant_id %in% part_select_blood,]

###remove people who first thought they had COVID >14 days before first positive swab test
part_select_suspected <- dat_full$participant_id[is.na(dat_full$time_since_suspected_covid) | dat_full$time_since_suspected_covid<14]
part_select_suspected <- part_select_suspected[!duplicated(part_select_suspected)]
dat_full <- dat_full[dat_full$participant_id %in% part_select_suspected,]

###remove people who had LC symptoms at any time before first positive swab test
part_select_lc <- dat_full$participant_id[is.na(dat_full$time_since_lc) | dat_full$time_since_lc<0]
part_select_lc <- part_select_lc[!duplicated(part_select_lc)]
dat_full <- dat_full[dat_full$participant_id %in% part_select_lc,]

### identify visits at least 12 weeks after first infection
dat_full <- dat_full %>% mutate(duration_after_infection = as.numeric(visit_date - first_pos_swab_date),
                                visit_84_after_infection = ifelse(duration_after_infection >= 84, 1, 0))

### derive vaccine manufacturer (first dose)
dat_full$vaccine_manufacturer <- ifelse(is.na(dat_full$covid_vaccine_type1), "None",
                                        ifelse(dat_full$covid_vaccine_type1==3, "Moderna",
                                               ifelse(dat_full$covid_vaccine_type1==4, "Oxford/AZ",
                                                      ifelse(dat_full$covid_vaccine_type1==5, "Pfizer/BioNTech",
                                                             "Other/unknown"))))

### derive vaccine type (first dose)
dat_full$vaccine_type <- ifelse(dat_full$vaccine_manufacturer=="Oxford/AZ", "Adenovirus vector",
                                ifelse(dat_full$vaccine_manufacturer %in% c("Pfizer/BioNTech", "Moderna"), "mRNA",
                                       ifelse(dat_full$vaccine_manufacturer=="None", "None", "Other/unknown")))

dat_full$vaccine_vector <- ifelse(dat_full$vaccine_type=="Adenovirus vector", 1, 0)

### flag if vaccinate manufacturer at second dose is different to at first dose
dat_full$vaccine_diff_manufacturer <- 0
dat_full$vaccine_diff_manufacturer[!is.na(dat_full$covid_vaccine_type1) &
                                     !is.na(dat_full$covid_vaccine_type2) &
                                     dat_full$covid_vaccine_type1!=dat_full$covid_vaccine_type2] <- 1

### exclude participants vaccinated with something other than AZ, Pfizer or Moderna,
### or those where the manufacturer at dose 2 was different to that at dose 1
dat_full <- dat_full %>% filter(vaccine_type!="Other/unknown" & vaccine_diff_manufacturer==0)

### aggregate to person level
dat_lc_status <- sqldf("
  select
    participant_id,
    max(long_covid_have_symptoms) as lc_ever,
    min(first_pos_swab_date) as infection_date
  from dat_full
  group by participant_id
")

dat_baseline_chars <- sqldf("
  select
    a.participant_id,
    age_at_visit,
    sex,
    imd_samp,
    ethnicityg,
    country,
    gor9d,
    health_conditions,
    health_conditions_impact,
    covid_vaccine_date1 as covid_vaccine_date1,
    covid_vaccine_date2 as covid_vaccine_date2,
    covid_vaccine_type1 as covid_vaccine_type1,
    covid_vaccine_type2 as covid_vaccine_type2,
    vaccine_manufacturer,
    vaccine_type,
    vaccine_vector
  from dat_full as a
  left join(
    select participant_id, min(visit_date) as enrolment_date
    from dat_full
    group by participant_id
  ) as b
  on a.participant_id = b.participant_id
  where visit_date = enrolment_date
")
dat_baseline_chars <- dat_baseline_chars[!duplicated(dat_baseline_chars$participant_id),]

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
    health_conditions as health_conditions_84,
    health_conditions_impact as health_conditions_impact_84
  from dat_full
  where visit_84_after_infection = 1 
  AND long_covid_have_symptoms is not null
  group by participant_id
")

dat_lc_baseline <- merge(x = dat_lc_status,
                         y = dat_baseline_chars,
                         by.x = "participant_id",
                         by.y = "participant_id",
                         all.x = TRUE,
                         all.y = FALSE)

dat_lc_status_date <- merge(x = dat_lc_baseline,
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

dat_infected_lc$imd_quintile <- as.factor(dat_infected_lc$imd_quintile)

### convert GOR to factor
dat_infected_lc$gor9d <- as.factor(dat_infected_lc$gor9d)
dat_infected_lc$gor9d <- relevel(dat_infected_lc$gor9d, ref="7")

### define white/non-white variable
dat_infected_lc$white <- ifelse(dat_infected_lc$ethnicityg==1, 1, 0)

### define health/disability status variable
dat_infected_lc$health_conditions[is.na(dat_infected_lc$health_conditions)] <- 0
dat_infected_lc$health_status <- 0
dat_infected_lc$health_status[dat_infected_lc$health_conditions==1 & dat_infected_lc$health_conditions_impact==0] <- 1
dat_infected_lc$health_status[dat_infected_lc$health_conditions==1 & dat_infected_lc$health_conditions_impact==1] <- 2
dat_infected_lc$health_status[dat_infected_lc$health_conditions==1 & dat_infected_lc$health_conditions_impact==2] <- 3
dat_infected_lc$health_status <- as.factor(dat_infected_lc$health_status)

### define activity limitation due to LC
dat_infected_lc$lc_activity_84_pooled <- ifelse(dat_infected_lc$lc_activity_84==4 | is.na(dat_infected_lc$lc_activity_84), 0, 1)

### define month/year of infection
dat_infected_lc$infection_month <- format(as.Date(dat_infected_lc$infection_date), "%Y-%m")

### define calendar time of infection (days since 24 Jan 2020)
dat_infected_lc$infection_day <- dat_infected_lc$infection_date - as.numeric(as.Date("2020-01-24"))

### derive infection period
dat_infected_lc$infection_period <- 1
dat_infected_lc$infection_period[as.Date(dat_infected_lc$infection_date) >= "2020-09-01"] <- 2
dat_infected_lc$infection_period[as.Date(dat_infected_lc$infection_date) >= "2021-05-01"] <- 3
dat_infected_lc$infection_period <- as.factor(dat_infected_lc$infection_period)

### derive Delta-dominant period
dat_infected_lc$delta_period <- 0
dat_infected_lc$delta_period[as.Date(dat_infected_lc$infection_date) >= "2021-05-17"] <- 1

### derive vaccinated before infection flags and remove people never vaccinated
dat_infected_lc <- dat_infected_lc %>% mutate(infection_date = as.Date(infection_date, origin = "1970-01-01"),
                                              first_response_84_date = as.Date(first_response_84_date, origin = "1970-01-01")) %>%
  #filter(!is.na(covid_vaccine_date1) & covid_vaccine_date1 <= cutoff_date) %>%
  filter(is.na(lc_onset_after_first_infection_date) | lc_onset_after_first_infection_date == 1)

dat_infected_lc <- dat_infected_lc %>%  mutate(diff_betweeen_vacc2_infection = as.Date(infection_date, origin = "1970-01-01") - covid_vaccine_date2,
                                               diff_betweeen_vacc1_infection = as.Date(infection_date, origin = "1970-01-01") - covid_vaccine_date1,
                                               vacc2_same_day_first_inf = ifelse(covid_vaccine_date2 == infection_date, 1, 0),
                                               duration_between_inf_response_84 = as.numeric(first_response_84_date - infection_date),
                                               first_vacc_between_infection_lc_response_84 = ifelse(covid_vaccine_date1 >= infection_date & covid_vaccine_date1 <= first_response_84_date, 1, 0),
                                               second_vacc_between_infection_lc_response_84 = ifelse(covid_vaccine_date2 >= infection_date & covid_vaccine_date2 <= first_response_84_date, 1, 0),
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
                                               lc_fu_time_84 = ifelse(diff_infection_lc_q_last_answered>=84, 1, 0),
                                               lc_yes_at_84_or_more = ifelse(diff_infection_lc_q_last_yes>=84, 1, 0))

### apply study sample inclusion criteria
dat_matching_84_2weeks_2nd_dose <- dat_infected_lc %>% filter(visit_84_after_infection == 1) %>%
  filter(age_at_visit >= 18 & age_at_visit < 70) %>%
  filter(is.na(vacc2_same_day_first_inf) | vacc2_same_day_first_inf == 0) %>%
  filter(first_vacc_between_infection_lc_response_84 == 0 | is.na(first_vacc_between_infection_lc_response_84)) %>%
  filter(second_vacc_between_infection_lc_response_84 == 0 | is.na(second_vacc_between_infection_lc_response_84)) %>%                                 
  filter(vacc2_1_to_13_days_before_infection == 0 | is.na(vacc2_1_to_13_days_before_infection)) %>%
  mutate(case = case_when(vacc2_14_days_before_infection == 1 ~ 1, 
                          vacc2_14_days_before_infection == 0 | is.na(vacc2_14_days_before_infection) ~ 0),
         duration_between_inf_response_84_grouped_7 = cut(duration_between_inf_response_84, seq(84, 238, 7), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_7 = ifelse(is.na(duration_between_inf_response_84_grouped_7), ">=238", duration_between_inf_response_84_grouped_7),
         duration_between_inf_response_84_grouped_14 = cut(duration_between_inf_response_84, seq(84, 238, 14), include.lowest = FALSE, right = FALSE),
         duration_between_inf_response_84_grouped_14 = ifelse(is.na(duration_between_inf_response_84_grouped_14), ">=238", duration_between_inf_response_84_grouped_14),
         health_conditions = ifelse(is.na(health_conditions), 0, health_conditions),
         age_group_10_year = case_when(age_at_visit >=18 & age_at_visit <= 29 ~ "18 to 29",
                                       age_at_visit >=30 & age_at_visit < 40 ~ "30 to 39",
                                       age_at_visit >=40 & age_at_visit < 50 ~ "40 to 49",
                                       age_at_visit >=50 & age_at_visit < 70 ~ "50 to 69"),
         age_40_plus = ifelse(age_at_visit>=40, 1, 0)) %>%
  filter(case == 1 | (case == 0 & vacc1_1_days_before_infection == 0)) %>%
  select(participant_id, case, age_at_visit, age_group_10_year, age_40_plus,
         imd_quintile, ethnicityg, white, health_conditions, health_status, sex, gor9d, country,
         infection_date, infection_day, first_response_84_date, duration_between_inf_response_84,
         duration_between_inf_response_84_grouped_7, duration_between_inf_response_84_grouped_14,
         covid_vaccine_date1, covid_vaccine_date2,
         lc_response_84, lc_activity_84_pooled,
         vacc1_1_days_before_infection, vacc2_1_days_before_infection,
         vacc1_21_days_before_infection, vacc2_21_days_before_infection,
         first_vacc_between_infection_lc_response_84, vacc1_1_to_20_days_before_infection,
         vacc2_1_to_20_days_before_infection, second_vacc_between_infection_lc_response_84,
         vaccine_manufacturer, vaccine_type, vaccine_vector)

### write out dataset
save(dat_matching_84_2weeks_2nd_dose, file=paste0(out_dir, "\\dataset.RData"))
