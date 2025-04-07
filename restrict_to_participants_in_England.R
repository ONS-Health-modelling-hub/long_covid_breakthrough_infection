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

infile = "filepath\\dataset.RData"
out_dir = "filepath"

### read in dataset
load(infile)

### matching
matched_84_2weeks_age_sex_2nd_dose <- matchit(case ~
                                                ns(duration_between_inf_response_84, df=2, Boundary.knots=quantile(duration_between_inf_response_84, c(.10, .90)))
                                              + ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(.10, .90)))
                                              + sex
                                              + white
                                              + gor9d
                                              + imd_quintile
                                              + health_status,
                                              data = dat_matching_84_2weeks_2nd_dose,
                                              method = "nearest",
                                              discard = "none",
                                              caliper = 0.1,
                                              std.caliper = FALSE)

matched_84_data_2weeks_age_sex_2nd_dose <- match.data(matched_84_2weeks_age_sex_2nd_dose)

### restrict sample to participants where both the matched case and control live in England
cases <- matched_84_data_2weeks_age_sex_2nd_dose[matched_84_data_2weeks_age_sex_2nd_dose$case==1, c("subclass", "country")]
controls <- matched_84_data_2weeks_age_sex_2nd_dose[matched_84_data_2weeks_age_sex_2nd_dose$case==0, c("subclass", "country")]
matched <- merge(x=cases, y=controls, by="subclass", all.x=TRUE, all.y=TRUE)
matched$both_eng <- ifelse(matched$country.x==0 & matched$country.y==0, 1, 0)
eng_id <- matched[matched$both_eng==1, "subclass"]
matched_84_data_2weeks_age_sex_2nd_dose <- matched_84_data_2weeks_age_sex_2nd_dose[matched_84_data_2weeks_age_sex_2nd_dose$subclass %in% eng_id,]
matched_84_data_2weeks_age_sex_2nd_dose$gor9d <- droplevels(matched_84_data_2weeks_age_sex_2nd_dose$gor9d)

### main estimates - long COVID of any severity

# calculate proportions
#proportions_84 <- matched_84_data_2weeks_age_sex_2nd_dose %>% group_by(case, lc_response_84) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
#write.csv(proportions_84, paste0(out_dir, "\\two_dose_proportions_84.csv"), row.names = F)

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

# unadjusted modelled estimates
mod_84_unadj_2_dose <- glm(lc_response_84 ~ case,
                           data = matched_84_data_2weeks_age_sex_2nd_dose,
                           family = binomial)

vcov_84_unadj_2_dose <- vcovCL(mod_84_unadj_2_dose, cluster = ~subclass)
output_84_unadj_2_dose <- coeftest(mod_84_unadj_2_dose, vcov = vcov_84_unadj_2_dose)

unadjust_84_results_2_dose <- data.frame(model = "unadjusted",
                                         variable = rownames(output_84_unadj_2_dose),
                                         OR = exp(output_84_unadj_2_dose[,1]),
                                         lower_ci = exp(output_84_unadj_2_dose[,1] - 1.96*output_84_unadj_2_dose[,2]),
                                         upper_ci = exp(output_84_unadj_2_dose[,1] + 1.96*output_84_unadj_2_dose[,2]),
                                         p = output_84_unadj_2_dose[,4],
                                         row.names = NULL,
                                         stringsAsFactors = FALSE)

# adjusted modelled estimates
mod_84_adj_2_dose <- glm(lc_response_84 ~ case
                         + ns(duration_between_inf_response_84, df=2, Boundary.knots=quantile(duration_between_inf_response_84, c(.10, .90)))
                         + ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(.10, .90)))
                         + sex
                         + white
                         + gor9d
                         + imd_quintile
                         + health_status,
                         data = matched_84_data_2weeks_age_sex_2nd_dose,
                         family = binomial)

vcov_84_adj_2_dose <- vcovCL(mod_84_adj_2_dose, cluster = ~subclass)
output_84_adj_2_dose <- coeftest(mod_84_adj_2_dose, vcov = vcov_84_adj_2_dose)

adjust_84_results_2_dose <- data.frame(model = "adjusted",
                                       variable = rownames(output_84_adj_2_dose),
                                       OR = exp(output_84_adj_2_dose[,1]),
                                       lower_ci = exp(output_84_adj_2_dose[,1] - 1.96*output_84_adj_2_dose[,2]),
                                       upper_ci = exp(output_84_adj_2_dose[,1] + 1.96*output_84_adj_2_dose[,2]),
                                       p = output_84_adj_2_dose[,4],
                                       row.names = NULL,
                                       stringsAsFactors = FALSE)

write.csv(rbind(unadjust_84_results_2_dose, adjust_84_results_2_dose), paste0(out_dir, "\\two_dose_OR_84.csv"), row.names = F) 

### main estimates - activity-limiting long COVID

# calculate proportions
#proportions_84 <- matched_84_data_2weeks_age_sex_2nd_dose %>% group_by(case, lc_activity_84_pooled) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
#write.csv(proportions_84, paste0(out_dir, "\\two_dose_proportions_84_lim.csv"), row.names = F)

proportions_84 <- matched_84_data_2weeks_age_sex_2nd_dose %>% group_by(case, lc_activity_84_pooled) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
case_84_counts <- proportions_84 %>% filter(case == 1) %>% select(n)
control_84_counts <- proportions_84 %>% filter(case == 0) %>% select(n)

case_84_ci_lc_yes <- scoreci(x = case_84_counts[2,2], n = sum(case_84_counts[,2]), conf.level = 0.95)
case_84_ci_lc_no <- scoreci(x = case_84_counts[1,2], n = sum(case_84_counts[,2]), conf.level = 0.95)

case_84_confints_lc_yes <- data.frame(case_84_ci_lc_yes$conf.int)
colnames(case_84_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
case_84_confints_lc_yes <- case_84_confints_lc_yes %>% mutate(case = 1, lc_activity_84_pooled = 1)

case_84_confints_lc_no <- data.frame(case_84_ci_lc_no$conf.int)
colnames(case_84_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
case_84_confints_lc_no <- case_84_confints_lc_no %>% mutate(case = 1, lc_activity_84_pooled = 0)

control_84_ci_lc_yes <- scoreci(x = control_84_counts[2,2], n = sum(control_84_counts[,2]), conf.level = 0.95)
control_84_ci_lc_no <- scoreci(x = control_84_counts[1,2], n = sum(control_84_counts[,2]), conf.level = 0.95)

control_84_confints_lc_yes <- data.frame(control_84_ci_lc_yes$conf.int)
colnames(control_84_confints_lc_yes) = c("lower_95_ci", "upper_95_ci")
control_84_confints_lc_yes <- control_84_confints_lc_yes %>% mutate(case = 0, lc_activity_84_pooled = 1)

control_84_confints_lc_no <- data.frame(control_84_ci_lc_no$conf.int)
colnames(control_84_confints_lc_no) = c("lower_95_ci", "upper_95_ci")
control_84_confints_lc_no <- control_84_confints_lc_no %>% mutate(case = 0, lc_activity_84_pooled = 0)

case_84_confints <- rbind(case_84_confints_lc_no, case_84_confints_lc_yes, control_84_confints_lc_yes, control_84_confints_lc_no)

results_84_fully_vacc <- left_join(proportions_84, case_84_confints, by = c("case", "lc_activity_84_pooled"))
write.csv(results_84_fully_vacc, paste0(out_dir, "\\two_dose_proportions_84_lim.csv"), row.names = F)

# unadjusted modelled estimates
mod_84_unadj_2_dose_lim <- glm(lc_activity_84_pooled ~ case,
                           data = matched_84_data_2weeks_age_sex_2nd_dose,
                           family = binomial)

vcov_84_unadj_2_dose_lim <- vcovCL(mod_84_unadj_2_dose_lim, cluster = ~subclass)
output_84_unadj_2_dose_lim <- coeftest(mod_84_unadj_2_dose_lim, vcov = vcov_84_unadj_2_dose_lim)

unadjust_84_results_2_dose_lim <- data.frame(model = "unadjusted",
                                         variable = rownames(output_84_unadj_2_dose_lim),
                                         OR = exp(output_84_unadj_2_dose_lim[,1]),
                                         lower_ci = exp(output_84_unadj_2_dose_lim[,1] - 1.96*output_84_unadj_2_dose_lim[,2]),
                                         upper_ci = exp(output_84_unadj_2_dose_lim[,1] + 1.96*output_84_unadj_2_dose_lim[,2]),
                                         p = output_84_unadj_2_dose_lim[,4],
                                         row.names = NULL,
                                         stringsAsFactors = FALSE)

# adjusted modelled estimates
mod_84_adj_2_dose_lim <- glm(lc_activity_84_pooled ~ case
                             + ns(duration_between_inf_response_84, df=2, Boundary.knots=quantile(duration_between_inf_response_84, c(.10, .90)))
                             + ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(.10, .90)))
                             + sex
                             + white
                             + gor9d
                             + imd_quintile
                             + health_status,
                             data = matched_84_data_2weeks_age_sex_2nd_dose,
                             family = binomial)

vcov_84_adj_2_dose_lim <- vcovCL(mod_84_adj_2_dose_lim, cluster = ~subclass)
output_84_adj_2_dose_lim <- coeftest(mod_84_adj_2_dose_lim, vcov = vcov_84_adj_2_dose_lim)

adjust_84_results_2_dose_lim <- data.frame(model = "adjusted",
                                       variable = rownames(output_84_adj_2_dose_lim),
                                       OR = exp(output_84_adj_2_dose_lim[,1]),
                                       lower_ci = exp(output_84_adj_2_dose_lim[,1] - 1.96*output_84_adj_2_dose_lim[,2]),
                                       upper_ci = exp(output_84_adj_2_dose_lim[,1] + 1.96*output_84_adj_2_dose_lim[,2]),
                                       p = output_84_adj_2_dose_lim[,4],
                                       row.names = NULL,
                                       stringsAsFactors = FALSE)

write.csv(rbind(unadjust_84_results_2_dose_lim, adjust_84_results_2_dose_lim), paste0(out_dir, "\\two_dose_OR_84_lim.csv"), row.names = F) 

### effect modification by vaccine type - long COVID of any severity

#unadjusted modelled estimates
mod_84_unadj_2_dose_int <- glm(lc_response_84 ~ case*vaccine_vector,
                               data = matched_84_data_2weeks_age_sex_2nd_dose,
                               family = binomial)

vcov_84_unadj_2_dose_int <- vcovCL(mod_84_unadj_2_dose_int, cluster = ~subclass)
output_84_unadj_2_dose_int <- coeftest(mod_84_unadj_2_dose_int, vcov = vcov_84_unadj_2_dose_int)

unadjust_84_results_2_dose_int <- data.frame(model = "unadjusted",
                                             variable = rownames(output_84_unadj_2_dose_int),
                                             OR = exp(output_84_unadj_2_dose_int[,1]),
                                             lower_ci = exp(output_84_unadj_2_dose_int[,1] - 1.96*output_84_unadj_2_dose_int[,2]),
                                             upper_ci = exp(output_84_unadj_2_dose_int[,1] + 1.96*output_84_unadj_2_dose_int[,2]),
                                             p = output_84_unadj_2_dose_int[,4],
                                             row.names = NULL,
                                             stringsAsFactors = FALSE)

#adjusted modelled estimates
mod_84_adj_2_dose_int <- glm(lc_response_84 ~ case*vaccine_vector 
                             + ns(duration_between_inf_response_84, df=2, Boundary.knots=quantile(duration_between_inf_response_84, c(.10, .90)))
                             + ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(.10, .90)))
                             + sex
                             + white
                             + gor9d
                             + imd_quintile
                             + health_status,
                             data = matched_84_data_2weeks_age_sex_2nd_dose,
                             family = binomial)

vcov_84_adj_2_dose_int <- vcovCL(mod_84_adj_2_dose_int, cluster = ~subclass)
output_84_adj_2_dose_int <- coeftest(mod_84_adj_2_dose_int, vcov = vcov_84_adj_2_dose_int)

adjust_84_results_2_dose_int <- data.frame(model = "adjusted",
                                           variable = rownames(output_84_adj_2_dose_int),
                                           OR = exp(output_84_adj_2_dose_int[,1]),
                                           lower_ci = exp(output_84_adj_2_dose_int[,1] - 1.96*output_84_adj_2_dose_int[,2]),
                                           upper_ci = exp(output_84_adj_2_dose_int[,1] + 1.96*output_84_adj_2_dose_int[,2]),
                                           p = output_84_adj_2_dose_int[,4],
                                           row.names = NULL,
                                           stringsAsFactors = FALSE)

write.csv(rbind(unadjust_84_results_2_dose_int, adjust_84_results_2_dose_int), paste0(out_dir, "\\two_dose_OR_84_interaction.csv"), row.names = F)

#data for unadjusted OR
coeff_comb_unadj <- output_84_unadj_2_dose["case",1]
coeff_mrna_unadj <- output_84_unadj_2_dose_int["case",1]
coeff_vector_unadj <- output_84_unadj_2_dose_int["case",1] + output_84_unadj_2_dose_int["case:vaccine_vector",1]

or_comb_unadj <- exp(coeff_comb_unadj)
or_mrna_unadj <- exp(coeff_mrna_unadj)
or_vector_unadj <- exp(coeff_vector_unadj)

se_comb_unadj <- sqrt(vcov_84_unadj_2_dose["case","case"])
se_mrna_unadj <- sqrt(vcov_84_unadj_2_dose_int["case","case"])
se_vector_unadj <- sqrt(vcov_84_unadj_2_dose_int["case","case"] +
                          vcov_84_unadj_2_dose_int["case:vaccine_vector","case:vaccine_vector"] +
                          2*vcov_84_unadj_2_dose_int["case","case:vaccine_vector"])

lcl_comb_unadj <- exp(coeff_comb_unadj - 1.96*se_comb_unadj)
lcl_mrna_unadj <- exp(coeff_mrna_unadj - 1.96*se_mrna_unadj)
lcl_vector_unadj <- exp(coeff_vector_unadj - 1.96*se_vector_unadj)

ucl_comb_unadj <- exp(coeff_comb_unadj + 1.96*se_comb_unadj)
ucl_mrna_unadj <- exp(coeff_mrna_unadj + 1.96*se_mrna_unadj)
ucl_vector_unadj <- exp(coeff_vector_unadj + 1.96*se_vector_unadj)

or_out_unadj <- data.frame(
  model = "Unadjusted",
  vaccine = c("Combined", "mRNA", "Adenovirus vector"),
  or = c(or_comb_unadj, or_mrna_unadj, or_vector_unadj),
  lcl = c(lcl_comb_unadj, lcl_mrna_unadj, lcl_vector_unadj),
  ucl = c(ucl_comb_unadj, ucl_mrna_unadj, ucl_vector_unadj)
)

#data for adjusted OR
coeff_comb_adj <- output_84_adj_2_dose["case",1]
coeff_mrna_adj <- output_84_adj_2_dose_int["case",1]
coeff_vector_adj <- output_84_adj_2_dose_int["case",1] + output_84_adj_2_dose_int["case:vaccine_vector",1]

or_comb_adj <- exp(coeff_comb_adj)
or_mrna_adj <- exp(coeff_mrna_adj)
or_vector_adj <- exp(coeff_vector_adj)

se_comb_adj <- sqrt(vcov_84_adj_2_dose["case","case"])
se_mrna_adj <- sqrt(vcov_84_adj_2_dose_int["case","case"])
se_vector_adj <- sqrt(vcov_84_adj_2_dose_int["case","case"] +
                        vcov_84_adj_2_dose_int["case:vaccine_vector","case:vaccine_vector"] +
                        2*vcov_84_adj_2_dose_int["case","case:vaccine_vector"])

lcl_comb_adj <- exp(coeff_comb_adj - 1.96*se_comb_adj)
lcl_mrna_adj <- exp(coeff_mrna_adj - 1.96*se_mrna_adj)
lcl_vector_adj <- exp(coeff_vector_adj - 1.96*se_vector_adj)

ucl_comb_adj <- exp(coeff_comb_adj + 1.96*se_comb_adj)
ucl_mrna_adj <- exp(coeff_mrna_adj + 1.96*se_mrna_adj)
ucl_vector_adj <- exp(coeff_vector_adj + 1.96*se_vector_adj)

or_out_adj <- data.frame(
  model = "Adjusted",
  vaccine = c("Combined", "mRNA", "Adenovirus vector"),
  or = c(or_comb_adj, or_mrna_adj, or_vector_adj),
  lcl = c(lcl_comb_adj, lcl_mrna_adj, lcl_vector_adj),
  ucl = c(ucl_comb_adj, ucl_mrna_adj, ucl_vector_adj)
)

#combined data
chart_data <- rbind(or_out_unadj, or_out_adj)
write.csv(chart_data, paste0(out_dir, "\\OR_chart_data.csv"), row.names = F) 

### effect modification by vaccine type - activity-limiting long COVID

#unadjusted modelled estimates
mod_84_unadj_2_dose_int <- glm(lc_activity_84_pooled ~ case*vaccine_vector,
                               data = matched_84_data_2weeks_age_sex_2nd_dose,
                               family = binomial)

vcov_84_unadj_2_dose_int <- vcovCL(mod_84_unadj_2_dose_int, cluster = ~subclass)
output_84_unadj_2_dose_int <- coeftest(mod_84_unadj_2_dose_int, vcov = vcov_84_unadj_2_dose_int)

unadjust_84_results_2_dose_int <- data.frame(model = "unadjusted",
                                             variable = rownames(output_84_unadj_2_dose_int),
                                             OR = exp(output_84_unadj_2_dose_int[,1]),
                                             lower_ci = exp(output_84_unadj_2_dose_int[,1] - 1.96*output_84_unadj_2_dose_int[,2]),
                                             upper_ci = exp(output_84_unadj_2_dose_int[,1] + 1.96*output_84_unadj_2_dose_int[,2]),
                                             p = output_84_unadj_2_dose_int[,4],
                                             row.names = NULL,
                                             stringsAsFactors = FALSE)

#adjusted modelled estimates
mod_84_adj_2_dose_int <- glm(lc_activity_84_pooled ~ case*vaccine_vector 
                             + ns(duration_between_inf_response_84, df=2, Boundary.knots=quantile(duration_between_inf_response_84, c(.10, .90)))
                             + ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(.10, .90)))
                             + sex
                             + white
                             + gor9d
                             + imd_quintile
                             + health_status,
                             data = matched_84_data_2weeks_age_sex_2nd_dose,
                             family = binomial)

vcov_84_adj_2_dose_int <- vcovCL(mod_84_adj_2_dose_int, cluster = ~subclass)
output_84_adj_2_dose_int <- coeftest(mod_84_adj_2_dose_int, vcov = vcov_84_adj_2_dose_int)

adjust_84_results_2_dose_int <- data.frame(model = "adjusted",
                                           variable = rownames(output_84_adj_2_dose_int),
                                           OR = exp(output_84_adj_2_dose_int[,1]),
                                           lower_ci = exp(output_84_adj_2_dose_int[,1] - 1.96*output_84_adj_2_dose_int[,2]),
                                           upper_ci = exp(output_84_adj_2_dose_int[,1] + 1.96*output_84_adj_2_dose_int[,2]),
                                           p = output_84_adj_2_dose_int[,4],
                                           row.names = NULL,
                                           stringsAsFactors = FALSE)

write.csv(rbind(unadjust_84_results_2_dose_int, adjust_84_results_2_dose_int), paste0(out_dir, "\\two_dose_OR_84_interaction_lim.csv"), row.names = F)

#data for unadjusted OR
coeff_comb_unadj <- output_84_unadj_2_dose_lim["case",1]
coeff_mrna_unadj <- output_84_unadj_2_dose_int["case",1]
coeff_vector_unadj <- output_84_unadj_2_dose_int["case",1] + output_84_unadj_2_dose_int["case:vaccine_vector",1]

or_comb_unadj <- exp(coeff_comb_unadj)
or_mrna_unadj <- exp(coeff_mrna_unadj)
or_vector_unadj <- exp(coeff_vector_unadj)

se_comb_unadj <- sqrt(vcov_84_unadj_2_dose_lim["case","case"])
se_mrna_unadj <- sqrt(vcov_84_unadj_2_dose_int["case","case"])
se_vector_unadj <- sqrt(vcov_84_unadj_2_dose_int["case","case"] +
                          vcov_84_unadj_2_dose_int["case:vaccine_vector","case:vaccine_vector"] +
                          2*vcov_84_unadj_2_dose_int["case","case:vaccine_vector"])

lcl_comb_unadj <- exp(coeff_comb_unadj - 1.96*se_comb_unadj)
lcl_mrna_unadj <- exp(coeff_mrna_unadj - 1.96*se_mrna_unadj)
lcl_vector_unadj <- exp(coeff_vector_unadj - 1.96*se_vector_unadj)

ucl_comb_unadj <- exp(coeff_comb_unadj + 1.96*se_comb_unadj)
ucl_mrna_unadj <- exp(coeff_mrna_unadj + 1.96*se_mrna_unadj)
ucl_vector_unadj <- exp(coeff_vector_unadj + 1.96*se_vector_unadj)

or_out_unadj <- data.frame(
  model = "Unadjusted",
  vaccine = c("Combined", "mRNA", "Adenovirus vector"),
  or = c(or_comb_unadj, or_mrna_unadj, or_vector_unadj),
  lcl = c(lcl_comb_unadj, lcl_mrna_unadj, lcl_vector_unadj),
  ucl = c(ucl_comb_unadj, ucl_mrna_unadj, ucl_vector_unadj)
)

#data for adjusted OR
coeff_comb_adj <- output_84_adj_2_dose_lim["case",1]
coeff_mrna_adj <- output_84_adj_2_dose_int["case",1]
coeff_vector_adj <- output_84_adj_2_dose_int["case",1] + output_84_adj_2_dose_int["case:vaccine_vector",1]

or_comb_adj <- exp(coeff_comb_adj)
or_mrna_adj <- exp(coeff_mrna_adj)
or_vector_adj <- exp(coeff_vector_adj)

se_comb_adj <- sqrt(vcov_84_adj_2_dose_lim["case","case"])
se_mrna_adj <- sqrt(vcov_84_adj_2_dose_int["case","case"])
se_vector_adj <- sqrt(vcov_84_adj_2_dose_int["case","case"] +
                        vcov_84_adj_2_dose_int["case:vaccine_vector","case:vaccine_vector"] +
                        2*vcov_84_adj_2_dose_int["case","case:vaccine_vector"])

lcl_comb_adj <- exp(coeff_comb_adj - 1.96*se_comb_adj)
lcl_mrna_adj <- exp(coeff_mrna_adj - 1.96*se_mrna_adj)
lcl_vector_adj <- exp(coeff_vector_adj - 1.96*se_vector_adj)

ucl_comb_adj <- exp(coeff_comb_adj + 1.96*se_comb_adj)
ucl_mrna_adj <- exp(coeff_mrna_adj + 1.96*se_mrna_adj)
ucl_vector_adj <- exp(coeff_vector_adj + 1.96*se_vector_adj)

or_out_adj <- data.frame(
  model = "Adjusted",
  vaccine = c("Combined", "mRNA", "Adenovirus vector"),
  or = c(or_comb_adj, or_mrna_adj, or_vector_adj),
  lcl = c(lcl_comb_adj, lcl_mrna_adj, lcl_vector_adj),
  ucl = c(ucl_comb_adj, ucl_mrna_adj, ucl_vector_adj)
)

#combined data
chart_data <- rbind(or_out_unadj, or_out_adj)
write.csv(chart_data, paste0(out_dir, "\\OR_chart_data_lim.csv"), row.names = F) 
