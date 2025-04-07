library(MatchIt)
library(splines)
library(car)

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

### define outcomes, exposure, covariates and effect-modifiers
outcomes = c(
  "lc_response_84",
  "lc_activity_84_pooled"
)

exposure = "case"

covariates = c(
  "ns(duration_between_inf_response_84, df=2, Boundary.knots=quantile(duration_between_inf_response_84, c(.10, .90)))",
  "ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(.10, .90)))",
  "sex",
  "white",
  "gor9d",
  "imd_quintile",
  "health_status"
)

modifiers <- c(
  "ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(.10, .90)))",
  "sex",
  "white",
  "imd_quintile",
  "health_status"
)

pvalue_list <- as.list(NULL)

### loop over outcomes
for(j in 1:length(outcomes)) {

  ### select outcome
  outcome <- outcomes[j]

  ### set up string for formula including main effects only
  formula_main <- paste0(outcome, " ~ ", exposure, " + ", paste(covariates, collapse=" + "))

  pvalues <- NULL

  ### loop over effect-modifiers
  for(i in 1:length(modifiers)) {

    ### select effect-modifier
    modifier <- modifiers[i]

    ### add effect-modifier to formula string
    formula_modifier <- paste0(formula_main, " + ", exposure, ":", modifier)

    ### fit model
    mod <- glm(as.formula(formula_modifier), data = matched_84_data_2weeks_age_sex_2nd_dose, family = binomial)

    ### extract p-value from likelihood ratio test of effect-modifier
    lr_test <- as.data.frame(Anova(mod))
    pvalues[i] <- lr_test[nrow(lr_test), 3]

  }

  pvalue_list[[j]] <- pvalues

}

### coerce listed p-values to data.frame
out_df <- data.frame(pvalue_list)
colnames(out_df) <- outcomes
rownames(out_df) <- modifiers

### write to working directory
write.csv(out_df, file=paste0(out_dir, "\\het_effects_pvalues.csv"))
