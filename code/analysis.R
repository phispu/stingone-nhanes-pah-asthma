########################################################################################################################################################
#
# Urinary metabolites of polycyclic aromatic hydrocarbons and short-acting beta agonist or systemic corticosteroid asthma medication use within NHANES
# Analysis
# Date: 4/27/2022
# Programmer: Stephen Uong
# Contributors: Jeanette Stingone
#
########################################################################################################################################################


####################################
########## LOAD LIBRARIES ########## 
####################################
library(rio)
library(tidyverse)
library(survey)
library(arsenal)
library(zoo) # rollmeans
library(here)
library(janitor)
library(broom)
#################################
########## IMPORT DATA ########## 
#################################
setwd(here("data/processed"))
pahas <- readRDS('nhanes-pah-asthma-analysis.RDS')
setwd(here('data'))

#############################
########## TABLE 1 ##########
#############################
### Weights / Svy Design
pahas_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_lab, nest = TRUE, data = pahas)
### Set up functions to generate freqs/percentages for table 1
n_decimal <- function(x, n) format(round(x, n), nsmall = n, big.mark=",")
# Overall
svy_cat <- function(var){
  f <- paste("~", var)
  perc <- do.call("svytable", list(as.formula(f), design = pahas_design)) %>% prop.table() %>% "*"(100) %>% n_decimal(1) 
  n <- pahas %>% janitor::tabyl(var, show_na = FALSE) %>% dplyr::select(n) %>% unlist() %>% as.vector() %>% formatC(format="d", big.mark=",")
  paste0(n, " (", trimws(perc), "%)") %>% return()
}
svy_cat_yn <- function(var){
  f <- paste("~", var)
  perc <- do.call("svytable", list(as.formula(f), design = pahas_design)) %>% prop.table() %>% "*"(100) %>% n_decimal(1) %>% magrittr::extract(2) 
  n <- pahas %>% janitor::tabyl(var, show_na = FALSE) %>% dplyr::select(n) %>% unlist() %>% as.vector() %>% magrittr::extract(2) %>% formatC(format="d", big.mark=",")
  paste0(n, " (", trimws(perc), "%)") %>% return()
}
svy_num <- function(var, ndec_mean, ndec_med){
  f <- paste("~", var)
  results_mean <- do.call("svymean", list(as.formula(f), design = pahas_design, na.rm = T)) %>% as.data.frame()
  results_var <- do.call("svyvar", list(as.formula(f), design = pahas_design, na.rm = T)) %>% as.data.frame()
  results_mean_clean <- results_mean %>% dplyr::transmute(mean = mean %>% n_decimal(ndec_mean)) 
  results_sd_clean <- results_var %>% dplyr::transmute(sd = variance %>% sqrt() %>% n_decimal(ndec_mean)) 
  mean_sd <- cbind(results_mean_clean, results_sd_clean) %>% dplyr::transmute(mean_sd = paste0(mean, " (", sd, ")"))
  mean_se <- paste0(results_mean[1] %>% n_decimal(ndec_mean), " (", results_mean[2] %>% n_decimal(ndec_mean), ")") 
  results_median <- do.call("svyquantile", list(as.formula(f), design = pahas_design, na.rm = T, quantile = c(0.25, 0.50, 0.75, ci = TRUE))) %>% as.data.frame()
  median_iqr <- paste0(results_median[2] %>% n_decimal(ndec_med), " (", results_median[1] %>% n_decimal(ndec_med), "-", results_median[3] %>% n_decimal(ndec_med), ")")
  cbind(mean_sd, median_iqr) %>% return()
}
# By PAH Quartile
svy_cat_q <- function(var){
  f <- paste("~", var, " + ratio_pahq")
  var_quo = enquo(var)
  perc <- do.call("svytable", list(as.formula(f), design = pahas_design)) %>% prop.table(2) %>% "*"(100) %>% n_decimal(1)
  n <- pahas %>% janitor::tabyl(!!sym(var), ratio_pahq, show_na = FALSE) %>% dplyr::mutate_if(is.numeric, format, big.mark = ",") %>% magrittr::extract(2:5) %>% as.matrix()
  matrix(paste0(n, " (", trimws(perc), "%)"), nrow = nrow(perc)) %>% return()
}
svy_cat_yn_q <- function(var){
  f <- paste("~", var, " + ratio_pahq")
  perc <- do.call("svytable", list(as.formula(f), design = pahas_design)) %>% prop.table(2) %>% tail(1) %>% "*"(100) %>% n_decimal(1)
  n <- pahas %>% janitor::tabyl(!!sym(var), ratio_pahq, show_na = FALSE) %>% dplyr::mutate_if(is.numeric, format, big.mark = ",") %>% magrittr::extract(2:5) %>% slice(2) %>% as.matrix()
  paste0(n, " (", trimws(perc), "%)") %>% return()
}
svy_num_q <- function(var, ndec_mean, ndec_med){
  f <- paste("~", var)
  var <- sym(var)
  results_mean <- do.call("svyby", list(as.formula(f), by = ~ ratio_pahq, design = pahas_design, na.rm = T, svymean)) %>% as.data.frame()
  results_var <- do.call("svyby", list(as.formula(f), by = ~ ratio_pahq, design = pahas_design, na.rm = T, svyvar)) %>% as.data.frame()
  results_mean_clean <- results_mean %>% dplyr::transmute(mean = !!var %>% n_decimal(ndec_mean)) 
  results_sd_clean <- results_var %>% dplyr::transmute(sd = !!var %>% sqrt() %>% n_decimal(ndec_mean)) 
  mean_sd <- cbind(results_mean_clean, results_sd_clean) %>% dplyr::transmute(mean_sd = paste0(mean, " (", sd, ")"))
  results_median <- do.call("svyby", list(as.formula(f), by = ~ ratio_pahq, design = pahas_design, na.rm = T, quantile = c(0.25, 0.50, 0.75), ci = TRUE, svyquantile)) %>% as.data.frame()
  median_iqr <- results_median %>% dplyr::transmute(median_iqr = paste0(`0.5` %>% n_decimal(ndec_med), " (", `0.25` %>% n_decimal(ndec_med), "-", `0.75` %>% n_decimal(ndec_med), ")"))
  cbind(mean_sd, median_iqr) %>% t() %>% as.data.frame() %>% dplyr::rename(Q1 = 1, Q2 = 2, Q3 = 3, Q4 = 4) %>% return()
}


### Table, overall
# Prep columnns of table
tabo_yr <- svy_cat("yr")
tabo_age <- svy_num("age", ndec_mean = 2, ndec_med = 0)
tabo_agecat <- svy_cat("agecat")
tabo_gender <- svy_cat("gender")
tabo_race <- svy_cat("race")
tabo_pir <- svy_num("pir", ndec_mean = 2, ndec_med = 2)
tabo_insurance <- svy_cat_yn("insurance")
tabo_cotcat <- svy_cat("cotcat")
tabo_cot <- svy_num("cot", ndec_mean = 2, ndec_med = 2)
tabo_pah <- svy_num("ratio_pah", ndec_mean = 2, ndec_med = 2)
tabo_asthma <- svy_cat_yn("asthma")
  # Meds
tabo_med_short <- svy_cat_yn("med_short")
tabo_med_steroid_no <- svy_cat_yn("med_short_steroid_no")
tabo_med_steroid <- svy_cat_yn("med_short_steroid")
tabo_med_alb <- svy_cat_yn("alb")
tabo_med_lalb <- svy_cat_yn("lalb")
tabo_med_pbut <- svy_cat_yn("pbut")
tabo_med_mprediso <- svy_cat_yn("mprediso")
tabo_med_prediso <- svy_cat_yn("prediso")
tabo_med_pred <- svy_cat_yn("pred")
tabo_med_dex <- svy_cat_yn("dex")
tabo_med_med_long <- svy_cat_yn("med_long")
tabo_med_becl <- svy_cat_yn("becl")
tabo_med_bud <- svy_cat_yn("bud")
tabo_med_flun <- svy_cat_yn("flun")
tabo_med_flut <- svy_cat_yn("flut")
tabo_med_tria <- svy_cat_yn("tria")
tabo_med_crom <- svy_cat_yn("crom") # none found
tabo_med_mont <- svy_cat_yn("mont")
tabo_med_zaf <- svy_cat_yn("zaf")
tabo_med_form <- svy_cat_yn("form")
tabo_med_sal <- svy_cat_yn("sal")
tabo_med_theo <- svy_cat_yn("theo")
tabo_med_mom <- svy_cat_yn("mom")
tabo_med_nedo <- svy_cat_yn("nedo")
tabo_med_oma <- svy_cat_yn("oma")

  # Create table df, Overall
tab1o <- c(tabo_yr, '','',
          tabo_agecat, 
          tabo_age, '','',
          tabo_gender, '','',
          tabo_race, '','',
          tabo_pir, '',
          tabo_insurance, '','',
          tabo_cotcat, '',
          tabo_asthma,
          tabo_med_short,
          tabo_med_steroid_no,
          tabo_med_steroid,
          tabo_med_med_long) %>% as.data.frame() %>% t() 
tabs3o <- c(tabo_med_short, 
            tabo_med_alb, 
            tabo_med_lalb, 
            tabo_med_pbut, 
            tabo_med_mprediso, 
            tabo_med_prediso, 
            tabo_med_pred, 
            tabo_med_dex, 
            '',
            tabo_med_med_long, 
            tabo_med_becl, 
            tabo_med_bud, 
            tabo_med_flun, 
            tabo_med_flut, 
            tabo_med_tria, 
            #tabo_med_crom, 
            tabo_med_mont, 
            tabo_med_zaf, 
            tabo_med_form, 
            tabo_med_sal, 
            tabo_med_theo,
            tabo_med_mom, 
            tabo_med_nedo, 
            tabo_med_oma
            ) %>% as.data.frame()

### Table, by PAH Quartile
# Prep
tabq_yr <- svy_cat_q("yr")
tabq_age <- svy_num_q("age", ndec_mean = 2, ndec_med = 0) %>% as.matrix()
tabq_agecat <- svy_cat_q("agecat")
tabq_gender <- svy_cat_q("gender")
tabq_race <- svy_cat_q("race")
tabq_pir <- svy_num_q("pir", ndec_mean = 2, ndec_med = 2) %>% as.matrix()
tabq_insurance <- svy_cat_yn_q("insurance")
tabq_cotcat <- svy_cat_q("cotcat")
tabq_cot <- svy_num_q("cot", ndec_mean = 2, ndec_med = 2) %>% as.matrix()
tabq_pah <- svy_num_q("ratio_pah", ndec_mean = 2, ndec_med = 2) %>% as.matrix()
tabq_asthma <- svy_cat_yn_q("asthma")
  # Meds
tabq_med_short <- svy_cat_yn_q("med_short")
tabq_med_steroid_no <- svy_cat_yn_q("med_short_steroid_no")
tabq_med_steroid <- svy_cat_yn_q("med_short_steroid")
tabq_med_alb <- svy_cat_yn_q("alb")
tabq_med_lalb <- svy_cat_yn_q("lalb")
tabq_med_pbut <- svy_cat_yn_q("pbut")
tabq_med_mprediso <- svy_cat_yn_q("mprediso")
tabq_med_prediso <- svy_cat_yn_q("prediso")
tabq_med_pred <- svy_cat_yn_q("pred")
tabq_med_dex <- svy_cat_yn_q("dex")
tabq_med_med_long <- svy_cat_yn_q("med_long")
tabq_med_becl <- svy_cat_yn_q("becl")
tabq_med_bud <- svy_cat_yn_q("bud")
tabq_med_flun <- svy_cat_yn_q("flun")
tabq_med_flut <- svy_cat_yn_q("flut")
tabq_med_tria <- svy_cat_yn_q("tria")
tabq_med_crom <- svy_cat_yn_q("crom")
tabq_med_mont <- svy_cat_yn_q("mont")
tabq_med_zaf <- svy_cat_yn_q("zaf")
tabq_med_form <- svy_cat_yn_q("form")
tabq_med_sal <- svy_cat_yn_q("sal")
tabq_med_theo <- svy_cat_yn_q("theo")
tabq_med_mom <- svy_cat_yn_q("mom")
tabq_med_nedo <- svy_cat_yn_q("nedo")
tabq_med_oma <- svy_cat_yn_q("oma")

  # Create table
tab1q <- rbind(tabq_yr, '', '',
              tabq_agecat, 
              tabq_age, '','',
              tabq_gender, '','',
              tabq_race, '','',
              tabq_pir, '',
              tabq_insurance, '','',
              tabq_cotcat, '',
              tabq_asthma,
              tabq_med_short,
              tabq_med_steroid_no,
              tabq_med_steroid,
              tabq_med_med_long)
tabs3q <- rbind(tabq_med_short, 
                tabq_med_alb, 
                tabq_med_lalb, 
                tabq_med_pbut, 
                tabq_med_mprediso, 
                tabq_med_prediso, 
                tabq_med_pred, 
                tabq_med_dex,
                '',
                tabq_med_med_long, 
                tabq_med_becl, 
                tabq_med_bud, 
                tabq_med_flun, 
                tabq_med_flut, 
                tabq_med_tria, 
                #tabq_med_crom, 
                tabq_med_mont, 
                tabq_med_zaf, 
                tabq_med_form, 
                tabq_med_sal, 
                tabq_med_theo,
                tabq_med_mom, 
                tabq_med_nedo, 
                tabq_med_oma) 


# View table results
cbind(tab1o, "", tab1q) %>% View()
cbind(tabs3o, "", tabs3q) %>% View()

# Just double check no capture
c(tabo_med_crom, "", tabq_med_crom) 
c(tabo_med_nedo, "", tabq_med_nedo) 
c(tabo_med_oma, "", tabq_med_oma) 
pahas %>% tabyl(crom, nedo, oma)


################################################## 
########## TABLE 2: REGRESSION ANALYSIS ########## 
##################################################
### Function to extract results
extract_mod_results <- function(mod){
  mod_suffix <- paste0("_", stringr::str_sub(deparse(substitute(mod)), -1))
  mod_result <- cbind(coef(mod), confint(mod)) %>% 
    exp() %>% 
    as_tibble(name_repair = "unique", rownames = "var") %>% 
    dplyr::rename_all(~ c("var", paste0("OR", mod_suffix), paste0("CI_low", mod_suffix), paste0("CI_high", mod_suffix)))
  return(mod_result)
}
########## Table 2 ##########
### Setup - general, variables
outcome_var <- "med_short"
pah_exp_1_cont <- "pahz"
pah_exp_2_cont <- "ratio_pahz"
pah_exp_3_cont <- "ratio_pahadjz"
pah_exp_4_cont <- "pahz + creat"
pah_exp_5_cont <- "pahz + creat_resid"
pah_exp_6_cont <- "ratio_pahz + creat"
pah_exp_7_cont <- "ratio_pahadjz + creat"
vars_a <- ""
vars_b <- "+ age + gender + race + pir + insurance + cot"

### OVERALL
  # Weights / Svy Design
pahas_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_lab, nest = TRUE, data = pahas)
  # Setup - overall
family_expr <- expr(quasipoisson(link = "log"))
  # Models
mod1a <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_1_cont,vars_a)), family = family_expr, design = pahas_design)
mod1b <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_1_cont,vars_b)), family = family_expr, design = pahas_design)
mod2a <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_2_cont,vars_a)), family = family_expr, design = pahas_design)
mod2b <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_2_cont,vars_b)), family = family_expr, design = pahas_design)
mod3a <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_3_cont,vars_a)), family = family_expr, design = pahas_design)
mod3b <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_3_cont,vars_b)), family = family_expr, design = pahas_design)
mod4a <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_4_cont,vars_a)), family = family_expr, design = pahas_design)
mod4b <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_4_cont,vars_b)), family = family_expr, design = pahas_design)
mod5a <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_5_cont,vars_a)), family = family_expr, design = pahas_design)
mod5b <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_5_cont,vars_b)), family = family_expr, design = pahas_design)
mod6a <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_6_cont,vars_a)), family = family_expr, design = pahas_design)
mod6b <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_6_cont,vars_b)), family = family_expr, design = pahas_design)
mod7a <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_cont,vars_a)), family = family_expr, design = pahas_design)
mod7b <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_cont,vars_b)), family = family_expr, design = pahas_design)

### ASTHMA
  # Subset
asthma_pahas_design <- subset(pahas_design, asthma == 1)
  # Setup
asthma_family_expr <- expr(quasipoisson(link = "log"))
  # Models
mod1a_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_1_cont,vars_a)), family = asthma_family_expr, design = asthma_pahas_design)
mod1b_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_1_cont,vars_b)), family = asthma_family_expr, design = asthma_pahas_design)
mod2a_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_2_cont,vars_a)), family = asthma_family_expr, design = asthma_pahas_design)
mod2b_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_2_cont,vars_b)), family = asthma_family_expr, design = asthma_pahas_design)
mod3a_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_3_cont,vars_a)), family = asthma_family_expr, design = asthma_pahas_design)
mod3b_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_3_cont,vars_b)), family = asthma_family_expr, design = asthma_pahas_design)
mod4a_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_4_cont,vars_a)), family = asthma_family_expr, design = asthma_pahas_design)
mod4b_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_4_cont,vars_b)), family = asthma_family_expr, design = asthma_pahas_design)
mod5a_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_5_cont,vars_a)), family = asthma_family_expr, design = asthma_pahas_design)
mod5b_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_5_cont,vars_b)), family = asthma_family_expr, design = asthma_pahas_design)
mod6a_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_6_cont,vars_a)), family = asthma_family_expr, design = asthma_pahas_design)
mod6b_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_6_cont,vars_b)), family = asthma_family_expr, design = asthma_pahas_design)
mod7a_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_cont,vars_a)), family = asthma_family_expr, design = asthma_pahas_design)
mod7b_asthma <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_cont,vars_b)), family = asthma_family_expr, design = asthma_pahas_design)

### Table 2 Results
# All variables, Model 7
extract_mod_results(mod7a) %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER = '') %>% 
  dplyr::right_join(extract_mod_results(mod7b), by = "var") %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER3 = '') %>% 
  dplyr::left_join(extract_mod_results(mod7a_asthma), by = "var") %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER5 = '') %>% 
  dplyr::right_join(extract_mod_results(mod7b_asthma), by = "var") %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER6 = '') %>% 
  View()

### Sensitivity analysis-- race adjustment in creatinine
  # results same
mod7b_race_creat <- svyglm(as.formula(paste0("med_short ~ ratio_pahadjz_race + creat + age + gender + race + pir + insurance + cot")), family = family_expr, design = pahas_design)
broom::tidy(mod7b_race_creat, exp = TRUE, conf.int = TRUE)
  # results same
mod7b_asthma_race_creat <- svyglm(as.formula(paste0("med_short ~ ratio_pahadjz_race + creat + age + gender + race + pir + insurance + cot")), family = family_expr, design = asthma_pahas_design)
broom::tidy(mod7b_asthma_race_creat, exp = TRUE, conf.int = TRUE)

### Sensitivity analysis-- race adjustment in creatinine
  # results same
mod7b_no_survey_creat <- svyglm(as.formula(paste0("med_short ~ ratio_pahadjz_no_survey + creat + age + gender + race + pir + insurance + cot")), family = family_expr, design = pahas_design)
broom::tidy(mod7b_no_survey_creat, exp = TRUE, conf.int = TRUE)
  # results same
mod7b_asthma_no_survey_creat <- svyglm(as.formula(paste0("med_short ~ ratio_pahadjz_no_survey + creat + age + gender + race + pir + insurance + cot")), family = family_expr, design = asthma_pahas_design)
broom::tidy(mod7b_asthma_no_survey_creat, exp = TRUE, conf.int = TRUE)

### Sensitivity analysis-- year adjustment (period effects)
  # results same
mod7b_creat_yr <- svyglm(as.formula(paste0("med_short ~ ratio_pahadjz + creat + age + gender + race + pir + insurance + cot + yr")), family = family_expr, design = pahas_design)
broom::tidy(mod7b_creat_yr, exp = TRUE, conf.int = TRUE)
  # results same
mod7b_asthma_creat_yr <- svyglm(as.formula(paste0("med_short ~ ratio_pahadjz + creat + age + gender + race + pir + insurance + cot + yr")), family = family_expr, design = asthma_pahas_design)
broom::tidy(mod7b_asthma_creat_yr, exp = TRUE, conf.int = TRUE)


### Sensitivity analysis-- BMI adjustment
  # results same
mod7b_creat_bmi <- svyglm(as.formula(paste0("med_short ~ ratio_pahadjz + creat + age + gender + race + pir + insurance + cot + bmi")), family = family_expr, design = pahas_design)
broom::tidy(mod7b_creat_bmi, exp = TRUE, conf.int = TRUE)
  # results same
mod7b_asthma_creat_bmi <- svyglm(as.formula(paste0("med_short ~ ratio_pahadjz + creat + age + gender + race + pir + insurance + cot + bmi")), family = family_expr, design = asthma_pahas_design)
broom::tidy(mod7b_asthma_creat_bmi, exp = TRUE, conf.int = TRUE)


################################################## 
########## TABLE 3 ########## 
##################################################
### Sensitivity analysis-- different rescue meds (steroid and non-steroid)
  # SABA (non-steroid)
mod7a_steroid_no <- svyglm(as.formula(paste0("med_short_steroid_no~",pah_exp_7_cont,vars_a)), family = family_expr, design = pahas_design)
mod7b_steroid_no <- svyglm(as.formula(paste0("med_short_steroid_no~",pah_exp_7_cont,vars_b)), family = family_expr, design = pahas_design)
mod7a_asthma_steroid_no <- svyglm(as.formula(paste0("med_short_steroid_no~",pah_exp_7_cont,vars_a)), family = family_expr, design = asthma_pahas_design)
mod7b_asthma_steroid_no <- svyglm(as.formula(paste0("med_short_steroid_no~",pah_exp_7_cont,vars_b)), family = family_expr, design = asthma_pahas_design)
  # steroid
mod7a_steroid <- svyglm(as.formula(paste0("med_short_steroid~",pah_exp_7_cont,vars_a)), family = family_expr, design = pahas_design)
mod7b_steroid <- svyglm(as.formula(paste0("med_short_steroid~",pah_exp_7_cont,vars_b)), family = family_expr, design = pahas_design)
mod7a_asthma_steroid <- svyglm(as.formula(paste0("med_short_steroid~",pah_exp_7_cont,vars_a)), family = family_expr, design = asthma_pahas_design)
mod7b_asthma_steroid <- svyglm(as.formula(paste0("med_short_steroid~",pah_exp_7_cont,vars_b)), family = family_expr, design = asthma_pahas_design)
  # View results - overall
extract_mod_results(mod7a_steroid_no) %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER = '') %>% 
  dplyr::right_join(extract_mod_results(mod7b_steroid_no), by = "var") %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER3 = '') %>% 
  dplyr::left_join(extract_mod_results(mod7a_steroid), by = "var") %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER5 = '') %>% 
  dplyr::right_join(extract_mod_results(mod7b_steroid), by = "var") %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER6 = '') %>% 
  View()
  # View results - asthma
extract_mod_results(mod7a_asthma_steroid_no) %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER = '') %>% 
  dplyr::right_join(extract_mod_results(mod7b_asthma_steroid_no), by = "var") %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER3 = '') %>% 
  dplyr::left_join(extract_mod_results(mod7a_asthma_steroid), by = "var") %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER5 = '') %>% 
  dplyr::right_join(extract_mod_results(mod7b_asthma_steroid), by = "var") %>% dplyr::filter(var == 'ratio_pahadjz') %>% add_column(SPACER6 = '') %>% 
  View()


broom::tidy(mod7b_steroid, exp = TRUE, conf.int = TRUE)
broom::tidy(mod7b_asthma_steroid, exp = TRUE, conf.int = TRUE)
broom::tidy(mod7b_steroid_no, exp = TRUE, conf.int = TRUE)
broom::tidy(mod7b_asthma_steroid, exp = TRUE, conf.int = TRUE)

################################################## 
########## SUPPLEMENTARY TABLE 2 ########## 
##################################################
cbind(
      rbind(
            cbind(extract_mod_results(mod1a) %>% dplyr::slice(2L), 
                  extract_mod_results(mod1b) %>% dplyr::slice(2L)),
            cbind(extract_mod_results(mod2a) %>% dplyr::slice(2L), 
                  extract_mod_results(mod2b) %>% dplyr::slice(2L)),
            cbind(extract_mod_results(mod3a) %>% dplyr::slice(2L), 
                  extract_mod_results(mod3b) %>% dplyr::slice(2L)),
            cbind(extract_mod_results(mod4a) %>% dplyr::slice(2L), 
                  extract_mod_results(mod4b) %>% dplyr::slice(2L)),
            cbind(extract_mod_results(mod5a) %>% dplyr::slice(2L), 
                  extract_mod_results(mod5b) %>% dplyr::slice(2L)),
            cbind(extract_mod_results(mod6a) %>% dplyr::slice(2L), 
                  extract_mod_results(mod6b) %>% dplyr::slice(2L)),
            cbind(extract_mod_results(mod7a) %>% dplyr::slice(2L), 
                  extract_mod_results(mod7b) %>% dplyr::slice(2L))),
      rbind(
            cbind(extract_mod_results(mod1a_asthma) %>% dplyr::slice(2L), 
                  extract_mod_results(mod1b_asthma) %>% dplyr::slice(2L)),
            cbind(extract_mod_results(mod2a_asthma) %>% dplyr::slice(2L), 
                  extract_mod_results(mod2b_asthma) %>% dplyr::slice(2L)), 
            cbind(extract_mod_results(mod3a_asthma) %>% dplyr::slice(2L), 
                  extract_mod_results(mod3b_asthma) %>% dplyr::slice(2L)),
            cbind(extract_mod_results(mod4a_asthma) %>% dplyr::slice(2L), 
                  extract_mod_results(mod4b_asthma) %>% dplyr::slice(2L)), 
            cbind(extract_mod_results(mod5a_asthma) %>% dplyr::slice(2L), 
                  extract_mod_results(mod5b_asthma) %>% dplyr::slice(2L)),
            cbind(extract_mod_results(mod6a_asthma) %>% dplyr::slice(2L), 
                  extract_mod_results(mod6b_asthma) %>% dplyr::slice(2L)),
            cbind(extract_mod_results(mod7a_asthma) %>% dplyr::slice(2L), 
                  extract_mod_results(mod7b_asthma) %>% dplyr::slice(2L)))
      ) %>% View()

################################################## 
########## SUPPLEMENTARY TABLE 1 ########## 
##################################################
# str_detect(term, "^ratio_pahadjz:")
pah_exp_7_q <- "ratio_pahadjq + creat"
mod7aq <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_q,vars_a)), family = family_expr, design = pahas_design)
mod7bq <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_q,vars_b)), family = family_expr, design = pahas_design)
mod7a_asthmaq <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_q,vars_a)), family = asthma_family_expr, design = asthma_pahas_design)
mod7b_asthmaq <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_q,vars_b)), family = asthma_family_expr, design = asthma_pahas_design)
# View Results
extract_mod_results(mod7aq) %>% dplyr::filter(str_detect(var, '^ratio_pahadjqQ')) %>% add_column(SPACER = '') %>% 
  dplyr::right_join(extract_mod_results(mod7bq), by = "var") %>% dplyr::filter(str_detect(var, '^ratio_pahadjqQ')) %>% add_column(SPACER3 = '') %>% 
  dplyr::left_join(extract_mod_results(mod7a_asthmaq), by = "var") %>% dplyr::filter(str_detect(var, '^ratio_pahadjqQ')) %>% add_column(SPACER5 = '') %>% 
  dplyr::right_join(extract_mod_results(mod7b_asthmaq), by = "var") %>% dplyr::filter(str_detect(var, '^ratio_pahadjqQ')) %>% add_column(SPACER6 = '') %>% 
  View()

################################################## 
########## FIGURE 1: INTERACTION PLOT (STRATIFIED BY AGE/CONTROLLER MED USE) ########## 
##################################################
## Model w/interaction
  # Controller meds
medlong_pah_exp_7 <- "ratio_pahadjz * med_long + creat"
vars_b_medlong <- "+ med_long + age + gender + race + pir + insurance + cot"
medlong_family_expr <- expr(quasipoisson(link = "log"))
mod7b_medlong <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_cont,vars_b_medlong)), family = medlong_family_expr,design = pahas_design)
mod7b_medlong_int <- svyglm(as.formula(paste0(outcome_var,"~",medlong_pah_exp_7,vars_b_medlong)), family = medlong_family_expr, design = pahas_design)
  # Test of interaction
anova(mod7b_medlong, mod7b_medlong_int) # p=0.72517
  # Age - 4 categories
vars_b_agecat4 <- "+ agecat4 + gender + race + pir + insurance + cot"
age_int_family_expr <- expr(quasipoisson(link = "log"))
mod7b_agecat4 <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_cont,vars_b_agecat4)), family = age_int_family_expr, data = data_expr, design = pahas_design)
mod7b_agecat4_int <- svyglm(as.formula(paste0(outcome_var,"~ ratio_pahadjz * agecat4 + creat",vars_b_agecat4)), family = age_int_family_expr, data = data_expr, design = pahas_design)
  # Test of interaction
anova(mod7b_agecat4, mod7b_agecat4_int) # p=0.22141
## Contrasts to extract stratified PRs + CIs
    # Age
mod7b_agecat4_pah_coef <- mod7b_agecat4_int %>% broom::tidy() %>% dplyr::filter(term == 'ratio_pahadjz') %>% pluck(2) # Grab coefficient for PAH, Alt: mod7b_agecat_int %>% coef() %>% pluck(2)
mod7b_agecat4_pah_se <- mod7b_agecat4_int %>% broom::tidy() %>% dplyr::filter(term == 'ratio_pahadjz') %>% pluck(3)
mod7b_agecat4_agecat_coef <- mod7b_agecat4_int %>% broom::tidy() %>% dplyr::filter(str_detect(term, "^ratio_pahadjz:")) %>% pluck(2)
contr_ageref4 <- tibble(nlcon = mod7b_agecat4_pah_coef, SE = mod7b_agecat4_pah_se)
contr_age19_4 <- svycontrast(mod7b_agecat4_int, quote(`ratio_pahadjz` + `ratio_pahadjz:agecat419-45`)) 
contr_age46_4 <- svycontrast(mod7b_agecat4_int, quote(`ratio_pahadjz` + `ratio_pahadjz:agecat446-70`)) 
contr_age70_4 <- svycontrast(mod7b_agecat4_int, quote(`ratio_pahadjz` + `ratio_pahadjz:agecat4>70`)) 
contr_age4 <- rbind(contr_ageref4, contr_age19_4 %>% as_tibble(), contr_age46_4 %>% as_tibble(), contr_age70_4 %>% as_tibble())
    # Controller meds
mod7b_medlong_pah_coef <- mod7b_medlong_int %>% broom::tidy() %>% dplyr::filter(term == 'ratio_pahadjz') %>% pluck(2) # Grab coefficient for PAH, Alt: mod7b_agecat_int %>% coef() %>% pluck(2)
mod7b_medlong_pah_se <- mod7b_medlong_int %>% broom::tidy() %>% dplyr::filter(term == 'ratio_pahadjz') %>% pluck(3)
contr_medlongref <- tibble(nlcon = mod7b_medlong_pah_coef, SE = mod7b_medlong_pah_se)
contr_medlongy <- svycontrast(mod7b_medlong_int, quote(`ratio_pahadjz` + `ratio_pahadjz:med_long`)) 
contr_medlong <- rbind(contr_medlongref, contr_medlongy %>% as_tibble())
## Prep data
  # Age
agecat4 <- c("18 or younger", "19-45", "46-70", "Older than 70")
agecat_int4 <- 
  tibble(var = "Age",
         cat = agecat4 %>% factor(ordered = TRUE, levels = agecat4),
         #est_pah = mod7b_agecat_pah_coef,
         est = contr_age4$nlcon,
         se = contr_age4$SE) 
  # Controller meds
medlongcat <- c("No Controller\nMedication Use", "Controller\nMedication Use")
medlong_int <- 
  tibble(var = "Controller Medication Use",
         cat = medlongcat %>% factor(ordered = TRUE, levels = medlongcat),
         #est_pah = mod7b_agecat_pah_coef,
         est = contr_medlong$nlcon,
         se = contr_medlong$SE)
  # combine
int_pr4 <- rbind(agecat_int4, medlong_int) %>% 
  dplyr::mutate(#est = est_pah + est_var,
    pr = est %>% exp(),
    pr_ln = log(pr),
    ci_low = (est - 1.96*se) %>% exp(),
    ci_high = (est + 1.96*se) %>% exp(),
    ci_low_ln = log(ci_low),
    ci_high_ln = log(ci_high))

##### Create plot
int_pr4 %>% 
  ggplot(aes(x = cat, y = pr)) +
  facet_grid(. ~ var, scales = "free_x") +
  geom_point() +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = 2) +
  theme_bw() + 
  scale_y_continuous(trans = 'log10', breaks = c(0.6, 0.8, 1.0, 1.2, 1.4)) +
  #scale_y_continuous(trans = 'log10', breaks = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)) +
  #scale_y_continuous(minor_breaks = seq(-0.6, 0.5, 0.1), breaks = seq(-0.4, 0.5, 0.2)) + # If want to show transformed PR instead
  #scale_y_continuous(expand = c(0,0), limits = c(0,1.5), minor_breaks = seq(0, 1.5, 0.1)) +
  labs(y = "Prevalence Ratio (log scale)") +
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 15), 
        axis.title.x = element_blank(), #element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 15),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        strip.text.x = element_text(size = 13))

## Export data
ggsave(paste0('exports/fig1_int_', Sys.Date() %>% format("%Y%m%d"), ".jpeg"), width = 11, height = 7, dpi = 600) 

################################################## 
########## Supplemental Figure 1. Age Interaction Plot Sensitivity Analysis ########## 
##################################################
##### 5-categorical age, Contrasts
  # Model
age_int_pah_exp_cont_7 <- "ratio_pahadjz * agecat + creat"
vars_b_agecat <- "+ agecat + gender + race + pir + insurance + cot"
mod7b_agecat <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_cont,vars_b_agecat)), family = age_int_family_expr, data = data_expr, design = pahas_design)
mod7b_agecat_int <- svyglm(as.formula(paste0(outcome_var,"~",age_int_pah_exp_cont_7,vars_b_agecat)), family = age_int_family_expr, data = data_expr, design = pahas_design)
  # Test
anova(mod7b_agecat, mod7b_agecat_int) # p=0.45163
  # Contrasts to extract stratified PRs + CIs
mod7b_agecat_int %>% broom::tidy(exponentiate = TRUE, conf.int = TRUE)
mod7b_agecat_int %>% broom::tidy(conf.int = TRUE)
mod7b_agecat_int %>% broom::tidy() %>% dplyr::mutate(low_ci = estimate - 1.96*std.error, high_ci = estimate + 1.96*std.error)
mod7b_agecat_pah_coef <- mod7b_agecat_int %>% broom::tidy() %>% dplyr::filter(term == 'ratio_pahadjz') %>% pluck(2) # Grab coefficient for PAH, Alt: mod7b_agecat_int %>% coef() %>% pluck(2)
mod7b_agecat_pah_se <- mod7b_agecat_int %>% broom::tidy() %>% dplyr::filter(term == 'ratio_pahadjz') %>% pluck(3)
mod7b_agecat_agecat_coef <- mod7b_agecat_int %>% broom::tidy() %>% dplyr::filter(str_detect(term, "^ratio_pahadjz:")) %>% pluck(2)
# mod7b_agecat_agecat_stderr <- mod7b_agecat_int %>% broom::tidy() %>% dplyr::filter(str_detect(term, "^ratio_pahadjz:")) %>% pluck(3)
contr_ageref <- tibble(nlcon = mod7b_agecat_pah_coef, SE = mod7b_agecat_pah_se)
contr_age12 <- svycontrast(mod7b_agecat_int, quote(`ratio_pahadjz` + `ratio_pahadjz:agecat12-18`)) 
contr_age19 <- svycontrast(mod7b_agecat_int, quote(`ratio_pahadjz` + `ratio_pahadjz:agecat19-45`)) 
contr_age46 <- svycontrast(mod7b_agecat_int, quote(`ratio_pahadjz` + `ratio_pahadjz:agecat46-70`)) 
contr_age70 <- svycontrast(mod7b_agecat_int, quote(`ratio_pahadjz` + `ratio_pahadjz:agecat>70`)) 
contr_age <- rbind(contr_ageref, contr_age12 %>% as_tibble(), contr_age19 %>% as_tibble(), contr_age46 %>% as_tibble(), contr_age70 %>% as_tibble())
  # Prepare data
agecat <- c("Younger\nthan 12", "12-18", "19-45", "46-70", "Older\nthan 70")
agecat_int <- 
  tibble(var = "Age (5 Categories)",
         cat = agecat %>% factor(ordered = TRUE, levels = agecat),
         #est_pah = mod7b_agecat_pah_coef,
         est = contr_age$nlcon,
         se = contr_age$SE) 

##### Dichotomous age
  # Model
vars_b_agecat2 <- "+ agecat2 + gender + race + pir + insurance + cot"
mod7b_agecat2 <- svyglm(as.formula(paste0(outcome_var,"~",pah_exp_7_cont,vars_b_agecat2)), family = age_int_family_expr, data = data_expr, design = pahas_design)
mod7b_agecat2_int <- svyglm(as.formula(paste0(outcome_var,"~ ratio_pahadjz * agecat2 + creat",vars_b_agecat2)), family = age_int_family_expr, data = data_expr, design = pahas_design)
  # Test
anova(mod7b_agecat2, mod7b_agecat2_int) # p=0.15599
  # Contrasts to extract stratified PRs + CIs
mod7b_agecat2_pah_coef <- mod7b_agecat2_int %>% broom::tidy() %>% dplyr::filter(term == 'ratio_pahadjz') %>% pluck(2) # Grab coefficient for PAH, Alt: mod7b_agecat_int %>% coef() %>% pluck(2)
mod7b_agecat2_pah_se <- mod7b_agecat2_int %>% broom::tidy() %>% dplyr::filter(term == 'ratio_pahadjz') %>% pluck(3)
mod7b_agecat2_agecat_coef <- mod7b_agecat2_int %>% broom::tidy() %>% dplyr::filter(str_detect(term, "^ratio_pahadjz:")) %>% pluck(2)
contr_ageref2 <- tibble(nlcon = mod7b_agecat2_pah_coef, SE = mod7b_agecat2_pah_se)
contr_age19_2 <- svycontrast(mod7b_agecat2_int, quote(`ratio_pahadjz` + `ratio_pahadjz:agecat2>19`)) 
contr_age2 <- rbind(contr_ageref2, contr_age19_2 %>% as_tibble())
  # Prepare data
agecat2 <- c("18 or younger", "Older than 18")
agecat_int2 <- 
  tibble(var = "Age (Dichotomous)",
         cat = agecat2 %>% factor(ordered = TRUE, levels = agecat2),
         #est_pah = mod7b_agecat_pah_coef,
         est = contr_age2$nlcon,
         se = contr_age2$SE) 


##### Concatenate data for two age categorizations
int_pr_s1 <- rbind(agecat_int2, agecat_int) %>% 
  dplyr::mutate(#est = est_pah + est_var,
    pr = est %>% exp(),
    pr_ln = log(pr),
    ci_low = (est - 1.96*se) %>% exp(),
    ci_high = (est + 1.96*se) %>% exp(),
    ci_low_ln = log(ci_low),
    ci_high_ln = log(ci_high))

##### Create plot
int_pr_s1 %>% 
  ggplot(aes(x = cat, y = pr)) +
  facet_grid(. ~ var, scales = "free_x") +
  geom_point() +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = 2) +
  theme_bw() + 
  scale_y_continuous(trans = 'log10', breaks = c(0.6, 0.8, 1.0, 1.2, 1.4)) +
  labs(y = "Prevalence Ratio (log scale)") +
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 15),
        axis.title.x = element_blank(), #element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        strip.text.x = element_text(size = 13))

##### Export
ggsave(paste0('exports/figs1_int_', Sys.Date() %>% format("%Y%m%d"), ".png"), width = 11, height = 7)
