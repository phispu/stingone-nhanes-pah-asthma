########################################################################################################################################################
#
# Urinary metabolites of polycyclic aromatic hydrocarbons and short-acting beta agonist or systemic corticosteroid asthma medication use within NHANES
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
# Data Import
setwd(here("data/other"))
rxq_drug <- import('RXQ_DRUG.XPT')
setwd(here("data/raw"))
  # 2005-06
demo06 <- import('2005_06/DEMO_D.XPT')
bmi06 <- import('2005_06/BMX_D.XPT')
albcr06 <- import('2005_06/ALB_CR_D.XPT')
pah06 <- import('2005_06/PAH_D.XPT')
cot06 <- import('2005_06/COT_D.XPT')
mcq06 <- import('2005_06/MCQ_D.XPT')
rx06 <- import('2005_06/RXQ_RX_D.XPT')
hiq06 <- import('2005_06/HIQ_D.XPT')
  # 2007-08
demo08 <- import('2007_08/DEMO_E.XPT')
bmi08 <- import('2007_08/BMX_E.XPT')
albcr08 <- import('2007_08/ALB_CR_E.XPT')
pah08 <- import('2007_08/PAH_E.XPT')
cot08 <- import('2007_08/COTNAL_E.XPT')
mcq08 <- import('2007_08/MCQ_E.XPT')
rx08 <- import('2007_08/RXQ_RX_E.XPT')
hiq08 <- import('2007_08/HIQ_E.XPT')
  # 2009-10
demo10 <- import('2009_10/DEMO_F.XPT')
bmi10 <- import('2009_10/BMX_F.XPT')
albcr10 <- import('2009_10/ALB_CR_F.XPT')
pah10 <- import('2009_10/PAH_F.XPT')
cot10 <- import('2009_10/COTNAL_F.XPT')
mcq10 <- import('2009_10/MCQ_F.XPT')
rx10 <- import('2009_10/RXQ_RX_F.XPT')
hiq10 <- import('2009_10/HIQ_F.XPT')
  # 2011-12
demo12 <- import('2011_12/DEMO_G.XPT')
bmi12 <- import('2011_12/BMX_G.XPT')
albcr12 <- import('2011_12/ALB_CR_G.XPT')
pah12 <- import('2011_12/PAH_G.XPT')
cot12 <- import('2011_12/COTNAL_G.XPT')
mcq12 <- import('2011_12/MCQ_G.XPT')
rx12 <- import('2011_12/RXQ_RX_G.XPT')
hiq12 <- import('2011_12/HIQ_G.XPT')
  # 2013-14
demo14 <- import('2013_14/DEMO_H.XPT')
bmi14 <- import('2013_14/BMX_H.XPT')
albcr14 <- import('2013_14/ALB_CR_H.XPT')
pah14 <- import('2013_14/PAH_H.XPT')
cot14 <- import('2013_14/COT_H.XPT')
mcq14 <- import('2013_14/MCQ_H.XPT')
rx14 <- import('2013_14/RXQ_RX_H.XPT')
hiq14 <- import('2013_14/HIQ_H.XPT')
  # 2015-16
demo16 <- import('2015_16/DEMO_I.XPT')
bmi16 <- import('2015_16/BMX_I.XPT')
albcr16 <- import('2015_16/ALB_CR_I.XPT')
pah16 <- import('2015_16/PAH_I.XPT')
cot16 <- import('2015_16/COT_I.XPT')
mcq16 <- import('2015_16/MCQ_I.XPT')
rx16 <- import('2015_16/RXQ_RX_I.XPT')
hiq16 <- import('2015_16/HIQ_I.XPT')

### Flag for lab data
pah06 <- pah06 %>% dplyr::mutate(flag_lab = 1)
pah08 <- pah08 %>% dplyr::mutate(flag_lab = 1)
pah10 <- pah10 %>% dplyr::mutate(flag_lab = 1)
pah12 <- pah12 %>% dplyr::mutate(flag_lab = 1)
pah14 <- pah14 %>% dplyr::mutate(flag_lab = 1)
pah16 <- pah16 %>% dplyr::mutate(flag_lab = 1)

################################
########## MERGE DATA ##########
################################
  # Variable List
varlist <- c('SEQN','RIAGENDR','RIDAGEYR','RIDAGEMN','RIDRETH1',
             'DMDCITZN','DMDYRSUS','DMDEDUC3','DMDEDUC2','DMDMARTL','FIALANG','DMDHHSIZ','DMDFMSIZ',
             'INDFMPIR',
             'WTINT2YR','WTMEC2YR',
             'SDMVPSU','SDMVSTRA',
             'URXCRS',
             'BMXBMI','BMXWT', 'BMXHT',
             'HIQ011', # Health insurance
             # PAH subsample weight: 05-10: WTSB2YR, 11-16: WTSA2YR
             'URXP02','URXP10',
             'MCQ010','MCQ025','MCQ035','MCQ040',
             'LBXCOT', 
             'flag_lab', 'yr') 
  # 05-06
    # Warning - keep a look out for errors: Column `SEQN` has different attributes on LHS and RHS of join --> fixed!
pahas_orig06 <- demo06 %>% 
  left_join(bmi06, by = "SEQN") %>% 
  left_join(albcr06, by = "SEQN") %>% 
  left_join(pah06, by = "SEQN") %>% 
  left_join(mcq06, by = "SEQN") %>% 
  left_join(cot06, by = "SEQN") %>%
  left_join(hiq06, by = "SEQN") %>% 
  dplyr::mutate(yr = '2005-2006') %>% 
  dplyr::select(one_of(c(varlist, 'WTSB2YR'))) %>% 
  dplyr::rename(wt_lab_orig = WTSB2YR) 
  # 07-08
pahas_orig08 <- demo08 %>% 
  left_join(bmi08, by = "SEQN") %>% 
  left_join(albcr08, by = "SEQN") %>% 
  left_join(pah08, by = "SEQN") %>% 
  left_join(mcq08, by = "SEQN") %>% 
  left_join(cot08, by = "SEQN") %>% 
  left_join(hiq08, by = "SEQN") %>% 
  dplyr::mutate(yr = '2007-2008') %>% 
  dplyr::select(one_of(c(varlist, 'WTSB2YR'))) %>% 
  dplyr::rename(wt_lab_orig = WTSB2YR)
  # 09-10
pahas_orig10 <- demo10 %>% 
  left_join(bmi10, by = "SEQN") %>% 
  left_join(albcr10, by = "SEQN") %>% 
  left_join(pah10, by = "SEQN") %>% 
  left_join(mcq10, by = "SEQN") %>% 
  left_join(cot10, by = "SEQN") %>% 
  left_join(hiq10, by = "SEQN") %>% 
  dplyr::mutate(yr = '2009-2010') %>% 
  dplyr::select(one_of(c(varlist, 'WTSB2YR'))) %>% 
  dplyr::rename(wt_lab_orig = WTSB2YR)
  # 11-12
pahas_orig12 <- demo12 %>% 
  left_join(bmi12, by = "SEQN") %>% 
  left_join(albcr12, by = "SEQN") %>% 
  left_join(pah12, by = "SEQN") %>% 
  left_join(mcq12, by = "SEQN") %>% 
  left_join(cot12, by = "SEQN") %>% 
  left_join(hiq12, by = "SEQN") %>% 
  dplyr::mutate(yr = '2011-2012') %>% 
  dplyr::select(one_of(c(varlist, 'WTSA2YR'))) %>% 
  dplyr::rename(wt_lab_orig = WTSA2YR)
  # 13-14
pahas_orig14 <- demo14 %>% 
  left_join(bmi14, by = "SEQN") %>% 
  left_join(albcr14, by = "SEQN") %>% 
  left_join(pah14, by = "SEQN") %>% 
  left_join(mcq14, by = "SEQN") %>% 
  left_join(cot14, by = "SEQN") %>% 
  left_join(hiq14, by = "SEQN") %>% 
  dplyr::mutate(yr = '2013-2014') %>% 
  dplyr::select(one_of(c(varlist, 'WTSA2YR'))) %>% 
  dplyr::rename(wt_lab_orig = WTSA2YR)
  # 15-16
pahas_orig16 <- demo16 %>%  
  left_join(bmi16, by = "SEQN") %>% 
  left_join(albcr16, by = "SEQN") %>% 
  left_join(pah16, by = "SEQN") %>% 
  left_join(mcq16, by = "SEQN") %>%  
  left_join(cot16, by = "SEQN") %>% 
  left_join(hiq16, by = "SEQN") %>% 
  dplyr::mutate(yr = '2015-2016') %>% 
  select(one_of(c(varlist, 'WTSA2YR'))) %>% 
  dplyr::rename(wt_lab_orig = WTSA2YR)
### Concatenate
pahas_orig <- rbind(pahas_orig06, pahas_orig08, pahas_orig10, pahas_orig12, pahas_orig14, pahas_orig16)
  # Check dupes
pahas_orig %>% nrow()
pahas_orig %>% distinct(SEQN) %>% nrow() # None found

# Side-- count missing BMI
is.na(pahas_orig$BMXBMI) %>% sum #191 missing BMI

############################################
########## COUNT N AND EXCLUSIONS ########## 
############################################
  # Initial: 60,936
pahas_orig %>% nrow() 
  # Exclusions
    # Not sampled 
pahas_orig %>% dplyr::filter(is.na(flag_lab)) %>% nrow() #44,133, 72.43%
(44133/60936)*100
pahas_orig <- pahas_orig %>% filter(flag_lab == 1) #16,803
pahas_orig %>% nrow() 
    # Missing or 0 lab weight
pahas_orig %>% filter(is.na(wt_lab_orig)|wt_lab_orig==0) %>% nrow() # 253
(253/16803)*100
pahas_orig <- pahas_orig %>% filter(!is.na(wt_lab_orig) & wt_lab_orig != 0)
pahas_orig %>% nrow() # 16,550


#################################
########## RECODE DATA ##########
#################################
  # Cutpoints
age_cutpt <- c(0, 11, 18, 45, 70, 100) # Max 85
age_labels <- c("<12", "12-18", "19-45", "46-70", ">70")

age_cutpt2 <- c(0, 18, 100) # Max 85
age_labels2 <- c("<18", ">19")

age_cutpt4 <- c(0, 18, 45, 70, 100) # Max 85
age_labels4 <- c("<18", "19-45", "46-70", ">70")

age_cutpt6 <- c(0, 4, 11, 18, 45, 70, 100) # Max 85
age_labels6 <- c("3-4", "5-11", "12-18", "19-45", "46-70", ">70")

### Recode 1: Variables
pahas_recode <- pahas_orig  %>% 
  dplyr::mutate(age = RIDAGEYR,
                agecat = cut(RIDAGEYR, include.lowest = TRUE, breaks = age_cutpt, labels = age_labels),
                agecat2 = cut(RIDAGEYR, include.lowest = TRUE, breaks = age_cutpt2, labels = age_labels2),
                agecat4 = cut(RIDAGEYR, include.lowest = TRUE, breaks = age_cutpt4, labels = age_labels4),
                agecat6 = cut(RIDAGEYR, include.lowest = TRUE, breaks = age_cutpt6, labels = age_labels6),
                gender = case_when(RIAGENDR == 1 ~ "M",
                                   RIAGENDR == 2 ~ "F") %>% as.factor(),
                # race = RIDRETH1 %>% as.factor, 
                race = case_when(RIDRETH1 == 1 ~ "Mexican American",
                                 RIDRETH1 == 2 ~ "Other Hispanic",
                                 RIDRETH1 == 3 ~ "Non-Hispanic White",
                                 RIDRETH1 == 4 ~ "Non-Hispanic Black", 
                                 RIDRETH1 == 5 ~ "Other Race (including multiracial)") %>% as.factor(),
                race = relevel(race, ref = "Non-Hispanic White"),
                pir = INDFMPIR,
                insurance = dplyr::recode(HIQ011, "1" = 1, "2" = 0, .default = as.double(NA)),    
                bmi = BMXBMI,
                cot = LBXCOT,
                cotcat = case_when(cot < 1 ~ 'Nonsmoker (<1)',
                                   cot >= 1  & cot <= 10 ~ 'Nonsmoker, heavy (1-10)',
                                   cot > 10 ~ 'Active Smoker (>10)') %>% factor(levels = c('Nonsmoker (<1)','Nonsmoker, heavy (1-10)','Active Smoker (>10)')),
                pah = URXP10,
                pahz = scale(pah),
                ratio_pah = (URXP10/URXCRS),
                ratio_pahz = scale(ratio_pah),
                creat = URXCRS,
                creat_log = log(creat),
                asthma = dplyr::recode(MCQ010, "1" = 1, "2" = 0, .default = as.double(NA)), # Other values, 7 and 9 become missing
                all = 1,
                wt_lab = wt_lab_orig/6
                ) 

### Generate adjusted PAH 
  # Survey design (pre), to get adjusted creatinine
pahas_design_pre <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_lab, nest = TRUE, data = pahas_recode)
  # Log quasibinomial
mod_creat_log <- svyglm(creat_log ~ age + factor(gender) + bmi,
                        data = pahas_recode, family = gaussian(link = "identity"), design = pahas_design_pre) # linear model
result_creat_log <- data.frame(SEQN = mod_creat_log$survey.design$variables$SEQN,
                               creat_p = (mod_creat_log %>% predict() %>% exp() %>% as.vector()))       
    # Sensitivity analysis with race
mod_creat_race_log <- svyglm(creat_log ~ age + factor(gender) + factor(race) + bmi,
                        data = pahas_recode, family = gaussian(link = "identity"), design = pahas_design_pre) # linear model
result_creat_race_log <- data.frame(SEQN = mod_creat_race_log$survey.design$variables$SEQN,
                                       creat_race_p = (mod_creat_race_log %>% predict() %>% exp() %>% as.vector()))       
    # Sensitivity analysis without complex survey design
mod_creat_no_survey_log <- glm(creat_log ~ age + factor(gender) + bmi,
                          data = pahas_recode, family = gaussian(link = "identity")) # linear model
result_creat_no_survey_log <- data.frame(SEQN = pahas_recode %>% dplyr::filter(!is.na(creat_log) & !is.na(age) & !is.na(gender) & !is.na(bmi)) %>% select(SEQN),
                                         creat_no_survey_p = (mod_creat_no_survey_log %>% predict() %>% exp() %>% as.vector()))       

### Residuals from model for Adjustment method 5
mod_creat_pah <- svyglm(creat ~ pah, family = gaussian(link = "identity"), data = pahas_recode, design = pahas_design_pre) # linear model
result_creat_pah <- data.frame(SEQN = mod_creat_pah$survey.design$variables$SEQN,
                               creat_resid = residuals(mod_creat_pah))     
### Add in Variables
pahas_recode <- pahas_recode %>% 
  dplyr::left_join(result_creat_log, by = "SEQN") %>% 
  dplyr::left_join(result_creat_race_log, by = "SEQN") %>% # For sensitivity analysis, without race
  dplyr::left_join(result_creat_no_survey_log, by = "SEQN") %>% # For sensitivity analysis, without survey design
  dplyr::left_join(result_creat_pah, by = "SEQN") %>% 
  dplyr::mutate(creat_diff = creat - creat_p,
                ratio_pahadj = (URXP10/(creat/creat_p)),
                ratio_pahadjz = scale(ratio_pahadj),
                ratio_pahadj_race = (URXP10/(creat/creat_race_p)),
                ratio_pahadjz_race = scale(ratio_pahadj_race),
                ratio_pahadj_no_survey = (URXP10/(creat/creat_no_survey_p)),
                ratio_pahadjz_no_survey = scale(ratio_pahadj_no_survey))                
pahas_design_pre <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_lab, nest = TRUE, data = pahas_recode)

### Recode 2: Quartiles of standardized ratio PAH/creatinine (+ adjusted creatinine) & Quartiles PIR
  # Also: mean value of each quartile
q_labels <- c("Q1","Q2","Q3","Q4")
tert_labels <- c("T1","T2","T3")
med_labels <- c("< Median",">= Median")
  # PAH
pah_cutpt <- svyquantile(~pahz, design = pahas_design_pre, quantile = c(0, 0.25, 0.50, 0.75, ci = TRUE), na.rm = TRUE)
pah_rollmeans <- rollmean(pah_cutpt %>% as.vector(), 2)
ratio_pah_cutpt <- svyquantile(~ratio_pahz, design = pahas_design_pre, quantile = c(0, 0.25, 0.50, 0.75, ci = TRUE), na.rm = TRUE)
ratio_pah_rollmeans <- rollmean(ratio_pah_cutpt %>% as.vector(), 2)
ratio_pahadj_cutpt <- svyquantile(~ratio_pahadjz, design = pahas_design_pre, quantile = c(0, 0.25, 0.50, 0.75, ci = TRUE), na.rm = TRUE)
ratio_pahadj_rollmeans <- rollmean(ratio_pahadj_cutpt %>% as.vector(), 2)
  # PIR
pir_cutpt <- svyquantile(~pir, design = pahas_design_pre, quantile = c(0, 0.25, 0.50, 0.75, ci = TRUE), na.rm = TRUE)
pir_cutpt_tert <- svyquantile(~pir, design = pahas_design_pre, quantile = c(0, 1/3, 2/3, ci = TRUE), na.rm = TRUE)
pir_cutpt_median <- svyquantile(~pir, design = pahas_design_pre, quantile = c(0, 0.5, ci = TRUE), na.rm = TRUE)
  # Recode
pahas_recode <- pahas_recode %>% 
  dplyr::mutate(# PAH
                pahq = cut(pahz, breaks = c(pah_cutpt$pahz[,1]), labels = q_labels),
                ratio_pahq = cut(ratio_pahz, breaks = c(ratio_pah_cutpt$ratio_pahz[,1]), labels = q_labels),
                ratio_pahadjq = cut(ratio_pahadjz, breaks = c(ratio_pahadj_cutpt$ratio_pahadjz[,1]), labels = q_labels),
                ratio_pahq_val = recode(ratio_pahq, "Q1" = ratio_pah_rollmeans[1], "Q2" = ratio_pah_rollmeans[2], "Q3" = ratio_pah_rollmeans[3], "Q4" = ratio_pah_rollmeans[4]),
                ratio_pahadjq_val = recode(ratio_pahadjq, "Q1" = ratio_pahadj_rollmeans[1], "Q2" = ratio_pahadj_rollmeans[2], "Q3" = ratio_pahadj_rollmeans[3], "Q4" = ratio_pahadj_rollmeans[4]),
                pahq_val = recode(pahq, "Q1" = pah_rollmeans[1], "Q2" = pah_rollmeans[2], "Q3" = pah_rollmeans[3], "Q4" = pah_rollmeans[4]),
                # PIR
                pirq = cut(pir, include.lowest = TRUE, breaks = pir_cutpt$pir[,1], labels = q_labels),
                pirtert = cut(pir, include.lowest = TRUE, breaks = pir_cutpt_tert$pir[,1], labels = tert_labels),
                pirmedian = cut(pir, include.lowest = TRUE, breaks = pir_cutpt_median$pir[,1], labels = med_labels),
                )
pahas_design_pre <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_lab, nest = TRUE, data = pahas_recode)

#######################################
########## CLEAN - MEDS DATA ########## 
#######################################
### Lookup meds
rx_comb %>% 
  dplyr::filter(str_detect(RXDDRUG,'PRED')) %>%
  group_by(RXDDRUG, RXDDRGID) %>% 
  tally()
rx_comb %>% 
  dplyr::filter(str_detect(RXDDRGID, c('d03112'))) %>% #grepl('00206', RXDDRGID)
  group_by(RXDDRUG, RXDDRGID) %>% 
  tally()
### Clean RX data
# Select variables
varlist_rx <- c("SEQN", "RXDUSE", "RXDDRUG", "RXDDRGID", "RXQSEEN", "RXDDAYS", "RXDCOUNT")
rx06_prep <- rx06 %>% dplyr::select(one_of(varlist_rx))
rx08_prep <- rx08 %>% dplyr::select(one_of(varlist_rx))
rx10_prep <- rx10 %>% dplyr::select(one_of(varlist_rx))
rx12_prep <- rx12 %>% dplyr::select(one_of(varlist_rx))
rx14_prep <- rx14 %>% dplyr::select(one_of(varlist_rx))
rx16_prep <- rx16 %>% dplyr::select(one_of(varlist_rx))
# Concatenate
rx_comb <- rbind(rx06_prep, rx08_prep, rx10_prep, rx12_prep, rx14_prep, rx16_prep)
# Flag meds, then summarize by study participant
rx_clean <- rx_comb %>% 
  dplyr::mutate(alb = ifelse(str_detect(RXDDRGID, 'd00749'), 1, 0),
                lalb = ifelse(str_detect(RXDDRGID, 'd04427'), 1, 0),
                pbut = ifelse(str_detect(RXDDRGID, 'd00755'), 1, 0),
                mprediso = ifelse(str_detect(RXDDRGID, 'd00293'), 1, 0),
                prediso = ifelse(str_detect(RXDDRGID, 'd00084'), 1, 0),
                pred = ifelse(str_detect(RXDDRGID, 'd00350'), 1, 0),
                dex = ifelse(str_detect(RXDDRGID, 'd00206'), 1, 0),
                becl = ifelse(str_detect(RXDDRGID, 'd00760'), 1, 0),
                bud = ifelse(str_detect(RXDDRGID, 'd04276'), 1, 0),
                bud = ifelse(str_detect(RXDDRGID, 'd04795'), 1, bud),
                flun = ifelse(str_detect(RXDDRGID, 'd00761'), 1, 0),
                flut = ifelse(str_detect(RXDDRGID, 'd01296'), 1, 0),
                flut = ifelse(str_detect(RXDDRGID, 'd04611'), 1, flut),
                flut = ifelse(str_detect(RXDDRGID, 'd08100'), 1, flut),
                tria = ifelse(str_detect(RXDDRGID, 'd00620'), 1, 0),
                crom = ifelse(str_detect(RXDDRGID, 'd00200'), 1, 0),
                mont = ifelse(str_detect(RXDDRGID, 'd04289'), 1, 0),
                zaf = ifelse(str_detect(RXDDRGID, 'd04053'), 1, 0),
                form = ifelse(str_detect(RXDDRGID, 'd04572'), 1, 0),
                form = ifelse(str_detect(RXDDRGID, 'd07660'), 1, form),
                form = ifelse(str_detect(RXDDRGID, 'd04795'), 1, form),
                sal = ifelse(str_detect(RXDDRGID, 'd03759'), 1, 0),
                sal = ifelse(str_detect(RXDDRGID, 'd04611'), 1, sal),
                theo = ifelse(str_detect(RXDDRGID, 'd00142'), 1, 0),
                mom = ifelse(str_detect(RXDDRGID, 'd05262'), 1, 0),
                mom = ifelse(str_detect(RXDDRGID, 'd07660'), 1, mom),
                nedo = ifelse(str_detect(RXDDRGID, 'd03112'), 1, 0),
                oma = ifelse(str_detect(RXDDRGID, 'd04881'), 1, 0)) %>% 
  group_by(SEQN) %>% 
  summarise(med_short = as.numeric(sum(alb, lalb, pbut, mprediso, prediso, pred, dex)>0),
            med_short_steroid = as.numeric(sum(mprediso, prediso, pred, dex)>0),
            med_short_steroid_no = as.numeric(sum(alb, lalb, pbut)>0),
            med_long = as.numeric(sum(becl, bud, flun, flut, tria, crom, mont, zaf, form, sal, theo, mom, nedo, oma)>0), # add in dex
            alb = as.numeric(sum(alb) > 0),
            lalb = as.numeric(sum(lalb) > 0),
            pbut = as.numeric(sum(pbut) > 0),
            mprediso = as.numeric(sum(mprediso) > 0),
            prediso = as.numeric(sum(prediso) > 0),
            pred = as.numeric(sum(pred) > 0),
            dex = as.numeric(sum(dex) > 0),
            becl = as.numeric(sum(becl) > 0),
            bud = as.numeric(sum(bud) > 0),
            flun = as.numeric(sum(flun) > 0),
            flut = as.numeric(sum(flut) > 0),
            tria = as.numeric(sum(tria) > 0),
            crom = as.numeric(sum(crom) > 0),
            mont = as.numeric(sum(mont) > 0),
            zaf = as.numeric(sum(zaf) > 0),
            form = as.numeric(sum(form) > 0),
            sal = as.numeric(sum(sal) > 0),
            theo = as.numeric(sum(theo) > 0),
            mom = as.numeric(sum(mom) > 0), 
            nedo = as.numeric(sum(nedo) > 0),
            oma = as.numeric(sum(oma) > 0)) # No category yet #
    # Double check coding
xtabs(~ med_short + alb + lalb + pbut + mprediso + prediso + pred + dex, data = rx_clean) %>% freqlist() %>% summary()
xtabs(~ med_short_steroid + alb + lalb + pbut + mprediso + prediso + pred + dex, data = rx_clean) %>% freqlist() %>% summary()
xtabs(~ med_short_steroid_no + alb + lalb + pbut + mprediso + prediso + pred + dex, data = rx_clean) %>% freqlist() %>% summary()
xtabs(~ med_long + becl + bud + flun + flut + tria + crom + mont + zaf + form + sal + theo, data = rx_clean) %>% freqlist() %>% summary()

### MERGE with main data
pahas <- pahas_recode %>% 
  left_join(rx_clean, by = "SEQN")


#########################################
########## EXPORT CLEANED DATA ##########
#########################################
setwd(here("data/processed"))
saveRDS(pahas,'nhanes-pah-asthma-analysis.RDS')


















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
