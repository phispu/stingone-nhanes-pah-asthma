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