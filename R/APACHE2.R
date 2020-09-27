# APACHE II Functions

################################################################################
# test data
# load("/Users/NewNasty/Box/COVID/Patients/Code/Data/covid_master.Rdata")
#
# covid1 <- covid %>%
#   filter(is.na(researcher_presentation)==F)%>%
#   dplyr::select(record_id,
#                 abg_done,
#                 abg_pao2,
#                 abg_paco2,
#                 abg_ph,
#                 ams,
#                 age,
#                 platelets,
#                 creatinine,
#                 na,
#                 k,
#                 hct,
#                 wbc,
#                 hr_max,
#                 temp_max,
#                 dopamine_sofa,
#                 ed_date,
#                 epinepherine_sofa,
#                 ett,
#                 gcs,
#                 intubation_date,
#                 min_map,
#                 max_map,
#                 niv,
#                 norepinephrine_sofa,
#                 rr_max,
#                 sao2_fio2,
#                 sao2_min,
#                 sbp_min,
#                 tbili)



#################################################################################
# Calculated variable functions

# create function to calculate Aa gradient with default atmospheric pressure
calc_A_a_gradient <- function(Pa_O2, FiO2, PaCO2, Age,
                              at_pressure=760, H20_pressure = 47){
  grad = (FiO2 * (at_pressure - H20_pressure)) - ((PaCO2/0.8) - Pa_O2)
}

# create a function to classify most extreme MAP values
MAP_extreme <- function(max_map, min_map){
  High_MAP_points <- case_when(
    max_map >=160                                   ~ 4,
    max_map >=130 & max_map <160        ~ 3,
    max_map >=110 & max_map <130        ~ 2,
    max_map >= 70 & max_map <110        ~ 0
  )
  Low_MAP_points <- case_when(
    min_map >= 70 & min_map <110        ~ 0,
    min_map >= 50 & min_map <70         ~ 2,
    min_map <  50                                   ~ 4
  )
  max_MAP_points <- ifelse(High_MAP_points > Low_MAP_points,
                           High_MAP_points,
                           Low_MAP_points)
  return(max_MAP_points)
}

# create a function to classify oxygenation
calc_oxygenation <- function(abg_pao2, sao2_fio2, aagrad){
  oxy_points <- ifelse(abg_pao2 > 50,
         case_when((sao2_fio2 < .5 & abg_pao2 > 70) |
                                   (aagrad <200 & sao2_fio2 > 50)   ~ 0,
                                 abg_pao2 >60 & abg_pao2 <=70       ~ 1,
                                 abg_pao2 >= 55 & abg_pao2 <= 60    ~ 3,
                                 abg_pao2 < 55                      ~ 4),
         case_when(aagrad > 499                  ~ 4,
                    aagrad >=350 & aagrad <= 499 ~ 3,
                    aagrad >=200 & aagrad <=349  ~ 2))
  return(oxy_points)
}



#################################################################################
# Components of APACHE II score pulled from MD-CALC
# https://www.mdcalc.com/apache-ii-score

# score is calculated by worst value from initial 24 hours in the ICU

######################################
# Calculated variables:
# 01 - PaO2/FiO2
# 02 - MAP value OR pressor dose
# 03 - Predicted mortality if initial
# 04 - Predicted mortality if highest
# 05 - Predicted mortality from mean SOFA
######################################

######################################
# Classified variables:
# 00 - Immunocompromise / chronic organ insufficiency
# Biopsy proven cirrhosis, documented portal hypertension,
# episodes of prior GI bleeding attributed to portal HTN,
# Prior episodes of hepatic failure / encephalopathy / coma
#     possible variable: cirrhosis
#     potentially missing: encephalopathy, GI bleeding, portal HTN

# NYHA Class IV Heart Failure
#     possible variable: chf
#     potentially missing: stage of heart failure

# Chronic restrictive, obstructive, or vascular disease,
# Documented chronic hypoxia, hypercapnia, secondary polycythemia,
# severe pulmonary hypertension (>40 mmHg), respirator dependency
#     possible variable: copd, other_lung

# Chronic dialysis
#     possible variable: chf

# Immunosuppression: chemotherapy, radiation, high-dose steroids
#     possible variable: transplant
#     potentially missing: chemo, radiation, steroid dose / cutoff

# these values get collapsed as binary yes/no, then stratified as post op or not

# 01 - Age
# 02 - Temperature
# 03 - Mean arterial pressure
# 04 - pH
# 05 - Heart rate
# 06 - Respiratory rate
# 07 - Na
# 08 - K
# 09 - Cr
# 10 - Acute Kidney Injury (Clinical judgment)
# 11 - Hematocrit
# 12 - White blood cell count
# 13 - Glasgow coma scale
# 14 - FiO2

# Potential data issues for calculation of APACHE II scores:
# - missingness of pH for ABG values?
# - AKI as clinically characterized
# - temperature assumed rectal temp?
# - immunocompromise or severe organ insufficiency?
#   stratified by: nonoperative or emergency post-op vs elective post-op vs not

######################################
# examine extreme aagradient values
# covid1 %>% filter(!is.na(aagrad) & aagrad>800) %>%
#   select(abg_paco2, abg_pao2, sao2_fio2, age) %>%
#   View()

######################################

#' Assign APACHE II points & calculate score
#'
#' This function takes a cleaned dataframe with critical care data
#' and assumes variable names as defined in the COVID project data
#' dictionary but can be modified with some elbow grease.
#'
#' The APACHEII_points function classifies data, assigns an APACHEII
#' score, and returns a data frame with these results.
#' @param sepsis_data is a dataframe with critical care data
#' @return dataframe with APACHE II score & component point values.
#' @export
APACHEII_points <- function(sepsis_data){

  # calculate A-a gradient as aagrad
  sepsis_data <- sepsis_data %>%
    mutate(
    aagrad = calc_A_a_gradient(
    Pa_O2 = abg_pao2,
    FiO2 = sao2_fio2,
    PaCO2 = abg_paco2,
    Age = age))

  apache_points <- sepsis_data %>% mutate(
  ######################################
  # Assign APACHEII points
  # age, temperature, MAP, pH, heart rate,
  # respiratory rate, Na, K, Cr, AKI,
  # Hct, WBC ct, GCS, FiO2

  age_points = case_when(
    age <=44                            ~ 0,
    age > 45 & age <55                  ~ 2,
    age >=55 & age <65                  ~ 3,
    age >=65 & age <75                  ~ 5,
    age >=75                            ~ 6
  ),
  temp_points = case_when(                     #apparently rectal temp presumed?
    temp_max >=41                           ~ 4,
    temp_max >=39 & temp_max <41                ~ 3,
    temp_max >=38.5 & temp_max <39              ~ 1,
    temp_max >=36 & temp_max <38.5              ~ 0,
    temp_max >=34 & temp_max <36                ~ 1,
    temp_max >=32 & temp_max <34                ~ 2,
    temp_max >=30 & temp_max <32                ~ 3,
    temp_max < 30                           ~ 4
  ),
  MAP_points = MAP_extreme(sepsis_data$max_map, sepsis_data$min_map),
  HR_points = case_when(
    hr_max >= 180                       ~ 4,
    hr_max >=140 & hr_max <180          ~ 3,
    hr_max >=110 & hr_max <140          ~ 2,
    hr_max >= 70 & hr_max <110          ~ 0,
    hr_max >= 55 & hr_max <70           ~ 2,
    hr_max >= 40 & hr_max <55           ~ 3,
    hr_max <  55                        ~ 4,
  ),
  RR_points = case_when(
    rr_max >= 50                        ~ 4,
    rr_max >= 35 & rr_max < 50          ~ 3,
    rr_max >= 25 & rr_max < 35          ~ 2,
    rr_max >= 12 & rr_max < 25          ~ 0,
    rr_max >= 10 & rr_max < 12          ~ 2,
    rr_max >=  6 & rr_max < 10          ~ 3,
    rr_max <   6                        ~ 4,
  ),
  Oxy_points = calc_oxygenation(abg_pao2 = abg_pao2,
                                sao2_fio2 = sao2_fio2,
                                aagrad = aagrad
  ),
  ArtpH_points = case_when(
    abg_ph >= 7.7                       ~ 4,
    abg_ph >= 7.6 & abg_ph <7.7         ~ 3,
    abg_ph >= 7.5 & abg_ph <7.6         ~ 1,
    abg_ph >= 7.33 & abg_ph <7.5        ~ 0,
    abg_ph >= 7.25 & abg_ph <7.33       ~ 2,
    abg_ph >= 7.15 & abg_ph <7.25       ~ 3,
    abg_ph < 7.15                       ~ 4
  ),
  Na_points = case_when(
    na >= 180                           ~ 4,
    na >= 160 & na <180                 ~ 3,
    na >= 155 & na <160                 ~ 2,
    na >= 130 & na <150                 ~ 0,
    na >= 120 & na <130                 ~ 2,
    na >= 111 & na <120                 ~ 3,
    na <  111                           ~ 4
  ),
  K_points = case_when(
    k >= 7.0                            ~ 4,
    k >= 6.0 & k <7.0                   ~ 3,
    k >= 5.5 & k <6.0                   ~ 2,
    k >= 3.5 & k <5.5                   ~ 0,
    k >= 3.0 & k <3.5                   ~ 2,
    k >= 2.5 & k <3.0                   ~ 3,
    k <  2.5                            ~ 4
  ),
  cr_points = case_when(
    creatinine >= 3.5                   ~ 4,
      # should be stratified by acuity
      # creatinine > 3.5 & acute = 1    ~ 8,
      # creatinine > 3.5 & acute = 0    ~ 4,

    creatinine >= 2.0 & creatinine <3.5 ~ 3,
      # should be stratified by acuity
      # creatinine > 2.0 & creatinine < 3.5 & acute = 1    ~ 6,
      # creatinine > 2.0 & creatinine < 3.5 & acute = 0    ~ 3,

    creatinine >=1.5 & creatinine < 2    ~ 2,
      # should be stratified by acuity
      # creatinine > 1.5 & creatinine < 2 & acute = 1    ~ 4,
      # creatinine > 1.5 & creatinine < 2 & acute = 0    ~ 2,

    creatinine >= 0.6 & creatinine <1.5  ~ 0,
    creatinine < 0.6                     ~ 2
  ),
  hct_points = case_when(
    hct >= 60                            ~ 4,
    hct >= 50 & hct < 60                 ~ 2,
    hct >= 46 & hct < 50                 ~ 1,
    hct >= 30 & hct < 46                 ~ 0,
    hct >= 20 & hct < 30                 ~ 2,
    hct < 20                             ~ 3
    ),
  wbc_points = case_when(
    wbc >= 40                            ~ 4,
    wbc >= 20 & wbc < 40                 ~ 2,
    wbc >= 15 & wbc < 20                 ~ 1,
    wbc >= 3  & wbc < 15                 ~ 0,
    wbc >= 1  & wbc < 3                  ~ 2,
    wbc < 1                              ~ 4
  ),
  gcs_points = (15 - gcs)
  )

  ## tally APACHE II score
  apache_points$APACHEII.score <- apache_points %>%
    select(age_points, temp_points,
           MAP_points, HR_points, RR_points,
           Oxy_points, ArtpH_points, Na_points,
           K_points, cr_points, hct_points,
           wbc_points, gcs_points) %>%
    rowSums(na.rm=FALSE)

  APACHEII_score <- apache_points %>%
    select(record_id, age_points, temp_points,
           MAP_points, HR_points, RR_points,
           Oxy_points, ArtpH_points, Na_points,
           K_points, cr_points, hct_points,
           wbc_points, gcs_points, APACHEII.score)

  return(APACHEII_score)
}


#' Estimate in-hospital mortality
#'
#' This function takes a APACHE II scores calculated with the
#' APACHEII_points function and optional parameter to indicate
#' operative status and estimates in-hospital mortality rates.
#'
#' The output is a dataframe with APACHE II scores, non-operative,
#' and operative mortality unless operative status is specified.
#'
#' @param APACHEII_score is a dataframe with APACHEII scores
#' @param op_status is a column name that indicates operative status,
#' which should be coded as 0 = non-operative, 1 = post-operative.
#' @return dataframe with APACHE II score from APACHEII_points function.
#' @export
APACHEII_mortality <- function(APACHEII_scores){

  mort <- APACHEII_scores %>%
    mutate(
      Nonop_mortality = case_when(
        APACHEII.score <= 4                       ~ 0.04,
        APACHEII.score >5  & APACHEII.score <= 9  ~ 0.08,
        APACHEII.score >9  & APACHEII.score <=14  ~ 0.15,
        APACHEII.score >14 & APACHEII.score <=19  ~ 0.25,
        APACHEII.score >19 & APACHEII.score <=24  ~ 0.40,
        APACHEII.score >24 & APACHEII.score <=29  ~ 0.55,
        APACHEII.score >29 & APACHEII.score <=34  ~ 0.73,
        APACHEII.score >34                        ~ 0.85
      ),
      Postop_mortality = case_when(
        APACHEII.score <= 4                       ~ 0.01,
        APACHEII.score >5  & APACHEII.score <= 9  ~ 0.03,
        APACHEII.score >9  & APACHEII.score <=14  ~ 0.07,
        APACHEII.score >14 & APACHEII.score <=19  ~ 0.12,
        APACHEII.score >19 & APACHEII.score <=24  ~ 0.30,
        APACHEII.score >24 & APACHEII.score <=29  ~ 0.35,
        APACHEII.score >29 & APACHEII.score <=34  ~ 0.73,
        APACHEII.score >34                        ~ 0.88
      )
    )

  # mort <- APACHEII_scores %>%
  #   mutate(
  #     Nonop_mortality = ifelse(is.null(op_status),
  #          case_when(
  #            APACHEII.score <= 4                       ~ 0.04,
  #            APACHEII.score >5  & APACHEII.score <= 9  ~ 0.08,
  #            APACHEII.score >9  & APACHEII.score <=14  ~ 0.15,
  #            APACHEII.score >14 & APACHEII.score <=19  ~ 0.25,
  #            APACHEII.score >19 & APACHEII.score <=24  ~ 0.40,
  #            APACHEII.score >24 & APACHEII.score <=29  ~ 0.55,
  #            APACHEII.score >29 & APACHEII.score <=34  ~ 0.73,
  #            APACHEII.score >34                        ~ 0.85
  #          ),
  #          case_when(
  #            op_status == 0 & APACHEII.score <= 4                       ~ 0.04,
  #            op_status == 0 & APACHEII.score >5  & APACHEII.score <= 9  ~ 0.08,
  #            op_status == 0 & APACHEII.score >9  & APACHEII.score <=14  ~ 0.15,
  #            op_status == 0 & APACHEII.score >14 & APACHEII.score <=19  ~ 0.25,
  #            op_status == 0 & APACHEII.score >19 & APACHEII.score <=24  ~ 0.40,
  #            op_status == 0 & APACHEII.score >24 & APACHEII.score <=29  ~ 0.55,
  #            op_status == 0 & APACHEII.score >29 & APACHEII.score <=34  ~ 0.73,
  #            op_status == 0 & APACHEII.score >34                        ~ 0.85
  #          )),
  #     Postop_mortality = ifelse(is.null(op_status),
  #           case_when(
  #             APACHEII.score <= 4                       ~ 0.01,
  #             APACHEII.score >5  & APACHEII.score <= 9  ~ 0.03,
  #             APACHEII.score >9  & APACHEII.score <=14  ~ 0.07,
  #             APACHEII.score >14 & APACHEII.score <=19  ~ 0.12,
  #             APACHEII.score >19 & APACHEII.score <=24  ~ 0.30,
  #             APACHEII.score >24 & APACHEII.score <=29  ~ 0.35,
  #             APACHEII.score >29 & APACHEII.score <=34  ~ 0.73,
  #             APACHEII.score >34                        ~ 0.88
  #           ),
  #           op_status == 1 & APACHEII.score <= 4                       ~ 0.01,
  #           op_status == 1 & APACHEII.score >5  & APACHEII.score <= 9  ~ 0.03,
  #           op_status == 1 & APACHEII.score >9  & APACHEII.score <=14  ~ 0.07,
  #           op_status == 1 & APACHEII.score >14 & APACHEII.score <=19  ~ 0.12,
  #           op_status == 1 & APACHEII.score >19 & APACHEII.score <=24  ~ 0.30,
  #           op_status == 1 & APACHEII.score >24 & APACHEII.score <=29  ~ 0.35,
  #           op_status == 1 & APACHEII.score >29 & APACHEII.score <=34  ~ 0.73,
  #           op_status == 1 & APACHEII.score >34                        ~ 0.88
  #     ))

  apache_mort <- mort %>%
    select(record_id,
           APACHEII.score,
           Nonop_mortality,
           Postop_mortality
           )

 return(apache_mort)
}

#' Summarize APACHEII missingness
#'
#' The APACHEII_missingness function summarizes missingness of
#' APACHEII points & scores that might contextualize estimates.
#' @param APACHEII_scores is a dataframe with APACHEII scores
#' @return dataframe with APACHEII score & component missingness
#' @export
APACHEII_missingness <- function(APACHEII_scores) {
  missingness.apache_points <- APACHEII_scores %>%
    summarise_all(~sum(is.na(.))) %>%
    pivot_longer(everything(), names_to = "APACHEII.point.var") %>%
    mutate(N.missing = value,
           percent.missing = round(100 * value/nrow(APACHEII_scores),
                                   digits=1)) %>%
    select(-value) %>%
    arrange(desc(percent.missing))
  return(missingness.apache_points)
}
