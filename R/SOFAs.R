# data:
# load("/Users/NewNasty/Box/COVID/Patients/Code/Data/covid_master.Rdata")

#################################################################################
# Components of SOFA score pulled from MD-CALC
# https://www.mdcalc.com/sequential-organ-failure-assessment-sofa-score

# score is calculated by worst value for each day, from initial admission or mean

######################################
# Calculated variables:
# 01 - PaO2/FiO2
# 02 - MAP value OR pressor dose
# 03 - Predicted mortality if initial
# 04 - Predicted mortality if highest
# 05 - Predicted mortality from mean SOFA
######################################

######################################
# Define functions for calculated variables

# 01a - Calculate PaO2/FiO2
PaO2_to_FiO2  <- function(pao2, fio2) {
  pf <- as_tibble(pao2 / fio2)                  # calculate pao2/fio2
  pf <- pf %>% rename(PaO2_to_FiO2 = value)
  return(pf)                                    # return value
}

# 01b - Assign SOFA points to PaO2/FiO2 values
PaO2_to_FiO2_points <- function(p2f, vent) {
  pf_pts <- p2f %>% mutate(P2F_points =         # add points by p/f and ETT values
                             case_when(p2f >= 400 ~ 0,
                                       (p2f >=300) & (p2f < 400) ~ 1,
                                       (p2f >=200) & (p2f < 300) ~ 2,
                                       (p2f >=200) & (p2f < 300) & (vent == 0) ~ 2,
                                       (p2f >=100) & (p2f < 200) & (vent == 1) ~ 3,
                                       (p2f >=0) & (p2f < 100) & (vent == 1) ~ 4))
  return(pf_pts)                                # return value
}

# 02 - Assign SOFA points to MAP/pressor values
# for purposes of function coding, will assume pressor doses > 0.1 mcg/kg/min
MAP_points <- function(min_map, dop_cat, epi_cat, norepi_cat) {
  MAP <- tibble(min_map, dop_cat, epi_cat, norepi_cat)

  MAP <- MAP %>% mutate(any_pressor = case_when(
    (dop_cat == 0 & epi_cat == 0 & norepi_cat == 0)  ~ 0,
    (dop_cat == 1 | epi_cat == 1 | norepi_cat == 1)  ~ 1))

  MAP <- MAP %>% mutate(MAPpressor_points = case_when(
    min_map >=70 & any_pressor == 0 ~ 0,
    min_map <70 & any_pressor == 0 ~ 1,
    #       dop_cat == 1 & epi_cat == 0 & norepi_cat ==0 ~ 2,  # dopamine <5 gets 2 points
    dop_cat == 1 & epi_cat == 0 & norepi_cat ==0 ~ 3,  # pressor dose of dopa gets 3 pts
    dop_cat == 1 & (epi_cat == 1 | norepi_cat ==1) ~ 4   # dopa >15 or epi or norepi get 4
  ))
  return(MAP)                                   # return value
}

######################################
# Assign SOFA, qSOFA, & modSOFA points to other variables

#' Calculate SOFA & qSOFA scores
#'
#' This function takes a cleaned dataframe with critical care data
#' and assumes variable names as defined in the COVID project data
#' dictionary but can be modified with some elbow grease.
#'
#' The SOFA_points function classifies data, assigns a SOFA
#' and qSOFA score, and returns a data frame with these results.
#' @param sepsis_data is a dataframe with critical care data
#' @return dataframe with SOFA & qSOFA scores & components
#' @export
SOFA_points <- function(sepsis_data){

sepsis_data <- sepsis_data %>% mutate(
  ######################################
  # Assign SOFA points
  # platelets, GCS, creatinine, total bili
  plt_points = case_when(
    platelets >= 150 ~ 0,
    platelets >= 100 & platelets < 150  ~ 1,
    platelets >= 50  & platelets < 100  ~ 2,
    platelets >= 20  & platelets < 50   ~ 3,
    platelets < 20                      ~ 4,
  ),
  gcs_points = case_when(
    gcs ==15                            ~ 0,
    gcs >=13 & gcs <15                  ~ 1,
    gcs >=10 & gcs <13                  ~ 2,
    gcs >= 6 & gcs <10                  ~ 3,
    gcs < 6                             ~ 4
  ),
  cr_points = case_when(
    creatinine < 1.2                    ~ 0,
    creatinine >=1.2 & creatinine <2    ~ 1,
    creatinine >=2   & creatinine <3.4  ~ 2,
    creatinine >=3.5 & creatinine <5    ~ 3,
    creatinine >= 5                     ~ 4
  ),
  bili_points = case_when(
    tbili < 1.2                         ~ 0,
    tbili >=1.2 & tbili < 2             ~ 1,
    tbili >= 2 & tbili < 6              ~ 2,
    tbili >= 6 & tbili < 12             ~ 3,
    tbili >= 12                         ~ 4
  ),
  ######################################
  # Assign qSOFA points
  # RR, GCS, SBP

  gcs_qsofa_points = case_when(
    gcs <15                            ~ 0,
    gcs ==15                           ~ 1
  ),
  rr_qsofa_points = case_when(
    rr_max < 22                        ~ 0,
    rr_max >= 22                       ~ 1
  ),
  sbp_qsofa_points = case_when(
    sbp_min > 100                      ~ 0,
    sbp_min <=100                      ~ 1
  )
)

######################################
# Calculate SOFA & qSOFA scores by tallying points

# calculate values with functions
pf <- PaO2_to_FiO2(sepsis_data$sao2_min, sepsis_data$sao2_fio2)
pf_points <- PaO2_to_FiO2_points(pf, sepsis_data$ett)

MAP_pts <- MAP_points(sepsis_data$min_map, sepsis_data$dopamine_sofa,
                      sepsis_data$epinepherine_sofa, sepsis_data$norepinephrine_sofa)

# bind & select point columns
covid_points <- bind_cols(sepsis_data, MAP_pts, pf_points)
covid_points <- covid_points %>%
  select(record_id, P2F_points, plt_points, gcs_points,
         bili_points, MAPpressor_points, cr_points,
         gcs_qsofa_points, rr_qsofa_points, sbp_qsofa_points)

# create SOFA & qSOFA score variables
covid_points <- covid_points %>%
  mutate(SOFA_Score = rowSums(.[2:7]),
         qSOFA_Score = rowSums(.[8:10]))

return(covid_points)
}

######################################
# summarize sofa scores
#summary(covid_points)

#' Summarize SOFA & qSOFA missingness
#'
#' The SOFA_missingness function summarizes missingness of SOFA &
#' qSOFA parameters that might contextualize variation in estimates.
#' @param sepsis_data is a dataframe with critical care data
#' @return dataframe with SOFA & qSOFA score & component missingness
#' @export
SOFA_missingness <- function(covid_points) {
missingness.covid_points <- covid_points %>%
  summarise_all(~sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "Sofa.point.var") %>%
  mutate(N.missing = value,
         percent.missing = round(100 * value/nrow(covid_points),
                                 digits=1)) %>%
  select(-value) %>%
  arrange(desc(percent.missing))
  return(missingness.covid_points)
}

# estimate mortality

#' Predict mortality risk with SOFA & qSOFA scores
#'
#' The SOFA_mortality function predicts mortality based on
#' SOFA and qSOFA scores calculated with the SOFA_points function.
#' @param SOFA_scores is a dataframe with scores from SOFA_points
#' @return mortality dataframe with SOFA & qSOFA scores and mortality risk
#' @export
SOFA_mortality <- function(SOFA_scores){
mortality <- SOFA_scores %>%
  mutate(
    sofa_mortality =
      case_when(SOFA_Score < 2                 ~ 0,
                SOFA_Score >=2 & SOFA_Score <4         ~ 0.064,
                SOFA_Score >=4 & SOFA_Score <6         ~ 0.202,
                SOFA_Score >=6 & SOFA_Score <8         ~ 0.215,
                SOFA_Score >=8 & SOFA_Score <10        ~ 0.333,
                SOFA_Score >=10 & SOFA_Score <12       ~ 0.5,
                SOFA_Score >=12 & SOFA_Score <=14      ~ 0.952,
                SOFA_Score >14                         ~ 0.952,
      ),
    qsofa_mortality =
      case_when(qSOFA_Score <2                 ~ "Not high risk",
                qSOFA_Score >= 2               ~ "High risk")
  )
mortality <- mortality %>%
  select(record_id, SOFA_Score, sofa_mortality,
         qSOFA_Score, qsofa_mortality)
  return(mortality)
}
