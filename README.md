# SepsisScores
R package to estimate mortality with SOFA, qSOFA, and APACHEII scores.

### SepsisScores functions:
- APACHEII_missingness
- APACHEII_mortality  
- APACHEII_points  
- SOFA_missingness  
- SOFA_mortality  
- SOFA_points

### Descriptions:
- APACHEII_points: takes a dataframe with critical care variables named by the Grady Memorial Hospital COVID chart review data dictionary and produces a dataframe with record_id and points for each component of the APACHEII score.

- APACHEII_missingness: takes APACHEII_points output as input dataframe and calculates number and percent missing by variable, including APACHEII score.

- APACHEII_mortality: takes APACHEII_points output table as input dataframe and calculates estimated in-hospital mortality by operative status.  Output is a dataframe with record_id, post-operative status, and non-operative status mortality predictions.

- SOFA_points: takes a dataframe with critical care variables named by the Grady Memorial Hospital COVID chart review data dictionary and produces a dataframe with record_id and points for each component of the SOFA & qSOFA scores.

- SOFA_missingness: takes SOFA_points output table as input dataframe and calculates number and percent missing by variable, including APACHEII score.

- SOFA_mortality: takes SOFA_points output table as input dataframe and calculates estimated in-hospital mortality by operative status.  Output is a dataframe with record_id, SOFA & qSOFA mortality predictions.
