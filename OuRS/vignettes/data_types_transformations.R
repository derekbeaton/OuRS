## ----setup, include=FALSE, warning=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(kableExtra)
library(ours)
library(magrittr)

## ----adni---------------------------------------------------------------------
data("synthetic_ADNI")

## ----categorical_head, echo = T, eval = F-------------------------------------
#  
#  # Defining and viewing the categorical data set
#  categorical_data <- synthetic_ADNI[, c("dx", "ptgender")]
#  head(categorical_data)
#  

## ----categorical_head_nice_table, echo = F------------------------------------

# Defining and viewing the categorical data set
categorical_data <- synthetic_ADNI[, c("dx", "ptgender")]
head(categorical_data) %>%
  kbl() %>%
  kable_minimal()


## ----categorical_table--------------------------------------------------------

# Examining groups within each categorical variable
apply(categorical_data, 2, table)


## ----categorical_transform, echo = T, eval = F--------------------------------
#  
#  categorical_data_disjunct <- disjunctive_coding(categorical_data)
#  head(categorical_data_disjunct)
#  

## ----categorical_transform_nice_table, echo = F-------------------------------

categorical_data_disjunct <- disjunctive_coding(categorical_data)
head(categorical_data_disjunct) %>%
  kbl() %>%
  kable_minimal() %>%
  # add_header_above(c(" ", "dx" = 3, "ptgender" = 2))
  column_spec(2:4, background = rgb(t(col2rgb("mediumorchid3")), alpha = 255, maxColorValue = 255)) %>%
  column_spec(5:6, background = rgb(t(col2rgb("olivedrab3")), alpha = 255, maxColorValue = 255))



## ----categorical_properties---------------------------------------------------

# Confirming each column sums to the observed frequency of each group
colSums_categorical <- colSums(categorical_data_disjunct)
colSums_categorical

# Confirming the sum of the column sums within each original variable equals the number of rows
sum(colSums_categorical[c("dx.MCI", "dx.Dementia", "dx.CN")])
sum(colSums_categorical[c("ptgender.Female", "ptgender.Male")])

# Confirming each row sums to the number of original variables
rowSums(head(categorical_data_disjunct))


## ----ordinal_head, echo = T, eval = F-----------------------------------------
#  
#  # Defining and viewing the ordinal data set
#  ordinal_data <- synthetic_ADNI[, c("gdtotal", "hmscore", "cdrsb")]
#  head(ordinal_data)
#  

## ----ordinal_head_nice_table, echo = F, eval = T------------------------------

# Defining and viewing the ordinal data set
ordinal_data <- synthetic_ADNI[, c("gdtotal", "hmscore", "cdrsb")]
head(ordinal_data) %>%
  kbl() %>%
  kable_minimal()


## ----ordinal_range------------------------------------------------------------

# Examining values within each ordinal variable
apply(ordinal_data, 2, function(i){sort(unique(i))})


## ----ordinal_transform, echo = T, eval = F------------------------------------
#  
#  # Ordinal data set after transformation
#  ordinal_data_poles <- thermometer_coding(ordinal_data)
#  head(round(ordinal_data_poles, 3))
#  

## ----ordinal_transform_nice_table, echo = F, eval = T-------------------------

# Ordinal data set after transformation
ordinal_data_poles <- thermometer_coding(ordinal_data)
head(round(ordinal_data_poles, 3)) %>%
  kbl() %>%
  kable_minimal() %>%
  column_spec(., column = 2:4, background = rgb(t(col2rgb("mediumorchid3")), alpha = 255, maxColorValue = 255)) %>%
  column_spec(., column = 5:6, background = rgb(t(col2rgb("olivedrab3")), alpha = 255, maxColorValue = 255))


## ----ordinal_properties-------------------------------------------------------

# Confirming the sum of each set of poles equals 1 for each participant
rowSums(ordinal_data_poles[, c("gdtotal+", "gdtotal-")])[1:3]
rowSums(ordinal_data_poles[, c("hmscore+", "hmscore-")])[1:3]
rowSums(ordinal_data_poles[, c("cdrsb+", "cdrsb-")])[1:3]

# Confirming each row sums to the number of original variables
rowSums(ordinal_data_poles)[1:3]

# Confirming the sum of column sums of each set of poles equals the number of rows
colSums_ordinal <- colSums(ordinal_data_poles)
sum(colSums_ordinal[c("gdtotal+", "gdtotal-")])
sum(colSums_ordinal[c("hmscore+", "hmscore-")])
sum(colSums_ordinal[c("cdrsb+", "cdrsb-")])


## ----continuous_mixed---------------------------------------------------------

# Defining and viewing the continuous data set
continuous_data <- synthetic_ADNI[, c("age", "mpacctrailsb","mpaccdigit", "hippocampus", "icv")]
round(continuous_data[1:3, ], 3)

# Defining and viewing the mixed data set
mixed_data <- cbind(categorical_data, ordinal_data, round(continuous_data, 3))
mixed_data[1:3, ]


## ----continuous_transform-----------------------------------------------------

# Continuous data set after transformation
continuous_data_escofier <- escofier_coding(continuous_data, center = TRUE, scale = TRUE)
round(continuous_data_escofier[1:3, ], 3)


## ----continuous_properties----------------------------------------------------

# Confirming the sum of each set of original variables equals 1 for each participant
rowSums(continuous_data_escofier [, c("age+", "age-")])[1:3]
rowSums(continuous_data_escofier [, c("mpacctrailsb+", "mpacctrailsb-")])[1:3]
rowSums(continuous_data_escofier [, c("mpaccdigit+", "mpaccdigit-")])[1:3]
rowSums(continuous_data_escofier [, c("hippocampus+", "hippocampus-")])[1:3]
rowSums(continuous_data_escofier [, c("icv+", "icv-")])[1:3]

# Confirming each row sums to the number of original variables
rowSums(continuous_data_escofier)[1:3]

# Confirming the sum of column sums of each set of poles equals the number of rows
# note that each column sum equals 311.5, which is the number of rows divided by 2
colSums_continuous <- colSums(continuous_data_escofier)
sum(colSums_continuous[c("age+", "age-")])
sum(colSums_continuous[c("mpacctrailsb+", "mpacctrailsb-")])
sum(colSums_continuous[c("mpaccdigit+", "mpaccdigit-")])
sum(colSums_continuous[c("hippocampus+", "hippocampus-")])
sum(colSums_continuous[c("icv+", "icv-")])


