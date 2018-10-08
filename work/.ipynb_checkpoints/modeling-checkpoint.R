##
'
Model Building & Evaluation ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
05 Oct 2018 - Finalized

This script...

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.5.1
'
##

setwd('..')

'
Load libraries & dataset
'
library(glmnet)

data_imp <- read.csv('data/final_dataset_cart_imputed.csv')

# structure list
str(data_imp, list.len=ncol(data_imp))

'
Test Model - All Variables
'
# copy data
data_copy <- data_imp

# extract observation weights
obs_weights <- data_copy$obs_weight
data_copy$obs_weight <- NULL

# define test model x & y
y <- as.numeric(data_copy$act_demented)
x <- model.matrix(~.-act_demented, data=data_copy)
x <- x[ , -1]

test_model <- glmnet(x, as.factor(y), family='binomial', weights=obs_weights, alpha=1, nlambda=100)

cv.test_model <- cv.glmnet(x, y, weights=obs_weights, family = 'binomial', type.measure = 'deviance')
plot(cv.test_model)

print(cv.test_model)

cv.test_model$lambda.min

coef(cv.test_model, s='lambda.min')

'
Drop Diagnosis Variables
'
# remove dsm_iv_clinical_diagnosis, nincds_arda_diagnosis, nia_reagan, braak, cerad
diagnosis_vars <- c('dsm_iv_clinical_diagnosis',
                    'nincds_arda_diagnosis',
                    'cerad',
                    'braak',
                    'nia_reagan')

data_no_diag <- data_copy[ , !names(data_copy) %in% diagnosis_vars]

y <- as.numeric(data_no_diag$act_demented)
x <- model.matrix(~.-act_demented, data=data_no_diag)
x <- x[ , -1]

cv.test_model_no_diag <- cv.glmnet(x, y, weights=obs_weights, family = 'binomial', lambda=seq(0,1,0.001), type.measure = 'deviance')
plot(cv.test_model_no_diag)

print(cv.test_model_no_diag)

cv.test_model_no_diag$lambda.min

coef(cv.test_model_no_diag, s='lambda.min')

# ElasticNet
elas.cv.test_model_no_diag <- cv.glmnet(x, y, weights=obs_weights, family = 'binomial', alpha=0.15, type.measure = 'deviance')
plot(elas.cv.test_model_no_diag)

print(elas.cv.test_model_no_diag)

elas.cv.test_model_no_diag$lambda.min

coef(elas.cv.test_model_no_diag, s='lambda.min')