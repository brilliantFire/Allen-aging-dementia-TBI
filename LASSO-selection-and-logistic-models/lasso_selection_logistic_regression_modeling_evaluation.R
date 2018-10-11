##
'
LASSO Variable Selection, Model Building, & Evaluation ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
05 Oct 2018 - Transfered notebook code to script
10 Oct 2018 - Final models & evaluation (confusion matrices & ROC curves)
11 Oct 2018 - Finalized

This script performs variables selection using LASSO and constructs two logistic regression
models with either all or a subset of the selected variables. Model evaluation with
confusion matrices and ROC curves is included.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.5.1
'
##

setwd('..')

'
Load libraries & dataset
'
library(glmnet)    # for LASSO
library(caret)     # confusion matrices
library(pROC)      # easy ROC curves

data_imp <- read.csv('data/final_dataset_cart_imputed.csv')

# structure list
str(data_imp, list.len=ncol(data_imp))

# remove X
data_imp$X <- NULL

'
Split into Train/Test Subsets
'
# RNG seed
set.seed(617)

# choose 27 random donors to hold out as test set
test <- data_imp[sample(nrow(data_imp), 27), ]
train <- data_imp[!rownames(data_imp) %in% rownames(test), ]

'
Prepare Model Matrices
'
# extract observation weights
train.weights <- as.numeric(train$obs_weight)
train$obs_weight <- NULL

test.weights <- as.numeric(test$obs_weight)
test$obs_weight <- NULL

# define model x & y for train & test sets
train.y <- as.numeric(train$act_demented)
train.x <- model.matrix(~.-act_demented, data=train)
train.x <- train.x[ , -1]

test.y <- as.numeric(test$act_demented)
test.x <- model.matrix(~.-act_demented, data=test)
test.x <- test.x[ , -1]

'
Model 1 - LASSO w/ All Variables
'
# 
cv.fit.model1.all <- cv.glmnet(train.x, train.y, weights=train.weights, 
                               nfolds = 4, family = 'binomial', 
                               lambda = rev(seq(0,1,0.001)), type.measure = 'mse')
plot(cv.fit.model1.all)
plot(cv.fit.model1.all$lambda, cv.fit.model1.all$cvm)

print(cv.fit.model1.all)

model1.lambda.min <- cv.fit.model1.all$lambda.min

model1.coef <- coef(cv.fit.model1.all, s='lambda.min')

'
Drop Diagnosis Variables
'
# RNG seed
set.seed(617)

# choose 27 random donors to hold out as test set (same 27 as above)
test <- data_imp[sample(nrow(data_imp), 27), ]
train <- data_imp[!rownames(data_imp) %in% rownames(test), ]

# remove dsm_iv_clinical_diagnosis, nincds_arda_diagnosis, nia_reagan, braak, cerad
diagnosis_vars <- c('dsm_iv_clinical_diagnosis',
                    'nincds_arda_diagnosis',
                    'cerad',
                    'braak',
                    'nia_reagan')

train_no_diag <- train[ , !names(train) %in% diagnosis_vars]
test_no_diag <- test[ , !names(test) %in% diagnosis_vars]

# extract observation weights
train.weights.no.diag <- as.numeric(train_no_diag$obs_weight)
train_no_diag$obs_weight <- NULL

test.weights.no.diag <- as.numeric(test_no_diag$obs_weight)
test_no_diag$obs_weight <- NULL

# define model x & y for train & test sets
train.y.no.diag <- as.numeric(train_no_diag$act_demented)
train.x.no.diag <- model.matrix(~.-act_demented, data=train_no_diag)
train.x.no.diag <- train.x.no.diag[ , -1]

test.y.no.diag <- as.numeric(test_no_diag$act_demented)
test.x.no.diag <- model.matrix(~.-act_demented, data=test_no_diag)
test.x.no.diag <- test.x.no.diag[ , -1]

'
Model 2 - LASSO w/o Diagnosis Variables
'
cv.fit.model2.no.diag <- cv.glmnet(train.x.no.diag, train.y.no.diag, weights=train.weights.no.diag, 
                                   nfolds = 4, family = 'binomial', alpha=1,
                                   lambda = rev(seq(0,1,0.001)), type.measure = 'mse')

plot(cv.fit.model2.no.diag)
plot(cv.fit.model2.no.diag$lambda, cv.fit.model2.no.diag$cvm)

model2.lambda.min <- cv.fit.model2.no.diag$lambda.min

model2.coef.min <- coef(cv.fit.model2.no.diag, s='lambda.min')

model2.yhat <- predict(cv.fit.model2.no.diag, s=cv.fit.model2.no.diag$lambda.1se, 
                       weights=test.weights.no.diag, newx=test.x.no.diag, type='response')

'
Model 3 - Weighted Logistic Regression with 14 Nonzero Features from Model 2
'
# subset original dataset, features selected by elasticnet only
nonzero.features <- c('ihc_ptdp_43_ffpe.TCx',
                      'ihc_a_syn.PCx',
                      'ihc_tau2_ffpe.HIP',
                      'ihc_at8_ffpe.PCx',
                      'ihc_a_beta_ffpe.HIP',
                      'ihc_gfap_ffpe.TCx',
                      'ihc_gfap_ffpe.PCx',
                      'apo_e4_allele',
                      'a_syn_pg_per_mg.PCx',
                      'mcp_1_pg_per_mg.TCx',
                      'gene_cluster02.HIP',
                      'ihc_iba1_ffpe.FWM',
                      'ihc_iba1_ffpe.PCx',
                      'obs_weight',
                      'act_demented')

train.nonzero <- train[ , names(train) %in% nonzero.features]
test.nonzero <- test[ , names(test) %in% nonzero.features]

# change act_demented to binary
train.nonzero$act_demented <- as.numeric(train.nonzero$act_demented)
train.nonzero$act_demented[train.nonzero$act_demented == 2] <- 0

test.nonzero$act_demented <- as.numeric(test.nonzero$act_demented)
test.nonzero$act_demented[test.nonzero$act_demented == 2] <- 0

model3.glm <- glm(act_demented ~.-obs_weight, data=train.nonzero, weights=obs_weight, family=binomial(link='logit'))

# make predictions on test set
model3.yhat <- round(predict(model3.glm, newdata=test.nonzero, type='response'), 4)

# make test prediction labels
model3.yhat.labels <- model3.yhat
model3.yhat.labels[model3.yhat.labels > 0.5] = 1
model3.yhat.labels[model3.yhat.labels <= 0.5] = 0

# make train prediction labels
model3.glm.fitted.labels <- model3.glm$fitted
model3.glm.fitted.labels[model3.glm.fitted.labels > 0.5] = 1
model3.glm.fitted.labels[model3.glm.fitted.labels <= 0.5] = 0

# test confusion matrix
confusionMatrix(as.factor(model3.yhat.labels), as.factor(test.nonzero$act_demented))

# train confusion matrix
confusionMatrix(as.factor(model3.glm.fitted.labels), as.factor(train.nonzero$act_demented))

# ROC curves - test
test.nonzero$prob <- model3.yhat
model3.roc.test <- roc(act_demented~prob, data=test.nonzero, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# train curve
train.nonzero$prob <- model3.glm$fitted
model3.roc.train <- roc(act_demented~prob, data=train.nonzero, plot=TRUE,  grid=TRUE, print.auc=TRUE)

'
Model 4 - Removed two lowest performing variables
'
# removed 3 least significant variables
best.features <- c('ihc_ptdp_43_ffpe.TCx',
                    'ihc_a_syn.PCx',
                    'ihc_tau2_ffpe.HIP',
                    'ihc_at8_ffpe.PCx',
                    'ihc_a_beta_ffpe.HIP',
                    'ihc_gfap_ffpe.TCx',
                    #'ihc_gfap_ffpe.PCx',
                    'apo_e4_allele',
                    'a_syn_pg_per_mg.PCx',
                    'mcp_1_pg_per_mg.TCx',
                    'gene_cluster02.HIP',
                    'ihc_iba1_ffpe.FWM',
                    #'ihc_iba1_ffpe.PCx',
                    'obs_weight',
                    'act_demented')

train.best <- train[ , names(train) %in% best.features]
test.best <- test[ , names(test) %in% best.features]

# change act_demented to binary
train.best$act_demented <- as.numeric(train.best$act_demented)
train.best$act_demented[train.best$act_demented == 2] <- 0

test.best$act_demented <- as.numeric(test.best$act_demented)
test.best$act_demented[test.best$act_demented == 2] <- 0

model4.best.glm <- glm(act_demented ~.-obs_weight, data=train.best, weights=obs_weight, family=binomial(link='logit'))

# make predictions on test set
model4.yhat <- round(predict(model4.best.glm, newdata=test.best, type='response'), 4)

# make test prediction labels
model4.yhat.labels <- model4.yhat
model4.yhat.labels[model4.yhat.labels > 0.5] = 1
model4.yhat.labels[model4.yhat.labels <= 0.5] = 0

# make train prediction labels
model4.glm.fitted.labels <- model4.best.glm$fitted
model4.glm.fitted.labels[model4.glm.fitted.labels > 0.5] = 1
model4.glm.fitted.labels[model4.glm.fitted.labels <= 0.5] = 0

# test confusion matrix
confusionMatrix(as.factor(model4.yhat.labels), as.factor(test.best$act_demented))

# train confusion matrix
confusionMatrix(as.factor(model4.glm.fitted.labels), as.factor(train.best$act_demented))

# ROC curves - test
test.best$prob <- model4.yhat
model4.roc.test <- roc(act_demented~prob, data=test.best, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# train curve
train.best$prob <- model4.best.glm$fitted
model4.roc.test <- roc(act_demented~prob, data=train.best, plot=TRUE,  grid=TRUE, print.auc=TRUE)
