##
'
Lasso Feature Selection ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
05 Oct 2018 - Transfered notebook code to script
10 Oct 2018 - Final models & evaluation (confusion matrices & ROC curves)
16 Oct 2018 - Split off logistic models into separate script

This script performs feature selection using Lasso.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.5.1
'
##

setwd('..')

'
Load libraries & dataset
'
library(glmnet)    # for LASSO

data_imp <- read.csv('data/final_dataset_cart_imputed.csv')

# structure list
str(data_imp, list.len=ncol(data_imp))

# remove X
data_imp$X <- NULL

'
Drop Diagnosis Variables & Make Matrices
'
set.seed(617)

# choose 27 random donors to hold out as test set (same 27 as above)
test.2 <- data_imp[sample(nrow(data_imp), 27), ]
train.2 <- data_imp[!rownames(data_imp) %in% rownames(test), ]

# remove dsm_iv_clinical_diagnosis, nincds_arda_diagnosis, nia_reagan, braak, cerad
diagnosis_vars <- c('dsm_iv_clinical_diagnosis',
                    'nincds_arda_diagnosis',
                    'cerad',
                    'braak',
                    'nia_reagan')

train.2 <- train.2[ , !names(train.2) %in% diagnosis_vars]
test.2 <- test.2[ , !names(test.2) %in% diagnosis_vars]

# extract observation weights
train.weights.2 <- as.numeric(train.2$obs_weight)
train.2$obs_weight <- NULL

test.weights.2 <- as.numeric(test.2$obs_weight)
test.2$obs_weight <- NULL

# define model x & y for train & test sets
train.y.2 <- as.numeric(train.2$act_demented)
train.x.2 <- model.matrix(~.-act_demented, data=train.2)
train.x.2 <- train.x.2[ , -1]

test.y.2 <- as.numeric(test.2$act_demented)
test.x.2 <- model.matrix(~.-act_demented, data=test.2)
test.x.2 <- test.x.2[ , -1]

'
LASSO w/o Diagnosis Variables
'
# cross-validated lasso with deviance measure
cv.fit.lasso2.dev <- cv.glmnet(train.x.2, train.y.2, weights=train.weights.2, 
                           nfolds = 4, family = 'binomial', alpha=1,
                           lambda = rev(seq(0,1,0.001)), type.measure = 'deviance')

# deviance versus log(lambda)
png('data/lasso2_dev_v_logLambda.png', height=500, width=1000)
par(cex=1.7)
plot(cv.fit.lasso2.dev, xlab=expression(log(lambda)), ylab='deviance',
     ylim=c(0,6))
dev.off()

# deviance versus lambda
plot(cv.fit.lasso2.dev$lambda, cv.fit.lasso2.dev$cvm, xlab=expression(lambda), ylab='deviance')

# print
print(cv.fit.lasso2.dev)

# lambda with the minimum cross-validation error
lasso2.lambda.min <- cv.fit.lasso2.dev$lambda.min

# coefficients of the model using minimum CV deviance lambda
coef(cv.fit.lasso2.dev, s='lambda.min')

# same CV model except with misclassification as measure
cv.fit.lasso2.misclass <- cv.glmnet(train.x.2, train.y.2, weights=train.weights.2, 
                                    nfolds = 4, family = 'binomial', alpha=1,
                                    lambda = rev(seq(0,1,0.001)), type.measure = 'class')

# deviance versus log(lambda)
png('data/lasso2_misclass_v_logLambda.png', height=500, width=1000)
par(cex=1.7)
plot(cv.fit.lasso2.misclass, xlab=expression(log(lambda)), ylab='misclassification error')
dev.off()

# print
print(cv.fit.lasso2.misclass)

# lowest CV misclassification error coefs
coef(cv.fit.lasso2.misclass, s='lambda.min')

# non-cv glmnet call for coefficient paths
lasso2.fit <- glmnet(train.x.2, train.y.2, weights=train.weights.2,
                     family = 'binomial', lambda = rev(seq(0,1,0.001)), alpha = 1)

# coefficient paths for lasso2, all variables
plot(lasso2.fit, xvar='lambda', xlab=expression(log(lambda)),
     ylim=c(-1000,1000))
