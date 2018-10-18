##
'
Lasso Feature Selection ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
05 Oct 2018 - Transfered notebook code to script
10 Oct 2018 - Final models & evaluation (confusion matrices & ROC curves)


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
library(caret)     # confusion matrices
library(pROC)      # easy ROC curves
library(car)       # VIFs
library(xtable)    # pretty LaTeX tables

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
train.2 <- data_imp[!rownames(data_imp) %in% rownames(test.2), ]

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
train.y.2[train.y.2 == 2] <- 0   # 'No Dementia' = 0
train.x.2 <- model.matrix(~.-act_demented, data=train.2)
train.x.2 <- train.x.2[ , -1]

test.y.2 <- as.numeric(test.2$act_demented)
test.y.2[test.y.2 == 2] <- 0   # 'No Dementia' = 0
test.x.2 <- model.matrix(~.-act_demented, data=test.2)
test.x.2 <- test.x.2[ , -1]

'
Model 1 - Minimum Deviance
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

# deviance versus lambda (not log)
png('data/lasso2_dev_v_lambda.png', height=400, width=900)
par(mar=c(4.5,5.5,1.5,1.5))
plot(cv.fit.lasso2.dev$lambda, cv.fit.lasso2.dev$cvm, xlab=expression(lambda), 
     ylab='deviance', pch=20, col='magenta4', xlim=c(0,0.2),
     ylim=range(c(cv.fit.lasso2.dev$cvm-cv.fit.lasso2.dev$cvsd, cv.fit.lasso2.dev$cvm+cv.fit.lasso2.dev$cvsd)),
     cex.axis=2, cex.lab=2.5, cex=2)
arrows(cv.fit.lasso2.dev$lambda, cv.fit.lasso2.dev$cvm-cv.fit.lasso2.dev$cvsd, 
       cv.fit.lasso2.dev$lambda, cv.fit.lasso2.dev$cvm+cv.fit.lasso2.dev$cvsd, 
       length=0.05, angle=90, code=3, col='gray50')
abline(v=cv.fit.lasso2.dev$lambda.min, lty='dashed', lwd=3, col='darkorange2')
text(0.063, 2, labels=expression(lambda[min] ~ '= 0.080'), cex=1.7)
dev.off()

# print
print(cv.fit.lasso2.dev)

# lambda with the minimum cross-validation error
lasso2.lambda.min.dev <- cv.fit.lasso2.dev$lambda.min

# minimum deviance
min.deviance <- cv.fit.lasso2.dev$cvm[which(cv.fit.lasso2.dev$lambda==lasso2.lambda.min.dev)]

# coefficients of the model using minimum CV deviance lambda
dev.min.coefs <- coef(cv.fit.lasso2.dev, s='lambda.min')

# make predicitions using labmda.min
lasso2.dev.yhat <- predict(cv.fit.lasso2.dev, newx=test.x.2, s='lambda.min')

'
Model 2 - Minimum Misclassification Error
'
# same CV model except with misclassification as measure
cv.fit.lasso2.misclass <- cv.glmnet(train.x.2, train.y.2, weights=train.weights.2, 
                                    nfolds = 4, family = 'binomial', alpha=1,
                                    lambda = rev(seq(0,1,0.001)), type.measure = 'class')

# deviance versus log(lambda)
png('data/lasso2_misclass_v_logLambda.png', height=500, width=1000)
par(cex=1.7)
plot(cv.fit.lasso2.misclass, xlab=expression(log(lambda)), ylab='misclassification error')
dev.off()

# deviance versus lambda (not log)
png('data/lasso2_misclass_v_lambda.png', height=400, width=900)
par(mar=c(4.5,5.5,1.5,1.5))
plot(cv.fit.lasso2.misclass$lambda, cv.fit.lasso2.misclass$cvm, xlab=expression(lambda), 
     ylab='misclassification error', pch=20, col='magenta4', xlim=c(0,0.2),
     ylim=range(c(cv.fit.lasso2.misclass$cvm-cv.fit.lasso2.misclass$cvsd, 
                  cv.fit.lasso2.misclass$cvm+cv.fit.lasso2.misclass$cvsd)),
     cex.axis=2, cex.lab=2.3, cex=2)
arrows(cv.fit.lasso2.misclass$lambda, cv.fit.lasso2.misclass$cvm-cv.fit.lasso2.misclass$cvsd, 
       cv.fit.lasso2.misclass$lambda, cv.fit.lasso2.misclass$cvm+cv.fit.lasso2.misclass$cvsd, 
       length=0.05, angle=90, code=3, col='gray50')
abline(v=cv.fit.lasso2.misclass$lambda.min, lty='dashed', lwd=3, col='darkorange2')
text(0.056, 0.4, labels=expression(lambda[min] ~ '= 0.040'), cex=1.7)
dev.off()

# lambda with the minimum cross-validation error
lasso2.lambda.min.misclass <- cv.fit.lasso2.misclass$lambda.min

# minimum deviance
min.misclass <- cv.fit.lasso2.misclass$cvm[which(cv.fit.lasso2.misclass$lambda==lasso2.lambda.min.misclass)]

# print
print(cv.fit.lasso2.misclass)

# lowest CV misclassification error coefs
misclass.min.coefs <- coef(cv.fit.lasso2.misclass, s='lambda.min')

# non-cv glmnet call for coefficient paths
lasso2.fit <- glmnet(train.x.2, train.y.2, weights=train.weights.2,
                     family = 'binomial', lambda = rev(seq(0,1,0.001)), alpha = 1)

# coefficient paths for lasso2, all variables
png('data/lasso2_coefficient_paths.png', height=400, width=900)
par(mar=c(4.5,5.5,1.0,1.5))
plot(lasso2.fit, xvar='lambda', xlab=expression(log(lambda)),
     ylim=c(-500,500),cex.axis=1.7, cex.lab=2, cex=2, lwd=2)
abline(v=log(cv.fit.lasso2.dev$lambda.min), lty='dashed', lwd=3, col='darkorange2')
abline(v=log(cv.fit.lasso2.misclass$lambda.min), lty='dashed', lwd=3, col='forestgreen')
text(-4.1, -400, labels=expression('misclassification' ~ lambda[min]), cex=1.7)
text(-1.9, -400, labels=expression('deviance' ~ lambda[min]), cex=1.7)
dev.off()

'
Predictions and Evaluation for the Two Models
'
# make predictions on models with lambda = 0.04 (misclass model) & lambda = 0.08 (deviance)
lasso2.yhats <- predict(lasso2.fit, newx=test.x.2, s=c(0.04, 0.08), type='response')

# predictions
lasso2.dev.yhat <- lasso2.yhats[ , 2]
lasso2.misclass.yhat <- lasso2.yhats[ , 1]

# labels for predictions - deviance model
lasso2.dev.yhat.labels <- lasso2.dev.yhat
lasso2.dev.yhat.labels[lasso2.dev.yhat.labels > 0.5000] = 1
lasso2.dev.yhat.labels[lasso2.dev.yhat.labels <= 0.5000] = 0

# labels for predictions - misclass model
lasso2.misclass.yhat.labels <- lasso2.misclass.yhat
lasso2.misclass.yhat.labels[lasso2.misclass.yhat.labels > 0.5000] = 1
lasso2.misclass.yhat.labels[lasso2.misclass.yhat.labels <= 0.5000] = 0

# confusion matrices
dev.test.conf.mat <- confusionMatrix(as.factor(lasso2.dev.yhat.labels), as.factor(test.y.2),
                                   positive='1')

misclass.test.conf.mat <- confusionMatrix(as.factor(lasso2.misclass.yhat.labels), as.factor(test.y.2),
                                          positive='1')

# roc
lasso2.dev.roc.test <- roc(lasso2.dev.yhat.labels, lasso2.dev.yhat, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       #ex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)