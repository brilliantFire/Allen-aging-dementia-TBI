##
'
Lasso Feature Selection ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
05 Oct 2018 - Transfered notebook code to script
10 Oct 2018 - Final models & evaluation (confusion matrices & ROC curves)
19 Oct 2018 - Documentation update
20 Oct 2018 - added low cutoff confusion matrices

This script uses cross-validation to determine optimal values of lambda
(the regularization hyperparameter) that minimize either A) the binomial
deviance or B) the misclassification error in lasso-regularized logistic 
regression models. Four values of lambda are used (e.g. four models are
constructed) using lambdas that minimize the CV error ("lambda.min") or 
the largest value of lambda giving a CV error that is still within 1 
standard error of the lowest value. Models are evaluated using confusion
matrics and ROC curves.


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
test <- data_imp[sample(nrow(data_imp), 27), ]
train <- data_imp[!rownames(data_imp) %in% rownames(test), ]

# extract observation weights
train.weights <- as.numeric(train$obs_weight)
train$obs_weight <- NULL

test.weights <- as.numeric(test$obs_weight)
test$obs_weight <- NULL

# define model x & y for train & test sets
train.y <- as.numeric(train$act_demented)
train.y[train.y == 2] <- 0   # 'No Dementia' = 0
train.x <- model.matrix(~.-act_demented, data=train)
train.x <- train.x[ , -1]

test.y <- as.numeric(test$act_demented)
test.y[test.y == 2] <- 0   # 'No Dementia' = 0
test.x <- model.matrix(~.-act_demented, data=test)
test.x <- test.x[ , -1]

'
Model 1 - Minimum Deviance
'
set.seed(617)

# cross-validated lasso with deviance measure
cv.fit.lasso1.dev <- cv.glmnet(train.x, train.y, weights=train.weights, 
                               nfolds = 4, family = 'binomial', alpha=1,
                               lambda = rev(seq(0.001,1,0.001)), type.measure = 'deviance')

# print
print(cv.fit.lasso1.dev)

# lambda with the minimum cross-validation error
lasso1.lambda.min.dev <- cv.fit.lasso1.dev$lambda.min

# minimum deviance
deviance.min <- cv.fit.lasso1.dev$cvm[which(cv.fit.lasso1.dev$lambda==lasso1.lambda.min.dev)]

# coefficients of the model using minimum CV deviance lambda
lasso1.coefs <- coef(cv.fit.lasso1.dev, s='lambda.min')

# obtain train fitted values
lasso1.train.yhat <- predict(cv.fit.lasso1.dev, newx=train.x, 
                             s='lambda.min', type='response')

# class assignments for train set
lasso1.train.yhat.labels <- predict(cv.fit.lasso1.dev, newx=train.x, 
                                    s='lambda.min', type='class')

# make predictions on test using labmda.min
lasso1.test.yhat <- predict(cv.fit.lasso1.dev, newx=test.x, 
                            s='lambda.min', type='response')

# predict class assignments at p=0.5
lasso1.test.yhat.labels <- predict(cv.fit.lasso1.dev, newx=test.x, 
                                   s='lambda.min', type='class')

# lower cutoff class assignments for train
lasso1.train.yhat.labels.lowCutoff <- lasso1.train.yhat
lasso1.train.yhat.labels.lowCutoff[lasso1.train.yhat.labels.lowCutoff > 0.25] <- 1
lasso1.train.yhat.labels.lowCutoff[lasso1.train.yhat.labels.lowCutoff <= 0.25] <- 0

# lower cutoff class assignments for test
lasso1.test.yhat.labels.lowCutoff <- lasso1.test.yhat
lasso1.test.yhat.labels.lowCutoff[lasso1.test.yhat.labels.lowCutoff > 0.25] <- 1
lasso1.test.yhat.labels.lowCutoff[lasso1.test.yhat.labels.lowCutoff <= 0.25] <- 0

# confusion matrix for train
lasso1.train.conf.mat <- confusionMatrix(as.factor(lasso1.train.yhat.labels), as.factor(train.y),
                                         positive='1')

# confusion matrix for test
lasso1.test.conf.mat <- confusionMatrix(as.factor(lasso1.test.yhat.labels), as.factor(test.y),
                                        positive='1')

# lower cutoff confusion matrix for train
lasso1.train.conf.mat.lowcut <- confusionMatrix(as.factor(lasso1.train.yhat.labels.lowCutoff), as.factor(train.y),
                                                positive='1')

# lower cutoff confusion matrix for test
lasso1.test.conf.mat.lowcut <- confusionMatrix(as.factor(lasso1.test.yhat.labels.lowCutoff), as.factor(test.y),
                                        positive='1')

##### train ROC curve #####
# make a copy of train set
train.copy <- train

# add probabilities
train.copy$prob <- lasso1.train.yhat

# make ROC curve
lasso1.roc.train <- roc(act_demented~as.numeric(prob), data=train.copy, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                        cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                        print.auc.col='black', identity.lwd=2.5)

##### test ROC curve #####
# make a copy of train set
test.copy <- test

# add probabilities
test.copy$prob <- lasso1.test.yhat

# make ROC curve
lasso1.roc.test <- roc(act_demented~as.numeric(prob), data=test.copy, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

'
Model 2 - Minimum Misclassification Error
'
set.seed(617)

# cross-validated lasso with missclassification measure
cv.fit.lasso2.misclass <- cv.glmnet(train.x, train.y, weights=train.weights, 
                                    nfolds = 4, family = 'binomial', alpha=1,
                                    lambda = rev(seq(0.001,1,0.001)), type.measure = 'class')

# print
print(cv.fit.lasso2.misclass)

# lambda with the minimum cross-validation error
lasso2.lambda.min.misclass <- cv.fit.lasso2.misclass$lambda.min

# minimum deviance
misclass.min <- cv.fit.lasso2.misclass$cvm[which(cv.fit.lasso2.misclass$lambda==lasso2.lambda.min.misclass)]

# coefficients of the model using minimum CV misclassification lambda
lasso2.coefs <- coef(cv.fit.lasso2.misclass, s='lambda.min')

# obtain train fitted values
lasso2.train.yhat <- predict(cv.fit.lasso2.misclass, newx=train.x, 
                             s='lambda.min', type='response')

# class assignments for train set
lasso2.train.yhat.labels <- predict(cv.fit.lasso2.misclass, newx=train.x,  
                                    s='lambda.min', type='class')

# make predictions on test using labmda.min
lasso2.test.yhat <- predict(cv.fit.lasso2.misclass, newx=test.x,  
                            s='lambda.min', type='response')

# predict class assignments
lasso2.test.yhat.labels <- predict(cv.fit.lasso2.misclass, newx=test.x,  
                                   s='lambda.min', type='class')

# lower cutoff predictions on train set
lasso2.train.yhat.labels.lowCutoff <- lasso2.train.yhat
lasso2.train.yhat.labels.lowCutoff[lasso2.train.yhat.labels.lowCutoff > 0.25] <- 1
lasso2.train.yhat.labels.lowCutoff[lasso2.train.yhat.labels.lowCutoff <= 0.25] <- 0

# lower cutoff predictions on test set
lasso2.test.yhat.labels.lowCutoff <- lasso2.test.yhat
lasso2.test.yhat.labels.lowCutoff[lasso2.test.yhat.labels.lowCutoff > 0.25] <- 1
lasso2.test.yhat.labels.lowCutoff[lasso2.test.yhat.labels.lowCutoff <= 0.25] <- 0

# confusion matrix for train
lasso2.train.conf.mat <- confusionMatrix(as.factor(lasso2.train.yhat.labels), as.factor(train.y),
                                         positive='1')

# confusion matrix for test
lasso2.test.conf.mat <- confusionMatrix(as.factor(lasso2.test.yhat.labels), as.factor(test.y),
                                        positive='1')

# confusion matrix for low cutoff train
lasso2.train.conf.mat.lowcut <- confusionMatrix(as.factor(lasso2.train.yhat.labels.lowCutoff), as.factor(train.y),
                                                positive='1')

# confusion matrix for low cutoff test
lasso2.test.conf.mat.lowcut <- confusionMatrix(as.factor(lasso2.test.yhat.labels.lowCutoff), as.factor(test.y),
                                               positive='1')

##### train ROC curve #####
# make a copy of train set
train.copy <- train

# add probabilities
train.copy$prob <- lasso2.train.yhat

# make ROC curve
lasso2.roc.train <- roc(act_demented~as.numeric(prob), data=train.copy, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                        cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                        print.auc.col='black', identity.lwd=2.5)

##### test ROC curve #####
# make a copy of train set
test.copy <- test

# add probabilities
test.copy$prob <- lasso2.test.yhat

# make ROC curve
lasso2.roc.test <- roc(act_demented~as.numeric(prob), data=test.copy, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

'
Model 3 - 1se Deviance
'
set.seed(617)

# cross-validated lasso with deviance measure (repeat of model 1 call)
cv.fit.lasso3.dev <- cv.glmnet(train.x, train.y, weights=train.weights, 
                               nfolds = 4, family = 'binomial', alpha=1,
                               lambda = rev(seq(0.001,1,0.001)), type.measure = 'deviance')

# print
print(cv.fit.lasso3.dev)

# lambda with the minimum cross-validation error
lasso3.lambda.1se.dev <- cv.fit.lasso3.dev$lambda.1se

# minimum deviance
deviance.1se <- cv.fit.lasso3.dev$cvm[which(cv.fit.lasso3.dev$lambda==lasso3.lambda.1se.dev)]

# coefficients of the model using 1se CV deviance lambda
lasso3.coefs <- coef(cv.fit.lasso3.dev, s='lambda.1se')

# obtain train fitted values
lasso3.train.yhat <- predict(cv.fit.lasso3.dev, newx=train.x,  
                             s='lambda.1se', type='response')

# class assignments for train set
lasso3.train.yhat.labels <- predict(cv.fit.lasso3.dev, newx=train.x,
                                    s='lambda.1se', type='class')

# make predictions on test using labmda.1se
lasso3.test.yhat <- predict(cv.fit.lasso3.dev, newx=test.x, 
                            s='lambda.1se', type='response')

# predict class assignments
lasso3.test.yhat.labels <- predict(cv.fit.lasso3.dev, newx=test.x,
                                   s='lambda.1se', type='class')

# lower cutoff labels for train
lasso3.train.yhat.labels.lowCutoff <- lasso3.train.yhat
lasso3.train.yhat.labels.lowCutoff[lasso3.train.yhat.labels.lowCutoff > 0.27] <- 1
lasso3.train.yhat.labels.lowCutoff[lasso3.train.yhat.labels.lowCutoff <= 0.27] <- 0

# lower cutoff labels for test
lasso3.test.yhat.labels.lowCutoff <- lasso3.test.yhat
lasso3.test.yhat.labels.lowCutoff[lasso3.test.yhat.labels.lowCutoff > 0.27] <- 1
lasso3.test.yhat.labels.lowCutoff[lasso3.test.yhat.labels.lowCutoff <= 0.27] <- 0

# confusion matrix for train
lasso3.train.conf.mat <- confusionMatrix(as.factor(lasso3.train.yhat.labels), as.factor(train.y),
                                         positive='1')

# confusion matrix for test
lasso3.test.conf.mat <- confusionMatrix(as.factor(lasso3.test.yhat.labels), as.factor(test.y),
                                        positive='1')

# confusion matrix for train, low cutoff
lasso3.train.conf.mat.lowcut <- confusionMatrix(as.factor(lasso3.train.yhat.labels.lowCutoff), as.factor(train.y),
                                                positive='1')

# confusion matrix for test, low cutoff
lasso3.test.conf.mat.lowcut <- confusionMatrix(as.factor(lasso3.test.yhat.labels.lowCutoff), as.factor(test.y),
                                                  positive='1')

##### train ROC curve #####
# make a copy of train set
train.copy <- train

# add probabilities
train.copy$prob <- lasso3.train.yhat

# make ROC curve
lasso3.roc.train <- roc(act_demented~as.numeric(prob), data=train.copy, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                        cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                        print.auc.col='black', identity.lwd=2.5)

##### test ROC curve #####
# make a copy of train set
test.copy <- test

# add probabilities
test.copy$prob <- lasso3.test.yhat

# make ROC curve
lasso3.roc.test <- roc(act_demented~as.numeric(prob), data=test.copy, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

'
Model 4 - 1se Misclassification Error
'
set.seed(617)

# cross-validated lasso with misclassification error measure
cv.fit.lasso4.misclass <- cv.glmnet(train.x, train.y, weights=train.weights, 
                                    nfolds = 4, family = 'binomial', alpha=1,
                                    lambda = rev(seq(0.001,1,0.001)), type.measure = 'class')

# print
print(cv.fit.lasso4.misclass)

# lambda with the minimum cross-validation error
lasso4.lambda.1se.misclass <- cv.fit.lasso4.misclass$lambda.1se

# minimum deviance
misclass.1se <- cv.fit.lasso4.misclass$cvm[which(cv.fit.lasso4.misclass$lambda==lasso4.lambda.1se.misclass)]

# coefficients of the model using 1se CV misclassification lambda
lasso4.coefs <- coef(cv.fit.lasso4.misclass, s='lambda.1se')

# obtain train fitted values
lasso4.train.yhat <- predict(cv.fit.lasso4.misclass, newx=train.x, 
                             s='lambda.1se', type='response')

# class assignments for train set
lasso4.train.yhat.labels <- predict(cv.fit.lasso4.misclass, newx=train.x, 
                                    s='lambda.1se', type='class')

# make predictions on test using labmda.1se
lasso4.test.yhat <- predict(cv.fit.lasso4.misclass, newx=test.x, 
                            s='lambda.1se', type='response')

# predict class assignments
lasso4.test.yhat.labels <- predict(cv.fit.lasso4.misclass, newx=test.x, 
                                   s='lambda.1se', type='class')

# lower cutoff labels for train
lasso4.train.yhat.labels.lowCutoff <- lasso4.train.yhat
lasso4.train.yhat.labels.lowCutoff[lasso4.train.yhat.labels.lowCutoff > 0.27] <- 1
lasso4.train.yhat.labels.lowCutoff[lasso4.train.yhat.labels.lowCutoff <= 0.27] <- 0

# lower cutoff labels for test
lasso4.test.yhat.labels.lowCutoff <- lasso4.test.yhat
lasso4.test.yhat.labels.lowCutoff[lasso4.test.yhat.labels.lowCutoff > 0.27] <- 1
lasso4.test.yhat.labels.lowCutoff[lasso4.test.yhat.labels.lowCutoff <= 0.27] <- 0

# confusion matrix for train
lasso4.train.conf.mat <- confusionMatrix(as.factor(lasso4.train.yhat.labels), as.factor(train.y),
                                         positive='1')

# confusion matrix for test
lasso4.test.conf.mat <- confusionMatrix(as.factor(lasso4.test.yhat.labels), as.factor(test.y),
                                        positive='1')

# confusion matrix for train, low cutoff
lasso4.train.conf.mat.lowcut <- confusionMatrix(as.factor(lasso4.train.yhat.labels.lowCutoff), as.factor(train.y),
                                                positive='1')

# confusion matrix for test, low cutoff
lasso4.test.conf.mat.lowcut <- confusionMatrix(as.factor(lasso4.test.yhat.labels.lowCutoff), as.factor(test.y),
                                               positive='1')


##### train ROC curve #####
# make a copy of train set
train.copy <- train

# add probabilities
train.copy$prob <- lasso4.train.yhat

# make ROC curve
lasso4.roc.train <- roc(act_demented~as.numeric(prob), data=train.copy, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                        cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                        print.auc.col='black', identity.lwd=2.5)

##### test ROC curve #####
# make a copy of train set
test.copy <- test

# add probabilities
test.copy$prob <- lasso4.test.yhat

# make ROC curve
lasso4.roc.test <- roc(act_demented~as.numeric(prob), data=test.copy, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

'
Plots
'
## ~*~*~*~* COEFFICIENT PATHS - ADD 1SE LINES (10.18) *~*~*~*~ ##
# non-cv glmnet call for coefficient paths
lasso.fit <- glmnet(train.x, train.y, weights=train.weights,
                    family = 'binomial', lambda = rev(seq(0.001,1,0.001)), alpha = 1)

# coefficient paths, all variables, lambda mins marked
png('data/lasso_coefficient_paths.png', height=400, width=900)
par(mar=c(4.5,5.5,1.0,5.5))
plot(lasso.fit, xvar='lambda', xlab=expression(log(lambda)),
     ylim=c(-500,500),cex.axis=1.7, cex.lab=2, cex=2, lwd=2)
abline(v=log(cv.fit.lasso1.dev$lambda.min), lty='dashed', lwd=3, col='darkorange2')
abline(v=log(cv.fit.lasso2.misclass$lambda.min), lty='dashed', lwd=3, col='hotpink2')
abline(v=log(cv.fit.lasso3.dev$lambda.1se), lty='dashed', lwd=3, col='blue3')
abline(v=log(cv.fit.lasso4.misclass$lambda.1se), lty='dashed', lwd=3, col='darkcyan')
legend('bottomleft', c(expression(lambda['min, dev'] ~ 'lasso1'),
                       expression(lambda['min, misclass'] ~ 'lasso2'),
                       expression(lambda['1se, dev'] ~ 'lasso3'),
                       expression(lambda['1se, misclass'] ~ 'lasso4')),
       col=c('darkorange2', 'hotpink2', 'blue3', 'darkcyan'),
       lty='dashed', lwd=3, cex=1.3)
dev.off()

## ~*~*~*~* DEVIANCE VERSUS LOG(LAMBDA) *~*~*~*~ ##
png('data/lasso1_deviance_v_logLambda.png', height=500, width=1000)
par(cex=1.7)
plot(cv.fit.lasso1.dev, xlab=expression(log(lambda)), ylab='deviance')
dev.off()

## ~*~*~*~* DEVIANCE VERSUS LOG(LAMBDA) *~*~*~*~ ##
png('data/lasso2_misclass_v_logLambda.png', height=500, width=1000)
par(cex=1.7)
plot(cv.fit.lasso2.misclass, xlab=expression(log(lambda)), ylab='misclassification error')
dev.off()

## ~*~*~*~* DEVIANCE VERSUS LAMBDA - ADD 1SE DASHED LINES (10.18) *~*~*~*~ ##
png('data/lasso1_deviance_v_lambda.png', height=400, width=900)
par(mar=c(4.5,5.5,1.5,1.5))
plot(cv.fit.lasso1.dev$lambda, cv.fit.lasso1.dev$cvm, xlab=expression(lambda), 
     ylab='deviance', pch=20, col='magenta4', xlim=c(0,0.2),
     ylim=range(c(cv.fit.lasso1.dev$cvm-cv.fit.lasso1.dev$cvsd, cv.fit.lasso1.dev$cvm+cv.fit.lasso1.dev$cvsd)),
     cex.axis=2, cex.lab=2.5, cex=2)
# error bars
arrows(cv.fit.lasso1.dev$lambda, cv.fit.lasso1.dev$cvm-cv.fit.lasso1.dev$cvsd, 
       cv.fit.lasso1.dev$lambda, cv.fit.lasso1.dev$cvm+cv.fit.lasso1.dev$cvsd, 
       length=0.05, angle=90, code=3, col='gray50')
abline(v=cv.fit.lasso1.dev$lambda.min, lty='dashed', lwd=3, col='darkorange2')
abline(v=cv.fit.lasso1.dev$lambda.1se, lty='dashed', lwd=3, col='blue3')
text(0.097, 3.0, labels=expression(lambda[min] ~ '= 0.113'), cex=1.7)
text(0.172, 3.0, labels=expression(lambda['1se'] ~ '= 0.189'), cex=1.7)
dev.off()

## ~*~*~*~* MISCLASS VERSUS LAMBDA - ADD 1SE DASHED LINES (10.18) *~*~*~*~ ##
png('data/lasso2_misclass_v_lambda.png', height=400, width=900)
par(mar=c(4.5,5.5,1.5,1.5))
plot(cv.fit.lasso2.misclass$lambda, cv.fit.lasso2.misclass$cvm, xlab=expression(lambda), 
     ylab='misclassification error', pch=20, col='magenta4', xlim=c(-0.01,0.2),
     ylim=range(c(cv.fit.lasso2.misclass$cvm-cv.fit.lasso2.misclass$cvsd, cv.fit.lasso2.misclass$cvm+cv.fit.lasso2.misclass$cvsd)),
     cex.axis=2, cex.lab=2.5, cex=2)
# error bars
arrows(cv.fit.lasso2.misclass$lambda, cv.fit.lasso2.misclass$cvm-cv.fit.lasso2.misclass$cvsd, 
       cv.fit.lasso2.misclass$lambda, cv.fit.lasso2.misclass$cvm+cv.fit.lasso2.misclass$cvsd, 
       length=0.05, angle=90, code=3, col='gray50', xlim=c(-0.001, 0.20))
abline(v=cv.fit.lasso2.misclass$lambda.min, lty='dashed', lwd=3, col='hotpink2')
abline(v=cv.fit.lasso2.misclass$lambda.1se, lty='dashed', lwd=3, col='darkcyan')
text(0.0055, 0.295, labels=expression(lambda[min] ~ '= 0.021'), cex=1.7)
text(0.065, 0.150, labels=expression(lambda['1se'] ~ '= 0.047'), cex=1.7)
dev.off()

'
Pretty tables
'
# confusion matrices for Model 1
lasso1.train.conf.mat

train.1.df <- data.frame(c(43,0),
                         c(32,5))
colnames(train.1.df) <- c('No Dementia', 'Dementia')
rownames(train.1.df) <- c('No Dementia', 'Dementia')

table07 <- xtable(train.1.df)
align(table07) <- 'r|c|c'
digits(table07) <- 0
print(table07, hline.after = c(0,1))

print.xtable(table07, type='latex', file='data/table07.tex')

lasso1.test.conf.mat
test.1.df <- data.frame(c(14,0),
                        c(11,2))
colnames(test.1.df) <- c('No Dementia', 'Dementia')
rownames(test.1.df) <- c('No Dementia', 'Dementia')

table08 <- xtable(test.1.df)
align(table08) <- 'r|c|c'
digits(table08) <- 0
print(table08, hline.after = c(0,1))

print.xtable(table08, type='latex', file='data/table08.tex')

### low cutoffs
lasso1.train.conf.mat.lowcut

train.1.df <- data.frame(c(30,13),
                         c(11,26))
colnames(train.1.df) <- c('No Dementia', 'Dementia')
rownames(train.1.df) <- c('No Dementia', 'Dementia')

table07.1 <- xtable(train.1.df)
align(table07.1) <- 'r|c|c'
digits(table07.1) <- 0
print(table07.1, hline.after = c(0,1))

print.xtable(table07.1, type='latex', file='data/table07.1.tex')

lasso1.test.conf.mat.lowcut
test.1.df <- data.frame(c(9,5),
                        c(3,10))
colnames(test.1.df) <- c('No Dementia', 'Dementia')
rownames(test.1.df) <- c('No Dementia', 'Dementia')

table08.1 <- xtable(test.1.df)
align(table08.1) <- 'r|c|c'
digits(table08.1) <- 0
print(table08.1, hline.after = c(0,1))

print.xtable(table08.1, type='latex', file='data/table08.1.tex')

# confusion matrices for Model 2
lasso2.train.conf.mat
train.2.df <- data.frame(c(42,1),
                         c(12,25))
colnames(train.2.df) <- c('No Dementia', 'Dementia')
rownames(train.2.df) <- c('No Dementia', 'Dementia')

table09 <- xtable(train.2.df)
align(table09) <- 'r|c|c'
digits(table09) <- 0
print(table09, hline.after = c(0,1))

print.xtable(table09, type='latex', file='data/table09.tex')

lasso2.test.conf.mat
test.2.df <- data.frame(c(12,2),
                        c(9,4))
colnames(test.2.df) <- c('No Dementia', 'Dementia')
rownames(test.2.df) <- c('No Dementia', 'Dementia')

table10 <- xtable(test.2.df)
align(table10) <- 'r|c|c'
digits(table10) <- 0
print(table10, hline.after = c(0,1))

print.xtable(table10, type='latex', file='data/table10.tex')

### low cutoffs
lasso2.train.conf.mat.lowcut
train.2.df <- data.frame(c(35,8),
                         c(8,29))
colnames(train.2.df) <- c('No Dementia', 'Dementia')
rownames(train.2.df) <- c('No Dementia', 'Dementia')

table09.1 <- xtable(train.2.df)
align(table09.1) <- 'r|c|c'
digits(table09.1) <- 0
print(table09.1, hline.after = c(0,1))

print.xtable(table09.1, type='latex', file='data/table09.1.tex')

lasso2.test.conf.mat.lowcut
test.2.df <- data.frame(c(7,7),
                        c(6,7))
colnames(test.2.df) <- c('No Dementia', 'Dementia')
rownames(test.2.df) <- c('No Dementia', 'Dementia')

table10.1 <- xtable(test.2.df)
align(table10.1) <- 'r|c|c'
digits(table10.1) <- 0
print(table10.1, hline.after = c(0,1))

print.xtable(table10.1, type='latex', file='data/table10.1.tex')


# confusion matrices for Model 3
lasso3.train.conf.mat
train.3.df <- data.frame(c(43,0),
                         c(37,0))
colnames(train.3.df) <- c('No Dementia', 'Dementia')
rownames(train.3.df) <- c('No Dementia', 'Dementia')

table11 <- xtable(train.3.df)
align(table11) <- 'r|c|c'
digits(table11) <- 0
print(table11, hline.after = c(0,1))

print.xtable(table11, type='latex', file='data/table11.tex')

lasso3.test.conf.mat
test.3.df <- data.frame(c(14,0),
                        c(13,0))
colnames(test.3.df) <- c('No Dementia', 'Dementia')
rownames(test.3.df) <- c('No Dementia', 'Dementia')

table12 <- xtable(test.3.df)
align(table12) <- 'r|c|c'
digits(table12) <- 0
print(table12, hline.after = c(0,1))

print.xtable(table12, type='latex', file='data/table12.tex')

### low cutoff confusion matrices, model 3
lasso3.train.conf.mat.lowcut
train.3.df <- data.frame(c(37,6),
                         c(14,23))
colnames(train.3.df) <- c('No Dementia', 'Dementia')
rownames(train.3.df) <- c('No Dementia', 'Dementia')

table11.1 <- xtable(train.3.df)
align(table11.1) <- 'r|c|c'
digits(table11.1) <- 0
print(table11.1, hline.after = c(0,1))

print.xtable(table11.1, type='latex', file='data/table11.1.tex')

lasso3.test.conf.mat.lowcut
test.3.df <- data.frame(c(12,2),
                        c(3,10))
colnames(test.3.df) <- c('No Dementia', 'Dementia')
rownames(test.3.df) <- c('No Dementia', 'Dementia')

table12.1 <- xtable(test.3.df)
align(table12.1) <- 'r|c|c'
digits(table12.1) <- 0
print(table12.1, hline.after = c(0,1))

print.xtable(table12.1, type='latex', file='data/table12.1.tex')

# confusion matrices for Model 4
lasso4.train.conf.mat
train.4.df <- data.frame(c(43,0),
                         c(19,18))
colnames(train.4.df) <- c('No Dementia', 'Dementia')
rownames(train.4.df) <- c('No Dementia', 'Dementia')

table13 <- xtable(train.4.df)
align(table13) <- 'r|c|c'
digits(table13) <- 0
print(table13, hline.after = c(0,1))

print.xtable(table13, type='latex', file='data/table13.tex')

lasso4.test.conf.mat
test.4.df <- data.frame(c(13,1),
                        c(9,4))
colnames(test.4.df) <- c('No Dementia', 'Dementia')
rownames(test.4.df) <- c('No Dementia', 'Dementia')

table14 <- xtable(test.4.df)
align(table14) <- 'r|c|c'
digits(table14) <- 0
print(table14, hline.after = c(0,1))

print.xtable(table14, type='latex', file='data/table14.tex')

### model 4, low cutoff confusion matrices
lasso4.train.conf.mat.lowcut
train.4.df <- data.frame(c(34,9),
                         c(9,27))
colnames(train.4.df) <- c('No Dementia', 'Dementia')
rownames(train.4.df) <- c('No Dementia', 'Dementia')

table13.1 <- xtable(train.4.df)
align(table13.1) <- 'r|c|c'
digits(table13.1) <- 0
print(table13.1, hline.after = c(0,1))

print.xtable(table13.1, type='latex', file='data/table13.1.tex')

lasso4.test.conf.mat.lowcut
test.4.df <- data.frame(c(9,5),
                        c(7,6))
colnames(test.4.df) <- c('No Dementia', 'Dementia')
rownames(test.4.df) <- c('No Dementia', 'Dementia')

table14.1 <- xtable(test.4.df)
align(table14.1) <- 'r|c|c'
digits(table14.1) <- 0
print(table14.1, hline.after = c(0,1))

print.xtable(table14.1, type='latex', file='data/table14.1.tex')
