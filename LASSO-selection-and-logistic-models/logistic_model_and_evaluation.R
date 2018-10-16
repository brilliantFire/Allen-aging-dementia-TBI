##
'
Logistic Models & Evaluation ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
16 Oct 2018 - Created fresh script to house logistic models & evaluations

This script constructs and evaluates five logistic regression models.

See http://aging.brain-map.org/overview/home for more details about the data.

R 3.5.1
'
##

setwd('..')

'
Load libraries & dataset
'
library(caret)     # confusion matrices
library(pROC)      # easy ROC curves
library(car)       # VIFs

data_imp <- read.csv('data/final_dataset_cart_imputed.csv')

# structure list
str(data_imp, list.len=ncol(data_imp))

# remove X
data_imp$X <- NULL

# RNG seed
set.seed(617)

# choose 27 random donors to hold out as test set
test <- data_imp[sample(nrow(data_imp), 27), ]
train <- data_imp[!rownames(data_imp) %in% rownames(test), ]

'
Model 1 - Weighted Logistic Regression with 10 Nonzero Features from Lasso 2 (lowest CV deviance)
'
features.1 <- c('ihc_ptdp_43_ffpe.TCx',
                'ihc_a_syn.PCx',
                'ihc_tau2_ffpe.HIP',
                'ihc_at8_ffpe.PCx',
                'ihc_a_syn.HIP',
                'ihc_gfap_ffpe.PCx',
                'a_syn_pg_per_mg.PCx',
                'apo_e4_allele',
                'mcp_1_pg_per_mg.TCx',
                'ihc_iba1_ffpe.PCx',
                'obs_weight',
                'act_demented')

train.1 <- train[ , names(train) %in% features.1]
test.1 <- test[ , names(test) %in% features.1]

# change act_demented to binary
train.1$act_demented <- as.numeric(train.1$act_demented)
train.1$act_demented[train.1$act_demented == 2] <- 0

test.1$act_demented <- as.numeric(test.1$act_demented)
test.1$act_demented[test.1$act_demented == 2] <- 0

# model 1 fit
model1.glm <- glm(act_demented ~.-obs_weight, data=train.1, weights=obs_weight, family=binomial(link='logit'))

# make predictions on test set
model1.yhat <- round(predict(model1.glm, newdata=test.1, type='response'), 4)

# make test prediction labels
model1.yhat.labels <- model1.yhat
model1.yhat.labels[model1.yhat.labels > 0.5] = 1
model1.yhat.labels[model1.yhat.labels <= 0.5] = 0

# make train prediction labels
model1.glm.fitted.labels <- model1.glm$fitted
model1.glm.fitted.labels[model1.glm.fitted.labels > 0.5] = 1
model1.glm.fitted.labels[model1.glm.fitted.labels <= 0.5] = 0

# summary
summary(model1.glm)

# train confusion matrix
confusionMatrix(as.factor(model1.glm.fitted.labels), as.factor(train.1$act_demented))

# test confusion matrix
confusionMatrix(as.factor(model1.yhat.labels), as.factor(test.1$act_demented))

# train ROC curve
train.1$prob <- model1.glm$fitted
model1.roc.train <- roc(act_demented~prob, data=train.1, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# test ROC curve
test.1$prob <- model1.yhat
model1.roc.test <- roc(act_demented~prob, data=test.1, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# VIFs
model1.vif <- vif(model1.glm)

'
Model 2 - Weighted Logistic Regression with 14 Nonzero Features from Lasso 2 (lowest CV missclassification)
'
features.2 <- c('ihc_ptdp_43_ffpe.TCx',
                'ihc_a_syn.PCx',
                'ihc_tau2_ffpe.HIP',
                'ihc_at8_ffpe.PCx',
                'ihc_gfap_ffpe.PCx',
                'ihc_a_beta_ffpe.HIP',
                'a_syn_pg_per_mg.PCx',
                'ihc_tau2_ffpe.TCx',
                'apo_e4_allele',
                'ever_tbi_w_loc',
                'mcp_1_pg_per_mg.TCx',
                'mcp_1_pg_per_mg.FWM',
                'ab42_pg_per_mg.FWM',
                'gene_cluster02.HIP',
                'ihc_iba1_ffpe.PCx',
                'ihc_iba1_ffpe.FWM',
                'obs_weight',
                'act_demented')

train.2 <- train[ , names(train) %in% features.2]
test.2 <- test[ , names(test) %in% features.2]

# change act_demented to binary
train.2$act_demented <- as.numeric(train.2$act_demented)
train.2$act_demented[train.2$act_demented == 2] <- 0

test.2$act_demented <- as.numeric(test.2$act_demented)
test.2$act_demented[test.2$act_demented == 2] <- 0

model2.glm <- glm(act_demented ~.-obs_weight, data=train.2, weights=obs_weight, family=binomial(link='logit'))

# summary
summary(model2.glm)

# make predictions on test set
model2.yhat <- round(predict(model2.glm, newdata=test.2, type='response'), 4)

# make test prediction labels
model2.yhat.labels <- model2.yhat
model2.yhat.labels[model2.yhat.labels > 0.5] = 1
model2.yhat.labels[model2.yhat.labels <= 0.5] = 0

# make train prediction labels
model2.glm.fitted.labels <- model2.glm$fitted
model2.glm.fitted.labels[model2.glm.fitted.labels > 0.5] = 1
model2.glm.fitted.labels[model2.glm.fitted.labels <= 0.5] = 0

# train confusion matrix
confusionMatrix(as.factor(model2.glm.fitted.labels), as.factor(train.2$act_demented))

# test confusion matrix
confusionMatrix(as.factor(model2.yhat.labels), as.factor(test.2$act_demented))

# train ROC curve
train.2$prob <- model2.glm$fitted
model2.roc.train <- roc(act_demented~prob, data=train.2, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# test ROC curve
test.2$prob <- model2.yhat
model2.roc.test <- roc(act_demented~prob, data=test.2, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# VIFs
model2.vif <- vif(model2.glm)

'
Model 3 - Remove ihc_a_syn.HIP from Model 3
'
features.3 <- c('ihc_ptdp_43_ffpe.TCx',
                'ihc_a_syn.PCx',
                'ihc_tau2_ffpe.HIP',
                'ihc_at8_ffpe.PCx',
                #'ihc_a_syn.HIP',
                'ihc_gfap_ffpe.PCx',
                'a_syn_pg_per_mg.PCx',
                'apo_e4_allele',
                'mcp_1_pg_per_mg.TCx',
                'ihc_iba1_ffpe.PCx',
                'obs_weight',
                'act_demented')

train.3 <- train[ , names(train) %in% features.3]
test.3 <- test[ , names(test) %in% features.3]

# change act_demented to binary
train.3$act_demented <- as.numeric(train.3$act_demented)
train.3$act_demented[train.3$act_demented == 2] <- 0

test.3$act_demented <- as.numeric(test.3$act_demented)
test.3$act_demented[test.3$act_demented == 2] <- 0

model3.glm <- glm(act_demented ~.-obs_weight, data=train.3, weights=obs_weight, family=binomial(link='logit'))

# summary
summary(model3.glm)

# make predictions on test set
model3.yhat <- round(predict(model3.glm, newdata=test.3, type='response'), 4)

# make test prediction labels
model3.yhat.labels <- model3.yhat
model3.yhat.labels[model3.yhat.labels > 0.5] = 1
model3.yhat.labels[model3.yhat.labels <= 0.5] = 0

# make train prediction labels
model3.glm.fitted.labels <- model3.glm$fitted
model3.glm.fitted.labels[model3.glm.fitted.labels > 0.5] = 1
model3.glm.fitted.labels[model3.glm.fitted.labels <= 0.5] = 0

# train confusion matrix
confusionMatrix(as.factor(model3.glm.fitted.labels), as.factor(train.3$act_demented))

# test confusion matrix
confusionMatrix(as.factor(model3.yhat.labels), as.factor(test.3$act_demented))

# train ROC curve
train.3$prob <- model3.glm$fitted
model3.roc.train <- roc(act_demented~prob, data=train.3, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# test ROC curve
test.3$prob <- model3.yhat
model3.roc.test <- roc(act_demented~prob, data=test.3, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# VIFs
model3.vif <- vif(model3.glm)

'
Model 4 - Removed Underperforming Features from Model 4
'
features.4 <- c('ihc_ptdp_43_ffpe.TCx',
                'ihc_a_syn.PCx',
                'ihc_tau2_ffpe.HIP',
                #'ihc_at8_ffpe.PCx',
                'ihc_gfap_ffpe.PCx',
                'ihc_a_beta_ffpe.HIP',
                'a_syn_pg_per_mg.PCx',
                #'ihc_tau2_ffpe.TCx',
                'apo_e4_allele',
                'ever_tbi_w_loc',
                'mcp_1_pg_per_mg.TCx',
                #'mcp_1_pg_per_mg.FWM',
                #'ab42_pg_per_mg.FWM',
                'gene_cluster02.HIP',
                #'ihc_iba1_ffpe.PCx',
                #'ihc_iba1_ffpe.FWM',
                'obs_weight',
                'act_demented')

train.4 <- train[ , names(train) %in% features.4]
test.4 <- test[ , names(test) %in% features.4]

# change act_demented to binary
train.4$act_demented <- as.numeric(train.4$act_demented)
train.4$act_demented[train.4$act_demented == 2] <- 0

test.4$act_demented <- as.numeric(test.4$act_demented)
test.4$act_demented[test.4$act_demented == 2] <- 0

model4.glm <- glm(act_demented ~.-obs_weight, data=train.4, weights=obs_weight, family=binomial(link='logit'))

# summary
summary(model4.glm)

# make predictions on test set
model4.yhat <- round(predict(model4.glm, newdata=test.4, type='response'), 4)

# make test prediction labels
model4.yhat.labels <- model4.yhat
model4.yhat.labels[model4.yhat.labels > 0.5] = 1
model4.yhat.labels[model4.yhat.labels <= 0.5] = 0

# make train prediction labels
model4.glm.fitted.labels <- model4.glm$fitted
model4.glm.fitted.labels[model4.glm.fitted.labels > 0.5] = 1
model4.glm.fitted.labels[model4.glm.fitted.labels <= 0.5] = 0

# train confusion matrix
confusionMatrix(as.factor(model4.glm.fitted.labels), as.factor(train.4$act_demented))

# test confusion matrix
confusionMatrix(as.factor(model4.yhat.labels), as.factor(test.4$act_demented))

# train ROC curve
train.4$prob <- model4.glm$fitted
model4.roc.train <- roc(act_demented~prob, data=train.4, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# test ROC curve
test.4$prob <- model4.yhat
model4.roc.test <- roc(act_demented~prob, data=test.4, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# VIFs
model4.vif <- vif(model4.glm)

'
Model 5 - Model 3 + HIP cluster 2 medoid, minus ihc_a_syn.HIP
'
features.5 <- c('ihc_ptdp_43_ffpe.TCx',
                'ihc_a_syn.PCx',
                'ihc_tau2_ffpe.HIP',
                'ihc_at8_ffpe.PCx',
                #'ihc_a_syn.HIP',
                'ihc_gfap_ffpe.PCx',
                'a_syn_pg_per_mg.PCx',
                'apo_e4_allele',
                'mcp_1_pg_per_mg.TCx',
                'ihc_iba1_ffpe.PCx',
                'gene_cluster02.HIP',
                'obs_weight',
                'act_demented')

train.5 <- train[ , names(train) %in% features.5]
test.5 <- test[ , names(test) %in% features.5]

# change act_demented to binary
train.5$act_demented <- as.numeric(train.5$act_demented)
train.5$act_demented[train.5$act_demented == 2] <- 0

test.5$act_demented <- as.numeric(test.5$act_demented)
test.5$act_demented[test.5$act_demented == 2] <- 0

# model 5 fit
model5.glm <- glm(act_demented ~.-obs_weight, data=train.5, weights=obs_weight, family=binomial(link='logit'))

# make predictions on test set
model5.yhat <- round(predict(model5.glm, newdata=test.5, type='response'), 4)

# make test prediction labels
model5.yhat.labels <- model5.yhat
model5.yhat.labels[model5.yhat.labels > 0.5] = 1
model5.yhat.labels[model5.yhat.labels <= 0.5] = 0

# make train prediction labels
model5.glm.fitted.labels <- model5.glm$fitted
model5.glm.fitted.labels[model5.glm.fitted.labels > 0.5] = 1
model5.glm.fitted.labels[model5.glm.fitted.labels <= 0.5] = 0

# summary
summary(model5.glm)

# train confusion matrix
confusionMatrix(as.factor(model5.glm.fitted.labels), as.factor(train.5$act_demented))

# test confusion matrix
confusionMatrix(as.factor(model5.yhat.labels), as.factor(test.5$act_demented))

# train ROC curve
train.5$prob <- model5.glm$fitted
model5.roc.train <- roc(act_demented~prob, data=train.5, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# test ROC curve
test.5$prob <- model5.yhat
model5.roc.test <- roc(act_demented~prob, data=test.5, plot=TRUE,  grid=TRUE, print.auc=TRUE)

# VIFs
model5.vif <- vif(model5.glm)
