##
'
Logistic Models & Evaluation ~ Allen Aging, Dementia, & TBI Data
Rebecca Vislay Wade
16 Oct 2018 - Created fresh script to house logistic models & evaluations
17 Oct 2018 - Pretty LaTeX confusion matrices

This script constructs and evaluates six logistic regression models.

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
library(xtable)    # pretty LaTeX tables

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
Model 1 - Weighted Logistic Regression with 10 Nonzero Features from lasso w/ lowest CV deviance
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
model1.yhat.labels[model1.yhat.labels > 0.5000] = 1
model1.yhat.labels[model1.yhat.labels <= 0.5000] = 0

# make train prediction labels
model1.glm.fitted.labels <- model1.glm$fitted
model1.glm.fitted.labels[model1.glm.fitted.labels > 0.5000] = 1
model1.glm.fitted.labels[model1.glm.fitted.labels <= 0.5000] = 0

# summary
summary(model1.glm)

# train confusion matrix
train.1.conf.mat <- confusionMatrix(as.factor(model1.glm.fitted.labels), as.factor(train.1$act_demented),
                                    positive='1')

# test confusion matrix
test.1.conf.mat <- confusionMatrix(as.factor(model1.yhat.labels), as.factor(test.1$act_demented),
                                   positive='1')

# train ROC curve
train.1$prob <- model1.glm$fitted
model1.roc.train <- roc(act_demented~prob, data=train.1, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                        cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                        print.auc.col='black', identity.lwd=2.5)

# test ROC curve
test.1$prob <- model1.yhat
model1.roc.test <- roc(act_demented~prob, data=test.1, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# VIFs
model1.vif <- vif(model1.glm)

'
Model 2 - Weighted Logistic Regression with 16 Nonzero Features from Lasso w/ lowest CV missclassification
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
model2.yhat.labels[model2.yhat.labels > 0.5000] = 1
model2.yhat.labels[model2.yhat.labels <= 0.5000] = 0

# make train prediction labels
model2.glm.fitted.labels <- model2.glm$fitted
model2.glm.fitted.labels[model2.glm.fitted.labels > 0.5000] = 1
model2.glm.fitted.labels[model2.glm.fitted.labels <= 0.5000] = 0

# train confusion matrix
train.2.conf.mat <- confusionMatrix(as.factor(model2.glm.fitted.labels), as.factor(train.2$act_demented),
                                    positive='1')

# test confusion matrix
test.2.conf.mat <- confusionMatrix(as.factor(model2.yhat.labels), as.factor(test.2$act_demented),
                                   positive='1')

# train ROC curve
train.2$prob <- model2.glm$fitted
model2.roc.train <- roc(act_demented~prob, data=train.2, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# test ROC curve
test.2$prob <- model2.yhat
model2.roc.test <- roc(act_demented~prob, data=test.2, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# VIFs
model2.vif <- vif(model2.glm)

'
Model 3 - Overlapping Variables from Deviance and Misclassification Lists Only
'
features.3 <- c('ihc_ptdp_43_ffpe.TCx',
                'ihc_a_syn.PCx',
                #'ihc_tau2_ffpe.HIP',
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

train.3 <- train[ , names(train) %in% features.3]
test.3 <- test[ , names(test) %in% features.3]

# change act_demented to binary
train.3$act_demented <- as.numeric(train.3$act_demented)
train.3$act_demented[train.3$act_demented == 2] <- 0

test.3$act_demented <- as.numeric(test.3$act_demented)
test.3$act_demented[test.3$act_demented == 2] <- 0

# model 3 fit
model3.glm <- glm(act_demented ~.-obs_weight, data=train.3, weights=obs_weight, family=binomial(link='logit'))

# make predictions on test set
model3.yhat <- round(predict(model3.glm, newdata=test.3, type='response'), 4)

# make test prediction labels
model3.yhat.labels <- model3.yhat
model3.yhat.labels[model3.yhat.labels > 0.5000] = 1
model3.yhat.labels[model3.yhat.labels <= 0.5000] = 0

# make train prediction labels
model3.glm.fitted.labels <- model3.glm$fitted
model3.glm.fitted.labels[model3.glm.fitted.labels > 0.5000] = 1
model3.glm.fitted.labels[model3.glm.fitted.labels <= 0.5000] = 0

# summary
summary(model3.glm)

# train confusion matrix
train.3.conf.mat <- confusionMatrix(as.factor(model3.glm.fitted.labels), as.factor(train.3$act_demented),
                                    positive='1')

# test confusion matrix
test.3.conf.mat <- confusionMatrix(as.factor(model3.yhat.labels), as.factor(test.3$act_demented),
                                   positive='1')

# train ROC curve
train.3$prob <- model3.glm$fitted
model3.roc.train <- roc(act_demented~prob, data=train.3, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# test ROC curve
test.3$prob <- model3.yhat
model3.roc.test <- roc(act_demented~prob, data=test.3, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# VIFs
model3.vif <- vif(model3.glm)

'
Model 4 - Overlapping list with hippocampus tau2 levels added back
'
features.4 <- c('ihc_ptdp_43_ffpe.TCx',
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
model4.yhat.labels[model4.yhat.labels > 0.5000] = 1
model4.yhat.labels[model4.yhat.labels <= 0.5000] = 0

# make train prediction labels
model4.glm.fitted.labels <- model4.glm$fitted
model4.glm.fitted.labels[model4.glm.fitted.labels > 0.5] = 1
model4.glm.fitted.labels[model4.glm.fitted.labels <= 0.5] = 0

# train confusion matrix
train.4.conf.mat <- confusionMatrix(as.factor(model4.glm.fitted.labels), as.factor(train.4$act_demented),
                                    positive='1')

# test confusion matrix
test.4.conf.mat <- confusionMatrix(as.factor(model4.yhat.labels), as.factor(test.4$act_demented),
                                   positive='1')

# train ROC curve
train.4$prob <- model4.glm$fitted
model4.roc.train <- roc(act_demented~prob, data=train.4, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# test ROC curve
test.4$prob <- model4.yhat
model4.roc.test <- roc(act_demented~prob, data=test.4, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# VIFs
model4.vif <- vif(model4.glm)

'
Model 5 - Removed Underperforming Features from Model 2 (misclassification error variable list)
'
features.5 <- c('ihc_ptdp_43_ffpe.TCx',
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

train.5 <- train[ , names(train) %in% features.5]
test.5 <- test[ , names(test) %in% features.5]

# change act_demented to binary
train.5$act_demented <- as.numeric(train.5$act_demented)
train.5$act_demented[train.5$act_demented == 2] <- 0

test.5$act_demented <- as.numeric(test.5$act_demented)
test.5$act_demented[test.5$act_demented == 2] <- 0

model5.glm <- glm(act_demented ~.-obs_weight, data=train.5, weights=obs_weight, family=binomial(link='logit'))

# summary
summary(model5.glm)

# make predictions on test set
model5.yhat <- round(predict(model5.glm, newdata=test.5, type='response'), 4)

# make test prediction labels
model5.yhat.labels <- model5.yhat
model5.yhat.labels[model5.yhat.labels > 0.5000] = 1
model5.yhat.labels[model5.yhat.labels <= 0.5000] = 0

# make train prediction labels
model5.glm.fitted.labels <- model5.glm$fitted
model5.glm.fitted.labels[model5.glm.fitted.labels > 0.5000] = 1
model5.glm.fitted.labels[model5.glm.fitted.labels <= 0.5000] = 0

# train confusion matrix
train.5.conf.mat <- confusionMatrix(as.factor(model5.glm.fitted.labels), as.factor(train.5$act_demented),
                                    positive='1')

# test confusion matrix
test.5.conf.mat <- confusionMatrix(as.factor(model5.yhat.labels), as.factor(test.5$act_demented),
                                    positive='1')

# train ROC curve
train.5$prob <- model5.glm$fitted
model5.roc.train <- roc(act_demented~prob, data=train.5, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# test ROC curve
test.5$prob <- model5.yhat
model5.roc.test <- roc(act_demented~prob, data=test.5, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# VIFs
model5.vif <- vif(model5.glm)

'
Model 6 - Model 4 + HIP cluster 2 medoid
'
features.6 <- c('ihc_ptdp_43_ffpe.TCx',
                'ihc_a_syn.PCx',
                'ihc_tau2_ffpe.HIP',
                'ihc_at8_ffpe.PCx',
                #'ihc_a_syn.HIP',
                'ihc_gfap_ffpe.PCx',
                'a_syn_pg_per_mg.PCx',
                'apo_e4_allele',
                'mcp_1_pg_per_mg.TCx',
                'ihc_iba1_ffpe.PCx',
                'gene_cluster02.HIP',  # added
                'obs_weight',
                'act_demented')

train.6 <- train[ , names(train) %in% features.6]
test.6 <- test[ , names(test) %in% features.6]

# change act_demented to binary
train.6$act_demented <- as.numeric(train.6$act_demented)
train.6$act_demented[train.6$act_demented == 2] <- 0

test.6$act_demented <- as.numeric(test.6$act_demented)
test.6$act_demented[test.6$act_demented == 2] <- 0

# model 6 fit
model6.glm <- glm(act_demented ~.-obs_weight, data=train.6, weights=obs_weight, family=binomial(link='logit'))

# make predictions on test set
model6.yhat <- round(predict(model6.glm, newdata=test.6, type='response'), 4)

# make test prediction labels
model6.yhat.labels <- model6.yhat
model6.yhat.labels[model6.yhat.labels > 0.5000] = 1
model6.yhat.labels[model6.yhat.labels <= 0.5000] = 0

# make train prediction labels
model6.glm.fitted.labels <- model6.glm$fitted
model6.glm.fitted.labels[model6.glm.fitted.labels > 0.5000] = 1
model6.glm.fitted.labels[model6.glm.fitted.labels <= 0.5000] = 0

# summary
summary(model6.glm)

# train confusion matrix
train.6.conf.mat <- confusionMatrix(as.factor(model6.glm.fitted.labels), as.factor(train.6$act_demented),
                                    positive='1')

# test confusion matrix
test.6.conf.mat <- confusionMatrix(as.factor(model6.yhat.labels), as.factor(test.6$act_demented),
                                   positive='1')

# train ROC curve
train.6$prob <- model6.glm$fitted
model6.roc.train <- roc(act_demented~prob, data=train.6, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# test ROC curve
test.6$prob <- model6.yhat
model6.roc.test <- roc(act_demented~prob, data=test.6, plot=TRUE,  grid=TRUE, print.auc=TRUE,
                       cex.axis=1.5, cex.lab=1.5, lwd=4, print.auc.cex=2, col='darkorchid3',
                       print.auc.col='black', identity.lwd=2.5)

# VIFs
model6.vif <- vif(model6.glm)

'
Pretty tables
'
# confusion matrices for Model 1
train.1.df <- data.frame(c(40,3),
                         c(14,23))
colnames(train.1.df) <- c('No Dementia', 'Dementia')
rownames(train.1.df) <- c('No Dementia', 'Dementia')

table07 <- xtable(train.1.df)
align(table07) <- 'r|c|c'
digits(table07) <- 0
print(table07, hline.after = c(0,1))

print.xtable(table07, type='latex', file='data/table07.tex')

test.1.df <- data.frame(c(12,2),
                        c(7,6))
colnames(test.1.df) <- c('No Dementia', 'Dementia')
rownames(test.1.df) <- c('No Dementia', 'Dementia')

table08 <- xtable(test.1.df)
align(table08) <- 'r|c|c'
digits(table08) <- 0
print(table08, hline.after = c(0,1))

print.xtable(table08, type='latex', file='data/table08.tex')

# confusion matrices for Model 2
train.2.df <- data.frame(c(41,2),
                         c(10,27))
colnames(train.2.df) <- c('No Dementia', 'Dementia')
rownames(train.2.df) <- c('No Dementia', 'Dementia')

table09 <- xtable(train.2.df)
align(table09) <- 'r|c|c'
digits(table09) <- 0
print(table09, hline.after = c(0,1))

print.xtable(table09, type='latex', file='data/table09.tex')

test.2.df <- data.frame(c(8,6),
                        c(8,5))
colnames(test.2.df) <- c('No Dementia', 'Dementia')
rownames(test.2.df) <- c('No Dementia', 'Dementia')

table10 <- xtable(test.2.df)
align(table10) <- 'r|c|c'
digits(table10) <- 0
print(table10, hline.after = c(0,1))

print.xtable(table10, type='latex', file='data/table10.tex')

# confusion matrices for Model 3
train.3.df <- data.frame(c(40,3),
                         c(16,21))
colnames(train.3.df) <- c('No Dementia', 'Dementia')
rownames(train.3.df) <- c('No Dementia', 'Dementia')

table11 <- xtable(train.3.df)
align(table11) <- 'r|c|c'
digits(table11) <- 0
print(table11, hline.after = c(0,1))

print.xtable(table11, type='latex', file='data/table11.tex')

test.3.df <- data.frame(c(11,3),
                        c(8,5))
colnames(test.3.df) <- c('No Dementia', 'Dementia')
rownames(test.3.df) <- c('No Dementia', 'Dementia')

table12 <- xtable(test.3.df)
align(table12) <- 'r|c|c'
digits(table12) <- 0
print(table12, hline.after = c(0,1))

print.xtable(table12, type='latex', file='data/table12.tex')

# confusion matrices for Model 4
train.4.df <- data.frame(c(40,3),
                         c(14,23))
colnames(train.4.df) <- c('No Dementia', 'Dementia')
rownames(train.4.df) <- c('No Dementia', 'Dementia')

table13 <- xtable(train.4.df)
align(table13) <- 'r|c|c'
digits(table13) <- 0
print(table13, hline.after = c(0,1))

print.xtable(table13, type='latex', file='data/table13.tex')

test.4.df <- data.frame(c(12,2),
                        c(7,6))
colnames(test.4.df) <- c('No Dementia', 'Dementia')
rownames(test.4.df) <- c('No Dementia', 'Dementia')

table14 <- xtable(test.4.df)
align(table14) <- 'r|c|c'
digits(table14) <- 0
print(table14, hline.after = c(0,1))

print.xtable(table14, type='latex', file='data/table14.tex')

# confusion matrices for Model 5
train.5.df <- data.frame(c(34,12),
                         c(12,25))
colnames(train.5.df) <- c('No Dementia', 'Dementia')
rownames(train.5.df) <- c('No Dementia', 'Dementia')

table15 <- xtable(train.5.df)
align(table15) <- 'r|c|c'
digits(table15) <- 0
print(table15, hline.after = c(0,1))

print.xtable(table15, type='latex', file='data/table15.tex')

test.5.df <- data.frame(c(8,6),
                        c(9,4))
colnames(test.5.df) <- c('No Dementia', 'Dementia')
rownames(test.5.df) <- c('No Dementia', 'Dementia')

table16 <- xtable(test.5.df)
align(table16) <- 'r|c|c'
digits(table16) <- 0
print(table16, hline.after = c(0,1))

print.xtable(table16, type='latex', file='data/table16.tex')

# confusion matrices for Model 6
train.6.df <- data.frame(c(41,2),
                         c(15,22))
colnames(train.6.df) <- c('No Dementia', 'Dementia')
rownames(train.6.df) <- c('No Dementia', 'Dementia')

table17 <- xtable(train.6.df)
align(table17) <- 'r|c|c'
digits(table17) <- 0
print(table17, hline.after = c(0,1))

print.xtable(table17, type='latex', file='data/table17.tex')

test.6.df <- data.frame(c(12,2),
                        c(8,5))
colnames(test.6.df) <- c('No Dementia', 'Dementia')
rownames(test.6.df) <- c('No Dementia', 'Dementia')

table18 <- xtable(test.6.df)
align(table18) <- 'r|c|c'
digits(table18) <- 0
print(table18, hline.after = c(0,1))

print.xtable(table18, type='latex', file='data/table18.tex')