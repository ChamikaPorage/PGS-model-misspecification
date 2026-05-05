#load required packages
library(devtools)
library(boot)
library(MASS) 
library(randomForest)
library(rpart)
library(xtable)
library(PSW)
library(iWeigReg)
library(cobalt)
library(WeightIt)
library(table1)
library(readxl)
library(caret)
library(xgboost)
library(dplyr)
library(ranger)
library(Matching)
#read data from dfmale.csv and for_table8.csv
datat <- read.csv(file = "C:/Users/chapo752/Dropbox/PhD work- Chamika Porage/Second paper/R codes/Emp_ana_new/df_allgenders.csv")

datat <- datat %>%
  mutate(
    gender = case_when(
      gender == 1 ~ 1,   # keep males as 1
      gender == 2 ~ 0,   # recode females as 0
      TRUE ~ NA_real_    # anything else becomes NA
    )
  )

#covariates
X.matrix <- model.matrix(smoking ~ married.living.with.partner + birth.country + edu  + race + income + army.service + c.age + c.age2 + c.family.size + gender, data = datat)

dat.X <- as.data.frame(X.matrix)
dat.X <- dat.X[,-1]
colnames(dat.X)
dat.X <- dat.X[, !names(dat.X) %in% "incomeOther"]
#names(dat.X) <- c(LETTERS[1:26])

dat <- cbind(datat$smoking, datat$lead, dat.X)

names(dat) <- c("smoking", "lead", LETTERS[1:26])


#PS MODELS

#true ps
mod.ps <- glm(smoking ~ A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z, family = binomial, data = dat)
ps <- fitted.values(mod.ps, type = "response") 
datat$ps <- ps

# Subsetting data
data0 <- subset(dat, smoking == 0)
data1 <- subset(dat, smoking == 1)

#Outcome Regression models(prognostic scores)

#TRUE OR
mod0 <- lm(lead ~ A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z, data = data0)
mu0 <- predict(mod0, newdata = dat, type = "response")

mod1 <- lm(lead ~ A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z, data = data1)
mu1 <- predict(mod1, newdata = dat, type = "response")



###Non-linear models###

folds <- createFolds(dat$lead, k = 5, list = TRUE)

# Initializing prediction vectors

N <- nrow(dat)
p0rf <- p1rf  <- rep(NA, nrow(dat))


for (k in 1:5) {
  # Training and validation indices for fold k
  train_id <- unlist(folds[-k])  # Use all folds except fold k for training
  test_id <- folds[[k]]          # Use fold k for testing
  
  # Separate training and validation data
  train_data <- dat[train_id, ]
  test_data <- dat[test_id, ]
  
  # Subset training data into treatment groups
  data0_train <- subset(train_data, smoking == 0)
  data1_train <- subset(train_data, smoking == 1)
  
  # True model
  
  # Train random forest on untreated group (smoking == 0)
  model_rf_0 <- ranger(lead ~ A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z, 
                       data = data0_train, 
                       num.trees = 300, 
                       mtry = 2)
  p0rf[test_id] <- predict(model_rf_0, data = test_data[, !(colnames(test_data) %in% c("smoking", "lead"))])$predictions
  
  # Train random forest on treated group (smoking == 1)
  model_rf_1 <- ranger(lead ~ A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z, 
                       data = data1_train, 
                       num.trees = 300, 
                       mtry = 2)
  p1rf[test_id] <- predict(model_rf_1, data = test_data[, !(colnames(test_data) %in% c("smoking", "lead"))])$predictions
  
} 


Dat<-cbind.data.frame(dat$lead,dat$smoking,ps, mu0,mu1,p0rf,p1rf)
colnames(Dat) <- c('lead','smoking','ps','mu0','mu1','p0rf','p1rf')

# Define Prognostic Score Matrices
x_pgs_t <- as.matrix(Dat[, c("mu0")])
x_pgs_nl_t <- as.matrix(Dat[, c("p0rf")])


#Outcome Regression models(full prognostic scores)

Dat1 <- subset(Dat, Dat$smoking == 1)
Dat0 <- subset(Dat, Dat$smoking == 0)

###FPGS estimators

# Linear model #

### True model ###
model1 <- lm(lead ~ mu1 + mu0, data = Dat1)
pred1 <- predict(model1, newdata = Dat)

model0 <- lm(lead ~ mu1 + mu0, data = Dat0)
pred0 <- predict(model0, newdata = Dat)

yp <- ifelse(Dat$smoking == 1, pred1, pred0)

# Non-linear model #

# Initializing prediction vectors

N <- nrow(dat)
pred1_rf <- pred0_rf <- pre1f_rf <- pre0f_rf <- rep(NA, N)

for (k in 1:5) {
  # Training and validation indices for fold k
  train_id <- unlist(folds[-k])
  test_id <- folds[[k]]
  
  # Separate training and validation data
  train_data <- Dat[train_id, ]
  test_data <- Dat[test_id, ]
  
  # Subset training data into treatment groups
  data0_train <- subset(train_data, smoking == 0)
  data1_train <- subset(train_data, smoking == 1)
  
  # True model 
  
  # Train Random Forest on untreated group
  model_rf_0 <- ranger(lead ~ p0rf +  p1rf, data = data0_train, num.tree = 300, mtry = 2)
  pred0_rf[test_id] <- predict(model_rf_0, data = test_data[,c("p0rf" ,"p1rf")])$predictions
  
  # Train Random Forest on treated group
  model_rf_1 <- ranger(lead ~ p0rf +  p1rf, data = data1_train, num.tree = 300, mtry = 2)
  pred1_rf[test_id] <- predict(model_rf_1, data = test_data[,c("p0rf" ,"p1rf")])$predictions
  
}

# propensity score for full prognostic scores
#mod.ps_g <- glm(smoking ~ mu1 + mu0 , family = binomial, data = dat)
#ps_g <- predict(mod.ps_g, type = "response")

#mod.psf_g <- glm(smoking ~ mu1f + mu0f , family = binomial(link = "logit"), data = dat)
#psf_g <- predict(mod.psf_g, type = "response")

Dat_t<-cbind.data.frame(dat$lead,dat$smoking,ps, mu0,mu1,p0rf,p1rf, pred1, pred0, pred1_rf, pred0_rf)
colnames(Dat_t) <- c('lead','smoking','ps','mu0','mu1','p0rf', 'p1rf', 'pred1' , 'pred0', 'pred1_rf', 'pred0_rf')

x_fpgs_t <- as.matrix(cbind(pred0, pred1))
x_fpgs_nl_t <- as.matrix(cbind(pred0_rf,pred1_rf))

### RI Estimators

###PGS estimators

#Linear models

#1.linear
ri_pgt <- mean(Dat_t$mu1) - mean(Dat_t$mu0)

#2.non-linear
rf_pgt <- mean(Dat_t$p1rf) - mean(Dat_t$p0rf)

###FPGS estimators

#3.linear
ri_fpgt <- mean(pred1) - mean(pred0)

#4. non-linear
rf_fpgt <- mean(pred1_rf) - mean(pred0_rf)


###. matching estimator

#5.-linear PGS
match_lin_pgs_t <- Match(Y=Dat_t$lead, Tr=Dat_t$smoking, X=x_pgs_t, estimand = "ATE", M = 1)
match_lin_pgs_t$est

#6.-non-linear PGS
match_nn_pgs_t <- Match(Y=Dat_t$lead, Tr=Dat_t$smoking, X=x_pgs_nl_t, estimand = "ATE", M = 1)
match_nn_pgs_t$est

#7.-linear FPGS
match_lin_fpgs_t <- Match(Y=Dat_t$lead, Tr=Dat_t$smoking, X=x_fpgs_t, estimand = "ATE", M = 1)
match_lin_fpgs_t$est

#8.-non-linear FPGS
match_nn_fpgs_t <- Match(Y=Dat_t$lead, Tr=Dat_t$smoking, X=x_fpgs_nl_t, estimand = "ATE", M = 1)
match_nn_fpgs_t$est


# Load the necessary library for bootstrap
library(boot)

# Set the number of bootstrap samples
num_bootstrap_samples <- 1000

### Define bootstrap functions for each of the 12 estimators

### 1. linear PGS (ri_pgt)
boot_ri_pgt <- function(data, indices) {
  data_s <- data[indices, ]
  ri_pgt <- mean(data_s$mu1) - mean(data_s$mu0)
  return(ri_pgt)
}


### 2. non-linear PGS (rf_pgt)
boot_rf_pgt <- function(data, indices) {
  data_s <- data[indices, ]
  rf_pgt <- mean(data_s$p1rf) - mean(data_s$p0rf)
  return(rf_pgt)
}


### 3. linear FPGS (ri_fpgt)
boot_ri_fpgt <- function(data, indices) {
  data_s <- data[indices, ]
  ri_fpgt <- mean(data_s$pred1) - mean(data_s$pred0)
  return(ri_fpgt)
}


### 4. non-linear FPGS (rf_fpgt)
boot_rf_fpgt <- function(data, indices) {
  data_s <- data[indices, ]
  rf_fpgt <- mean(data_s$pred1_rf) - mean(data_s$pred0_rf)
  return(rf_fpgt)
}


###5. Matching - PGS-parametric
boot_mat_lin_pgs_t <- function(data, indices) {
  data_s <- data[indices, ]  # resample the data
  mat_lin_pgs_t <- Match(Y = data_s$lead, Tr = data_s$smoking, X = x_pgs_t[indices, ], estimand = "ATE", M = 1)
  return(mat_lin_pgs_t$est)
}


###6. Matching -PGS - nonparametric
boot_mat_nn_pgs_t <- function(data, indices) {
  data_s <- data[indices, ]  # resample the data
  mat_nn_pgs_t <- Match(Y = data_s$lead, Tr = data_s$smoking, X = x_pgs_nl_t[indices, ], estimand = "ATE", M = 1)
  return(mat_nn_pgs_t$est)
}


###7. Matching - FPGS-parametric
boot_mat_lin_fpgs_t <- function(data, indices) {
  data_s <- data[indices, ]  # resample the data
  mat_lin_fpgs_t <- Match(Y = data_s$lead, Tr = data_s$smoking, X = x_fpgs_t[indices, ], estimand = "ATE", M = 1)
  return(mat_lin_fpgs_t$est)
}

###8. Matching -PGS - nonparametric
boot_mat_nn_fpgs_t <- function(data, indices) {
  data_s <- data[indices, ]  # resample the data
  mat_nn_fpgs_t <- Match(Y = data_s$lead, Tr = data_s$smoking, X = x_fpgs_nl_t[indices, ], estimand = "ATE", M = 1)
  return(mat_nn_fpgs_t$est)
}


### Perform Bootstrap for Each Estimator
# Apply bootstrapping to calculate standard errors and confidence intervals for all 12 estimators.

# 1. True-linear PGS (ri_pgt)
res_ri_pgt <- boot(data = Dat_t, statistic = boot_ri_pgt, R = num_bootstrap_samples)
ri_pgt_se <- sd(res_ri_pgt$t)
conf_int_ri_pgt <- boot.ci(res_ri_pgt, type = "perc")


# 2. True-non-linear PGS (rf_pgt)
res_rf_pgt <- boot(data = Dat_t, statistic = boot_rf_pgt, R = num_bootstrap_samples)
rf_pgt_se <- sd(res_rf_pgt$t)
conf_int_rf_pgt <- boot.ci(res_rf_pgt, type = "perc")


# 3. True-linear FPGS (ri_fpgt)
res_ri_fpgt <- boot(data = Dat_t, statistic = boot_ri_fpgt, R = num_bootstrap_samples)
ri_fpgt_se <- sd(res_ri_fpgt$t)
conf_int_ri_fpgt <- boot.ci(res_ri_fpgt, type = "perc")

# 4. True-non-linear FPGS (rf_fpgt)
res_rf_fpgt <- boot(data = Dat_t, statistic = boot_rf_fpgt, R = num_bootstrap_samples)
rf_fpgt_se <- sd(res_rf_fpgt$t)
conf_int_rf_fpgt <- boot.ci(res_rf_fpgt, type = "perc")

# 5. Matching - True PGS-parametric
res_mat_lin_pgs_t <- boot(data = Dat_t, statistic = boot_mat_lin_pgs_t, R = num_bootstrap_samples)
mat_lin_pgs_t_se <- sd(res_mat_lin_pgs_t$t)
conf_int_mat_lin_pgs_t <- boot.ci(res_mat_lin_pgs_t, type = "perc")

# 6. Matching - True PGS - non-parametric
res_mat_nn_pgs_t <- boot(data = Dat_t, statistic = boot_mat_nn_pgs_t, R = num_bootstrap_samples)
mat_nn_pgs_t_se <- sd(res_mat_nn_pgs_t$t)
conf_int_mat_nn_pgs_t<- boot.ci(res_mat_nn_pgs_t, type = "perc")

# 7. Matching - True FPGS-parametric
res_mat_lin_fpgs_t <- boot(data = Dat_t, statistic = boot_mat_lin_fpgs_t, R = num_bootstrap_samples)
mat_lin_fpgs_t_se <- sd(res_mat_lin_fpgs_t$t)
conf_int_mat_lin_fpgs_t <- boot.ci(res_mat_lin_fpgs_t, type = "perc")

# 8. Matching - True FPGS - non-parametric
res_mat_nn_fpgs_t <- boot(data = Dat_t, statistic = boot_mat_nn_fpgs_t, R = num_bootstrap_samples)
mat_nn_fpgs_t_se <- sd(res_mat_nn_fpgs_t$t)
conf_int_mat_nn_fpgs_t<- boot.ci(res_mat_nn_fpgs_t, type = "perc")

### Print Results for Standard Error and 95% Confidence Intervals

# 1. linear PGS (ri_pgt)
print(paste("Bootstrap SE for ri_pgt:", ri_pgt_se))
print(paste("95% CI for ri_pgt:", conf_int_ri_pgt$percent[4], "-", conf_int_ri_pgt$percent[5]))

# 2. non-linear PGS (rf_pgt)
print(paste("Bootstrap SE for rf_pgt:", rf_pgt_se))
print(paste("95% CI for rf_pgt:", conf_int_rf_pgt$percent[4], "-", conf_int_rf_pgt$percent[5]))


# 3. linear FPGS (ri_fpgt)
print(paste("Bootstrap SE for ri_fpgt:", ri_fpgt_se))
print(paste("95% CI for ri_fpgt:", conf_int_ri_fpgt$percent[4], "-", conf_int_ri_fpgt$percent[5]))


# 4. non-linear FPGS (rf_fpgt)
print(paste("Bootstrap SE for rf_fpgt:", rf_fpgt_se))
print(paste("95% CI for rf_fpgt:", conf_int_rf_fpgt$percent[4], "-", conf_int_rf_fpgt$percent[5]))


# 5. Matching: linear PGS 
print(paste("Bootstrap SE for mat_pgt:", mat_lin_pgs_t_se))
print(paste("95% CI for mat_pgt:", conf_int_mat_lin_pgs_t$percent[4], "-", conf_int_mat_lin_pgs_t$percent[5]))

# 6. Matching:non-linear PGS 
print(paste("Bootstrap SE for mat_rf_pgt:", mat_nn_pgs_t_se))
print(paste("95% CI for mat_rf_pgt:", conf_int_mat_nn_pgs_t$percent[4], "-", conf_int_mat_nn_pgs_t$percent[5]))

# 7. Matching:linear FPGS 
print(paste("Bootstrap SE for mat_fpgt:", mat_lin_fpgs_t_se))
print(paste("95% CI for mat_fpgt:", conf_int_mat_lin_fpgs_t$percent[4], "-", conf_int_mat_lin_fpgs_t$percent[5]))

# 8. Matching:non-linear FPGS 
print(paste("Bootstrap SE for mat_rf_fpgt:", mat_nn_fpgs_t_se))
print(paste("95% CI for mat_rf_fpgt:", conf_int_mat_nn_fpgs_t$percent[4], "-", conf_int_mat_nn_fpgs_t$percent[5]))

