library(here)
library(tidyverse)
library(viridis)
library(glmnet)
data.dir <- here('data', 'cleandata')
#Loading the data
enif2021_subset <- read_csv(here(data.dir, 'risk_assesment_nsfi2021.csv'),
                            show_col_types = F)
enif2021_subset <- enif2021_subset %>%
  mutate(across(where(is.character), as.factor))
enif2021_subset <-
  enif2021_subset %>%
  mutate(
    school_level =
      fct_relevel(school_level,
                  c('No school', 'Preschool','Primary','Secundary',
                    'Technical-secundary','Basic Normal school',
                    'High school', 'Technical-high school',
                    'Bachelors', 'Graduate')),
    marital_status = fct_relevel(marital_status,
                                 c('Single' , 'Not married', 'Married',
                                   'Apart', 'Divorced', 'Widow'))
  )
#Splitting the data in K-folds (K = 18)
K_folds <- 3; set.seed(2023-13-11) #A seed for reproducibility
fold_id <-
  #Create a vector that indicates the fold number for each sample
  rep(1:K_folds, each = nrow(enif2021_subset)/K_folds) %>%
  #Compute a random permutation of the vector that contains
  #the fold number for each sample. In this way we try to
  #avoid any particular clustering of the data in a given fold
  sample()
enif2021_subset <-
  enif2021_subset %>% mutate('fold_class' = as.factor(fold_id))
#Splitting the data in three small data sets for the
#different analysis to be performed:
#a) Data for model selection
model_selection.data <- enif2021_subset %>% filter(fold_class==1) %>%
  select(!c('FOLIO', 'fold_class'))
#b) Data for model fitting and inference
fitting_inference.data <- enif2021_subset %>% filter(fold_class==2) %>%
  select(!c('FOLIO', 'fold_class'))
#c) Data for model prediction and performance assessment
model_assesment.data <- enif2021_subset %>% filter(fold_class==3) %>%
  select(!c('FOLIO', 'fold_class'))
#We will fit a first model that includes only all the main effects
model1_formula <- earned_enough ~
  .
model1.matrix <- model.matrix(model1_formula, data = model_selection.data)
model1.matrix <- model1.matrix[,-1] #Removing intercept for glmnet
#glmnet variable selection
set.seed(2023-16-12) #A seed for reproducibility
variable_selection.canonical <-
  cv.glmnet(
    y = model_selection.data$earned_enough, x = model1.matrix,
    family = 'binomial', nfolds = 18)
coefs.canonical <- coef(variable_selection.canonical)
set.seed(2023-16-12) #A seed for reproducibility
variable_selection.probit <-
  cv.glmnet(
    y = model_selection.data$earned_enough, x = model1.matrix,
    family = binomial(link = 'probit'), nfolds = 18)
coefs.probit <- coef(variable_selection.probit)
set.seed(2023-16-12) #A seed for reproducibility
variable_selection.cloglog <-
  cv.glmnet(
    y = model_selection.data$earned_enough, x = model1.matrix,
    family = binomial(link = 'cloglog'), nfolds = 18
  )
coefs.cloglog <- coef(variable_selection.cloglog)
cbind(
  rownames(coefs.canonical)[as.matrix(coefs.canonical) != 0],
  rownames(coefs.probit)[as.matrix(coefs.probit) != 0],
  rownames(coefs.cloglog)[as.matrix(coefs.cloglog) != 0]
)

#Formulas for selected models
model_canonical <- earned_enough ~ school_level+ speak_indigenouslang +
  any_smartphone + earning_expenses_history + automatic_monthly_charges +
  savings_account + num_bathrooms + any_car + internet_access
model_probit <- earned_enough ~ school_level+ speak_indigenouslang +
  any_smartphone + earning_expenses_history + automatic_monthly_charges +
  savings_account + num_bathrooms + any_car + internet_access
model_cloglog <- earned_enough ~ school_level+ speak_indigenouslang +
  any_smartphone + earning_expenses_history + automatic_monthly_charges +
  savings_account + num_bathrooms + any_car + internet_access

##### Model fitting steps
inference_canonical <- glm(model_canonical, family = binomial(),
                           data = fitting_inference.data)
inference_probit <- glm(model_probit, family = binomial(link = 'probit'),
                        data = fitting_inference.data)
inference_cloglog <- glm(model_cloglog, family = binomial(link = 'cloglog'),
                         data = fitting_inference.data)
inference_canonical %>% summary(): inference_probit %>% summary()
inference_cloglog %>% summary()

######### Model diagnostics
#Checking the systematic component
qresids_canonical <- statmod::qresid(inference_canonical)
qresids_probit <- statmod::qresid(inference_probit)
qresids_cloglog <- statmod::qresid(inference_cloglog)
par(mfrow = c(1,3))
scatter.smooth(qresids_canonical ~ fitted(inference_canonical),
               ylab = 'Quantile residuals', xlab = 'Fitted probabilities',
               main = 'Logit link', lpars = list(col = 'red', lwd = 3))
abline(h = 0, col = 'blue', lwd = 2)
scatter.smooth(qresids_probit ~ fitted(inference_canonical),
               ylab = 'Quantile residuals', xlab = 'Fitted probabilities',
               main = 'Probit link', lpars = list(col = 'red', lwd = 3))
abline(h = 0, col = 'blue', lwd = 2)
scatter.smooth(qresids_cloglog ~ fitted(inference_canonical),
               ylab = 'Quantile residuals', xlab = 'Fitted probabilities',
               main = 'Complementary Log-Log link',
               lpars = list(col = 'red', lwd = 3))
abline(h = 0, col = 'blue', lwd = 2); par(mfrow = c(1,1))
#Cheking the random component
par(mfrow = c(1,3))
qqnorm(qresids_canonical, main = 'Normal Q-Q Plot: Logit link')
qqline(qresids_canonical)
qqnorm(qresids_probit, main = 'Normal Q-Q Plot: Probit link')
qqline(qresids_probit)
qqnorm(qresids_cloglog, main = 'Normal Q-Q Plot: Complementary Log-Log link')
qqline(qresids_cloglog)
par(mfrow = c(1,1))
#Cheking for outliers
rstud_canonical <- rstudent(inference_canonical)
rstud_probit <- rstudent(inference_probit)
rstud_cloglog <- rstudent(inference_cloglog)
par(mfrow = c(1,3))
plot(abs(rstud_canonical), ylab = 'abs(Studentized deviance residuals)',
     main = 'Logit link', ylim = c(0, 3.5),)
abline(h = 2, col = 'orange', lwd = 2); abline(h = 3, col = 'red', lwd = 2)
plot(abs(rstud_probit), ylab = ''
     , main = 'Probit link', ylim = c(0, 3.5),)
abline(h = 2, col = 'orange', lwd = 2); abline(h = 3, col = 'red', lwd = 2)
plot(abs(rstud_cloglog), ylab = ''
     , main = 'Complementary Log-Log link',
     ylim = c(0, 3.5))
abline(h = 2, col = 'orange', lwd = 2); abline(h = 3, col = 'red', lwd = 2)
par(mfrow = c(1,1))
#Cheking for influential onbservations
cooksD_canonical <- cooks.distance(inference_canonical)
cooksD_probit <- cooks.distance(inference_probit)
cooksD_cloglog <- cooks.distance(inference_cloglog)
flag_val <- 4/nrow(fitting_inference.data)
par(mfrow = c(1,3))
plot(cooksD_canonical, type = 'h', main = 'Logit link',
     ylab = "Cook's distance")
abline(h = flag_val, col = 'red', lwd = 2)
plot(cooksD_probit, type = 'h',
     main = 'Probit link', ylab = "Cook's distance")
abline(h = flag_val, col = 'red', lwd = 2)
plot(cooksD_cloglog, type = 'h',
     main = 'Complementary Log-Log link', ylab = "Cook's distance")
abline(h = flag_val, col = 'red', lwd = 2); par(mfrow = c(1,1))
#Removing potential influential point to see what happens
inference_canonical.2 <- update(inference_canonical,
                                subset = (-which.max(cooksD_canonical)))
inference_probit.2 <- update(inference_probit,
                             subset = (-which.max(cooksD_probit)))
inference_cloglog.2 <- update(inference_cloglog,
                              subset = (-which.max(cooksD_cloglog)))
#Comparing the coefficients to see if the results change
cbind(inference_canonical$coefficients,inference_canonical.2$coefficients)
cbind(inference_probit$coefficients, inference_probit.2$coefficients)
cbind(inference_cloglog$coefficients, inference_cloglog.2$coefficients)
# ROC and AUC for trainning data
logit_ROC <- pROC::roc(inference_canonical.2$model$earned_enough ~
                         inference_canonical.2$fitted.values)
probit_ROC <- pROC::roc(inference_probit.2$model$earned_enough ~
                          inference_canonical.2$fitted.values)
cloglog_ROC <- pROC::roc(inference_cloglog.2$model$earned_enough ~
                           inference_cloglog.2$fitted.values)
logit_ROC %>% plot( xlim = c(1,0), ylim = c(0,1))
probit_ROC %>% plot(add = T, col = 'red')
cloglog_ROC %>% plot(add = T, col = 'blue')
#Assesing predictive performance
mask.assesment <- sample(1:nrow(model_assesment.data),
                         size = nrow(model_assesment.data)/2)
model_assesment.data1 <- model_assesment.data[mask.assesment,]
model_assesment.data2 <- model_assesment.data[-mask.assesment,]
#ROC for selecting the final model
logit.assesment_predictions <- predict(inference_canonical.2,
                                       newdata = model_assesment.data1,
                                       type = 'response')
probit.assesment_predictions <- predict(inference_probit.2,
                                        newdata = model_assesment.data1,
                                        type = 'response')
cloglog.assesment_predictions <- predict(inference_cloglog.2,
                                         newdata = model_assesment.data1,
                                         type = 'response')
logit_ROC.assesment <- pROC::roc(model_assesment.data1$earned_enough ~
                                   logit.assesment_predictions)
probit_ROC.assesment <- pROC::roc(model_assesment.data1$earned_enough ~
                                    probit.assesment_predictions)
cloglog_ROC.assesment <- pROC::roc(model_assesment.data1$earned_enough ~
                                     cloglog.assesment_predictions)
auc_byModel <-
  data.frame(
    'Model' = c('Logit', 'Probit', 'CLogLog'),
    'ROC_fit' = c(logit_ROC$auc, probit_ROC$auc, cloglog_ROC$auc),
    'ROC_assesment' = c(logit_ROC.assesment$auc, probit_ROC.assesment$auc,
                        cloglog_ROC.assesment$auc)
  )
#Final predictions
logit.final_preds <- predict(inference_canonical.2,
                             newdata = model_assesment.data2,
                             type = 'response')
logit.finalROC <- pROC::roc(model_assesment.data2$earned_enough ~
                              logit.final_preds)
probit.final_preds <- predict(inference_probit.2,
                              newdata = model_assesment.data2,
                              type = 'response')
probit.finalROC <- pROC::roc(model_assesment.data2$earned_enough ~
                               probit.final_preds)
cloglog.final_preds <- predict(inference_cloglog.2,
                               newdata = model_assesment.data2,
                               type = 'response')
cloglog.finalROC <- pROC::roc(model_assesment.data2$earned_enough ~
                                cloglog.final_preds)
auc_byModel$'Real_auc' <- c( logit.finalROC$auc, probit.finalROC$auc,
                             cloglog.finalROC$auc)
auc_byModel