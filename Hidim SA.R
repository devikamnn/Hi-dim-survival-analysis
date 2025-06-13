library(readxl)
data=read_excel("survival1.xlsx")
str(data)
#preprocessing 
data=na.omit(data)
sum(is.na(data))
data <- data[data$PFI.time != 0, ]
dim(data)
subset=data[,7:6006]
str(subset)
subset<- as.data.frame(subset)
subset[] <- lapply(subset, function(col) { 
  if (is.character(col)) {
    as.numeric(as.character(col))  
  } else {
    col  
  }
})

library(survival)
library(glmnet)
str(subset)
sum(is.na(subset))
Surv_obj <- Surv(data$PFI.time,data$PFI)
num_censored <- sum(Surv_obj[, "status"] == 0)
total <- length(Surv_obj[, "status"])
percent_censored <- (num_censored / total) * 100
subset=as.matrix(subset)

#alpha = 0 for ridge regression
fit_ridge <- glmnet(subset, Surv_obj, family = "cox", alpha = 0)
cv_fit <- cv.glmnet(subset, Surv_obj, family = "cox", alpha = 0) #find optimal lambda
best_lambda <- cv_fit$lambda.min
final_ridge_model <- glmnet(subset, Surv_obj, family = "cox", alpha = 0, lambda = best_lambda)
ridge_coefs <- coef(final_ridge_model, s = "lambda.min")
selected_features <- which(abs(ridge_coefs) > 0.01) 
coef(final_ridge_model)
risk_scores <- predict(final_ridge_model, subset, type = "link")
head(risk_scores)
plot(final_ridge_model , xvar = "lambda", label = TRUE, main = "Ridge: Coefficient Paths")
######################################################################################
#lasso model
fit_lasso <- glmnet(subset, Surv_obj, family = "cox", alpha = 1)  # Lasso (alpha = 1)
cv_fit <- cv.glmnet(subset, Surv_obj, family = "cox", alpha = 1)
best_lambda <- cv_fit$lambda.min  # Best lambda values
lambda_1se <- cv_fit$lambda.1se   # Simpler model within 1 SE of min
final_lasso_model <- glmnet(subset, Surv_obj, family = "cox", alpha = 1, lambda = best_lambda)
selected_features <- which(coef(final_lasso_model) != 0)
selected_features
risk_scores1 <- predict(final_lasso_model, subset, type = "link")
head(risk_scores1)
plot(final_lasso_model, xvar = "lambda", label = TRUE, main = "Lasso: Coefficient Paths",col="blue")
lasso_coefs <- coef(final_lasso_model)

nonzero_lasso <- lasso_coefs[lasso_coefs != 0]
barplot(nonzero_lasso, las = 2)
###################################################################################
##elastic model
cv_elastic <- cv.glmnet(subset, Surv_obj, family = "cox", alpha = 0.5)
best_lambda <- cv_elastic$lambda.min  # Minimizes cross-validation error
lambda_1se <- cv_elastic$lambda.1se 
elastic_coefs <- coef(cv_elastic, s = "lambda.min")
selected_features <- which(elastic_coefs != 0)
selected_features
final_model <- glmnet(subset, Surv_obj, family = "cox", alpha = 0.5, lambda = best_lambda)
risk_scores2 <- predict(final_model, subset, type = "link")
head(risk_scores2)
#cindex 
concordance_ridge <- survConcordance(Surv_obj ~ risk_scores)$concordance
concordance_lasso <- survConcordance(Surv_obj ~ risk_scores1)$concordance
concordance_elastic <- survConcordance(Surv_obj ~ risk_scores2)$concordance

################################################################################
subset$PFI=data$PFI
subset$time=data$PFI.time
z=Surv(subset$time,subset$PFI)
library(rpart)
fitsurv=rpart(Surv(time,PFI)~.,data=subset,method="exp")
fitsurv
pred1=predict(fitsurv,subset,method="survival")
pred1
concordance_r<- survConcordance(z ~ pred1)$concordance
plot(fitsurv)
text(fitsurv,use.n=TRUE,cex=0.75,all=TRUE)
path.rpart(fitsurv,node=2)
plotcp(fitsurv)
fitsurv$where
km=survfit(z~fitsurv$where)
plot(km,col=c("green","blue","red","purple","brown","black","orange"))
legend("topright",legend=c("node 4","node 10","node 22","node 46","node 47","node 6","node 7"),col=c("green","blue","red","purple","brown","black","orange"),lty=1)
########################################################################
X_selected <- subset[, selected_features, drop = FALSE]  # Keep only selected columns
dim(X_selected)  # Check dimensions of the new dataset
new_data=cbind(data$PFI.time,data$PFI,X_selected)
dim(new_data)
#########################################################################################
cox_model=coxph(Surv(time,PFI)~.,data=subset)
summary(cox_model)
cox_fit=survfit(cox_model)
plot(cox_fit)
pred=predict(cox_model,data=data,type="survival")
concordance_cox <- survConcordance(Surv_obj ~ pred)$concordance
###########################################################################################
install.packages("randomForestSRC")
library(randomForestSRC)
library(survival)
rsf_model <- rfsrc(Surv(time,PFI)~., data =subset, ntree = 200,block.size = 1 ,importance = TRUE)  # Calculate variable importance
pred=predict(rsf_model,data=data,type="survival")
concordance_cox <- survConcordance(z ~ pred$predicted)$concordance
print(rsf_model)
plot(rsf_model)
library(ggRandomForests)
library(ggplot2)
plot(gg_vimp(rsf_model), nvar=50,main = "Variable Importance (RSF)",col="blue")
