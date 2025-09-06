library(survival)
library(data.table)
library(boot)

# Read in and clean data
zhouyanfeng <- read_sas("C:/Users/Documents/Kailuan/pwv_pro_base.sas7bdat")
zhouyanfeng$t_ckd <- as.numeric(zhouyanfeng$t_ckd)

# define model
formula1 <- Surv(t_ckd, ckd == 1) ~ age+sex+hypertension+diabetes+anemia+BMI+smoke+egfr
formula2 <- Surv(t_ckd, ckd == 1) ~ age+sex+hypertension+diabetes+anemia+BMI+smoke+egfr+pwv5
formula3 <- Surv(t_ckd, ckd == 1) ~  age+sex+hypertension+diabetes+anemia+BMI+smoke+egfr+pwv5+pwv_progression5
formula4 <- Surv(t_ckd, ckd == 1) ~ age+sex+hypertension+diabetes+BMI+smoke+egfr+anemia+alcohol+crp+physical+TG+LDL+HDL
formula5 <- Surv(t_ckd, ckd == 1) ~ age+sex+hypertension+diabetes+BMI+smoke+egfr+anemia+alcohol+crp+physical+TG+LDL+HDL+pwv5
formula6 <- Surv(t_ckd, ckd == 1) ~ age+sex+hypertension+diabetes+BMI+smoke+egfr+anemia+alcohol+crp+physical+TG+LDL+HDL+pwv5+pwv_progression5
# store
formulas <- list(formula1, formula2, formula3, formula4, formula5, formula6)
n_models <- length(formulas)

# fit model and calculate C-index
models <- list()
original_c <- numeric(n_models)
for (i in 1:n_models) {
  models[[i]] <- coxph(formulas[[i]], data = zhouyanfeng)
  c_stat <- concordance(models[[i]])
  original_c[i] <- c_stat$concordance
  cat(paste0("model", i, " original C-index: ", round(original_c[i], 4), "\n"))
}

# define comparision model
comparisons <- list(
  c(1, 2),  # model1 vs model2
  c(1, 3),  # model1 vs model3
  c(4, 5),  # model4 vs model5
  c(4, 6)   # model4 vs model6
)

# calcalate difference in original C-index
original_diff <- list()
for (comp in comparisons) {
  i <- comp[1]
  j <- comp[2]
  diff_val <- original_c[j] - original_c[i]
  original_diff[[paste0("M", i, "_vs_M", j)]] <- diff_val
  cat(paste0("model", j, " vs model", i, " original C-index difference: ", round(diff_val, 4), "\n"))
}
cat("\n")

# define Bootstrap，calculate C-statistic and optimism
calc_boot_stats <- function(data, indices) {
  boot_data <- data[indices, ] 
  
  # Fit all the models on the Bootstrap sample
  boot_models <- list()
  boot_c <- numeric(n_models)
  test_c <- numeric(n_models)
  optimism <- numeric(n_models)
  adj_c <- numeric(n_models)
  
  for (i in 1:n_models) {
    boot_models[[i]] <- tryCatch(
      coxph(formulas[[i]], data = boot_data),
      error = function(e) NULL
    )
    
    if (is.null(boot_models[[i]])) {
      boot_c[i] <- NA
      test_c[i] <- NA
      optimism[i] <- NA
      adj_c[i] <- NA
      next
    }
    
    # Calculate the C-index on the Bootstrap sample
    boot_c[i] <- tryCatch(
      concordance(boot_models[[i]])$concordance,
      error = function(e) NA
    )
    
    # Calculate the performance on the original sample
    lp_orig <- predict(boot_models[[i]], newdata = zhouyanfeng, type = "lp")
    
    test_c[i] <- tryCatch(
      concordance(Surv(t_ckd, ckd) ~ lp_orig, 
                  data = zhouyanfeng, reverse = TRUE)$concordance,
      error = function(e) NA
    )
    
    # calculate optimism
    optimism[i] <- boot_c[i] - test_c[i]
    
    # Calculate the optimism-adjusted C-index
    adj_c[i] <- original_c[i] - optimism[i]
  }
  
  # calculate the difference
  comp_diffs <- numeric(length(comparisons))
  comp_adj_diffs <- numeric(length(comparisons))
  
  for (k in 1:length(comparisons)) {
    comp <- comparisons[[k]]
    i <- comp[1]
    j <- comp[2]
    
    comp_diffs[k] <- boot_c[j] - boot_c[i]
    comp_adj_diffs[k] <- adj_c[j] - adj_c[i]
  }
  
  # return
  result <- c(
    boot_c,        # 1-6: Bootstrap sample C-index of each model
    test_c,        # 7-12: Test C-index for each model
    optimism,      # 13-18: optimism for each model
    adj_c,         # 19-24: optimism-adjusted C-index
    comp_diffs,    # 25-28: the original C-index difference
    comp_adj_diffs # 29-32: adjusted C-index difference
  )
  
  return(result)
}

# run Bootstrap
set.seed(123)
boot_results <- boot(data = zhouyanfeng,
                     statistic = calc_boot_stats,
                     R = 1000) 

# extract Bootstrap results
n_stats <- ncol(boot_results$t)
boot_c <- boot_results$t[, 1:n_models]
test_c <- boot_results$t[, (n_models+1):(2*n_models)]
optimism <- boot_results$t[, (2*n_models+1):(3*n_models)]
adj_c <- boot_results$t[, (3*n_models+1):(4*n_models)]
comp_diffs <- boot_results$t[, (4*n_models+1):(4*n_models+length(comparisons))]
comp_adj_diffs <- boot_results$t[, (4*n_models+length(comparisons)+1):n_stats]

# Calculate the average optimistic bias and the adjusted C-index
mean_optimism <- colMeans(optimism, na.rm = TRUE)
adjusted_c <- original_c - mean_optimism

# calculate estimate 
adjusted_diff <- list()
for (k in 1:length(comparisons)) {
  comp <- comparisons[[k]]
  i <- comp[1]
  j <- comp[2]
  diff_val <- adjusted_c[j] - adjusted_c[i]
  adjusted_diff[[paste0("M", i, "_vs_M", j)]] <- diff_val
}

# output results
cat("=== Model performance ===\n\n")

for (i in 1:n_models) {
  cat(paste0("model", i, " original C-index: ", round(original_c[i], 4), "\n"))
}
cat("\n")

for (k in 1:length(comparisons)) {
  comp <- comparisons[[k]]
  i <- comp[1]
  j <- comp[2]
  cat(paste0("model", j, " vs model", i, " original C-index difference: ", 
             round(original_diff[[paste0("M", i, "_vs_M", j)]], 4), "\n"))
}
cat("\n")

for (i in 1:n_models) {
  cat(paste0("model", i, " mean_optimism: ", round(mean_optimism[i], 4), "\n"))
}
cat("\n")

for (i in 1:n_models) {
  cat(paste0("model", i, " optimism-adjusted C-index: ", round(adjusted_c[i], 4), "\n"))
}
cat("\n")

for (k in 1:length(comparisons)) {
  comp <- comparisons[[k]]
  i <- comp[1]
  j <- comp[2]
  cat(paste0("model", j, " vs model", i, " optimism-adjusted C-index difference: ", 
             round(adjusted_diff[[paste0("M", i, "_vs_M", j)]], 4), "\n"))
}
cat("\n")

# calculate boot.ci
cat("=== 95% CI ===\n\n")

# 95% CI for the original C-index
for (i in 1:n_models) {
  ci <- boot.ci(boot_results, conf = 0.95, type = "basic", index = i)$basic[4:5]
  cat(paste0("model", i, " Bootstrap 95% CI for original C-index: [", 
             round(ci[1], 4), ", ", round(ci[2], 4), "]\n"))
}
cat("\n")

# 95% CI for the optimism
for (i in 1:n_models) {
  ci <- boot.ci(boot_results, conf = 0.95, type = "basic", index = 2*n_models + i)$basic[4:5]
  cat(paste0("model", i, " Bootstrap 95% CI for optimism: [", 
             round(ci[1], 4), ", ", round(ci[2], 4), "]\n"))
}
cat("\n")

# 95% CI for optimism-adjusted C-index
for (i in 1:n_models) {
  ci <- boot.ci(boot_results, conf = 0.95, type = "basic", index = 3*n_models + i)$basic[4:5]
  cat(paste0("model", i, " Bootstrap 95%CI for optimism-adjusted C-index: [", 
             round(ci[1], 4), ", ", round(ci[2], 4), "]\n"))
}
cat("\n")

# 95% CI for C-index change
for (k in 1:length(comparisons)) {
  comp <- comparisons[[k]]
  i <- comp[1]
  j <- comp[2]
  ci <- boot.ci(boot_results, conf = 0.95, type = "basic", index = 4*n_models + k)$basic[4:5]
  cat(paste0("model", j, " vs model", i, " Bootstrap 95% CI for original C-index: [", 
             round(ci[1], 4), ", ", round(ci[2], 4), "]\n"))
}
cat("\n")

# 95% CI for optimism-adjusted C-index change
for (k in 1:length(comparisons)) {
  comp <- comparisons[[k]]
  i <- comp[1]
  j <- comp[2]
  ci <- boot.ci(boot_results, conf = 0.95, type = "basic", index = 4*n_models + length(comparisons) + k)$basic[4:5]
  cat(paste0("model", j, " vs model", i, " Bootstrap 95% CI for optimism-adjusted C-index: [", 
             round(ci[1], 4), ", ", round(ci[2], 4), "]\n"))
}
cat("\n")



#calculate IDI

#model 1 and model 2
#create model
a <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd)))
#model 1
indata0 <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd, age, sex, hypertension, diabetes, anemia, BMI, smoke, egfr)))
#model 2
Indata1 <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd, age, sex, hypertension, diabetes, anemia, BMI, smoke, egfr,  pwv5)))

#The predictor variable matrix in model
covs0 <- indata0[, -c(1,2)]
covs1 <- indata1[, -c(1,2)]

library(survIDINRI)
library(survC1)
#Set a random seed
set.seed(123)

#calculate IDI and NRI at year 5
x_5 <- IDI.INF(a, covs0, covs1, t0=5, npert = 200)
result_5year <- IDI.INF.OUT(x_5)


#model 1 and model3
a <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd)))
indata0 <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd, age, sex, hypertension, diabetes, anemia, BMI, smoke, egfr)))
Indata1 <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd, age, sex, hypertension, diabetes, anemia, BMI, smoke, egfr,  pwv5, pwv_progression5)))

covs0 <- indata0[, -c(1,2)]
covs1 <- indata1[, -c(1,2)]

library(survIDINRI)
library(survC1)
set.seed(123)

#calculate IDI and NRI at year 5
x_5 <- IDI.INF(a, covs0, covs1, t0=5, npert = 200)
result_5year <- IDI.INF.OUT(x_5)


#model 4 and model5
a <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd)))
indata0 <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd, age, sex, hypertension, diabetes, BMI, smoke, egfr, anemia, alcohol, crp, physical,  TG, LDL, HDL)))
Indata1 <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd, age, sex, hypertension, diabetes, BMI, smoke, egfr, anemia, alcohol, crp, physical,  TG, LDL, HDL, pwv5)))

covs0 <- indata0[, -c(1,2)]
covs1 <- indata1[, -c(1,2)]

library(survIDINRI)
library(survC1)
set.seed(123)

#calculate IDI and NRI at year 5
x_5 <- IDI.INF(a, covs0, covs1, t0=5, npert = 200)
result_5year <- IDI.INF.OUT(x_5)


#model 4 and model 6
a <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd)))
indata0 <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd, age, sex, hypertension, diabetes, BMI, smoke, egfr, anemia, alcohol, crp, physical,  TG, LDL, HDL)))
Indata1 <- with(zhouyanfeng, as.matrix(cbind(t_ckd, ckd, age, sex, hypertension, diabetes, BMI, smoke, egfr, anemia, alcohol, crp, physical,  TG, LDL, HDL, pwv5, pwv_progression5)))

covs0 <- indata0[, -c(1,2)]
covs1 <- indata1[, -c(1,2)]

library(survIDINRI)
library(survC1)
set.seed(123)

#calculate IDI and NRI at year 5
x_5 <- IDI.INF(a, covs0, covs1, t0=5, npert = 200)
result_5year <- IDI.INF.OUT(x_5)


#proportional hazards test
library(survival)
library(haven)
library(survminer)

# Read in and clean data
zhouyanfeng <- read_sas("C:/Users/Documents/Kailuan/pwv_pro.sas7bdat")

#model 1
y <- with(zhouyanfeng, Surv(t_ckd,ckd==1))
mod <- coxph(y ~ cluster(ID)+pwv_progression5,data=zhouyanfeng)
summary(mod)
phtest <- cox.zph(mod)
phtest
#phtest <- coxp.zph(mod, transform=rank)
#ggcoxzph(phtest)
Plot(phtest, se=T, var="pwv_progression5")

#model 3
y <- with(zhouyanfeng, Surv(t_ckd,ckd==1))
mod <- coxph(y ~  cluster(ID)+age+ sex+ hypertension +diabetes+BMI +smoke +egfr + anemia + alcohol + crp +physical +TG+LDL+HDL+pwv5+pwv_progression5, data=zhouyanfeng)
summary(mod)
phtest <- cox.zph(mod)
phtest
#phtest <- coxp.zph(mod, transform=rank)
#ggcoxzph(phtest)
Plot(phtest, se=T, var="pwv_progression5")

