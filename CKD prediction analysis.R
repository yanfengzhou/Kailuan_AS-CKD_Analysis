library(survival)
library(haven)

# Read in and clean data
zhouyanfeng <- read_sas("C:/Users/Documents/Kailuan/pwv_pro.sas7bdat")
zhouyanfeng$t_ckd <- as.numeric(zhouyanfeng$t_ckd)

# calculate follow-up year
mean(zhouyanfeng$t_ckd/365.25, na.rm = T)
sum(is.na(zhouyanfeng$t_ckd))

### C-index ###
# Model 1
s1 <- coxph(Surv(t_ckd, ckd == 1) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + cluster(ID), 
            data = zhouyanfeng)
summary(s1)
0.650  # C-index
0.650 - 1.96 * 0.01  # 95% CI
0.650 + 1.96 * 0.01  # 95% CI


# Model 2
s2 <- coxph(Surv(t_ckd, ckd == 1) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + pwv5 + cluster(ID), 
            data = zhouyanfeng)
summary(s2)
0.674  # C-index
0.674 - 1.96 * 0.01  # 95% CI
0.674 + 1.96 * 0.01  # 95% CI


# Model 3
s3 <- coxph(Surv(t_ckd, ckd == 1) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + pwv5 +
               pwv_progression5 + cluster(ID), 
            data = zhouyanfeng)
summary(s3)
0.690  # C-index
0.690 - 1.96 * 0.01  # 95% CI
0.690 + 1.96 * 0.01  # 95% CI

#Model 4
s4 <- coxph(Surv(t_ckd, ckd == 1) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + 
              alcohol +  physical + tg + ldl + hdl + crp + cluster(ID), 
            data = zhouyanfeng)
summary(s4)
0.674  # C-index
0.674 - 1.96 * 0.01  # 95% CI
0.674 + 1.96 * 0.01  # 95% CI


# Model 5
s5 <- coxph(Surv(t_ckd, ckd == 1) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + 
              alcohol +  physical + tg + ldl + hdl + crp + pwv5 + cluster(ID), 
            data = zhouyanfeng)
summary(s5)
0.691  # C-index
0.691 - 1.96 * 0.01  # 95% CI
0.691 + 1.96 * 0.01  # 95% CI


# Model 6
s6 <- coxph(Surv(t_ckd, ckd == 1) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + 
              alcohol +  physical + tg + ldl + hdl + crp + pwv5 + pwv_progression5 + cluster(ID), 
            data = zhouyanfeng)
summary(s6)
0.705  # C-index
0.705 - 1.96 * 0.01  # 95% CI
0.705 + 1.96 * 0.01  # 95% CI



### C-index change ###
library(survival)
library(boot)


# C-index change between model 1 and model 2
calc_cindex_diff <- function(data, indices) {
  d <- data[indices, ]
  fit1 <- coxph(Surv(t_ckd, ckd) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + cluster(ID), data = d)
  fit2 <- coxph(Surv(t_ckd, ckd) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + pwv5 + cluster(ID), data = d)
  
  c1 <- summary(fit1)$concordance[1]
  c2 <- summary(fit2)$concordance[1]
  return(c2 - c1)
}

# Bootstrap（performed 1000 times）
set.seed(123)
boot_results <- boot(zhouyanfeng, calc_cindex_diff, R = 1000)

# calculation
delta_c <- boot_results$t0
ci_percent <- boot.ci(boot_results, type = "perc")
p_value <- mean(boot_results$t < 0) * 2  # two-sided test

# output results
cat("C-index change (new - old):", round(delta_c, 3), "\n")
cat("95% CI: [", round(ci_percent$percent[4], 3), ",", round(ci_percent$percent[5], 3), "]\n")
cat("p-value:", format.pval(p_value, digits = 4), "\n")


# C-index change between model 1 and model 3
calc_cindex_diff 2<- function(data, indices) {
  d <- data[indices, ]
  fit1 <- coxph(Surv(t_ckd, ckd) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + cluster(ID), data = d)
  fit3 <- coxph(Surv(t_ckd, ckd) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + pwv5 +
                  pwv_progression5 + cluster(ID), data = d)
  
  c1 <- summary(fit1)$concordance[1]
  c3 <- summary(fit3)$concordance[1]
  return(c3 - c1)
}

# Bootstrap（performed 1000 times）
set.seed(123)
boot_results2 <- boot(zhouyanfeng, calc_cindex_diff2, R = 1000)

# calculation
delta_c2 <- boot_results2$t0
ci_percent2 <- boot.ci(boot_results2, type = "perc")
p_value2 <- mean(boot_results2$t < 0) * 2  # two-sided test

# output results
cat("C-index change (new - old):", round(delta_c2, 3), "\n")
cat("95% CI: [", round(ci_percent2$percent[4], 3), ",", round(ci_percent2$percent[5], 3), "]\n")
cat("p-value:", format.pval(p_value2, digits = 4), "\n")



# C-index change between model 4 and model 5
calc_cindex_diff3 <- function(data, indices) {
  d <- data[indices, ]
  fit4 <- coxph(Surv(t_ckd, ckd) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + 
                  alcohol +  physical + tg + ldl + hdl + crp + cluster(ID), data = d)
  fit5 <- coxph(Surv(t_ckd, ckd) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + 
                  alcohol +  physical + tg + ldl + hdl + crp + pwv5 + cluster(ID), data = d)
  
  c4 <- summary(fit4)$concordance[1]
  c5 <- summary(fit5)$concordance[1]
  return(c5 - c4)
}

# Bootstrap（performed 1000 times）
set.seed(123)
boot_results3 <- boot(zhouyanfeng, calc_cindex_diff3, R = 1000)

# calculation
delta_c3 <- boot_results3$t0
ci_percent3 <- boot.ci(boot_results3, type = "perc")
p_value3 <- mean(boot_results3$t < 0) * 2  # two-sided test

# output results
cat("C-index change (new - old):", round(delta_c3, 3), "\n")
cat("95% CI: [", round(ci_percent3$percent[4], 3), ",", round(ci_percent3$percent[5], 3), "]\n")
cat("p-value:", format.pval(p_value3, digits = 4), "\n")


# C-index change between model 4 and model 6
calc_cindex_diff4 <- function(data, indices) {
  d <- data[indices, ]
  fit4 <- coxph(Surv(t_ckd, ckd) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + 
                  alcohol +  physical + tg + ldl + hdl + crp + cluster(ID), data = d)
  fit6 <- coxph(Surv(t_ckd, ckd) ~ age + sex + hypertension + diabetes + anemia + BMI + smoke + egfr + 
                  alcohol +  physical + tg + ldl + hdl + crp + pwv5 + pwv_progression5 + cluster(ID), data = d)
  
  c4 <- summary(fit4)$concordance[1]
  c5 <- summary(fit6)$concordance[1]
  return(c6 - c4)
}

# Bootstrap（performed 1000 times）
set.seed(123)
boot_results4 <- boot(zhouyanfeng, calc_cindex_diff4, R = 1000)

# calculation
delta_c4 <- boot_results4$t0
ci_percent4 <- boot.ci(boot_results4, type = "perc")
p_value4 <- mean(boot_results4$t < 0) * 2  # two-sided test

# output results
cat("C-index change (new - old):", round(delta_c4, 3), "\n")
cat("95% CI: [", round(ci_percent4$percent[4], 3), ",", round(ci_percent4$percent[5], 3), "]\n")
cat("p-value:", format.pval(p_value4, digits = 4), "\n")