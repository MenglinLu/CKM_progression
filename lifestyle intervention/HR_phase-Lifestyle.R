library("survival")
library("survminer")
library(epiDisplay) 
library(ggsurvfit)
library(survival)
library(survminer)
library(KMunicate)
library(stringr)

data_all <- read.csv('result_sustain_lifestyle.csv')
data <- data_all[data_all$Flag_disease==1,]

data$phase <- ifelse(data$ml_stage >= 0 & data$ml_stage <= 8, 1,
                     ifelse(data$ml_stage >= 9 & data$ml_stage <= 12, 2,
                            ifelse(data$ml_stage >= 13 & data$ml_stage <= 20, 3, NA)))

lifestyle_task <- c('lifestyle_class'
)
result_df <- data.frame()
disease_i <- 'cvd'
data_disease <- data[data[str_c('label_',disease_i)]>=0,]
data_disease <- data_disease[data_disease[str_c('survial_day_',disease_i)]>0,]

for(subtype_i in c(1,2,3,4)){
  data_subtype_i <- data_disease[data_disease$ml_subtype==subtype_i-1,]
  
  for (phase_i in c(1,2,3)){
    data_subtype_stage_i <- data_subtype_i[data_subtype_i$phase==phase_i,]
    for(task in lifestyle_task){
      data_subtype_reference_exp <- data_subtype_stage_i[!is.na(data_subtype_stage_i[task]),]
      data_subtype_reference_exp$subgroup <- as.factor(data_subtype_reference_exp[,task])
      
      cox_res <- coxph(Surv(get(str_c('survial_day_',disease_i)), get(str_c('label_',disease_i)))~subgroup+Age+Sex+Creatinine+CysC+eGFR+Phosphate+Urate+Urea+BMI+Waist+Glucose+HbA1c+Albumin+CRP+CHOL+HDLc+LDLc+TG+SBP+DBP,
                       data=data_subtype_reference_exp)
      cox_summary<- summary(cox_res)
      
      HR_stage <- cox_summary$conf.int['subgroupFavorable','exp(coef)']
      HR_stage_upper <- cox_summary$conf.int['subgroupFavorable','upper .95']
      HR_stage_lower <- cox_summary$conf.int['subgroupFavorable','lower .95']
      HR_p <- cox_summary$coefficients['subgroupFavorable','Pr(>|z|)']
      
      result_df <- rbind(result_df, c(str_c('Subtype ',as.character(subtype_i)), str_c('Phase ',as.character(phase_i)), disease_i, str_c(task,'Favorable'), HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
          }
  }
}

