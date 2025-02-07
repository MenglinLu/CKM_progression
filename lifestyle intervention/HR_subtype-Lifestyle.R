library("survival")
library("survminer")
library(epiDisplay) 
library(ggsurvfit)
library(survival)
library(survminer)
library(KMunicate)
library(ComparisonSurv)
library(stringr)

data_all <- read.csv('result_sustain_lifestyle.csv')
data <- data_all[data_all$Flag_disease==1,]

lifestyle_task <- c('smoking_healthy',  'alcohol_healthy', 'diet_healthy',
                    'activity_healthy','sedentary_healthy',
                    'sleep_healthy','social_healthy'
)
result_df <- data.frame()
for(subtype_i in c(1,2,3,4)){
  data_subtype_i <- data[data$ml_subtype==subtype_i-1,]

  for(task in lifestyle_task){
    data_subtype_reference_exp <- data_subtype_i[!is.na(data_subtype_i[task]),]
    data_subtype_reference_exp$subgroup <- as.factor(data_subtype_reference_exp[,task])
    if(task == 'smoking_healthy'){
    cox_res <- coxph(Surv(get(str_c('survial_day_cvd')), get(str_c('label_cvd')))~subgroup+Age+Sex+Creatinine+CysC+eGFR+Phosphate+Urate+Urea+BMI+Waist+Glucose+HbA1c+Albumin+CRP+CHOL+HDLc+LDLc+TG+SBP+DBP,
                     data=data_subtype_reference_exp)}
    if(task == 'alcohol_healthy'){
      cox_res <- coxph(Surv(get(str_c('survial_day_cvd')), get(str_c('label_cvd')))~subgroup+Age+Sex+Creatinine+CysC+eGFR+Phosphate+Urate+Urea+BMI+Waist+Glucose+HbA1c+Albumin+CRP+CHOL+HDLc+LDLc+TG+SBP+DBP,
                       data=data_subtype_reference_exp)}
    if(task == 'diet_healthy'){
      cox_res <- coxph(Surv(get(str_c('survial_day_cvd')), get(str_c('label_cvd')))~subgroup+Age+Sex+Creatinine+CysC+eGFR+Phosphate+Urate+Urea+BMI+Waist+Glucose+HbA1c+Albumin+CRP+CHOL+HDLc+LDLc+TG+SBP+DBP,
                       data=data_subtype_reference_exp)}
    if(task == 'activity_healthy'){
      cox_res <- coxph(Surv(get(str_c('survial_day_cvd')), get(str_c('label_cvd')))~subgroup+Age+Sex+Creatinine+CysC+eGFR+Phosphate+Urate+Urea+BMI+Waist+Glucose+HbA1c+Albumin+CRP+CHOL+HDLc+LDLc+TG+SBP+DBP,
                       data=data_subtype_reference_exp)}
    if(task == 'sedentary_healthy'){
      cox_res <- coxph(Surv(get(str_c('survial_day_cvd')), get(str_c('label_cvd')))~subgroup+Age+Sex+Creatinine+CysC+eGFR+Phosphate+Urate+Urea+BMI+Waist+Glucose+HbA1c+Albumin+CRP+CHOL+HDLc+LDLc+TG+SBP+DBP,
                       data=data_subtype_reference_exp)}
    if(task == 'sleep_healthy'){
      cox_res <- coxph(Surv(get(str_c('survial_day_cvd')), get(str_c('label_cvd')))~subgroup+Age+Sex+Creatinine+CysC+eGFR+Phosphate+Urate+Urea+BMI+Waist+Glucose+HbA1c+Albumin+CRP+CHOL+HDLc+LDLc+TG+SBP+DBP,
                       data=data_subtype_reference_exp)}
    if(task == 'social_healthy'){
      cox_res <- coxph(Surv(get(str_c('survial_day_cvd')), get(str_c('label_cvd')))~subgroup+Age+Sex+Creatinine+CysC+eGFR+Phosphate+Urate+Urea+BMI+Waist+Glucose+HbA1c+Albumin+CRP+CHOL+HDLc+LDLc+TG+SBP+DBP,
                       data=data_subtype_reference_exp)}
    
    cox_summary<- summary(cox_res)
    
    HR_stage <- cox_summary$conf.int['subgroup1','exp(coef)']
    HR_stage_upper <- cox_summary$conf.int['subgroup1','upper .95']
    HR_stage_lower <- cox_summary$conf.int['subgroup1','lower .95']
    HR_p <- cox_summary$coefficients['subgroup1','Pr(>|z|)']
    
    result_df <- rbind(result_df, c(str_c('Subtype ',as.character(subtype_i)), task, HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  }
}

for(subtype_i in c(1,2,3,4)){
  data_subtype_i <- data[data$ml_subtype==subtype_i-1,]
  
  for(task in c('lifestyle_class')){
    data_subtype_reference_exp <- data_subtype_i[!is.na(data_subtype_i[task]),]
    data_subtype_reference_exp$subgroup <- as.factor(data_subtype_reference_exp[,task])
    
    cox_res <- coxph(Surv(get(str_c('survial_day_cvd')), get(str_c('label_cvd')))~subgroup+Age+Sex+Creatinine+CysC+eGFR+Phosphate+Urate+Urea+BMI+Waist+Glucose+HbA1c+Albumin+CRP+CHOL+HDLc+LDLc+TG+SBP+DBP,
                     data=data_subtype_reference_exp)
    cox_summary<- summary(cox_res)
    
    HR_stage <- cox_summary$conf.int['subgroupFavorable','exp(coef)']
    HR_stage_upper <- cox_summary$conf.int['subgroupFavorable','upper .95']
    HR_stage_lower <- cox_summary$conf.int['subgroupFavorable','lower .95']
    HR_p <- cox_summary$coefficients['subgroupFavorable','Pr(>|z|)']
    
    result_df <- rbind(result_df, c(str_c('Subtype ',as.character(subtype_i)), str_c(task,'Favorable'), HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  }
  print(str_c('End task ',task))
}


for(subtype_i in c(1,2,3,4)){
  data_subtype_i <- data[data$ml_subtype==subtype_i-1,]
  
  for(task in c('lifestyle_class')){
    data_subtype_reference_exp <- data_subtype_i[!is.na(data_subtype_i[task]),]
    data_subtype_reference_exp$subgroup <- as.factor(data_subtype_reference_exp[,task])
    
    cox_res <- coxph(Surv(get(str_c('survial_day_cvd')), get(str_c('label_cvd')))~subgroup+Age+Sex+Creatinine+CysC+eGFR+Phosphate+Urate+Urea+BMI+Waist+Glucose+HbA1c+Albumin+CRP+CHOL+HDLc+LDLc+TG+SBP+DBP,
                     data=data_subtype_reference_exp)
    cox_summary<- summary(cox_res)
    
    HR_stage <- cox_summary$conf.int['subgroupIntermediate','exp(coef)']
    HR_stage_upper <- cox_summary$conf.int['subgroupIntermediate','upper .95']
    HR_stage_lower <- cox_summary$conf.int['subgroupIntermediate','lower .95']
    HR_p <- cox_summary$coefficients['subgroupIntermediate','Pr(>|z|)']
    
    result_df <- rbind(result_df, c(str_c('Subtype ',as.character(subtype_i)), str_c(task,'Favorable'), HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  }
  print(str_c('End task ',task))
}

colnames(result_df) <- c('Subtype','Lifestyle','HR','HR_lower','HR_upper','pval')

