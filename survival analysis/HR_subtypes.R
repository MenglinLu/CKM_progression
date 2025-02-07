library("survival")
library("survminer")
library(epiDisplay) 
library(ggsurvfit)
library(survival)
library(survminer)
library(KMunicate)
library(ComparisonSurv)
library(stringr)

data_all <- read.csv('./result_sustain.csv')

data <- data_all[data_all$Flag_disease==1,]

data_control <- data_all[data_all$Flag_disease==0,]

data_surv <- read.csv('data_disease_surv.csv')

data <- merge(data,data_surv,by='Eid')
data_control <- merge(data_control,data_surv,by='Eid')

title_dict <- c(
                'hf'='Heart failure',
                'af'='Atrial fibrillation',
                'chd'='Coronary Heart Disease',
                'cvd'='Cardiovascular disease',
                'stroke_haemorr'='Cerebral hemorrhage',
                'stroke_infarc'='Cerebral infarction',
                'pad'='Peripheral artery disease',
                'death'='Mortality')

result_df <- data.frame()

for(task in c('hf','af','chd','pad','stroke_haemorr','stroke_infarc',
              'cvd','death')){
  
  data_temp_task <- data[data[str_c('label_',task)]>=0,]
  data_temp_task <- data_temp_task[data_temp_task[str_c('survial_day_',task)]>0,]
  
  data_temp_task_control <- data_control[data_control[str_c('label_',task)]>=0,]
  data_temp_task_control <- data_temp_task_control[data_temp_task_control[str_c('survial_day_',task)]>0,]
  
  data_subtype1_task <- data_temp_task[data_temp_task$ml_subtype==0,]
  data_subtype2_task <- data_temp_task[data_temp_task$ml_subtype==1,]
  data_subtype3_task <- data_temp_task[data_temp_task$ml_subtype==2,]
  data_subtype4_task <- data_temp_task[data_temp_task$ml_subtype==3,]
  
  data_subtype_reference <- data_temp_task_control
  
  for(subtype_i in c(1,2,3,4)){
    data_subtype_exp <- get(str_c('data_subtype',as.character(subtype_i),'_task'))
    data_subtype_reference$subtype_group <- 0
    data_subtype_exp$subtype_group <- 1
    
    data_subtype_reference_exp <- rbind(data_subtype_reference,data_subtype_exp)
    data_subtype_reference_exp$subtype_group <- as.factor(data_subtype_reference_exp$subtype_group)
    
    cox_res <- coxph(Surv(get(str_c('survial_day_',task)), get(str_c('label_',task)))~subtype_group,data=data_subtype_reference_exp)
    cox_summary<- summary(cox_res)
    
    HR_stage <- cox_summary$conf.int['subtype_group1','exp(coef)']
    HR_stage_upper <- cox_summary$conf.int['subtype_group1','upper .95']
    HR_stage_lower <- cox_summary$conf.int['subtype_group1','lower .95']
    HR_p <- cox_summary$coefficients['subtype_group1','Pr(>|z|)']
    
    result_df <- rbind(result_df, c(title_dict[task], str_c('Subtype ',as.character(subtype_i)), HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  }
  print(str_c('End task ',task))
  }
  
colnames(result_df) <- c('Disease','Subtype','HR','HR_lower','HR_upper','pval')