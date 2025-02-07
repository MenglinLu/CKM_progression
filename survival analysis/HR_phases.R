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

result_df <- data.frame()

for(subtype_i in c(0,1,2,3)){
  data_subtype5 <- data[data$ml_subtype==subtype_i,]
  data_subtype5_stage1 <- data_subtype5[data_subtype5$ml_stage>=0 & data_subtype5$ml_stage<=8,]
  data_subtype5_stage2 <- data_subtype5[data_subtype5$ml_stage>=9 & data_subtype5$ml_stage<=12,]
  data_subtype5_stage3 <- data_subtype5[data_subtype5$ml_stage>=13 & data_subtype5$ml_stage<=20,]
  
  data_subtype5_stage1$ml_stage_new <- 1
  data_subtype5_stage2$ml_stage_new <- 2
  data_subtype5_stage3$ml_stage_new <- 3
  
  data_subtype5_new <- rbind(data_subtype5_stage1,
                             data_subtype5_stage2,data_subtype5_stage3)

  title_dict <- c(
    'hf'='Heart failure',
    'af'='Atrial fibrillation',
    'chd'='Coronary Heart Disease',
    'cvd'='Cardiovascular disease',
    'stroke_haemorr'='Cerebral hemorrhage',
    'stroke_infarc'='Cerebral infarction',
    'pad'='Peripheral artery disease',
    'death'='Mortality')
  
  for(task in c('hf','af','chd',
                'stroke_haemorr','stroke_infarc','pad','death')){
    data_temp <- data_subtype5_new
    
    data_temp_task <- data_temp[data_temp[str_c('label_',task)]>=0,]
    data_temp_task <- data_temp_task[data_temp_task[str_c('survial_day_',task)]>0,]
    
    data_temp_task_control <- data_control[data_control[str_c('label_',task)]>=0,]
    data_temp_task_control <- data_temp_task_control[data_temp_task_control[str_c('survial_day_',task)]>0,]
    
    data_temp_task_reference <- data_temp_task_control
    data_temp_task_stage1 <- data_temp_task[data_temp_task$ml_stage_new==1,]
    data_temp_task_stage2 <- data_temp_task[data_temp_task$ml_stage_new==2,]
    data_temp_task_stage3 <- data_temp_task[data_temp_task$ml_stage_new==3,]
    
    for(stage_i in c(1,2,3)){
      data_temp_task_reference$stage_group <- 0
      data_temp_task_exp <- get(str_c('data_temp_task_stage',as.character(stage_i)))
      data_temp_task_exp$stage_group <- 1
      data_temp_task_reference_exp <- rbind(data_temp_task_reference[,c('Age','Sex','Smoker1','stage_group',str_c('survial_day_',task),str_c('label_',task))],data_temp_task_exp[,c('Age','Sex','Smoker1','stage_group',str_c('survial_day_',task),str_c('label_',task))])
      data_temp_task_reference_exp$stage_group <- as.factor(data_temp_task_reference_exp$stage_group)
      
      cox_res <- coxph(Surv(get(str_c('survial_day_',task)), get(str_c('label_',task)))~stage_group,data=data_temp_task_reference_exp)
      cox_summary<- summary(cox_res)
      
      HR_stage <- cox_summary$conf.int['stage_group1','exp(coef)']
      HR_stage_upper <- cox_summary$conf.int['stage_group1','upper .95']
      HR_stage_lower <- cox_summary$conf.int['stage_group1','lower .95']
      HR_p <- cox_summary$coefficients['stage_group1','Pr(>|z|)']
      
      result_df <- rbind(result_df, c(str_c('Subtype ',as.character(subtype_i)), title_dict[task], stage_i, HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
    }
    print(str_c('End task ',task))
  }}

colnames(result_df) <- c('Subtype','Disease','Stage','HR','HR_lower','HR_upper','pval')