library(epiDisplay) 
library(ggsurvfit)
library(survival)
library(survminer)
library(KMunicate)
library(ComparisonSurv)
library(stringr)
library(condSURV)
library(RColorBrewer)
library(tibble)
library(ggpp)

options(stringsAsFactors = F)

data <- read.csv('data_disease_surv.csv')

result_k41 <- read.csv('result_sustain.csv'))

result_k4 <- result_k41[result_k41$Flag_disease==1,]

data_merge <- merge(result_k4, data, by='Eid')

data_merge[data_merge$ml_subtype==0,'ml_subtype'] = 'Subtype 1'
data_merge[data_merge$ml_subtype==1,'ml_subtype'] = 'Subtype 2'
data_merge[data_merge$ml_subtype==2,'ml_subtype'] = 'Subtype 3'
data_merge[data_merge$ml_subtype==3,'ml_subtype'] = 'Subtype 4'

title_dict <- c(
                'hf'='Heart failure',
                'af'='Atrial fibrillation',
                'chd'='Coronary Heart Disease',
                'cvd'='Cardiovascular disease',
                'stroke_haemorr'='Cerebral hemorrhage',
                'stroke_infarc'='Cerebral infarction',
                'pad'='Peripheral artery disease',
                'death'='Mortality')

task <- 'death'
data_temp_task <- data_merge[data_merge[str_c('label_',task)]>=0,]
data_temp_task <- data_temp_task[data_temp_task[str_c('survial_day_',task)]>0,]

data_subtype1_task <- data_temp_task[data_temp_task$ml_subtype=='Subtype 1',]
data_subtype2_task <- data_temp_task[data_temp_task$ml_subtype=='Subtype 2',]
data_subtype3_task <- data_temp_task[data_temp_task$ml_subtype=='Subtype 3',]
data_subtype4_task <- data_temp_task[data_temp_task$ml_subtype=='Subtype 4',]

res_subtype <- data.frame()
data_subtype_reference <- data_subtype1_task 
for(subtype_i in c(2,3,4)){
  data_subtype_exp <- get(str_c('data_subtype',as.character(subtype_i),'_task'))
  data_subtype_reference$subtype_group <- 0
  data_subtype_exp$subtype_group <- 1
  
  data_subtype_reference_exp <- rbind(data_subtype_reference,data_subtype_exp)
  data_subtype_reference_exp$subtype_group <- as.factor(data_subtype_reference_exp$subtype_group)
  
  cox_res <- coxph(Surv(get(str_c('survial_day_',task)), get(str_c('label_',task)))~subtype_group,data=data_subtype_reference_exp)
  cox_summary<- summary(cox_res)
  
  HR_stage <- round(cox_summary$conf.int['subtype_group1','exp(coef)'],2)
  HR_stage_upper <- round(cox_summary$conf.int['subtype_group1','upper .95'],2)
  HR_stage_lower <- round(cox_summary$conf.int['subtype_group1','lower .95'],2)
  HR_p <- cox_summary$coefficients['subtype_group1','Pr(>|z|)']
  HR_show <- paste0('HR = ',HR_stage, '(',HR_stage_lower,'-',HR_stage_upper,')')
  
  res_subtype <- rbind(res_subtype, c(1, subtype_i, HR_show, HR_p))
}

data_subtype_reference <- data_subtype2_task 
for(subtype_i in c(3,4)){
  data_subtype_exp <- get(str_c('data_subtype',as.character(subtype_i),'_task'))
  data_subtype_reference$subtype_group <- 0
  data_subtype_exp$subtype_group <- 1
  
  data_subtype_reference_exp <- rbind(data_subtype_reference,data_subtype_exp)
  data_subtype_reference_exp$subtype_group <- as.factor(data_subtype_reference_exp$subtype_group)
  
  cox_res <- coxph(Surv(get(str_c('survial_day_',task)), get(str_c('label_',task)))~subtype_group,data=data_subtype_reference_exp)
  cox_summary<- summary(cox_res)
  
  HR_stage <- round(cox_summary$conf.int['subtype_group1','exp(coef)'],2)
  HR_stage_upper <- round(cox_summary$conf.int['subtype_group1','upper .95'],2)
  HR_stage_lower <- round(cox_summary$conf.int['subtype_group1','lower .95'],2)
  HR_p <- cox_summary$coefficients['subtype_group1','Pr(>|z|)']
  HR_show <- paste0('HR = ',HR_stage, '(',HR_stage_lower,'-',HR_stage_upper,')')
  
  res_subtype <- rbind(res_subtype, c(2, subtype_i, HR_show, HR_p))
}

data_subtype_reference <- data_subtype3_task 
for(subtype_i in c(4)){
  data_subtype_exp <- get(str_c('data_subtype',as.character(subtype_i),'_task'))
  data_subtype_reference$subtype_group <- 0
  data_subtype_exp$subtype_group <- 1
  
  data_subtype_reference_exp <- rbind(data_subtype_reference,data_subtype_exp)
  data_subtype_reference_exp$subtype_group <- as.factor(data_subtype_reference_exp$subtype_group)
  
  cox_res <- coxph(Surv(get(str_c('survial_day_',task)), get(str_c('label_',task)))~subtype_group,data=data_subtype_reference_exp)
  cox_summary<- summary(cox_res)
  
  HR_stage <- round(cox_summary$conf.int['subtype_group1','exp(coef)'],2)
  HR_stage_upper <- round(cox_summary$conf.int['subtype_group1','upper .95'],2)
  HR_stage_lower <- round(cox_summary$conf.int['subtype_group1','lower .95'],2)
  HR_p <- cox_summary$coefficients['subtype_group1','Pr(>|z|)']
  HR_show <- paste0('HR = ',HR_stage, '(',HR_stage_lower,'-',HR_stage_upper,')')
  
  res_subtype <- rbind(res_subtype, c(3, subtype_i, HR_show, HR_p))
}
colnames(res_subtype) <- c('Subtype_x','Subtype_y','HR','P')

res_subtype$P <- p.adjust(res_subtype$P, method = "fdr")

###survival curves
fit <- survfit(Surv(get(str_c('survial_day_',task)), get(str_c('label_',task))) ~ ml_subtype, data = data_temp_task)
pp = surv_pvalue(fit)$pval
ppval1 = paste0("Log-rank test P", ifelse(pp < 0.001, " < 0.001", # 若P值<0.001则标记为“<0.001”
                           paste0(" = ",round(pp, 3))))

col <- c('#f4cddd','#bfd4d7','#ced9f5','#ebdaca')
p <- survfit2(Surv(get(str_c('survial_day_',task)), get(str_c('label_',task))) ~ ml_subtype, data = data_temp_task) %>% 
  ggsurvfit(linewidth = 0.8, type = "survival") +
  add_confidence_interval() + # 添加置信区间
  theme_classic(base_size = 12) + theme_bw()+
  theme(legend.position = 'bottom',text=element_text(size=9,family="Arial"),
        plot.title = element_text(hjust = 0.5,size = 10.5,family="Arial"),
        panel.grid = element_blank())+
  add_risktable(size = 2.8, risktable_height = 0.3,
                theme=theme_test()+
                  theme(plot.title = element_text(color="black",size=9),
                    axis.title = element_blank(),
                        axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.text.y = element_text(color="black",size=8))) +
  labs(title = title_dict[task], x="Follow up time, days",y = "Survival probability")+
  scale_fill_manual(values = col)

label_i <- paste0(ppval1,'\n',
  '1-2: ',res_subtype[1,3],', ',res_subtype[1,4], '\n',
                  '1-3: ',res_subtype[2,3],', ',res_subtype[2,4], '\n',
                  '1-4: ',res_subtype[3,3],', ',res_subtype[3,4], '\n',
                  '2-3: ',res_subtype[4,3],', ',res_subtype[4,4], '\n',
                  '2-4: ',res_subtype[5,3],', ',res_subtype[5,4], '\n',
                  '3-4: ',res_subtype[6,3],', ',res_subtype[6,4])
p <- p + annotate("text", x = 100, y = 0.79, label = label_i,size = 3,hjust=0,vjust=0,family="Arial")



