library(tidyverse)
library(ggalluvial)
library(data.table)
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(reshape)

##Subtype 1
Rdata1 <- read.csv('result_subtype1_disease_proportion.csv')
p_tra1 <- Rdata1 %>% 
  ggplot(aes(x=Stage, y=value, fill=variable, stratum=variable, alluvium=variable)) +
  geom_flow(width = 0.33, curve_type = "linear", alpha=0.8) +
  geom_stratum(width = 0.33, alpha=0.9, color=NA) +
  geom_alluvium(width = 0.33,curve_type = "linear", fill=NA, color="#f4e2de") +
  scale_y_continuous(
    name = NULL,limits = c(0,300),breaks=c(0,50,100,150,200,250)) +
  labs(title = '',x='',y='Proportion',fill='') +
  theme_minimal() +
  theme(plot.background = element_rect(fill='white', color='white'),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        legend.text=element_blank())+
  scale_fill_manual(values = c('#B02323', '#E8503F', '#EDB39A', '#E6F1F0', '#97BFD3', '#21689F', '#CDCBD2'))+
  guides(
    color = guide_legend(nrow = 3, byrow = TRUE))

##Subtype 2
Rdata2 <- read.csv('result_subtype2_disease_proportion.csv')
p_tra2 <- Rdata2 %>% 
  ggplot(aes(x=Stage, y=value, fill=variable, stratum=variable, alluvium=variable)) +
  geom_flow(width = 0.33, curve_type = "linear", alpha=0.8) +
  geom_stratum(width = 0.33, alpha=0.9, color=NA) +
  geom_alluvium(width = 0.33,curve_type = "linear", fill=NA, color="#f4e2de") +
  scale_y_continuous(
    name = NULL,limits = c(0,300),breaks=c(0,50,100,150,200,250)) +
  labs(title = '',x='',y='Proportion',fill='') +
  theme_minimal() +
  theme(plot.background = element_rect(fill='white', color='white'),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        legend.text=element_blank())+
  scale_fill_manual(values = c('#B02323', '#E8503F', '#EDB39A', '#E6F1F0', '#97BFD3', '#21689F', '#CDCBD2'))+
  guides(
    color = guide_legend(nrow = 3, byrow = TRUE))

##Subtype 3
Rdata3 <- read.csv('result_subtype3_disease_proportion.csv')
p_tra3 <- Rdata3 %>% 
  ggplot(aes(x=Stage, y=value, fill=variable, stratum=variable, alluvium=variable)) +
  geom_flow(width = 0.33, curve_type = "linear", alpha=0.8) +
  geom_stratum(width = 0.33, alpha=0.9, color=NA) +
  geom_alluvium(width = 0.33,curve_type = "linear", fill=NA, color="#f4e2de") +
  scale_y_continuous(
    name = NULL,limits = c(0,300),breaks=c(0,50,100,150,200,250)) +
  labs(title = '',x='',y='Proportion',fill='') +
  theme_minimal() +
  theme(plot.background = element_rect(fill='white', color='white'),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        legend.text=element_blank())+
  scale_fill_manual(values = c('#B02323', '#E8503F', '#EDB39A', '#E6F1F0', '#97BFD3', '#21689F', '#CDCBD2'))+
  guides(
    color = guide_legend(nrow = 3, byrow = TRUE))

##Subtype 4
Rdata4 <- read.csv('result_subtype4_disease_proportion.csv')
p_tra4 <- Rdata4 %>% 
  ggplot(aes(x=Stage, y=value, fill=variable, stratum=variable, alluvium=variable)) +
  geom_flow(width = 0.33, curve_type = "linear", alpha=0.8) +
  geom_stratum(width = 0.33, alpha=0.9, color=NA) +
  geom_alluvium(width = 0.33,curve_type = "linear", fill=NA, color="#f4e2de") +
  scale_y_continuous(
    name = NULL,limits = c(0,300),breaks=c(0,50,100,150,200,250)) +
  labs(title = '',x='',y='Proportion',fill='') +
  theme_minimal() +
  theme(plot.background = element_rect(fill='white', color='white'),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        legend.text=element_blank())+
  scale_fill_manual(values = c('#B02323', '#E8503F', '#EDB39A', '#E6F1F0', '#97BFD3', '#21689F', '#CDCBD2'))+
  guides(
    color = guide_legend(nrow = 3, byrow = TRUE))

ggsave('disease_progression_s1.png',p_tra1, width = 6, height = 4, dpi = 600, bg='white')
ggsave('disease_progression_s2.png',p_tra2, width = 6, height = 4, dpi = 600, bg='white')
ggsave('disease_progression_s3.png',p_tra3, width = 6, height = 4, dpi = 600, bg='white')
ggsave('disease_progression_s4.png',p_tra4, width = 6, height = 4, dpi = 600, bg='white')