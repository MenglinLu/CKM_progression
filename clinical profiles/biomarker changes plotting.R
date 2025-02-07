library(ggiraphExtra)
library(ggplot2)

data <- read.csv('result_s3_biomarker_mean.csv')
group_colors <- c('#E8503F','#A8CBDF','#72A3C0')
p1 <- ggRadar(data=data[1:3,2:18],aes(group=phase),rescale=FALSE, size=2,
              fill = TRUE,
              fill.alpha = 0.1)+ 
  scale_fill_manual(values=group_colors) +
  scale_color_manual(values=group_colors) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    #axis.line = element_line(color = "gray"),
    panel.grid.major = element_line(color = "#dfe2e4"),
    panel.grid.minor = element_line(color = "#dfe2e4"))+labs(title='') + scale_y_continuous(limits = c(-0.2, max(data[,3:18])))
p2 <- ggRadar(data=data[4:6,2:18],aes(group=phase),rescale=FALSE, size=2,
              fill = TRUE,
              fill.alpha = 0.1)+ 
  scale_fill_manual(values=group_colors) +
  scale_color_manual(values=group_colors) +
  theme(
    legend.position = "none",
    
    panel.background = element_rect(fill = "white"),
    #axis.line = element_line(color = "gray"),
    panel.grid.major = element_line(color = "#dfe2e4"),
    panel.grid.minor = element_line(color = "#dfe2e4"))+labs(title='') + scale_y_continuous(limits = c(-0.2, max(data[,3:18])))
p3 <- ggRadar(data=data[7:9,2:18],aes(group=phase),rescale=FALSE, size=2,
              fill = TRUE,
              fill.alpha = 0.1)+ 
  scale_fill_manual(values=group_colors) +
  scale_color_manual(values=group_colors) +
  theme(
    legend.position = "none",
    
    panel.background = element_rect(fill = "white"),
    #axis.line = element_line(color = "gray"),
    panel.grid.major = element_line(color = "#dfe2e4"),
    panel.grid.minor = element_line(color = "#dfe2e4"))+labs(title='') + scale_y_continuous(limits = c(-0.2, max(data[,3:18])))
p4 <- ggRadar(data=data[10:12,2:18],aes(group=phase),rescale=FALSE, size=2,
              fill = TRUE,
              fill.alpha = 0.1)+ 
  scale_fill_manual(values=group_colors) +
  scale_color_manual(values=group_colors) +
  theme(
    legend.position = "none",
    
    panel.background = element_rect(fill = "white"),
    #axis.line = element_line(color = "gray"),
    panel.grid.major = element_line(color = "#dfe2e4"),
    panel.grid.minor = element_line(color = "#dfe2e4"))+labs(title='') + scale_y_continuous(limits = c(-0.2, max(data[,3:18])))

ggsave('biomarker_changes_s1.png',p1, width = 4, height = 4, dpi = 600, bg='white')
ggsave('biomarker_changes_s2.png',p2, width = 4, height = 4, dpi = 600, bg='white')
ggsave('biomarker_changes_s3.png',p3, width = 4, height = 4, dpi = 600, bg='white')
ggsave('biomarker_changes_s4.png',p4, width = 4, height = 4, dpi = 600, bg='white')
