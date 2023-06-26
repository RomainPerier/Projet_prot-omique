library(stringr)
library(ggplot2)
load('save/pval.RData')
load('save/res_ordered.RData')

## Data_frame ------------------------------------------------------
m=nrow(res_ordered)

vct_cond=res_ordered$Species
vct_cond[vct_cond=='levure']<-'H0'
vct_cond[vct_cond=='ups']<-'H1'

df=data.frame(Pvalue=res_ordered$Pvalue,Hypothesis=vct_cond)
df <- df[line_sorted_by_pval,]

ggplot(df,aes(x=Pvalue,color=Hypothesis))+
  stat_ecdf(geom='step',lwd=1)+
  xlab('')+ylab('')+
  ggtitle("Cdf under null/alternative hypothesis")+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),show.legend = NA,lwd=1,linetype='dashed',color='black')

rm(list = ls())
