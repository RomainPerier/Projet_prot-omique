library(stringr)
library(ggplot2)
#load('ARATH/save/pval.RData')
#load('ARATH/save/res_ordered.RData')
load('ARATH/save/pval2.RData')
load('ARATH/save/res_ordered2.RData')
source("00_Theme.R")
theme_set(theme_ben())

## Data_frame ------------------------------------------------------
m=nrow(res_ordered)

vct_cond=res_ordered$Species
vct_cond[vct_cond=='ARATH']<-'H0'
vct_cond[vct_cond=='UPS1']<-'H1'

df=data.frame(Pvalue=res_ordered$Pvalue,Condition=vct_cond)
df <- df[line_sorted_by_pval,]

ggplot(df,aes(x=Pvalue,color=Condition))+
  stat_ecdf(geom='step',lwd=1)+
  xlab('')+ylab('')+
  ggtitle('Fonction de réparition empirique des p-valeurs selon la condition')+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),show.legend = NA,lwd=1,linetype='dashed',color='black')



