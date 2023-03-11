library(ggplot2)
library(stringr)
source("00_Theme.R")
theme_set(theme_ben())

bdd_avant=read.csv("ARATH/peptides_ARATHUPS_MBR.txt",header=TRUE,sep="",blank.lines.skip=TRUE)
bdd_apres=read.csv("ARATH/proteom_treated_arath.txt",header=TRUE,sep="",blank.lines.skip=TRUE)

## Tableau des mesures ---------------------------------------------------------
# Avant Imputation :
tab_data_av_low=bdd_avant[,3:11]
tab_data_av_high=bdd_avant[,12:20]

data_av_low = unlist(tab_data_av_low,use.names=FALSE)
data_av_high = unlist(tab_data_av_high,use.names=FALSE)

data_av=data.frame(Méthode='no_imputation',
                   rbind(data.frame(Protéine= rep(bdd_avant$Leading_razor_protein,9)[data_av_low!=0],
                                    Condition="Low Spike",
                                    Value=data_av_low[data_av_low!=0])
                         ,data.frame(Protéine= rep(bdd_avant$Leading_razor_protein,9)[data_av_high!=0],
                                     Condition="High Spike",
                                     Value=data_av_high[data_av_high!=0])))

#Après Imputation :
tab_data_ap_low=bdd_apres[,3:11]
tab_data_ap_high=bdd_apres[,12:20]

data_ap_low = data.frame(Condition="Low Spike",Value=unlist(tab_data_ap_low,use.names=FALSE))
data_ap_high = data.frame(Condition="High Spike",Value=unlist(tab_data_ap_high,use.names=FALSE))


data_ap=data.frame(Méthode='imputation',
                   Protéine=rep(bdd_apres$Leading_razor_protein,18),
                   rbind(data_ap_low,data_ap_high))

## Graphique -------------------------------------------------------------------

data = rbind(data_ap,data_av)

ggplot(data,aes(x=Value,fill=Méthode))+
  geom_density(alpha=0.2)+
  xlim(c(0,10000000))+
  facet_grid(Condition ~ .)