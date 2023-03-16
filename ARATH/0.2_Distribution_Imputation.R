library(ggplot2)
library(stringr)
source("00_Theme.R")
source('Fonction/pretraitement.sets.R')
source('Fonction/no.pretraitement.R')
theme_set(theme_ben())

bdd=read.csv("ARATH/peptides_ARATHUPS_MBR.txt",header=TRUE,sep="",blank.lines.skip=TRUE)

data <- rbind(no.pretraitement(bdd),pretraitement.sets(bdd,1))

## Graphique -------------------------------------------------------------------

ggplot(data,aes(x=Value,fill=Méthode))+
  geom_density(alpha=0.2)+
  xlim(c(0,10000000))+
  facet_grid(Condition ~ .)
