library(ggplot2)
library(stringr)
library(philentropy)
##source("00_Theme.R")
##theme_set(theme_ben())

## Méthode ---------------------------------------------------------------------

source('Fonction/pretraitement.seuil.R')
source('Fonction/pretraitement.sets.R')
source('Fonction/no.pretraitement.R')

bdd=read.csv("Fichier/peptides_YEASTUPS.txt",header=TRUE,sep="",blank.lines.skip=TRUE)

listseuil=c(0.4,0.5,0.6)
listp=c(1,2)
data=data.frame()

# Data_avant : sans prétraitement
data <- rbind(data,no.pretraitement(bdd))
  
# Data_apres : avec imputation 
  for (seuil in listseuil){
    data=rbind(data,pretraitement.seuil(bdd,seuil))
  }
  
  for (p in listp){
    data = rbind(data,pretraitement.sets(bdd,p))
  }
  

## Graphique -------------------------------------------------------------------


ggplot(data,aes(x=Value,fill=Méthode))+
  geom_density(alpha=0.2,lwd=0.5)+
  xlim(c(0,30000000))+
  facet_grid(Condition ~ .)



meth=unique(data$Méthode)
keep = meth[c(1,5)]


ggplot(data[data$Méthode==keep,],aes(x=Value,fill=Méthode))+
  geom_density(alpha=0.2,lwd=0.5)+
  xlim(c(0,30000000))+
  facet_grid(Condition ~ .)

