#Il s’agit d’une expérience similaire que précédemment, sauf qu’on a utilisé un fond constant de 800ng protéines de levure au lieu des 400ng de protéines d’arabidopsis.
#Il y a 6 conditions: 0,5fmol, 1fmol, 2,5fmol, 5fmol, 10fmol et 25fmol de protéines UPS. 
#Pour distinguer les protéines UPS des protéines de levure cette fois: si le nom de la protéine termine par « ups » c’est que c’est une protéine UPS, sinon, c’est que c’est une protéine de levure. 
#Je vous propose à nouveau de construire deux groupes en agrégeant les 3 plus petites conditions d’une part et les 3 plus grandes d’autre part. 


# Code identique que le précédent car tableau de données à la même forme 


#On importe les données que l'on souhaite pré-traiter
bdd=read.csv("peptides_YEASTUPS.txt",header=TRUE,sep="",blank.lines.skip=TRUE)

## Filtration -------------------------------------------------------
#On retire des données toutes les lignes possédant 
#un échantillon complet manquant :

Count_Intensity_Sets=cbind(rowSums(bdd[,(3:5)]>0),#On compte le nombre 
                           rowSums(bdd[,(6:8)]>0),     #de valeurs mesurée 
                           rowSums(bdd[,(9:11)]>0),    #pour chaque set :
                           rowSums(bdd[,(12:14)]>0),   #S'il manque un sets
                           rowSums(bdd[,(15:17)]>0),   #complet, une de ces
                           rowSums(bdd[,(18:20)]>0))   #valeurs est à 0.

Count_Intensity = apply(Count_Intensity_Sets,1,min)   
  # Permet de repérer s'il y a des sets manquants 

bdd=bdd[Count_Intensity!=0,]
  # --> on retire les lignes ayant plus d'un set entier manquant 

rm(Count_Intensity)
rm(Count_Intensity_Sets)

Prop_Intensity = rowSums(bdd[,(3:20)]>0)/18*100
  # Proportion de Na par ligne 

bdd_save = bdd 

bdd = bdd_save[Prop_Intensity>=70,]

## Imputation -------------------------------------------------------
#Pour chaque set d'échantillons, on remplace la valeur manquante par la
#moyenne de celles mesurées.

###################################################
####                     Sets 1                 ### 
###################################################
i=3
j=i+2

bddaux=bdd[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

###################################################
####                     Sets 2                 ### 
###################################################
i=6
j=i+2

bddaux=bdd[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

###################################################
####                     Sets 3                 ### 
###################################################
i=9
j=i+2

bddaux=bdd[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

###################################################
####                     Sets 5                 ### 
###################################################
i=12
j=i+2

bddaux=bdd[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

###################################################
####                     Sets 6                 ### 
###################################################

i=15
j=i+2

bddaux=bdd[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

###################################################
####                     Sets 7                 ### 
###################################################

i=18
j=i+2

bddaux=bdd[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

#############################
## Prétraitement terminé ---------------------------------------------


proteom_treated=bdd
rm(bdd)



save(proteom_treated,file="proteom_treated_seuil_70.RData")
write.table(proteom_treated, "proteom_treated_seuil_70.txt",row.names=FALSE,quote=FALSE,sep="\t")
