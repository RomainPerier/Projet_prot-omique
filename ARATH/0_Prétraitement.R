#On importe les données que l'on souhaite pré-traiter
bdd=read.csv("peptides_ARATHUPS_MBR.txt",header=TRUE,sep="",
             blank.lines.skip=TRUE)

## Filtration -------------------------------------------------------
#On retire des données toutes les lignes possédant 
#un échantillon complet manquant :

Count_Intensity_Sets=cbind(rowSums(bdd[,(3:5)]>0),#On compte le nombre 
                      rowSums(bdd[,(6:8)]>0),     #de valeurs mesurée 
                      rowSums(bdd[,(9:11)]>0),    #pour chaque set :
                      rowSums(bdd[,(12:14)]>0),   #S'il manque un sets
                      rowSums(bdd[,(15:17)]>0),   #complet, une de ces
                      rowSums(bdd[,(18:20)]>0))   #valeurs est à 0.

Count_Intensity=apply(Count_Intensity_Sets,1,min)

bdd_no_na=bdd[Count_Intensity!=0,]

rm(Count_Intensity,Count_Intensity_Sets)

## Imputatuin -------------------------------------------------------
#Pour chaque set d'échantillons, on remplace la valeur manquante par la
#moyenne de celles mesurées.

###################################################
####                     Sets 1                 ### 
###################################################
i=3
j=i+2

bddaux=bdd_no_na[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd_no_na[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

###################################################
####                     Sets 2                 ### 
###################################################
i=6
j=i+2

bddaux=bdd_no_na[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd_no_na[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

###################################################
####                     Sets 3                 ### 
###################################################
i=9
j=i+2

bddaux=bdd_no_na[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd_no_na[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

###################################################
####                     Sets 5                 ### 
###################################################
i=12
j=i+2

bddaux=bdd_no_na[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd_no_na[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

###################################################
####                     Sets 6                 ### 
###################################################

i=15
j=i+2

bddaux=bdd_no_na[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd_no_na[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

###################################################
####                     Sets 7                 ### 
###################################################

i=18
j=i+2

bddaux=bdd_no_na[,(i:j)]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd_no_na[,(i:j)]<-bddaux
rm(bddaux,i,j,ind)

#############################

## Prétraitement terminé ---------------------------------------------


bdd_fin=bdd_no_na
rm(bdd_no_na)

write.table(bdd_fin, "proteom_treated.txt",row.names=FALSE,quote=FALSE,sep="\t")
