#Il y a 6 conditions: 0,5fmol, 1fmol, 2,5fmol, 5fmol, 10fmol et 25fmol de protéines UPS
# deux groupes en agrégeant les 3 plus petites conditions d’une part et les 3 plus grandes d’autre part. 

bdd=read.csv("Fichier/peptides_YEASTUPS.txt",header=TRUE,sep="",blank.lines.skip=TRUE)

## Filtration -------------------------------------------------------

## code pour récupérer les infos sur les colonnes de manière auto
vectcond=sapply(colnames(bdd), function(i) strsplit(i,"[_]")[[1]][3])
vectcond=="0.5fmol"
vectcond=="10fmol"
# etc... 


#On retire des données toutes les lignes possédant 
#un échantillon complet manquant :

Count_Intensity_Sets=cbind(rowSums(bdd[,(3:5)]>0),     # On compte le nombre 
                           rowSums(bdd[,(6:8)]>0),     # de valeurs mesurées 
                           rowSums(bdd[,(9:11)]>0),    # pour chaque set :
                           rowSums(bdd[,(12:14)]>0),   # S'il manque un sets
                           rowSums(bdd[,(15:17)]>0),   # complet, une de ces
                           rowSums(bdd[,(18:20)]>0))   # valeurs est à 0.

Count_Intensity=apply(Count_Intensity_Sets,1,min)

bdd=bdd[Count_Intensity!=0,]

rm(Count_Intensity,Count_Intensity_Sets)

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


proteom=bdd
rm(bdd)
# -> to load in Prostar : 
write.table(proteom, "Fichier/proteom_treated.txt",row.names=FALSE,quote=FALSE,sep="\t")

# -> to sort later : 
row.names(proteom) <- 1:nrow(proteom)
save(proteom,file='save/proteom.RData')


##ça push ou pas?--------