#Il y a 6 conditions: 0,5fmol, 1fmol, 2,5fmol, 5fmol, 10fmol et 25fmol de protéines UPS
# deux groupes en agrégeant les 3 plus petites conditions d’une part et les 3 plus grandes d’autre part. 

bdd=read.csv("Fichier/peptides_YEASTUPS.txt",header=TRUE,sep="",blank.lines.skip=TRUE)

library("sanssouci")
## Filtration -------------------------------------------------------


## code pour récupérer les infos sur les colonnes de manière auto
vectcond=lapply(colnames(bdd), function(i) strsplit(i,"[_]")[[1]][3])
vectcond[is.na(vectcond)]<-FALSE
set1 <- vectcond=="0.5fmol"
set2 <- vectcond=="1fmol"
set3 <- vectcond=="2.5fmol"
set4 <- vectcond=="5fmol"
set5 <- vectcond=="10fmol"
set6 <- vectcond=="25fmol"



#On retire des données toutes les lignes possédant 
#un échantillon complet manquant :

Count_Intensity_Sets=cbind(rowSums(bdd[,set1]>0),     # On compte le nombre 
                           rowSums(bdd[,set2]>0),     # de valeurs mesurées 
                           rowSums(bdd[,set3]>0),    # pour chaque set :
                           rowSums(bdd[,set4]>0),   # S'il manque un sets
                           rowSums(bdd[,set5]>0),   # complet, une de ces
                           rowSums(bdd[,set6]>0))   # valeurs est à 0.

Count_Intensity=apply(Count_Intensity_Sets,1,min)

bdd=bdd[Count_Intensity==3,]

rm(Count_Intensity,Count_Intensity_Sets)

## Pas d'imputation

bdd <- cbind(bdd[,(1:2)],rowMeans(bdd[,set1]),rowMeans(bdd[,set2]),rowMeans(bdd[,set3]),rowMeans(bdd[,set4]),rowMeans(bdd[,set5]),rowMeans(bdd[,set6]))

# Calcul des p-valeurs

categ <- rep(c(0,1), times = c(3,3))

prot_matr <- data.matrix(bdd[,(3:8)])

tests<- rowWelchTests(prot_matr, categ, alternative = "greater")

Pvalue <- tests$p.value

res <- cbind(proteom[,(1:2)],Pvalue)

save(res,file="save/res_test.Rdata")


rm(list = ls())

