library(stringr)

bdd=read.csv("Data/Fichier/peptides_DBT_MgSO4.txt",header=TRUE,sep="",blank.lines.skip=TRUE)

## Filtration -------------------------------------------------------


Count_Intensity_Sets=cbind(rowSums(bdd[,(3:6)]>0),     # On compte le nombre 
                           rowSums(bdd[,(7:10)]>0),     # de valeurs mesurées 
                           rowSums(bdd[,(11:14)]>0),    # pour chaque set :
                           rowSums(bdd[,(15:18)]>0),   # S'il manque un sets
                           rowSums(bdd[,(19:22)]>0),   # complet, une de ces
                           rowSums(bdd[,(23:26)]>0),
                           rowSums(bdd[,(27:30)]>0),
                           rowSums(bdd[,(31:34)]>0))   # valeurs est à 0.

Count_Intensity=apply(Count_Intensity_Sets,1,min)

bdd=bdd[Count_Intensity!=0,]

rm(Count_Intensity,Count_Intensity_Sets)

## Imputation -------------------------------------------------------
cond=c('Intensity_DBT_E',
       'Intensity_DBT_L',
       'Intensity_DBT_M',
       'Intensity_DBT_S',
       'Intensity_MgSO4_E',
       'Intensity_MgSO4_L',
       'Intensity_MgSO4_M',
       'Intensity_MgSO4_S')
###################################################
####                     Sets 1                 ### 
###################################################
line=str_detect(colnames(bdd),cond[1])

bddaux=bdd[,line]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,line]<-bddaux
rm(bddaux,line,ind)

###################################################
####                     Sets 2                 ### 
###################################################
line=str_detect(colnames(bdd),cond[2])

bddaux=bdd[,line]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,line]<-bddaux
rm(bddaux,line,ind)

###################################################
####                     Sets 3                 ### 
###################################################
line=str_detect(colnames(bdd),cond[3])

bddaux=bdd[,line]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,line]<-bddaux
rm(bddaux,line,ind)

###################################################
####                     Sets 4                 ### 
###################################################
line=str_detect(colnames(bdd),cond[4])

bddaux=bdd[,line]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,line]<-bddaux
rm(bddaux,line,ind)

###################################################
####                     Sets 5                 ### 
###################################################
line=str_detect(colnames(bdd),cond[5])

bddaux=bdd[,line]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,line]<-bddaux
rm(bddaux,line,ind)

###################################################
####                     Sets 6                 ### 
###################################################
line=str_detect(colnames(bdd),cond[6])

bddaux=bdd[,line]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,line]<-bddaux
rm(bddaux,line,ind)

###################################################
####                     Sets 7                 ### 
###################################################
line=str_detect(colnames(bdd),cond[7])

bddaux=bdd[,line]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,line]<-bddaux
rm(bddaux,line,ind)

###################################################
####                     Sets 8                 ### 
###################################################
line=str_detect(colnames(bdd),cond[8])

bddaux=bdd[,line]
bddaux[bddaux==0]<-NA

ind=which(is.na(bddaux),arr.ind=TRUE)
bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]

bdd[,line]<-bddaux
rm(bddaux,line,ind)

#############################
## Prétraitement terminé ---------------------------------------------


proteom=bdd
rm(bdd)
# -> to load in Prostar : 
write.table(proteom, "Data/Fichier/proteom_treated.txt",row.names=FALSE,quote=FALSE,sep="\t")

# -> to sort later : 
row.names(proteom) <- 1:nrow(proteom)
save(proteom,file='Data/save/proteom.RData')


