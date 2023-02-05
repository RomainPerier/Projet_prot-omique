library(sanssouci)
library(doBy)
library(stringr)

## Prétraitement  -------------------------------------------------------
bdd=read.csv("peptides_YEASTUPS.txt",header=TRUE,sep="",blank.lines.skip=TRUE)


Count_Intensity_Sets=cbind(rowSums(bdd[,(3:5)]>0),#On compte le nombre 
                           rowSums(bdd[,(6:8)]>0),     #de valeurs mesurée 
                           rowSums(bdd[,(9:11)]>0),    #pour chaque set :
                           rowSums(bdd[,(12:14)]>0),   #S'il manque un sets
                           rowSums(bdd[,(15:17)]>0),   #complet, une de ces
                           rowSums(bdd[,(18:20)]>0))   #valeurs est à 0.

Count_Intensity=apply(Count_Intensity_Sets,1,sum)

bdd_no_na=bdd[Count_Intensity==18,]

rm(Count_Intensity,Count_Intensity_Sets)


write.table(bdd_no_na, "proteom2_treated_na.txt",row.names=FALSE,quote=FALSE,sep="\t")



## Groupe ------------------------------------------------------------------------------------

res=read.csv("proteom2_treated_na_test.txt",header=TRUE,sep="\t")[,c('id','Pvalue','Leading_razor_protein')]

Species=!str_detect(res[,"Leading_razor_protein"],"ups") 
Species[Species==TRUE]<-"levure"
Species[Species==FALSE]<-"ups"
res=cbind(res,Species)
rm(Species)

res_ordered=orderBy(~ Species + Leading_razor_protein,res)

save(res_ordered,file='res_ordered_na.RData')


split_data=splitBy(formula = ~ Species + Leading_razor_protein,res_ordered)[]
save(split_data, file = "split_data_na.RData")


## Hiérarchie -------------------------------------------------------------------------------

## Les fichiers de données -------------------------------------------------------------

load("split_data_na.RData")
load("res_ordered_na.RData") 


n_leaf=length(split_data)

taille_leaf=rep(0,n_leaf) 
for (i in 1:n_leaf){
  taille_leaf[i]<-dim(split_data[[i]])[1]
}

taille_leaf_cum=cumsum(taille_leaf) 

leaf_list=list(1:(taille_leaf_cum[1]))   

for (i in 2:(n_leaf)){
  j=i-1
  leaf_list=append(leaf_list,list((taille_leaf_cum[j]+1):(taille_leaf_cum[i])))
}

aux=!str_detect(names(split_data),"ups")  
n_levure=length(aux[aux])  
n_ups=length(aux[!aux])   
n_species_cum=c(n_levure,n_levure+n_ups)   

rm(aux) 

C3=list(c(1,1))   

for (i in 2:n_leaf){
  C3<-append(C3,list(c(i,i)))
}


C=list(list(c(1,n_leaf)),list(c(1,n_levure),c(n_levure+1,n_leaf)),C3)

rm(C3)

final_tree=list(leaf_list,C)

save(final_tree,file="final_tree_na.RData")


