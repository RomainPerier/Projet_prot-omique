library(doBy)
library(stringr)

## Les fichiers de données -------------------------------------------------------------

load("ARATH/save/proteom.RData")
load("ARATH/save/res_ordered.RData")
load("ARATH/save/split_data.RData")


## Création des feuilles au bon format--------------------------------------------------

m=nrow(res_ordered)
n_leaf=length(split_data)


taille_leaf=rep(0,n_leaf) 
for (i in 1:n_leaf){
  taille_leaf[i]<-nrow(split_data[[i]])
}
taille_leaf_cum=cumsum(taille_leaf) 


leaf_list=list(1:(taille_leaf_cum[1]))
for (i in 2:n_leaf){
  j=i-1
  leaf_list=append(leaf_list,list((taille_leaf_cum[j]+1):(taille_leaf_cum[i])))
}
rm(taille_leaf,taille_leaf_cum,i,j)

## Création de C -----------------------------------------------------------------------

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
rm(leaf_list,C,m,n_levure,n_species_cum,n_ups,n_leaf,i)
save(final_tree,file="ARATH/save/final_tree.RData")
