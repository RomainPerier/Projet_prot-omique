library(doBy)
library(stringr)

## Les fichiers de données -------------------------------------------------------------

load("save/proteom.RData")
load("save/res_ordered.RData")
load("save/split_data.RData")


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


C=list(C3)  
final_tree=list(leaf_list,C)
save(final_tree,file="save/final_tree.RData")


## Création de l'arbre amélioré regroupant tous les positifs -----
final_tree_1 = list(leaf_list,list(list(c(n_levure+1,n_leaf)),C3))
save(final_tree_1,file="save/final_tree_1.RData")

## Création de l'arbre amélioré regroupant tous les positifs puis tous les peptides -----
C=list(list(c(1,n_leaf)),list(c(n_levure+1,n_leaf)),C3)
final_tree_2=list(leaf_list,C)
save(final_tree_2,file="save/final_tree_2.RData")


rm(list = ls())
