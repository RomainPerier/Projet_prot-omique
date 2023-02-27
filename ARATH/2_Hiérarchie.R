library(doBy)
library(stringr)

## Les fichiers de données -------------------------------------------------------------

load("split_data.RData")   #Les peptides et les protéines 
load("res_ordered.RData")  #Toutes les données mais tri


## Création des feuilles au bon format--------------------------------------------------
# leaf_list = les petides de chaques protéines = aux data (faut récupérer leur ordre....)
# leaf_list=L1,L2,L3,.........
#  Chaque feuille L correspond aux bornes des lignes qu'elle contient 

n_leaf=length(split_data)    #nombre de feuilles = nombres de protéines 

taille_leaf=rep(0,n_leaf)    #tailles des feuilles = nbr de peptides de chaque protéines 
for (i in 1:n_leaf){
  taille_leaf[i]<-dim(split_data[[i]])[1]
}

taille_leaf_cum=cumsum(taille_leaf) #tailles cumulées => retrouver les bonnes lignes de chaque peptides  
#La feuille i se termine ligne l[i] 
leaf_list=list(1:(taille_leaf_cum[1]))   

for (i in 2:(n_leaf)){
  j=i-1
  leaf_list=append(leaf_list,list((taille_leaf_cum[j]+1):(taille_leaf_cum[i])))
}



## Création de C -----------------------------------------------------------------------
# C a trois niveaux 
#                 TOUT                      C1
#      ARATH             UPS1               C2 
# prot prot prot      prot prot prot prot   C3
# On a n_leaf protéines, on repère celle associé à ARATH ou non 



aux=str_detect(names(split_data),"ARATH")
n_ARATH=length(aux[aux])    #nombre arath = nombre de protéines d'arabidopsis 
n_UPS1=length(aux[!aux])    #      UPS1 = prot humaines
n_species_cum=c(n_ARATH,n_ARATH+n_UPS1)   #cumulée pour retrouver les bonnes lignes de chaques espèces 

rm(aux) #adios amigos 

C3=list(c(1,1))   #A AMELIORER (voir pour créer C3 d'un coup ?????) u osef 

for (i in 2:n_leaf){
  C3<-append(C3,list(c(i,i)))
}


C=list(list(c(1,n_leaf)),list(c(1,n_ARATH),c(n_ARATH+1,n_leaf)),C3)

rm(C3)

final_tree=list(leaf_list,C)

save(final_tree,file="final_tree.RData")

