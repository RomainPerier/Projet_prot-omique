library(doBy)
library(stringr)

## Les fichiers de données -------------------------------------------------------------

load("split_data_70.RData")   #Les peptides et les protéines 
load("res_ordered_70.RData")  #Toutes les données mais tri


## Création des feuilles au bon format--------------------------------------------------
# leaf_list = les petides de chaques protéines = aux data (faut récupérer leur ordre....)
# leaf_list=L1,L2,L3,.........
#  Chaque feuille L correspond aux bornes des lignes qu'elle contient 
n=dim(res_ordered)[1]

n_leaf=length(split_data)    #nombre de feuilles = nombres de protéines 

taille_leaf=rep(0,n_leaf)    #tailles des feuilles = nbr de peptides de chaque protéines 
for (i in 1:n_leaf){
  taille_leaf[i]<-dim(split_data[[i]])[1]
}


# On veut créer une liste qui indique la première dernière ligne de chaque prot : 
taille_leaf_cum=cumsum(taille_leaf) #tailles cumulées => retrouver les bonnes lignes de chaque peptides  


#taille_leaf_cum[i]=dernière ligne de la ie protéine 

#maintenant, on va créer la list des feuilles 

leaf_list=list(1:(taille_leaf_cum[1])) #première feuille initialiser 

for (i in 2:n_leaf){
  j=i-1
  leaf_list=append(leaf_list,list((taille_leaf_cum[j]+1):(taille_leaf_cum[i])))
}


## Création de C -----------------------------------------------------------------------
# C a trois niveaux 
#                 TOUT                      C1
#      levure             ups               C2 
# prot prot prot      prot prot prot prot   C3
# On a n_leaf protéines, on repère celle associé à levure ou non 



aux=!str_detect(names(split_data),"ups")   #on repère les protéines de levure
n_levure=length(aux[aux])    #nombre levure = nombre de protéines de levure 
n_ups=length(aux[!aux])    #      ups = prot humaines
n_species_cum=c(n_levure,n_levure+n_ups)   #cumulée pour retrouver les bonnes lignes de chaques espèces 

rm(aux) #adios amigos 

C3=list(c(1,1))   #A AMELIORER (voir pour créer C3 d'un coup ?????) u osef 

for (i in 2:n_leaf){
  C3<-append(C3,list(c(i,i)))
}


C=list(list(c(1,n_leaf)),list(c(1,n_levure),c(n_levure+1,n_leaf)),C3)    
#list(c(1,n_leaf)) -> C1
#list(c(1,n_levure) , c(n_levure+1,n_leaf))-> C2

rm(C3)

final_tree=list(leaf_list,C)

save(final_tree,file="final_tree_70.RData")

