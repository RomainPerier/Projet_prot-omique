library(doBy)
library(cgwtools)
library(stringr)

## Les fichiers de données -------------------------------------------------------------

load("save/data.RData")


## Création des feuilles au bon format--------------------------------------------------
# leaf_list = les petides de chaques protéines = aux data (faut récupérer leur ordre....)
# leaf_list=L1,L2,L3,.........
# tailles des feuilles = nombre de peptides de chaque protéines 
# nbr de feuilles = nbr de protéines 
# taille_leaf_cum[i]=dernière ligne de la ie protéine 

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
# C a trois niveaux 
#                 TOUT                      C1
#      levure             ups               C2 
# prot prot prot      prot prot prot prot   C3
# On a n_leaf protéines, on repère celle associé à levure ou non 



aux=!str_detect(names(split_data),"ups")
n_levure=length(aux[aux])    
n_ups=length(aux[!aux])    
n_species_cum=c(n_levure,n_levure+n_ups)
rm(aux) 

C3=list(c(1,1))   
for (i in 2:n_leaf){
  C3<-append(C3,list(c(i,i)))
}

# C1 <- list(c(1,n_leaf))
# C2 <-list(c(1,n_levure) , c(n_levure+1,n_leaf))
C=list(list(c(1,n_leaf)),list(c(1,n_levure),c(n_levure+1,n_leaf)),C3)    
rm(C3)

final_tree=list(leaf_list,C)
rm(leaf_list,C,m,n_levure,n_species_cum,n_ups,n_leaf,i)
resave(final_tree,file="save/data.RData")

## Useful object --------------------------------------------------------------

pval = res_ordered[,'Pvalue']
line_sorted_by_pval = order(pval)
pval_sorted = sort(pval)
proteom = cbind(proteom,pval)
resave(proteom,pval,line_sorted_by_pval,pval_sorted,file='save/data.RData')
