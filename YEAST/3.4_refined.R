#Enfin, si les protéines sont assez petites en terme de nombre de peptides, on peut peut-être enfin essayer ma méthode de raffinement des zeta HB (zetas.tree.refined).


library(sanssouci)
library(cherry)
library(doBy)

load('data.RData')


######################### LIGNE 28 ###############################################"



# Vstar.no.id ----------------------------------------------------------------
# V.star(S =,C = ,ZL = ,leaf_list = )
# zetas.tree(C = ,leaf_list = ,method = ,pvalues = ,alpha = )
pval = res_ordered[,'Pvalue']

C=final_tree[[2]]
leaf_list=final_tree[[1]]
alpha=0.05

ZL_refined=zetas.tree.refined(C,leaf_list,zeta.DKWM,pval,alpha)


iter_pointeur = 0
aux = 1869   # m//10
Vstar_nrefined=rep(0,aux)


for (i in 1:aux){
  S=line_sorted_by_pValue[1:i*10]
  Vstar_refined[i]<-V.star(S,C,ZL_DKWM,leaf_list)
  iter_pointeur=iter_pointeur+1
  print(iter_pointeur)
}
save(Vstar_refined,file='Test_refined.RData')