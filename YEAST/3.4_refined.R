library(sanssouci)
library(doBy)

load('data.RData')


## But -------------------------------------------------------------------------
# On veut tester la fonction de raffinement zetas.tree.refined, comparer les
# V.star obtenu avec zeta_HB et zeta_refined 
#
# Attention : 
# || Error in nb.elements(C) : could not find function "nb.elements"
#  >>>> nv.elements indisponible dans sanssouci
#  >>>> mais sanssouci:::nb.elements marche 
#
#
#
#
## Vstar -----------------------------------------------------------------------
pval = res_ordered[,'Pvalue']
m=length(pval)

C=final_tree[[2]]
leaf_list=final_tree[[1]]
alpha=0.05

ZL_refined=zetas.tree.refined(C,leaf_list,zeta.HB,pval,alpha)
ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval,alpha)

iter_pointeur = 0
# On va tester sur les 300 dernières 
Vstar_refined=rep(0,300)
Vstar_aux=rep(0,300)
Vsimes_aux=rep(0,300)

thr=alpha/m*(1:m)


for (i in 1301:1600){
  S=line_sorted_by_pValue[1:i*10]
  Vstar_refined[i]<-V.star(S,C,ZL_refined,leaf_list)
  
  Vstar_aux[i]<-V.star(S,C,ZL_HB,leaf_list)
  
  Vsimes_aux[i]<-sanssouci:::curveMaxFP(pval[S],thr[S]) 
  iter_pointeur=iter_pointeur+1
  
  print(iter_pointeur)
}


save(Vstar_aux, Vsimes_aux, Vstar_refined,file='Test_refined.RData')
