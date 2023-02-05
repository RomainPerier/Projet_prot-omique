library(sanssouci)
library(stringr)
library(doBy)

## Données  -----------------------------------------------------------------

load('res_ordered.RData')
load('final_tree.RData')
#load('proteom_treated.RData')
load('proteom_reorganised.RData')

#On veut faire correspondre les lignes de proteom_treated et res_ordered 

proteom<-proteom[as.numeric(rownames(res_ordered)),]

#> all(proteom[,'Sequence']==res_ordered[,'id'])
#[1] TRUE

pval = res_ordered[,'Pvalue']

proteom=cbind(proteom,pval)

# Trie par p-valeurs -----------------------------------------------------------

n=length(pval)

aux <- res_ordered

rownames(aux)<-1:n  #Le nom de la ligne dans aux correspondra à sa position dans res_ordered et ON NE TOUCHE JAMAIS A RES_ORDERED 

aux=orderBy(~ Pvalue,aux)  #En triant, on récupère les lignes de aux trié par p-val mais donc les positions dans res_ordered dans l'ordre 
# typiquement res_ordered[as.numeric(rownames(aux)),] -> trié par p-val

pval_sorted=aux[,'Pvalue']
line_sorted_by_pValue=as.numeric(rownames(aux)) 

rm(aux) 

# On veut regarder les valeurs des 100 plus petites p-valeurs : 

name=proteom[line_sorted_by_pValue[1:100],'Sequence']
small_before_imputation = rep(0,100)
i=1

for (seq in name){
  small_before_imputation[i]<- as.numeric(rownames(bdd[bdd$Sequence==seq,]))
  i=i+1
}

bdd_small_pval = bdd[small_before_imputation,]

bdd_small_pval[bdd_small_pval==0]<-NA
summary(bdd_small_pval)

save(bdd_small_pval,file='bdd_small_pval.RData')


## True Positif -------------------------------------------------------
# On veut ici retrouver les ups 

TP=cumsum(aux$Species=='ups')

plot(TP)
lines(c(0,1),c(0,1),col='red',lwd=2.0)

FP = 1:n - TP


## V.star  ----------------------------------------------------------------------
# V.star(S =,C = ,ZL = ,leaf_list = )
# zetas.tree(C = ,leaf_list = ,method = ,pvalues = ,alpha = )


C=final_tree[[2]]
leaf_list=final_tree[[1]]
alpha=0.05

ZL=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha)


iter_pointeur = 0      #où tout les 10
Vstar=rep(0,n)


for (i in 1:n){
  S=line_sorted_by_pValue[1:i]
  Vstar[i]<-V.star(S,C,ZL,leaf_list)
  iter_pointeur=iter_pointeur+1
  print(iter_pointeur)
}

save(Vstar,file='Vstar_save.RData')

## BySimes -------------------------------------------------------------------
# post hoc bound obtained from Simes' inequality
#
# Lower bound on the number of correct rejections (TP) using Simes' reference
# family   -> Vsimes
#


### Graph : 

plt.plot()



#curve Max fp  pour l'utiliser sanssouci:::curveMaxFP(pval,thr)   -> Simes


# hybride ? 
# faire BF BH et voir combien de fp sont estimer par V* 
# tester BH sur les p-val restreintes m= avec les plus petites
# Sk={pvaleur rejete par BH par les k plus petites valeurs} -> pas de garantis statistiques 
# mais tester V* qui a des garantis 
#  VDKWM, V HB zetak 
#     zeta -> dekwm au dessus puis hb sur feuilles 
#  hybride Vdkwm et Vhb 
# V*DKHM(zeta.dkwm(alpha)), et  V*HB(zeta.HB(alpha))
# plusieurs manières d'hybdrider 
# 1. Vhybri = min (VDKWM(S),VHB(S))   -> mais perds garantis stats sauf si avec poids 
#     i.e. Vhybrid = min (VDKWM(zeta(gamma alpha)),VHB(zeta(1-gammaalpha))) 
# 2. faire un min avec les zeta ? 
#   car même R_k ! 
#     Vhybdri' = V*(min(zeta.DKWM(gamma alpha),zeta.BH(1- gamma alpha))) tjrs inférieur à la premièreb manièer -> garantis stats 
#    Egalité stricte ??? -> si oui mieux ! 
#    boucle for qui parcourt les zeta et 
#
#  varier gqamma aussi 
#
#  en théorie gamma petit n'empêche pas la perforcmance de v* 
#
# papier dbnr gamma=0.2

