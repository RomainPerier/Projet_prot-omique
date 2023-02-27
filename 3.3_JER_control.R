library(sanssouci)
library(ggplot2)
library(cgwtools)
library(doBy)
load('save/data.RData')
# JER : voir dans quel Rk JER violé dans notre régionnement 
# => trouver quelles zeta_k foire le JER !
# voir ce qu'il y a a conservé, qu'est ce qui fonctionne fonctionne pas 
# on va regarder quel Rk viole le JER 
# zeta contrôle le JER sur RK ie IP(FP dans Rk <= zeta_k)>= 1-alpha 

## Contrôle du JER sur les feuilles  -------------------------------------------
alpha = 0.05
m = nrow(proteom)
H = length(C)

df_JER_control = data.frame()
ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha)

ZL_leaf = ZL_DKWM[[length(final_tree[[2]])]]
leaf_list=final_tree[[1]]
n_leaf = length(leaf_list)

for (k in 1:n_leaf){
  Rk=leaf_list[[k]]
  Zk=ZL_leaf[k]
  FP = sum(!str_detect(proteom[Rk,'Leading_razor_protein'],"ups"))
  df_JER_control <- rbind(df_JER_control,data.frame(Index=k,Region=0,Zeta=Zk,FP=FP,is_violated=(FP>Zk)))
  df_JER_control[k,]$Region <- list(Rk)
}

View(df_JER_control[df_JER_control$is_violated,])
Rk_violated=unlist(df_JER_control[df_JER_control$is_violated,]$Region)
line=df_JER_control[df_JER_control$is_violated,]$Index
View(proteom[Rk_violated,])
View(proteom[line,])

name_prot_violated=names(splitBy(formula = ~ Leading_razor_protein,proteom[Rk_violated,]))
prot_violated=split_data[line]

len_prop_violated=unlist(lapply(splitBy(formula = ~ Leading_razor_protein,proteom[Rk_violated,]),function(x){dim(x[1])}))
len_prot=unlist(lapply(prot_violated,function(x){dim(x[1])}))


