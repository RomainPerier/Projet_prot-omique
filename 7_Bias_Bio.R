library(sanssouci)
library(cgwtools)
library(ggplot2)
library(reshape2)
load('save/proteom.RData')
load('save/res_ordered.RData')
load('save/final_tree.RData')
load('save/pval.RData')
load('save/final_tree.RData')


## Procédure de rejet ----------------------------------------------------------

BH_procedure <- function(pval,alpha=0.05){
  pval_sorted = sort(pval)
  m = length(pval)
  thr =alpha/m*(1:m)
  df = data.frame(Index=1:m,res=(pval_sorted<=thr))
  khat = max(df[df$res==TRUE,]$Index)
  
  df2 = data.frame(Index=1:m,res2=(pval<=alpha*khat/m))
  R = df2[df2$res2==TRUE,]$Index
  return(R)
}

BF_procedure <- function(pval,alpha=0.05){
  m = length(pval)
  df = data.frame(Index=1:m,res=(pval<=alpha/m))
  R = df[df$res==TRUE,]$Index
  return(R)
  
}

BY_procedure <- function(pval,alpha=0.05){
  pval_sorted = sort(pval)
  m = length(pval)
  Hm = sum(1/(1:m))
  thr =alpha/(m*Hm)*(1:m)
  df = data.frame(Index=1:m,res=(pval_sorted<=thr))
  khat = max(df[df$res==TRUE,]$Index)
  
  df2 = data.frame(Index=1:m,res2=(pval<=alpha*khat/(m*Hm)))
  R = df2[df2$res2==TRUE,]$Index
  return(R)
}

## Biais Biol ------------------------------------------------------------------
# Sk={pvaleur rejete par BH par les k plus petites valeurs} -> pas de garantis statistiques 
# mais tester V* qui a des garantis 

listk = 500*(1:10)
df_res = data.frame(Méthode=c('BH','BF','BY'))
alpha=0.05
ZL = zetas.tree(C = final_tree[[2]],leaf_list = final_tree[[1]],zeta.DKWM,pval,alpha)

for (k in listk){
  pvalk = sort(pval)[1:k]
  Sk1 = BH_procedure(pvalk)
  Sk2 = BF_procedure(pvalk)
  Sk3 = BY_procedure(pvalk)
  V1 = V.star(Sk1,C = final_tree[[2]], leaf_list = final_tree[[1]], ZL)
  V2 = V.star(Sk2,C = final_tree[[2]], leaf_list = final_tree[[1]], ZL)
  V3 = V.star(Sk3,C = final_tree[[2]], leaf_list = final_tree[[1]], ZL)
  
  df_res = cbind(df_res,col=c(V1,V2,V3))
  colnames(df_res)[ncol(df_res)] <- paste('Vstar_k',k,sep='_')
}