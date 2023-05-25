library(sanssouci)
library(cgwtools)
library(stringr)
library(ggplot2)
library(reshape2)
load('Low_pvalue_pull_out/save/proteom.RData')
load('Low_pvalue_pull_out/save/res_ordered.RData')
load('Low_pvalue_pull_out/save/final_tree.RData')
load('Low_pvalue_pull_out/save/pval.RData')
load('Low_pvalue_pull_out/save/final_tree.RData')


## Procédure de rejet ----------------------------------------------------------

BH_procedure <- function(pval,alpha=0.05){#Benjamini-Hochberg
  pval_sorted = sort(pval)
  m = length(pval)
  thr =alpha/m*(1:m)
  df = data.frame(Index=1:m,res=(pval_sorted<=thr))
  khat = max(df[df$res==TRUE,]$Index)
  
  df2 = data.frame(Index=1:m,res2=(pval<=alpha*khat/m))
  R = df2[df2$res2==TRUE,]$Index
  return(R)
}

BF_procedure <- function(pval,alpha=0.05){#Bonferroni
  m = length(pval)
  df = data.frame(Index=1:m,res=(pval<=alpha/m))
  R = df[df$res==TRUE,]$Index
  return(R)
  
}

HB_procedure <- function(pval,alpha=0.05){#Holm-Bonferroni
  m=length(pval)
  pval_order=order(pval)
  thr = alpha/(m-pval_order+1)
  
  df = data.frame(Index=1:m,res=(pval<=thr))
  R = df[df$res==TRUE,]$Index
  return(R)
}

## Biais Biol ------------------------------------------------------------------
# Sk={pvaleur rejete par BH par les k plus petites valeurs} -> pas de garantis statistiques 
# mais tester V* qui a des garantis 

listk = c(500*(1:10),nrow(proteom))
df_res = data.frame()
alpha=0.05
ZL = zetas.tree(C = final_tree[[2]],leaf_list = final_tree[[1]],zeta.DKWM,pval,alpha)

for (k in listk){
  pvalk = sort(pval)[1:k]
  line=order(pval)[1:k]
  Sk1 = BH_procedure(pvalk)
  FP1=sum(res_ordered[line[Sk1],]$Species=='levure')
  Sk2 = BF_procedure(pvalk)
  FP2=sum(res_ordered[line[Sk2],]$Species=='levure')
  Sk3 = HB_procedure(pvalk)
  FP3=sum(res_ordered[line[Sk3],]$Species=='levure')
  V1 = V.star(Sk1,C = final_tree[[2]], leaf_list = final_tree[[1]], ZL)
  V2 = V.star(Sk2,C = final_tree[[2]], leaf_list = final_tree[[1]], ZL)
  V3 = V.star(Sk3,C = final_tree[[2]], leaf_list = final_tree[[1]], ZL)
  
  
  Dlength=rep(c(length(Sk1),length(Sk2),length(Sk3)),2)
  Dlength[Dlength==0]<- 1 # FDP = V(S)/min(|S|,1)
  FDPk = c(V1,V2,V3,FP1,FP2,FP3)/Dlength
  
  
  df_res=rbind(df_res,data.frame(k=k,
               Méthode=c('BH','BF','HB','BH_oracle','BF_oracle','HB_oracle'),
               Cardinal=rep(c(length(Sk1),length(Sk2),length(Sk3)),2),
               Vstar=c(V1,V2,V3,FP1,FP2,FP3),
               FDP = FDPk))
}