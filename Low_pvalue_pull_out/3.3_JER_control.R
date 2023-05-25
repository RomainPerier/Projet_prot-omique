library(sanssouci)
library(ggplot2)
library(cgwtools)
library(stringr)
library(reshape2)
library(doBy)
load('Low_pvalue_pull_out/save/final_tree.RData')
load('Low_pvalue_pull_out/save/proteom.RData')
load('Low_pvalue_pull_out/save/pval.RData')
load('Low_pvalue_pull_out/save/split_data.RData')
load('Low_pvalue_pull_out/save/res_ordered.RData')

## Contrôle du JER sur les feuilles  -------------------------------------------
alpha = 0.05
C=final_tree[[2]]
leaf_list=final_tree[[1]]
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

Rk_violated=unlist(df_JER_control[df_JER_control$is_violated,]$Region)
line=df_JER_control[df_JER_control$is_violated,]$Index

name_prot_violated=names(splitBy(formula = ~ Leading_razor_protein,proteom[Rk_violated,]))
prot_violated=split_data[line]

len_prop_violated=unlist(lapply(splitBy(formula = ~ Leading_razor_protein,proteom[Rk_violated,]),function(x){dim(x[1])}))
len_prot=unlist(lapply(prot_violated,function(x){dim(x[1])}))

## Suppression des peptides problèmatiques -------------------------------------

proteom_jer_proof = proteom[-Rk_violated,]
write.table(proteom_jer_proof, "Fichier/proteom_jer.txt",row.names=FALSE,quote=FALSE,sep="\t")
res=read.csv("Fichier/proteom_jer_test.txt",header=TRUE,sep="",blank.lines.skip=TRUE)


Species=!str_detect(res[,"Leading_razor_protein"],"ups")
Species[Species==TRUE]<-"levure"
Species[Species==FALSE]<-"ups"
res=cbind(res,Species)
rm(Species)

res_ordered_jer_proof=orderBy(~ Species + Leading_razor_protein,res)
proteom_jer_proof <- proteom_jer_proof[as.numeric(row.names(res_ordered)),]
rm(res)

split_data=splitBy(formula = ~ Species + Leading_razor_protein,res_ordered_jer_proof)[]

# Structure 
  
  split_data_jer_proof=splitBy(formula = ~ Species + Leading_razor_protein,res_ordered_jer_proof)[]
  
  
  m=nrow(res_ordered_jer_proof)
  n_leaf=length(split_data_jer_proof)
  
  
  taille_leaf=rep(0,n_leaf) 
  for (i in 1:n_leaf){
    taille_leaf[i]<-nrow(split_data_jer_proof[[i]])
  }
  taille_leaf_cum=cumsum(taille_leaf) 
  
  
  leaf_list=list(1:(taille_leaf_cum[1]))
  for (i in 2:n_leaf){
    j=i-1
    leaf_list=append(leaf_list,list((taille_leaf_cum[j]+1):(taille_leaf_cum[i])))
  }
  rm(taille_leaf,taille_leaf_cum,i,j)
  aux=!str_detect(names(split_data_jer_proof),"ups")
  n_levure=length(aux[aux])    
  n_ups=length(aux[!aux])    
  n_species_cum=c(n_levure,n_levure+n_ups)
  rm(aux) 
  
  C3=list(c(1,1))   
  for (i in 2:n_leaf){
    C3<-append(C3,list(c(i,i)))
  }
  C=list(list(c(1,n_leaf)),list(c(1,n_levure),c(n_levure+1,n_leaf)),C3)    
  rm(C3)
  
  final_tree_jer_proof=list(leaf_list,C)
  
## Bornes ---------------------------------------------------------------------- 
C=final_tree_jer_proof[[2]]
leaf_list=final_tree_jer_proof[[1]]
m=nrow(proteom_jer_proof)
pval=res_ordered_jer_proof$Pvalue

# TP 
TP=cumsum(res_ordered_jer_proof[order(pval),]$Species=='ups')
df_bornes = data.frame(Index=1:m)

df_bornes <- cbind(df_bornes,TP)

alpha=0.05
K = 200 # Calculer tout les K lignes
ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha)
ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval,alpha)
source("fonction/zetas.tree.refined.R")
ZL_refined=zetas.tree.refined(C,leaf_list,zeta.DKWM,pval,alpha)

Vstar_DKWM=rep(0,m%/%K)
Vstar_HB=rep(0,m%/%K)
Vstar_refined=rep(0,m%/%K)
line_sorted_by_pval_jer_proof = order(pval)
for (i in 1:(m%/%K)){
  S=line_sorted_by_pval_jer_proof[1:(K*i)]
  Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
  Vstar_HB[i]<-V.star(S,C,ZL_HB,leaf_list)
  Vstar_refined[i] <- V.star(S,C,ZL_refined,leaf_list)
  print(i)
}
# Vsimes : 
thr=alpha/m*(1:m)
Vsimes=sanssouci:::curveMaxFP(pval,thr)[1:(m%/%K)*K]

df <- cbind(df_bornes[1:(m%/%K)*K,],Vstar_DKWM,Vstar_HB,Vsimes,Vstar_refined)

#Plot : 
df_plot =  cbind(df[,1:2],TP_Vsimes=df[,'Index']-df[,'Vsimes']
                 ,TP_Vstar_DKWM=df[,'Index']-df[,'Vstar_DKWM']
                 ,TP_Vstar_HB=df[,'Index']-df[,'Vstar_HB']
                ,TP_Vstar_refined=df[,'Index']-df[,'Vstar_refined'])

df_plot = melt(df_plot, id.vars = "Index")


ggplot(df_plot,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ylim(c(0,400))+
  ggtitle('Lower Bound on True Positive in Yeast Data')  
