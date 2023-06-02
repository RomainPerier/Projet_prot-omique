library(sanssouci)
library(ggplot2)
library(ggthemes)
library(stringr)
library(cgwtools)
library(reshape2)
load('save/proteom.RData')
load('save/res_ordered.RData')
load('save/final_tree.RData')
load('save/pval.RData')

##source("00_Theme.R")

##theme_set(theme_ben())

m=length(pval)
df_bornes = data.frame(Index=1:m)

## Bornes ----------------------------------------------------------------------

C=final_tree[[2]]
leaf_list=final_tree[[1]]
m=nrow(proteom)

## True Positive ---------------------------------------------------------------

TP=cumsum(res_ordered[line_sorted_by_pval,]$Species=='ups')
df_bornes <- cbind(df_bornes,TP)

## Nb de peptides par feuilles ------------------
nb_pept <- c()
for (leaf in leaf_list){nb_pept <- cbind(nb_pept,length(leaf))}

## alpha = 0.05-------------------------------------

alpha=0.05
K = 10 # Calculer tout les K lignes
ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha) 
ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval,alpha)
source("Fonction/zetas.tree.refined.R")
ZL_refined=zetas.tree.refined(C,leaf_list,zeta.DKWM,pval,alpha)

Vstar_DKWM=rep(0,m%/%K)
Vstar_HB=rep(0,m%/%K)
Vstar_refined=rep(0,m%/%K)

"for (i in 1:(m%/%K)){
  S=line_sorted_by_pval[1:(K*i)]
  Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
  Vstar_HB[i]<-V.star(S,C,ZL_HB,leaf_list)
  Vstar_refined[i] <- V.star(S,C,ZL_refined,leaf_list)
  print(i)
}"

# Vsimes : 
thr=alpha/m*(1:m)
Vsimes=sanssouci:::curveMaxFP(pval,thr)[1:(m%/%K)*K]
##Sauvegarde des bornes -----
"save(Vstar_HB,file=paste('save/Vstar_HB',as.character(alpha),'.RData', sep = ''))
save(Vstar_DKWM,file=paste('save/Vstar_DKWM',as.character(alpha),'.RData', sep = ''))
save(Vsimes,file=paste('save/Vsimes',as.character(alpha),'.RData', sep = ''))
save(Vstar_refined,file=paste('save/Vstar_refined',as.character(alpha),'.RData', sep = ''))"
     
## Calcul de la borne DKWM avec l'arbre amélioré --------
load('save/final_tree_1.RData')
C <- final_tree_1[[2]]
ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha) 
for (i in 1:(m%/%K)){
  S=line_sorted_by_pval[1:(K*i)]
  Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
  print(i)
}
save(Vstar_DKWM,file=paste('save/Vstar_DKWM',as.character(alpha),'_improved_tree.RData', sep = ''))

#Plot : 
"df <- cbind(df_bornes[1:(m%/%K)*K,],Vstar_DKWM,Vstar_HB,Vsimes,Vstar_refined)
df_plot =  cbind(df[,1:2]
                 ,TP_Vsimes=df[,'Index']-df[,'Vsimes']
                 ,TP_Vstar_DKWM=df[,'Index']-df[,'Vstar_DKWM']
                 ,TP_Vstar_HB=df[,'Index']-df[,'Vstar_HB']
                  ,TP_Vstar_refined=df[,'Index']-df[,'Vstar_refined'])

df_plot = melt(df_plot, id.vars = "Index")

save(df_plot,file=paste('save/df_plot',as.character(alpha),'.RData', sep = ''))


ggplot(df_plot,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ylim(c(0,200))+
  ggtitle(paste('Lower Bound on True Positive in Yeast Data with alpha = ', as.character(alpha), sep = ''))
"

#Plot DKWM : 
df <- cbind(df_bornes[1:(m%/%K)*K,],Vstar_DKWM)
df_plot =  cbind(df[,1:2]
                 ,TP_Vstar_DKWM=df[,'Index']-df[,'Vstar_DKWM'])

df_plot = melt(df_plot, id.vars = "Index")




ggplot(df_plot,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ylim(c(0,200))+
  ggtitle('Lower Bound on True Positive in Yeast Data')

rm(list = ls())

