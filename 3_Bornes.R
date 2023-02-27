library(sanssouci)
library(ggplot2)
library(ggthemes)
library(stringr)
library(reshape2)
library(cgwtools)
load('save/data.RData')
#theme_set(theme_ben())

m=length(pval)
df_bornes = data.frame(Index=1:m)

## True Positive ---------------------------------------------------------------

TP=cumsum(res_ordered[line_sorted_by_pval,]$Species=='ups')

df_bornes <- cbind(df_bornes,TP)

## Vstar -----------------------------------------------------------------------
C=final_tree[[2]]
leaf_list=final_tree[[1]]
alpha=0.05

pas = 1 # -> calculer tout les pas p-values

ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha)
ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval,alpha)
ZL_multiple = zeta.multiple(C,leaf_list,zeta.DKWM,zeta.HB,pval,0.05)
ZL_refined = zetas.tree.refined(C,leaf_list,zeta.HB,pval,0.05)

iter_pointeur = 0     
Vstar_DKWM=rep(0,m%/%pas)
Vstar_HB=rep(0,m%/%pas)
Vstar_mlt=rep(0,m%/%pas)
Vstar_refined=rep(0,m%/%pas)


for (i in 4056:m%/%pas){
  S=line_sorted_by_pval[1:i*pas]
  Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
  Vstar_HB[i]<-V.star(S,C,ZL_HB,leaf_list)
  Vstar_mlt[i]<-V.star(S,C,ZL_multiple,leaf_list)
  Vstar_refined[i] <- V.star(S,C,ZL_refined,leaf_list)
  iter_pointeur=iter_pointeur+1*pas
  print(iter_pointeur)
}

df_bornes <- cbind(df_bornes,Vstar_DKWM,Vstar_HB,Vstar_mlt,Vstar_refined)

resave(df_bornes,file='save/data.RData')

## Vsimes ----------------------------------------------------------------------
thr=alpha/m*(1:m)
Vsimes=sanssouci:::curveMaxFP(pval,thr)

df_bornes <- cbind(df_bornes,Vsimes)

## Plot  -----------------------------------------------------------------------
df_plot =  cbind(df_bornes[,1:2],TP_Vsimes=df_bornes[,'Index']-df_bornes[,'Vsimes']
                 ,TP_Vstar_DKWM=df_bornes[,'Index']-df_bornes[,'Vstar_DKWM']
                 ,TP_Vstar_HB=df_bornes[,'Index']-df_bornes[,'Vstar_HB']
                 ,TP_Vstar_mlt=df_bornes[,'Index']-df_bornes[,'Vstar_mlt']
                 ,TP_Vstar_refined=df_bornes[,'Index']-df_bornes[,'Vstar_refined'])

df_plot = melt(df_plot, id.vars = "Index")

resave(df_plot,file='save/data.RData')

ggplot(df_plot,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ylim(c(0,200))+
  ggtitle('Lower Bound on True Positive in Yeast Data')
