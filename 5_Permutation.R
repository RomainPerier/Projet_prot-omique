library(sanssouci)
library(cgwtools)
library(ggplot2)
library(reshape2)
load('save/proteom.RData')
load('save/final_tree.RData')
load('save/pval.RData')
load('save/res_ordered.RData')

## Permutation -----------------------------------------------------------------
#sample(1:n,n,replace=FALSE)

# On veut permuter les valeurs des expressions dans les différentes conditions sous Ho

m_yeast = nrow(proteom[proteom$Species=='levure',])
m = nrow(proteom)


permut= sample(4:21,18,replace=FALSE)
proteom_permut=proteom
for (i in 1:m_yeast){
    proteom_permut[i,4:21] <- proteom[i,permut]
}

save(proteom_permut,file='save/permut.RData')

## Pvaleur ---------------------------------------------------------------------

Y <- as.matrix(proteom_permut[,4:21]) 

groups <- c(rep(0,9),rep(1,9))

truth <- as.numeric(proteom_permut[,'Species']=='ups')

SS_obj <- SansSouci(Y,groups,truth)

alpha = 0.05
n = 18

cal <- fit(SS_obj,alpha=0.05,B=0,family='Simes')
pval_permut = pValues(cal)

plot(ecdf(pval_permut))

plot(ecdf(pval_permut[proteom_permut$Species=='ups']))
plot(ecdf(pval_permut[proteom_permut$Species=='levure']))

line_sorted_by_pval_permut = order(pval_permut)
resave(pval_permut,line_sorted_by_pval_permut,file='save/permut.RData')

res_ordered_permut = cbind(res_ordered[,-2],Pvalue=pval)

resave(res_ordered_permut,file='save/permut.RData')

## Bornes ----------------------------------------------------------------------

load('save/permut.RData')

C=final_tree[[2]]
leaf_list=final_tree[[1]]
m=nrow(proteom_permut)

# TP :

TP=cumsum(res_ordered_permut[line_sorted_by_pval_permut,]$Species=='ups')
df_bornes = data.frame(Index=1:m)

df_bornes <- cbind(df_bornes,TP)

alpha=0.05
K = 200 # Calculer tout les K lignes
ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval_permut,alpha)
ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval_permut,alpha)
ZL_refined=zetas.tree.refined(C,leaf_list,zeta.DKWM,pval_permut,alpha)

Vstar_DKWM=rep(0,m%/%K)
Vstar_HB=rep(0,m%/%K)
Vstar_refined=rep(0,m%/%K)

for (i in 1:(m%/%K)){
  S=line_sorted_by_pval_permut[1:(K*i)]
  Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
  Vstar_HB[i]<-V.star(S,C,ZL_HB,leaf_list)
  Vstar_refined[i] <- V.star(S,C,ZL_refined,leaf_list)
  print(i)
}
# Vsimes : 
thr=alpha/m*(1:m)
Vsimes=sanssouci:::curveMaxFP(pval_permut,thr)[1:(m%/%K)*K]

df <- cbind(df_bornes[1:(m%/%K)*K,],Vstar_DKWM,Vstar_HB,Vstar_refined,Vsimes)

resave(df,file='save/permut.RData')

#Plot : 
df_plot =  cbind(df[,1:2],TP_Vsimes=df[,'Index']-df[,'Vsimes']
                     ,TP_Vstar_DKWM=df[,'Index']-df[,'Vstar_DKWM']
                     ,TP_Vstar_HB=df[,'Index']-df[,'Vstar_HB'])
                     #,TP_Vstar_refined=df[,'Index']-df[,'Vstar_refined'])

df_plot = melt(df_plot, id.vars = "Index")


ggplot(df_plot,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ggtitle('Lower Bound on True Positive in Yeast Data')