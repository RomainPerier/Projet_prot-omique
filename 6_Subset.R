library(sanssouci)
library(cgwtools)
library(stringr)
library(ggplot2)
library(reshape2)
load('save/final_tree.RData')
load('save/res_ordered.RData')
load('save/proteom.RData')
load('save/pval.RData')

C=final_tree[[2]]
leaf_list=final_tree[[1]]
m=length(leaf_list)

## Arbre des Levure ------------------------------------------------------------
m_yeast = length(res_ordered[!str_detect(res_ordered$Leading_razor_protein,'ups'),'Pvalue'])
first_yeast=C[[2]][[1]][1]
last_yeast =C[[2]][[1]][2]
leaf_list_yeast= leaf_list[first_yeast:last_yeast]

C3=list()
for (i in 1:length(leaf_list_yeast)){
  C3 <- append(C3,list(c(i,i)))
}

C_yeast = list(list(c(first_yeast,last_yeast)),C3)
tree_yeast = list(leaf_list_yeast,C_yeast)

res_ordered_yeast = res_ordered[!str_detect(res_ordered$Leading_razor_protein,'ups'),]
proteom_yeast = res_ordered[!str_detect(res_ordered$Leading_razor_protein,'ups'),]

save(res_ordered_yeast,proteom_yeast,m_yeast,tree_yeast,file='save/save_yeast.RData')

## Arbre UPS -------------------------------------------------------------------
m_ups = length(res_ordered[str_detect(res_ordered$Leading_razor_protein,'ups'),'Pvalue'])
# ATTENTION IL Y A UN DECALLAGE 

first_ups=C[[2]][[2]][1]
last_ups =C[[2]][[2]][2]
leaf_list_ups= leaf_list[first_ups:last_ups]

C3=list()
for (i in 1:length(leaf_list_ups)){
  leaf_list_ups[[i]] <- leaf_list_ups[[i]]-first_ups
  C3 <- append(C3,list(c(i,i)))
}

C_ups = list(list(c(first_ups-first_ups+1,last_ups-first_ups+1)),C3)
tree_ups = list(leaf_list_ups,C_ups)

rm(C3)

res_ordered_ups = res_ordered[str_detect(res_ordered$Leading_razor_protein,'ups'),]
proteom_ups = proteom[str_detect(res_ordered$Leading_razor_protein,'ups'),]

save(res_ordered_ups,proteom_ups,m_ups,tree_ups,file='save/save_ups.RData')

## Bornes sur tree_yeast -------------------------------------------------------

pval_yeast = res_ordered_yeast[,'Pvalue']
df_yeast = data.frame(Index=1:m_yeast)

# Plot : 
  plot(ecdf(pval_yeast))
  lines(c(0,1),c(0,1),col='red')

# TP :
  TP=1:m_yeast
  df_yeast <- cbind(df_yeast,TP)
# Vstar : 
    alpha=0.05
    K = 100 # Calculer tout les K lignes
    ZL_DKWM=zetas.tree(C_yeast,leaf_list_yeast,zeta.DKWM,pval,alpha)
    ZL_HB = zetas.tree(C_yeast,leaf_list_yeast,zeta.HB,pval_yeast,alpha)
    ZL_refined=zetas.tree.refined(C_yeast,leaf_list_yeast,zeta.DKWM,pval_yeast,alpha)
      
    Vstar_DKWM=rep(0,m_yeast%/%K)
    Vstar_HB=rep(0,m_yeast%/%K)
    Vstar_refined=rep(0,m_yeast%/%K)
    
    for (i in 1:m_yeast%/%K){
      S=order(pval_yeast)[1:(K*i)]
      Vstar_DKWM[i]<-V.star(S,C_yeast,ZL_DKWM,leaf_list_yeast)
      Vstar_HB[i]<-V.star(S,C_yeast,ZL_HB,leaf_list_yeast)
      Vstar_refined[i] <- V.star(S,C_yeast,ZL_refined,leaf_list_yeast)
    }
# Vsimes : 
    thr=alpha/m_yeast*(1:m_yeast)
    Vsimes=sanssouci:::curveMaxFP(pval_yeast,thr)[1:(m_yeast%/%K)*K]
    
df_yeast <- cbind(df_yeast[1:(m_yeast%/%K)*K,],Vstar_DKWM,Vstar_HB,Vstar_refined,Vsimes)

resave(df_yeast,file='save/save_yeast.RData')

#Plot : 
df_plot_yeast =  cbind(df_yeast[,1:2],TP_Vsimes=df_yeast[,'Index']-df_yeast[,'Vsimes']
                 ,TP_Vstar_DKWM=df_yeast[,'Index']-df_yeast[,'Vstar_DKWM']
                 ,TP_Vstar_HB=df_yeast[,'Index']-df_yeast[,'Vstar_HB']
                 ,TP_Vstar_refined=df_yeast[,'Index']-df_yeast[,'Vstar_refined'])

df_plot_yeast = melt(df_plot_yeast, id.vars = "Index")

resave(df_plot_yeast,file='save/save_yeast.RData')

ggplot(df_plot_yeast,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ggtitle('Lower Bound on True Positive in Yeast Data (only positive data)')

## Bornes sur tree_ups ---------------------------------------------------------

pval_ups = res_ordered_ups[,'Pvalue']
df_ups = data.frame(Index=1:m_ups)

# Plot : 
  plot(ecdf(pval_ups))
  lines(c(0,1),c(0,1),col='red')

# TP :
  TP=rep(0,m_ups)
  df_ups <- cbind(df_ups,TP)
# Vstar : 
  alpha=0.05
  K = 1 # Calculer tout les K lignes
  ZL_DKWM=zetas.tree(C_ups,leaf_list_ups,zeta.DKWM,pval,alpha)
  ZL_HB = zetas.tree(C_ups,leaf_list_ups,zeta.HB,pval_ups,alpha)
  ZL_refined=zetas.tree.refined(C_ups,leaf_list_ups,zeta.DKWM,pval_ups,alpha)

  Vstar_DKWM=rep(0,m_ups%/%K)
  Vstar_HB=rep(0,m_ups%/%K)
  Vstar_refined=rep(0,m_ups%/%K)

  for (i in 1:m_ups%/%K){
    S=order(pval_ups)[1:(K*i)]
    Vstar_DKWM[i]<-V.star(S,C_ups,ZL_DKWM,leaf_list_ups)
    Vstar_HB[i]<-V.star(S,C_ups,ZL_HB,leaf_list_ups)
    Vstar_refined[i] <- V.star(S,C_ups,ZL_refined,leaf_list_ups)
  }
# Vsimes : 
  thr=alpha/m_ups*(1:m_ups)
  Vsimes=sanssouci:::curveMaxFP(pval_ups,thr)[1:(m_ups%/%K)*K]

df_ups <- cbind(df_ups[1:(m_ups%/%K)*K,],Vstar_DKWM,Vstar_HB,Vstar_refined,Vsimes)

resave(df_ups,file='save/save_ups.RData')

#Plot : 
df_plot_ups =  cbind(df_ups[,1:2],TP_Vsimes=df_ups[,'Index']-df_ups[,'Vsimes']
                     ,TP_Vstar_DKWM=df_ups[,'Index']-df_ups[,'Vstar_DKWM']
                     ,TP_Vstar_HB=df_ups[,'Index']-df_ups[,'Vstar_HB']
                     ,TP_Vstar_refined=df_ups[,'Index']-df_ups[,'Vstar_refined'])

df_plot_ups = melt(df_plot_ups, id.vars = "Index")

resave(df_plot_ups,file='save/save_ups.RData')

ggplot(df_plot_ups,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ggtitle('Lower Bound on True Positive in Yeast Data (only false data)')