library(sanssouci)
library(cgwtools)
library(ggplot2)
load('save/df_plot.RData')
load('save/res_ordered.RData')
load('save/proteom.RData')
load('save/pval.RData')
load('save/final_tree.RData')

## Hyper paramètres ------------------------------------------------------------

list_gamma = (1:9)*0.1
alpha = 0.05

## Première hybridation --------------------------------------------------------
m=nrow(proteom)
C=final_tree[[2]]
leaf_list=final_tree[[1]]
df_hybrid=data.frame()

for (gamma in list_gamma){
  # Calcul tout les K
  K = 50
  max=m%/%K
  V_hybrid = rep(0,max)
  TP_V_hybrid = rep(0,max)
  ZL_1 = zetas.tree(C,leaf_list,zeta.DKWM,pval,gamma*alpha)
  ZL_2 = zetas.tree(C,leaf_list,zeta.HB,pval,(1-gamma)*alpha)
  for (i in 1:max){
    S = line_sorted_by_pval[1:(i*K)]
    V_hybrid[i] <- min(V.star(S,C,ZL_1,leaf_list),V.star(S,C,ZL_2,leaf_list))
    TP_V_hybrid[i] <- i*K - V_hybrid[i]
  }
  df_aux <- data.frame(Index=(1:max)*K,Gamma=as.character(gamma),V_hybrid=V_hybrid,TP_V_hybrid=TP_V_hybrid)
  df_hybrid <- rbind(df_hybrid,df_aux)
}


## Deuxième hybridation --------------------------------------------------------

m=nrow(proteom)
C=final_tree[[2]]
leaf_list=final_tree[[1]]
df_hybrid2=data.frame()

for (gamma in list_gamma){
  df_aux=data.frame()
  # Calcul tout les K
  K = 50
  max=m%/%K
  V_hybrid = rep(0,max)
  TP_V_hybrid = rep(0,max)
  ZL_1 = zetas.tree(C,leaf_list,zeta.DKWM,pval,gamma*alpha)
  ZL_2 = zetas.tree(C,leaf_list,zeta.HB,pval,(1-gamma)*alpha)
  ZL=list()
  H = length(C)
  for (h in 1:H){
    ZL[[h]] <- apply(data.frame(ZL_1[[h]],ZL_2[[h]]),1,FUN=min)
  }
  for (i in 1:max){
    S = line_sorted_by_pval[1:(i*K)]
    auxV=V.star(S,C,ZL,leaf_list)
    TP_V_hybrid[i] <- i*K - auxV
    V_hybrid[i] <- auxV
  }
  df_aux <- data.frame(Index=(1:max)*K,Gamma=as.character(gamma),V_hybrid=V_hybrid,TP_V_hybrid=TP_V_hybrid)
  df_hybrid2 <- rbind(df_hybrid2,df_aux)
}

save(df_hybrid,df_hybrid2,file='save/hybrid.RData')
## Plot ------------------------------------------------------------------------

load('save/hybrid.RData')

ggplot(df_hybrid2,aes(x=Index,y=TP_V_hybrid,color=Gamma))+
  geom_line(lwd=1) +
  ylim(c(100,190)) +
  xlim(c(1000,15000))+
  ggtitle('Hybrid 2 - Yeast Data')

ggplot(df_hybrid,aes(x=Index,y=TP_V_hybrid,color=Gamma))+
  geom_line(lwd=1) +  
  ylim(c(100,190)) +
  xlim(c(1000,15000))+
  ggtitle('Hybrid 1 - Yeast Data')

df <- df_hybrid
df[,'TP_V_hybrid'] <- df_hybrid2[,'TP_V_hybrid'] - df_hybrid[,'TP_V_hybrid']

ggplot(df,aes(x=Index,y=TP_V_hybrid,color=Gamma))+
  geom_line(lwd=1) +  
  ylim(c(-1,2))+
  ylab('')+
  ggtitle('Différence entre les deux bornes hybrides')


# la 2e mieux dans certain cas 

ggplot(df,aes(x=Index,y=TP_V_hybrid))+
  geom_line(lwd=1,color='mediumblue')+
  ylab('')+
  xlim(c(10000,15000))+
  ggtitle('Comparaison des deux bornes hybrides')+
  facet_grid(Gamma~ .)
