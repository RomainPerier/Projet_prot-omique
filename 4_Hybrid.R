library(sanssouci)
library(cgwtools)
library(ggplot2)
load('save/bornes.RData')
load('save/res_ordered.RData')
load('save/proteom.RData')
load('save/pval.RData')
#  hybride Vdkwm et Vhb 
# V*DKHM(zeta.dkwm(alpha)), et  V*HB(zeta.HB(alpha))
# plusieurs manières d'hybdrider 
#
#  varier gqamma aussi 
#
#  en théorie gamma petit n'empêche pas la perforcmance de v* 
#
# papier dbnr gamma=0.2

## Hyper paramètres ------------------------------------------------------------

list_gamma = 0.1*(1:9)
alpha = 0.05
hg
## Première hybridation --------------------------------------------------------
# 1. Vhybri = min (VDKWM(S),VHB(S))   -> mais perds garantis stats sauf si avec poids 
#     i.e. Vhybrid = min (VDKWM(zeta(gamma alpha)),VHB(zeta(1-gammaalpha))) 
m=nrow(proteom)
C=final_tree[[2]]
leaf_list=final_tree[[1]]
df_hybrid=data.frame()

for (gamma in list_gamma){
  # Calcul tout les K
  K = 100 
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
  df_aux <- data.frame(Index=(1:max),Gamma=as.character(gamma),V_hybrid=V_hybrid,TP_V_hybrid=TP_V_hybrid)
  df_hybrid <- rbind(df_hybrid,df_aux)
}


## Deuxième hybridation --------------------------------------------------------
# 2. faire un min avec les zeta ? 
#   car même R_k ! 
#     Vhybdri' = V*(min(zeta.DKWM(gamma alpha),zeta.BH(1- gamma alpha))) tjrs inférieur à la premièreb manièer -> garantis stats 
#    Egalité stricte ??? -> si oui mieux ! 
#    boucle for qui parcourt les zeta et 

m=nrow(proteom)
C=final_tree[[2]]
leaf_list=final_tree[[1]]
df_hybrid2=data.frame()

for (gamma in list_gamma){
  df_aux=data.frame()
  # Calcul tout les K
  K = 100 
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
  df_aux <- data.frame(Index=(1:max),Gamma=as.character(gamma),V_hybrid=V_hybrid,TP_V_hybrid=TP_V_hybrid)
  df_hybrid2 <- rbind(df_hybrid2,df_aux)
}

save(df_hybrid,df_hybrid2,file='save/hybrid.RData')
## Plot ------------------------------------------------------------------------

ggplot(df_hybrid2,aes(x=Index,y=TP_V_hybrid,color=Gamma))+
  geom_line(lwd=1) +  
  ylim(c(0,200))+
  ggtitle('Hybrid - Yeast Data')

ggplot(df_hybrid,aes(x=Index,y=TP_V_hybrid,color=Gamma))+
  geom_line(lwd=1) +  
  ylim(c(150,175))+
  xlim(c(0,140))
  ggtitle('Hybrid - Yeast Data')

df <- df_hybrid
df[,'TP_V_hybrid'] <- df_hybrid2[,'TP_V_hybrid'] - df_hybrid[,'TP_V_hybrid']

ggplot(df,aes(x=Index,y=TP_V_hybrid,color=Gamma))+
  geom_line(lwd=1) +  
  ylim(c(-2,2))+
  ggtitle('Diff Hybrid - Yeast Data')

# la 2e mieux dans certain cas 





