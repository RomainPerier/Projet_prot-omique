library(sanssouci)
library(cgwtools)
library(stringr)
library(ggplot2)
library(reshape2)
load('save/final_tree.RData')
load('save/res_ordered.RData')
load('save/proteom.RData')
load('save/pval.RData')

## Sous-tableau de 500 de tailles ----------------------------------------------
N = 500

line_yeast = row.names(res_ordered[res_ordered$Species=='levure',])
n_ups = nrow(res_ordered[res_ordered$Species=='ups',])

tirage = sample(line_yeast,(N-n_ups),replace=FALSE)

res_sub = rbind(res_ordered[tirage,],res_ordered[res_ordered$Species=='ups',])

proteom_sub = proteom[row.names(res_sub),]

# Plot 
  plot(ecdf(res_sub$Pvalue))

# Structure 
  
  split_data_sub=splitBy(formula = ~ Species + Leading_razor_protein,res_sub)[]
  
  
  m=nrow(res_sub)
  n_leaf=length(split_data_sub)
  
  
  taille_leaf=rep(0,n_leaf) 
  for (i in 1:n_leaf){
    taille_leaf[i]<-nrow(split_data_sub[[i]])
  }
  taille_leaf_cum=cumsum(taille_leaf) 
  
  
  leaf_list=list(1:(taille_leaf_cum[1]))
  for (i in 2:n_leaf){
    j=i-1
    leaf_list=append(leaf_list,list((taille_leaf_cum[j]+1):(taille_leaf_cum[i])))
  }
  rm(taille_leaf,taille_leaf_cum,i,j)
  aux=!str_detect(names(split_data_sub),"ups")
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
  
  final_tree_sub=list(leaf_list,C)
  
# Bornes  
  C=final_tree_sub[[2]]
  leaf_list=final_tree_sub[[1]]
  m=nrow(proteom_sub)
  
# TP 
  pval=res_sub$Pvalue
  TP=cumsum(res_sub[order(pval),]$Species=='ups')
  df_bornes = data.frame(Index=1:m)
  
  df_bornes <- cbind(df_bornes,TP)
  
  alpha=0.05
  K = 20 # Calculer tout les K lignes
  ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha)
  ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval,alpha)
  #ZL_refined=zetas.tree.refined(C,leaf_list,zeta.DKWM,pval,alpha)
  
  Vstar_DKWM=rep(0,m%/%K)
  Vstar_HB=rep(0,m%/%K)
  #Vstar_refined=rep(0,m%/%K)
  line_sorted_by_pval_sub = order(pval)
  for (i in 1:(m%/%K)){
    S=line_sorted_by_pval_sub[1:(K*i)]
    Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
    Vstar_HB[i]<-V.star(S,C,ZL_HB,leaf_list)
    #Vstar_refined[i] <- V.star(S,C,ZL_refined,leaf_list)
    print(i)
  }
# Vsimes : 
  thr=alpha/m*(1:m)
  Vsimes=sanssouci:::curveMaxFP(pval,thr)[1:(m%/%K)*K]
  
  df <- cbind(df_bornes[1:(m%/%K)*K,],Vstar_DKWM,Vstar_HB,Vsimes)
  
#Plot : 
  df_plot =  cbind(df[,1:2],TP_Vsimes=df[,'Index']-df[,'Vsimes']
                   ,TP_Vstar_DKWM=df[,'Index']-df[,'Vstar_DKWM']
                   ,TP_Vstar_HB=df[,'Index']-df[,'Vstar_HB'])
  #,TP_Vstar_refined=df[,'Index']-df[,'Vstar_refined'])
  
  df_plot = melt(df_plot, id.vars = "Index")
  
  
  ggplot(df_plot,aes(x=Index,y=value,color=variable))+
    geom_line(lwd=1) +  
    ggtitle('Lower Bound on True Positive in Yeast Data')  
  
  # still bad 
## Sous-tableau de 2 000 de tailles --------------------------------------------
  N = 2000
  
  line_yeast = row.names(res_ordered[res_ordered$Species=='levure',])
  n_ups = nrow(res_ordered[res_ordered$Species=='ups',])
  
  tirage = sample(line_yeast,(N-n_ups),replace=FALSE)
  
  res_sub = rbind(res_ordered[tirage,],res_ordered[res_ordered$Species=='ups',])
  
  proteom_sub = proteom[row.names(res_sub),]
  
  # Plot 
  plot(ecdf(res_sub$Pvalue))
  
  # Structure 
  
  split_data_sub=splitBy(formula = ~ Species + Leading_razor_protein,res_sub)[]
  
  
  m=nrow(res_sub)
  n_leaf=length(split_data_sub)
  
  
  taille_leaf=rep(0,n_leaf) 
  for (i in 1:n_leaf){
    taille_leaf[i]<-nrow(split_data_sub[[i]])
  }
  taille_leaf_cum=cumsum(taille_leaf) 
  
  
  leaf_list=list(1:(taille_leaf_cum[1]))
  for (i in 2:n_leaf){
    j=i-1
    leaf_list=append(leaf_list,list((taille_leaf_cum[j]+1):(taille_leaf_cum[i])))
  }
  rm(taille_leaf,taille_leaf_cum,i,j)
  aux=!str_detect(names(split_data_sub),"ups")
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
  
  final_tree_sub=list(leaf_list,C)
  
# Bornes  
  C=final_tree_sub[[2]]
  leaf_list=final_tree_sub[[1]]
  m=nrow(proteom_sub)
  
# TP 
  pval=res_sub$Pvalue
  TP=cumsum(res_sub[order(pval),]$Species=='ups')
  df_bornes = data.frame(Index=1:m)
  
  df_bornes <- cbind(df_bornes,TP)
  
  alpha=0.05
  K = 20 # Calculer tout les K lignes
  ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha)
  ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval,alpha)
  #ZL_refined=zetas.tree.refined(C,leaf_list,zeta.DKWM,pval,alpha)
  
  Vstar_DKWM=rep(0,m%/%K)
  Vstar_HB=rep(0,m%/%K)
  #Vstar_refined=rep(0,m%/%K)
  line_sorted_by_pval_sub = order(pval)
  for (i in 1:(m%/%K)){
    S=line_sorted_by_pval_sub[1:(K*i)]
    Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
    Vstar_HB[i]<-V.star(S,C,ZL_HB,leaf_list)
    #Vstar_refined[i] <- V.star(S,C,ZL_refined,leaf_list)
    print(i)
  }
# Vsimes : 
  thr=alpha/m*(1:m)
  Vsimes=sanssouci:::curveMaxFP(pval,thr)[1:(m%/%K)*K]
  
  df <- cbind(df_bornes[1:(m%/%K)*K,],Vstar_DKWM,Vstar_HB,Vsimes)
  
#Plot : 
  df_plot =  cbind(df[,1:2],TP_Vsimes=df[,'Index']-df[,'Vsimes']
                   ,TP_Vstar_DKWM=df[,'Index']-df[,'Vstar_DKWM']
                   ,TP_Vstar_HB=df[,'Index']-df[,'Vstar_HB'])
  #,TP_Vstar_refined=df[,'Index']-df[,'Vstar_refined'])
  
  df_plot = melt(df_plot, id.vars = "Index")
  
  
  ggplot(df_plot,aes(x=Index,y=value,color=variable))+
    geom_line(lwd=1) +  
    ggtitle('Lower Bound on True Positive in Yeast Data')  
  
  # still bad 

## Sous-tableau de 10 000 de tailles -------------------------------------------

N = 1000
  
line_yeast = row.names(res_ordered[res_ordered$Species=='levure',])
n_ups = nrow(res_ordered[res_ordered$Species=='ups',])
  
tirage = sample(line_yeast,(N-n_ups),replace=FALSE)
  
res_sub = rbind(res_ordered[tirage,],res_ordered[res_ordered$Species=='ups',])
  
proteom_sub = proteom[row.names(res_sub),]
  
# Plot 
  plot(ecdf(res_sub$Pvalue))
  
# Structure 
  
  split_data_sub=splitBy(formula = ~ Species + Leading_razor_protein,res_sub)[]
  
  
  m=nrow(res_sub)
  n_leaf=length(split_data_sub)
  
  
  taille_leaf=rep(0,n_leaf) 
  for (i in 1:n_leaf){
    taille_leaf[i]<-nrow(split_data_sub[[i]])
  }
  taille_leaf_cum=cumsum(taille_leaf) 
  
  
  leaf_list=list(1:(taille_leaf_cum[1]))
  for (i in 2:n_leaf){
    j=i-1
    leaf_list=append(leaf_list,list((taille_leaf_cum[j]+1):(taille_leaf_cum[i])))
  }
  rm(taille_leaf,taille_leaf_cum,i,j)
  aux=!str_detect(names(split_data_sub),"ups")
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
  
  final_tree_sub=list(leaf_list,C)
  
# Bornes  
  C=final_tree_sub[[2]]
  leaf_list=final_tree_sub[[1]]
  m=nrow(proteom_sub)
  
# TP 
  pval=res_sub$Pvalue
  TP=cumsum(res_sub[order(pval),]$Species=='ups')
  df_bornes = data.frame(Index=1:m)
  
  df_bornes <- cbind(df_bornes,TP)
  
  alpha=0.05
  K = 100 # Calculer tout les K lignes
  ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha)
  ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval,alpha)
  #ZL_refined=zetas.tree.refined(C,leaf_list,zeta.DKWM,pval,alpha)
  
  Vstar_DKWM=rep(0,m%/%K)
  Vstar_HB=rep(0,m%/%K)
  #Vstar_refined=rep(0,m%/%K)
  line_sorted_by_pval_sub = order(pval)
  for (i in 1:(m%/%K)){
    S=line_sorted_by_pval_sub[1:(K*i)]
    Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
    Vstar_HB[i]<-V.star(S,C,ZL_HB,leaf_list)
    #Vstar_refined[i] <- V.star(S,C,ZL_refined,leaf_list)
    print(i)
  }
# Vsimes : 
  thr=alpha/m*(1:m)
  Vsimes=sanssouci:::curveMaxFP(pval,thr)[1:(m%/%K)*K]
  
  df <- cbind(df_bornes[1:(m%/%K)*K,],Vstar_DKWM,Vstar_HB,Vsimes)
  
#Plot : 
  df_plot =  cbind(df[,1:2],TP_Vsimes=df[,'Index']-df[,'Vsimes']
                   ,TP_Vstar_DKWM=df[,'Index']-df[,'Vstar_DKWM']
                   ,TP_Vstar_HB=df[,'Index']-df[,'Vstar_HB'])
  #,TP_Vstar_refined=df[,'Index']-df[,'Vstar_refined'])
  
  df_plot = melt(df_plot, id.vars = "Index")
  
  
  ggplot(df_plot,aes(x=Index,y=value,color=variable))+
    geom_line(lwd=1) +  
    ggtitle('Lower Bound on True Positive in Yeast Data')  
  
  # still bad 
