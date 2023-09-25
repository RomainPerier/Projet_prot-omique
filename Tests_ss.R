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

## Dev sanssouci ---- 
forest.completion <- function(C, ZL, leaf_list) {
  H <- length(C)
  
  leaves.to.place <- 1:length(leaf_list)
  len.to.place <- length(leaves.to.place)
  to.delete <- numeric(0)
  
  for (h in 1:H) {
    
    j <- 1
    l <- 1
    
    while (j <= len.to.place) {
      expected_leaf <- leaves.to.place[j]
      end_of_line <- l > length(C[[h]])
      if (! end_of_line) {Chl <- C[[h]][[l]]}
      if ((! end_of_line) && expected_leaf == Chl[[1]]) {
        if (expected_leaf == Chl[[2]]) {
          to.delete <- c(to.delete, j)
        }
        j <- j + Chl[2] - Chl[1] +1
      } else {
        C[[h]] <- append(C[[h]], list(c(expected_leaf, expected_leaf)), l - 1)
        ZL[[h]] <- append(ZL[[h]], length(leaf_list[[expected_leaf]]), l - 1)
        to.delete <- c(to.delete, j)
        j <- j + 1
      }
      l <- l + 1
    }
    
    leaves.to.place <- leaves.to.place[-to.delete]
    len.to.place <- length(leaves.to.place)
    to.delete <- numeric(0)
    
  }
  
  if (len.to.place > 0) {
    h <- H + 1
    C[[h]] <- list()
    ZL[[h]] <- numeric(0)
    for (expected_leaf in leaves.to.place) {
      C[[h]] <- append(C[[h]], list(c(expected_leaf, expected_leaf)))
      ZL[[h]] <- append(ZL[[h]], length(leaf_list[[expected_leaf]]))
    }
  }
  
  return(list(C = C, ZL = ZL))
}
V.star.no.extension <- function(S, C, ZL, leaf_list) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  Vec <- numeric(nb_leaves) 
  for (i in 1:nb_leaves) {
    Vec[i] <- sum(S %in% leaf_list[[i]])
  }
  # the initialization term for each atom P_i
  # is equivalent to completing the family if it isn't,
  # assuming that leaf_list does indeed contain all leaves
  # and some were just eventually missing in C and ZL
  for (h in H:1) {
    nb_regions <- length(C[[h]])
    if (nb_regions>0) {
      for (j in 1:nb_regions) {
        Chj <- C[[h]][[j]]
        if (Chj[1]==Chj[2]) { # means this is an atom, no need to compute 
          # len_inter given that we already did it during initialization,
          # furthermore there are no successors
          len_inter <- Vec[Chj[1]]
          res <- min(ZL[[h]][j], len_inter)
        } else {
          region_vector <- unlist(leaf_list[Chj[1]:Chj[2]])
          len_inter <- sum(S %in% region_vector)
          sum_succ <- sum(Vec[Chj[1]:Chj[2]]) 
          res <- min(ZL[[h]][j], len_inter, sum_succ)
        }
        Vec[Chj[1]:Chj[2]] <- 0
        Vec[Chj[1]] <- res
      }
    }
  }
  return(sum(Vec))
}


curve.V.star.forest.fast <- function(perm, C, ZL, leaf_list, is.pruned = FALSE, is.complete = FALSE, pruning = FALSE, delete.gaps = FALSE){
  
  vstars <- numeric(length(perm))
  
  if (! is.pruned) {
    if (! is.complete) {
      # the fast version needs a proper completion of the
      # forest structure, and for the same reason
      # it must not use super pruning
      completed <- forest.completion(C, ZL, leaf_list)
      C <- completed$C
      ZL <- completed$ZL
    }
    
    if (pruning) {
      is.pruned <- TRUE
      pruned <- pruning(C, ZL, leaf_list, super.prune = FALSE, delete.gaps = delete.gaps)
      C <- pruned$C
      ZL <- pruned$ZL
      m <- length(unlist(leaf_list))
      
      if (length(perm) == m) {
        # means that length(perm) = m,
        # but the pruning already computed
        # V^*({1, ..., m}) as a by-product so we
        # might as well use it:
        vstars[m] <- pruned$VstarNm
        perm <- perm[-m]
      }
      
    }
  }
  
  
  H <- length(C)
  
  etas <- ZL
  K.minus <- list()
  for (h in 1:H){
    etas[[h]] <- rep(0, length(ZL[[h]]))
    K.minus[[h]] <- list()
    if (length(ZL[[h]]) > 0){
      for (j in 1:length(ZL[[h]])){
        if (ZL[[h]][j] == 0){
          K.minus[[h]][[j]] <- C[[h]][[j]]
        }
      }
    }
  }
  
  for (t in 1:length(perm)) {
    
    i.t <- perm[t]
    if (t > 1) {
      previous.vstar <- vstars[t - 1]
    } else {
      previous.vstar <- 0
    }
    
    ################################
    # SEARCHING IF i_t IS IN K MINUS
    # if so, go.next == TRUE
    # and we just go next to step t+1
    go.next <- FALSE
    for (h in 1:H) {
      if (go.next) {
        break
      }
      for (couple in K.minus[[h]]) {
        if (! is.null(couple)) {
          lower_leaf <- leaf_list[[couple[1]]]
          lower_hyp <- lower_leaf[1]
          upper_leaf <- leaf_list[[couple[2]]]
          upper_hyp <- upper_leaf[length(upper_leaf)]
          if ((i.t >= lower_hyp) && (i.t <= upper_hyp)) {
            go.next <- TRUE
            # print(paste0(i.t, " is in K minus"))
            break
          }
        }
      }
    }
    # print(paste0(i.t, " isn't in K minus"))
    #########################################
    
    # COMPUTING V.STAR AND UPDATING K.MINUS AND ETAS
    ################################################
    if (go.next) {
      vstars[t] <- previous.vstar
    } else {
      # Here, i_t isn't in K minus
      for (h in 1:H) {
        nb_regions <- length(C[[h]])
        if(nb_regions > 0){
          is.found <- FALSE
          for (j in 1:nb_regions) {
            couple <- C[[h]][[j]]
            lower_leaf <- leaf_list[[couple[1]]]
            lower_hyp <- lower_leaf[1]
            upper_leaf <- leaf_list[[couple[2]]]
            upper_hyp <- upper_leaf[length(upper_leaf)]
            if((i.t >= lower_hyp) && (i.t <= upper_hyp)){
              # we found k^{(t,h)}
              is.found <- TRUE
              break
            }
          }
          if (! is.found) {
            next
          }
          etas[[h]][[j]] <- etas[[h]][[j]] + 1
          if(etas[[h]][[j]] < ZL[[h]][[j]]){
            # pass
          } else {
            K.minus[[h]][[j]] <- C[[h]][[j]]
            break
          }
        }
      }
      vstars[t] <- previous.vstar + 1
    }
    ################################################
    
  }
  return(vstars)
}

m=length(pval)
df_bornes = data.frame(Index=1:m)


## Bornes ----------------------------------------------------------------------

C=final_tree[[2]]
leaf_list=final_tree[[1]]
m=nrow(proteom)

## True Positive ---------------------------------------------------------------

TP=cumsum(res_ordered[line_sorted_by_pval,]$Species=='ups')

df_bornes <- cbind(df_bornes,TP)

alpha=0.05
K = 10 # Calculer tout les K lignes
ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha) 

ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval,alpha)
source("Fonction/zetas.tree.refined.R")
ZL_refined=zetas.tree.refined(C,leaf_list,zeta.DKWM,pval,alpha)

Vstar_DKWM=rep(0,m%/%K)
Vstar_HB=rep(0,m%/%K)
Vstar_refined=rep(0,m%/%K)

for (i in 1:(m%/%K)){
  S=line_sorted_by_pval[1:(K*i)]
  Vstar_DKWM[i]<-V.star.no.extension(S,C,ZL_DKWM,leaf_list)
  Vstar_HB[i]<-V.star.no.extension(S,C,ZL_HB,leaf_list)
  Vstar_refined[i] <- V.star.no.extension(S,C,ZL_refined,leaf_list)
  print(i)
}
# Vsimes : 
thr=alpha/m*(1:m)
Vsimes=sanssouci:::curveMaxFP(pval,thr)[1:(m%/%K)*K]

df <- cbind(df_bornes[1:(m%/%K)*K,],Vstar_DKWM,Vstar_HB,Vsimes,Vstar_refined)

#Plot : 
df_plot =  cbind(df[,1:2]
                 ,TP_Vsimes=df[,'Index']-df[,'Vsimes']
                 ,TP_Vstar_DKWM=df[,'Index']-df[,'Vstar_DKWM']
                 ,TP_Vstar_HB=df[,'Index']-df[,'Vstar_HB'])

df_plot = melt(df_plot, id.vars = "Index")

save(df_plot,file='save/df_plot_ss_verif.RData')


ggplot(df_plot,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ylim(c(0,200))+
  ggtitle('Lower Bound on True Positive')

ggplot(df_plot,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  xlim(c(0,1000))+
  ylim(c(0,200))+
  ggtitle('Lower Bound on True Positive')

Vstar_DKWM1=rep(0,m%/%K)
Vstar_HB1=rep(0,m%/%K)
Vstar_refined1=rep(0,m%/%K)

for (i in 1:(m%/%K)){
  S=line_sorted_by_pval[1:(K*i)]
  Vstar_DKWM1[i]<-V.star(S,C,ZL_DKWM,leaf_list)
  Vstar_HB1[i]<-V.star(S,C,ZL_HB,leaf_list)
  Vstar_refined1[i] <- V.star(S,C,ZL_refined,leaf_list)
  print(i)
}
nb_err1 = sum(Vstar_DKWM1!=Vstar_DKWM)
nb_err2 = sum(Vstar_HB1!=Vstar_HB)
nb_err3 = sum(Vstar_refined1!=Vstar_refined)

Vstar_DKWM2=curve.V.star.forest.fast(line_sorted_by_pval,C,ZL_DKWM,leaf_list)
Vstar_HB2=curve.V.star.forest.fast(line_sorted_by_pval,C,ZL_HB,leaf_list)
Vstar_refined2=curve.V.star.forest.fast(line_sorted_by_pval,C,ZL_refined,leaf_list)

nb_err1_bis <- 0
nb_err2_bis <- 0
nb_err3_bis <- 0

for (i in 1:(m%/%K)){
  if (Vstar_DKWM2[i*K]!=Vstar_DKWM[i]){
    nb_err1_bis <- nb_err1_bis + 1
  }
  if (Vstar_HB2[i*K]!=Vstar_HB[i]){
    nb_err2_bis <- nb_err2_bis + 1
  }
  if (Vstar_refined2[i*K]!=Vstar_refined[i]){
    nb_err3_bis <- nb_err3_bis + 1
  }
} 

curve.boucle.DKWM <- function(mlim,pval_sorted,C,leaf_list) {
  Vstar_DKWM <- rep(0,mlim)
  for (i in 1:mlim){
    Vstar_DKWM[i]<-V.star(pval_sorted[1:i],C,ZL_DKWM,leaf_list)
  }
}

t0 <- system.time()
t1 <- system.time()
tdif <- difftime(t0,t1)
