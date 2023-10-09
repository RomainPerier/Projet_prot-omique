bdd=read.csv("Fichier/peptides_YEASTUPS.txt",header=TRUE,sep="",blank.lines.skip=TRUE)
library(sanssouci)
library(ggplot2)
library(ggthemes)
library(stringr)
library(reshape2)
library(DAPAR)
library(cherry)

## Dev sanssouci ---- 
zeta.HB.no.extension <- function(pval, lambda) {
  m <- length(pval)
  sorted.pval <- sort(pval)
  
  thresholds <- lambda / (m:1)
  v <- sorted.pval - thresholds
  indexes <- which(v > 0)
  if (! length(indexes)) {
    return(0)
  }
  else{
    return(m - indexes[1] + 1)
  }
}
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
zetas.tree.no.extension <- function(C, leaf_list, method, pvalues, alpha, refine=FALSE, verbose=FALSE) {
  H <- length(C)
  K <- nb.elements.no.extension(C)
  ZL <- list()
  new_K <- K
  continue <- TRUE
  nb_loop <- 0
  while (continue) {
    usage_K <- new_K
    new_K <- K
    for (h in H:1) {
      Ch <- C[[h]]
      len <- length(Ch)
      zeta_inter <- numeric(len)
      for (j in 1:len) {
        Chj <- Ch[[j]]
        pvals <- pvalues[unlist(leaf_list[Chj[1]:Chj[2]])]
        if (typeof(method) == "list") {
          if (typeof(method[[h]]) == "list") {
            zeta_method <- method[[h]][[j]]
          } else {
            zeta_method <- method[[h]]
          }
        } else {
          zeta_method <- method
        }
        zeta_inter[j] <- zeta_method(pvals, alpha / usage_K)
        if (refine && (zeta_inter[j] == 0) )
          new_K <- new_K - 1
      }
      ZL[[h]] <- zeta_inter
    }
    if (verbose) {
      nb_loop <- nb_loop + 1
      print(paste0("loop number=", nb_loop,", usage_K=",usage_K,", new_K=",new_K))
    }
    continue <- refine && (new_K < usage_K)
  }
  return(ZL)
}
nb.elements.no.extension <- function(C) {
  H <- length(C)
  count <- 0
  for (h in H:1) {
    count <- count + length(C[[h]])
  }
  return(count)
}

## Functions proteomic --------------

preprocessing <- function(bdd,type){
  vectcond=lapply(colnames(bdd), function(i) strsplit(i,"[_]")[[1]][3])
  vectcond[is.na(vectcond)]<-FALSE
  set1 <- vectcond=="0.5fmol"
  set2 <- vectcond=="1fmol"
  set3 <- vectcond=="2.5fmol"
  set4 <- vectcond=="5fmol"
  set5 <- vectcond=="10fmol"
  set6 <- vectcond=="25fmol"
  Count_Intensity_Sets=cbind(rowSums(bdd[,set1]>0),    
                             rowSums(bdd[,set2]>0),    
                             rowSums(bdd[,set3]>0),   
                             rowSums(bdd[,set4]>0),   
                             rowSums(bdd[,set5]>0),   
                             rowSums(bdd[,set6]>0))   
  nb_col <- sum(set1,set2,set3,set4,set5,set6)
  Count_Intensity=apply(Count_Intensity_Sets,1,min)
  Count_Intensity_bis=apply(Count_Intensity_Sets,1,max)
  if (type == "strict"){
    bdd<-bdd[rowSums(Count_Intensity_Sets) == nb_col]
    return(bdd)
  }
  else if (type == "wild"){
    bdd<-bdd[Count_Intensity!=0,]
  }
  else if (type == "nwild"){
    bdd<-bdd[Count_Intensity!=0 & Count_Intensity_bis>1,]
  }
  bddaux<-bdd[,set1]
  bddaux[bddaux==0]<-NA
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  bdd[,set1]<-bddaux
  
  bddaux<-bdd[,set2]
  bddaux[bddaux==0]<-NA
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  bdd[,set2]<-bddaux
  
  bddaux<-bdd[,set3]
  bddaux[bddaux==0]<-NA
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  bdd[,set3]<-bddaux
  
  bddaux<-bdd[,set4]
  bddaux[bddaux==0]<-NA
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  bdd[,set4]<-bddaux
  
  bddaux<-bdd[,set5]
  bddaux[bddaux==0]<-NA
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  bdd[,set5]<-bddaux
  
  bddaux<-bdd[,set6]
  bddaux[bddaux==0]<-NA
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  bdd[,set6]<-bddaux
  
  return(bdd)
}

get_p_values_ss <- function(bdd,alternative = "two.sided"){
  library(sanssouci)
  vectcond=lapply(colnames(bdd), function(i) strsplit(i,"[_]")[[1]][3])
  vectcond[is.na(vectcond)]<-FALSE
  set1 <- vectcond=="0.5fmol"
  set2 <- vectcond=="1fmol"
  set3 <- vectcond=="2.5fmol"
  set4 <- vectcond=="5fmol"
  set5 <- vectcond=="10fmol"
  set6 <- vectcond=="25fmol"
  categ <- rep(c(0,1), times = c(3,3))
  prot_matr <- data.matrix(cbind(rowMeans(cbind(bdd[,set1][1],bdd[,set2][1],bdd[,set3][1])),
                           rowMeans(cbind(bdd[,set1][2],bdd[,set2][2],bdd[,set3][2])),
                           rowMeans(cbind(bdd[,set1][3],bdd[,set2][3],bdd[,set3][3])),
                           rowMeans(cbind(bdd[,set4][1],bdd[,set5][1],bdd[,set6][1])),
                           rowMeans(cbind(bdd[,set4][2],bdd[,set5][2],bdd[,set6][2])),
                           rowMeans(cbind(bdd[,set4][3],bdd[,set5][3],bdd[,set6][3]))))
  tests<- rowWelchTests(prot_matr, categ, alternative = alternative)
  Pvalue <- tests$p.value
  return(cbind(bdd[,(1:2)],Pvalue))
}

get_p_values_DAPAR_mean <- function(bdd,comp.type){
  vectcond=lapply(colnames(bdd), function(i) strsplit(i,"[_]")[[1]][3])
  vectcond[is.na(vectcond)]<-FALSE
  set1 <- vectcond=="0.5fmol"
  set2 <- vectcond=="1fmol"
  set3 <- vectcond=="2.5fmol"
  set4 <- vectcond=="5fmol"
  set5 <- vectcond=="10fmol"
  set6 <- vectcond=="25fmol"
  prot_matr <- data.matrix(cbind(Low_spike1 = rowMeans(cbind(bdd[,set1][1],bdd[,set2][1],bdd[,set3][1])),
                                 Low_spike2 = rowMeans(cbind(bdd[,set1][2],bdd[,set2][2],bdd[,set3][2])),
                                 Low_spike3 = rowMeans(cbind(bdd[,set1][3],bdd[,set2][3],bdd[,set3][3])),
                                 High_spike1 = rowMeans(cbind(bdd[,set4][1],bdd[,set5][1],bdd[,set6][1])),
                                 High_spike1 = rowMeans(cbind(bdd[,set4][2],bdd[,set5][2],bdd[,set6][2])),
                                 High_spike1 = rowMeans(cbind(bdd[,set4][3],bdd[,set5][3],bdd[,set6][3]))))
  Tab_prot <- data.frame(cbind(Sample.name = c("Low_spike1","Low_spike2","Low_spike3","High_spike1","High_spike2","High_spike3"),Condition = c("Low_spike","Low_spike","Low_spike","High_spike","High_spike","High_spike"),Bio.Rep = 1:6))
  limma <- limmaCompleteTest(prot_matr,Tab_prot,comp.type)
  Pvalue <- limma[[2]][,1]
  return(cbind(bdd[,1:2],Pvalue))
}

get_forest_and_sorted_pvalues <- function(pvalues){
  library(stringr)
  Species <- !str_detect(pvalues[,"Leading_razor_protein"],"ups")
  Species[Species==TRUE] <-"yeast"
  Species[Species==FALSE] <-"ups"
  pvalues <- cbind(pvalues,Species)
  res_ordered <- pvalues[order(pvalues[,colnames(pvalues)=="Species"]),]
  n_peptid_ups <- sum(Species == "ups") 
  ord <- c(order(res_ordered[res_ordered$Species=="ups","Leading_razor_protein"]), n_peptid_ups + order(res_ordered[res_ordered$Species=="yeast","Leading_razor_protein"]))
  res_ordered <- res_ordered[ord,]
  ordered_species <- res_ordered$Species
  m <- nrow(res_ordered)
  i <- 1
  leaf_list <- list()
  n_ups <- 0
  n_leaves <- 0
  while (i<length(Species)){
    j <- i
    while(res_ordered[j,"Leading_razor_protein"] == res_ordered[j+1,"Leading_razor_protein"] && j < length(Species)){
      j <- j+1
    }
    leaf_list <- append(leaf_list,list(i:j))
    i<-j+1
    n_leaves <- n_leaves + 1
    if (ordered_species[j]=="ups"){
      n_ups <- n_ups + 1
    }
  }
  C3=list()   
  for (i in 1:n_leaves){
    C3 <- append(C3,list(c(i,i)))
  }
  C <- list(list(c(1,n_leaves)),list(c(1,n_ups),c(n_ups+1,n_leaves)),C3)    
  return(list(leaf_list,C,res_ordered))
}

get_p_values_DAPAR_nomean <- function(bdd,comp.type){
  vectcond=lapply(colnames(bdd), function(i) strsplit(i,"[_]")[[1]][3])
  vectcond[is.na(vectcond)]<-FALSE
  set1 <- vectcond=="0.5fmol"
  set2 <- vectcond=="1fmol"
  set3 <- vectcond=="2.5fmol"
  set4 <- vectcond=="5fmol"
  set5 <- vectcond=="10fmol"
  set6 <- vectcond=="25fmol"
  prot_matr <- data.matrix(cbind(bdd[,set1],bdd[,set2],bdd[,set3],bdd[,set4],bdd[,set5],bdd[,set6]))
  Tab_prot <- data.frame(cbind(Sample.name = c(rep("0.5fmol",3),rep("1fmol",3),rep("2.5fmol",3),rep("5fmol",3),rep("10fmol",3),rep("25fmol",3)),Condition = c(rep(1,9),rep(2,9)),Bio.Rep = 1:18))
  limma <- limmaCompleteTest(prot_matr,Tab_prot,comp.type)
  Pvalue <- limma[[2]][,1]
  return(cbind(bdd[,1:2],Pvalue))
}

cdf_p_values <- function(pvalues,bdd){
  p_values <- pvalues[,3]
  m=nrow(p_values)
  Species <- !str_detect(new_bdd[,"Leading_razor_protein"],"ups")
  Species[Species==TRUE] <- "H0"
  Species[Species==FALSE] <- "H1"
  df=data.frame(Pvalue =p_values,Condition=Species)
  line_sorted_by_pval <- order(p_values)
  df <- df[line_sorted_by_pval,]
  ggplot(df,aes(x=Pvalue,color=Condition))+
    stat_ecdf(geom='step',lwd=1)+
    xlab('')+ylab('')+
    ggtitle('Fonction de réparition empirique des p-valeurs selon la condition')+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),show.legend = NA,lwd=1,linetype='dashed',color='black')
}

get_envelop <- function(p_values,alpha){
  L <- get_forest_and_sorted_pvalues(p_values)
  sorted_p_values <- L[[3]]
  pval <- sorted_p_values[,3]
  leaf_list <- L[[1]]
  forest <- L[[2]]
  line_sorted_by_pval <- order(pval)
  m <- length(line_sorted_by_pval)


  ZL_DKWM <- zetas.tree.no.extension(forest,leaf_list,zeta.DKWM,pval,alpha/2) 
  ZL_HB <- zetas.tree.no.extension(forest,leaf_list,zeta.HB.no.extension,pval,alpha/2)
  ZL_DKWM_bis <- zetas.tree.no.extension(forest,leaf_list,zeta.DKWM,1-pval,alpha/2) 
  ZL_HB_bis <- zetas.tree.no.extension(forest,leaf_list,zeta.HB.no.extension,1-pval,alpha/2)
  Vstar_DKWM <- curve.V.star.forest.fast(line_sorted_by_pval,forest,ZL_DKWM,leaf_list)
  print("DKWM done")
  Vstar_HB <- curve.V.star.forest.fast(line_sorted_by_pval,forest,ZL_HB,leaf_list)
  print("HB done")
  Vstar_DKWM_bis <- curve.V.star.forest.fast(line_sorted_by_pval,forest,ZL_DKWM_bis,leaf_list)
  print("DKWM_bis done")
  Vstar_HB_bis <- curve.V.star.forest.fast(line_sorted_by_pval,forest,ZL_HB_bis,leaf_list)
  print("HB_bis done")
  
  thr <- alpha/(m*2)*(1:m)
  Vsimes <- sanssouci:::curveMaxFP(pval,thr)
  Vsimes_bis <- sanssouci:::curveMaxFP(1-pval,thr)
  print("Simes done")
  
  TP <- cumsum(sorted_p_values[line_sorted_by_pval,"Species"]=='ups')
  df_bornes <- data.frame(Index=1:m)
  df_bornes <- cbind(df_bornes,TP)
  df <- cbind(df_bornes,Vstar_DKWM,Vstar_DKWM_bis,Vstar_HB,Vstar_HB_bis,Vsimes,Vsimes_bis)
  df_plot =  cbind(df[,1:2]
                   ,TP_Vsimes=df[,'Index']-df[,'Vsimes']
                   ,TP_Vstar_DKWM=df[,'Index']-df[,'Vstar_DKWM']
                   ,TP_Vstar_HB=df[,'Index']-df[,'Vstar_HB']                   
                   ,TP_Vsimes_sup=df[,'Vsimes_bis']
                   ,TP_Vstar_DKWM_sup=df[,'Vstar_DKWM_bis']
                   ,TP_Vstar_HB_sup=df[,'Vstar_HB_bis'])
  
  df_plot = melt(df_plot, id.vars = "Index")
  ggplot(df_plot,aes(x=Index,y=value,color=variable))+
    geom_line(lwd=1) +  
    ylim(c(0,5500))+
    ggtitle('Lower Bound on True Positive in Yeast Data')
  
}

get_bounds <- function(p_values,alpha,methods){
  L <- get_forest_and_sorted_pvalues(p_values)
  sorted_p_values <- L[[3]]
  pval <- sorted_p_values[,3]
  leaf_list <- L[[1]]
  forest <- L[[2]]
  line_sorted_by_pval <- order(pval)
  m <- length(line_sorted_by_pval)
  Oracle <-cumsum(sorted_p_values[line_sorted_by_pval,"Species"]=='ups')
  df <- data.frame(Index=1:m)
  df <- cbind.data.frame(df,Oracle)
  for (method in methods){
    print(method)
    if (method == "simes"){
      thr <- alpha/m*(1:m)
      Vsimes <- sanssouci:::curveMaxFP(pval,thr)
      df <- cbind(df,TP_simes=df[,'Index']-Vsimes)
    }
    else if (method == "cherry"){
      hom <- hommelFast(pval,simes=TRUE)
      TP_cherry <- sapply(1:m, FUN=function(i) pickSimes(hom, select= 1:i))
      df <- cbind(df,TP_cherry)
    }
    else {
      fun <- get(method)
      ZL <- zetas.tree.no.extension(forest,leaf_list,fun,pval,alpha) 
      Vstar <- curve.V.star.forest.fast(line_sorted_by_pval,forest,ZL,leaf_list)
      df <- cbind(df,df[,'Index'] - Vstar)
      colnames(df)[length(df)]<- paste0("TP_",strsplit(method,"[.]")[[1]][2])
    }
  }
  print("ok")
  return(df)
}

plot_bounds<- function(df_bounds,xbound,ybound){
  df_plot = melt(df_bounds, id.vars = "Index")
  ggplot(df_plot,aes(x=Index,y=value,color=variable))+
    geom_line(lwd=1) +  
    ylim(c(0,ybound))+
    ggtitle('Lower Bound on True Positive')+
    ylab("Lower confidence envelope on the number of true positives") + 
    xlab("Hypotheses sorted by p-value") + 
    xlim(c(0,xbound))+
    labs(color = "Method")
}

##Appli -----------------

new_bdd <- preprocessing(bdd,"nwild")

p_values <- get_p_values_ss(new_bdd,alternative = "greater")
cdf_p_values(p_values,new_bdd)

bounds <- get_bounds(p_values,0.05,list("cherry","simes","zeta.DKWM"))
plot_bounds(bounds,500,200)
plot_bounds(bounds,length(bounds[,1]),200)



#rm(list = ls())
