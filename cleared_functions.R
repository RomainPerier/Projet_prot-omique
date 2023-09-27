bdd=read.csv("Fichier/peptides_YEASTUPS.txt",header=TRUE,sep="",blank.lines.skip=TRUE)
library(ggplot2)

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

get_p_values_DAPAR <- function(bdd,comp.type){
  library(DAPAR)
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
  limmaCompleteTest(prot_matr,Tab_prot,comp.type)
}

get_forest <- function(pvalues){
  library(stringr)
  Species <- !str_detect(pvalues[,"Leading_razor_protein"],"ups")
  Species[Species==TRUE] <-"yeast"
  Species[Species==FALSE] <-"ups"
  n_yeast <- sum(Species == "yeast")    
  n_ups <- sum(Species == "ups")    
  pvalues <- cbind(pvalues,Species)
  res_ordered <- pvalues[order(pvalues[,colnames(pvalues)=="Species"]),]
  ord <- c(order(res_ordered[res_ordered$Species=="ups","Leading_razor_protein"]), n_ups + order(res_ordered[res_ordered$Species=="yeast","Leading_razor_protein"]))
  res_ordered <- res_ordered[ord,]
  m <- nrow(res_ordered)
  i <- 1
  leaf_list <- list()
  while (i<length(Species)){
    j <- i
    print(i)
    while(res_ordered[j,"Leading_razor_protein"] == res_ordered[j+1,"Leading_razor_protein"] && j < length(Species)){
      j <- j+1
    }
    leaf_list <- append(leaf_list,list(c(i,j)))
    i<-j+1
  }
  n_leaf <- n_yeast + n_ups
  C3=list()   
  for (i in 1:n_leaf){
    C3 <- append(C3,list(c(i,i)))
  }
  C <- list(list(c(1,n_leaf)),list(c(1,n_ups),c(n_ups+1,n_yeast + n_ups)),C3)    
  return(list(leaf_list,C))
}

cdf_p_values <- function(p_values,bdd){
  m=nrow(res_ordered)
  Species <- !str_detect(new_bdd[,"Leading_razor_protein"],"ups")
  Species[Species==TRUE] <- "H0"
  Species[Species==FALSE] <- "H1"
  df=data.frame(Pvalue=p_values,Condition=Species)
  df <- df[line_sorted_by_pval,]
  ggplot(df,aes(x=Pvalue,color=Condition))+
    stat_ecdf(geom='step',lwd=1)+
    xlab('')+ylab('')+
    ggtitle('Fonction de réparition empirique des p-valeurs selon la condition')+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),show.legend = NA,lwd=1,linetype='dashed',color='black')
}



new_bdd <- preprocessing(bdd,"wild")
limma <- get_p_values_DAPAR(new_bdd,'OnevsAll')
p_values <- limma[[2]]
p_values <- p_values[,1]
cdf_p_values(p_values,new_bdd)


p_values <- get_p_values_ss(new_bdd,alternative = "two.sided")$Pvalue
cdf_p_values(p_values,new_bdd)
