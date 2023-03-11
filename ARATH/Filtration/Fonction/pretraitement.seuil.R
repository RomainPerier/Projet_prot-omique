# Une fonction pour récupérer la liste des valeurs de notre bdd pour un prétraitement avec filtration a un seuil s fixé, si la ligne présente moins de s proportion de valeurs mesurées, alors on ne la garde pas (En plus d'aucun sets complet )

pretraitement.seuil <- function(bdd,s){
    ## Filtration -------------------------------------------------------
  
  Count_Intensity_Sets=cbind(rowSums(bdd[,(3:5)]>0),     # On compte le nombre 
                             rowSums(bdd[,(6:8)]>0),     # de valeurs mesurées 
                             rowSums(bdd[,(9:11)]>0),    # pour chaque set :
                             rowSums(bdd[,(12:14)]>0),   # S'il manque un sets
                             rowSums(bdd[,(15:17)]>0),   # complet, une de ces
                             rowSums(bdd[,(18:20)]>0))   # valeurs est à 0.
  n=nrow(bdd)
  keep = rep(TRUE,n)
  seuil = s
  for (i in 1:n){
    if(sum(bdd[i,3:20]==0)<=seuil*18 | any(Count_Intensity_Sets[i,]==0) ){keep[i]=FALSE}
  }
  
  
  bdd=bdd[keep,]
  
  ## Imputation -------------------------------------------------------
  #Pour chaque set d'échantillons, on remplace la valeur manquante par la
  #moyenne de celles mesurées.
  
  ###################################################
  ####                     Sets 1                 ### 
  ###################################################
  i=3
  j=i+2
  
  bddaux=bdd[,(i:j)]
  bddaux[bddaux==0]<-NA
  
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  
  bdd[,(i:j)]<-bddaux
  rm(bddaux,i,j,ind)
  
  ###################################################
  ####                     Sets 2                 ### 
  ###################################################
  i=6
  j=i+2
  
  bddaux=bdd[,(i:j)]
  bddaux[bddaux==0]<-NA
  
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  
  bdd[,(i:j)]<-bddaux
  rm(bddaux,i,j,ind)
  
  ###################################################
  ####                     Sets 3                 ### 
  ###################################################
  i=9
  j=i+2
  
  bddaux=bdd[,(i:j)]
  bddaux[bddaux==0]<-NA
  
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  
  bdd[,(i:j)]<-bddaux
  rm(bddaux,i,j,ind)
  
  ###################################################
  ####                     Sets 5                 ### 
  ###################################################
  i=12
  j=i+2
  
  bddaux=bdd[,(i:j)]
  bddaux[bddaux==0]<-NA
  
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  
  bdd[,(i:j)]<-bddaux
  rm(bddaux,i,j,ind)
  
  ###################################################
  ####                     Sets 6                 ### 
  ###################################################
  
  i=15
  j=i+2
  
  bddaux=bdd[,(i:j)]
  bddaux[bddaux==0]<-NA
  
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  
  bdd[,(i:j)]<-bddaux
  rm(bddaux,i,j,ind)
  
  ###################################################
  ####                     Sets 7                 ### 
  ###################################################
  
  i=18
  j=i+2
  
  bddaux=bdd[,(i:j)]
  bddaux[bddaux==0]<-NA
  
  ind=which(is.na(bddaux),arr.ind=TRUE)
  bddaux[ind] <- rowMeans(bddaux, na.rm=TRUE)[ind[,1]]
  
  bdd[,(i:j)]<-bddaux
  rm(bddaux,i,j,ind)
  
  #############################
  
  tab_data_ap_low=bdd[,3:11]
  tab_data_ap_high=bdd[,12:20]
  
  data_ap_low = data.frame(Condition="Low Spike",Value=unlist(tab_data_ap_low,use.names=FALSE))
  data_ap_high = data.frame(Condition="High Spike",Value=unlist(tab_data_ap_high,use.names=FALSE))
  
  
  data_ap=data.frame(Méthode=paste('imputation avec seuil :',s),
                           rbind(data_ap_low,data_ap_high))
  
  return(data_ap)
}