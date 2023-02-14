library(sanssouci)
library(doBy)
library(stringr)

## Prétraitement  -------------------------------------------------------
bdd=read.csv("peptides_YEASTUPS.txt",header=TRUE,sep="",blank.lines.skip=TRUE)


Count_Intensity_Sets=cbind(rowSums(bdd[,(3:5)]>0),#On compte le nombre 
                           rowSums(bdd[,(6:8)]>0),     #de valeurs mesurée 
                           rowSums(bdd[,(9:11)]>0),    #pour chaque set :
                           rowSums(bdd[,(12:14)]>0),   #S'il manque un sets
                           rowSums(bdd[,(15:17)]>0),   #complet, une de ces
                           rowSums(bdd[,(18:20)]>0))   #valeurs est à 0.

Count_Intensity=apply(Count_Intensity_Sets,1,sum)

bdd_no_na=bdd[Count_Intensity==18,]

rm(Count_Intensity,Count_Intensity_Sets)


write.table(bdd_no_na, "proteom2_treated_na.txt",row.names=FALSE,quote=FALSE,sep="\t")



## Groupe ------------------------------------------------------------------------------------

res=read.csv("proteom2_treated_na_test.txt",header=TRUE,sep="\t")[,c('id','Pvalue','Leading_razor_protein')]

Species=!str_detect(res[,"Leading_razor_protein"],"ups") 
Species[Species==TRUE]<-"levure"
Species[Species==FALSE]<-"ups"
res=cbind(res,Species)
rm(Species)

res_ordered=orderBy(~ Species + Leading_razor_protein,res)

save(res_ordered,file='res_ordered_na.RData')


split_data=splitBy(formula = ~ Species + Leading_razor_protein,res_ordered)[]
save(split_data, file = "split_data_na.RData")


## Hiérarchie -------------------------------------------------------------------------------

## Les fichiers de données -------------------------------------------------------------

load("split_data_na.RData")
load("res_ordered_na.RData") 


n_leaf=length(split_data)

taille_leaf=rep(0,n_leaf) 
for (i in 1:n_leaf){
  taille_leaf[i]<-dim(split_data[[i]])[1]
}

taille_leaf_cum=cumsum(taille_leaf) 

leaf_list=list(1:(taille_leaf_cum[1]))   

for (i in 2:(n_leaf)){
  j=i-1
  leaf_list=append(leaf_list,list((taille_leaf_cum[j]+1):(taille_leaf_cum[i])))
}

aux=!str_detect(names(split_data),"ups")  
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

final_tree=list(leaf_list,C)

save(final_tree,file="final_tree_na.RData")

## Vstar test -------------------------------------------------------------------------

# Trie par p-valeurs -----------------------------------------------------------
m=length(pval)

aux <- res_ordered

rownames(aux)<-1:m  #Le nom de la ligne dans aux correspondra à sa position dans res_ordered et ON NE TOUCHE JAMAIS A RES_ORDERED 

aux=orderBy(~ Pvalue,aux)  #En triant, on récupère les lignes de aux trié par p-val mais donc les positions dans res_ordered dans l'ordre 
# typiquement res_ordered[as.numeric(rownames(aux)),] -> trié par p-val

pval_sorted=aux[,'Pvalue']
line_sorted_by_pValue=as.numeric(rownames(aux)) 


# True Positive -------------------------------------------------------

TP=cumsum(aux$Species=='ups')

plot(TP,
     type='l',col='cadetblue',lwd=3,
     main='True Positif in Yeast Data',
     xlab='pval[1:n]',ylab='TP')

# Vstar ----------------------------------------------------------------
# V.star(S =,C = ,ZL = ,leaf_list = )
# zetas.tree(C = ,leaf_list = ,method = ,pvalues = ,alpha = )


C=final_tree[[2]]
leaf_list=final_tree[[1]]
alpha=0.05

# DKWM 

ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha)


iter_pointeur = 0     
Vstar_DKWM=rep(0,m)


for (i in 1:m){
  S=line_sorted_by_pValue[1:i]
  Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
  iter_pointeur=iter_pointeur+1
  print(iter_pointeur)
}

m=dim(res_ordered)[1]
plot(seq(m)-Vstar_DKWM,
     type='l',col='chocolate',lwd=3,
     main='V*',
     xlab='pval[1:n]',ylab='V.star')

# Simes -------------------------------------------------------------------------------
thr=alpha/m*(1:m)

Vsimes=sanssouci:::curveMaxFP(pval,thr)


plot(1:m-Vsimes,
     type='l',col='darkorchid',lwd=3,
     main='Vsimes',
     xlab='pval[1:n]',ylab='Vsimes')


# All ---------------------------------------------------------------------------------

plot(TP,
     type='l',col='cadetblue',lwd=3,
     main='Yeast Data no Na',
     xlab='pval[1:n]')
lines(1:m-Vstar_DKWM,
      type='l',col='chocolate',lwd=3)
lines(1:m-Vsimes,
      type='l',col='darkorchid',lwd=3,)
legend(x="topright", legend=c("Oracle","k-Vstar","k-Vsimes"), col=c("cadetblue",'chocolate',"chocolate4","darkorchid"), lty=c(1,1,1,1))

save(Vstar_DKWM,Vsimes,TP,file='Borne_na.RData')