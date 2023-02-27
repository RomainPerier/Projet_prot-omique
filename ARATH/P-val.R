library(sanssouci)
library(doBy)
library(stringr)

res=read.csv("proteom_treated_test.txt",header=TRUE,sep="\t")[,c('id','Pvalue','Leading_razor_protein')]
res_na=read.csv("proteom_treated_na_test.txt",header=TRUE,sep="\t")[,c('id','Pvalue','Leading_razor_protein')]
## Ajout de la colonne espèce----------
# (j'ai tester et l'ordre des sequence est bien identique ))

Species=str_detect(res[,"Leading_razor_protein"],"ARATH") 
Species[Species==TRUE]<-"ARATH"
Species[Species==FALSE]<-"UPS1"
res=cbind(res,Species)

Species=str_detect(res_na[,"Leading_razor_protein"],"ARATH") 
Species[Species==TRUE]<-"ARATH"
Species[Species==FALSE]<-"UPS1"
res_na=cbind(res_na,Species)
rm(Species)



# P-val ne sont pas exactement égale (all(res_limma[,'Pvalue']==res_premier[,'Pvalue'])  False !)

## tracer des p-val --------------

pval=res[,'Pvalue']
pval_na=res_na[,'Pvalue']


plot(ecdf(pval),lwd=2.0,main="Fonction de Répartition Empirique des p-valeurs (Student) ")
lines(c(0,1),c(0,1),col='red',lwd=2.0)


plot(ecdf(pval_na),lwd=2.0,main="Fonction de Répartition Empirique des p-valeurs (Student, NA) ")
lines(c(0,1),c(0,1),col='red',lwd=2.0)


## BH student #on fait que sur un comme les résultat sont sensiblement pareil 

#On notres vecteur de p-val pval_student 

pval_sorted=sort(pval)
pval_na_sorted=sort(pval_na)

#on détermine le khat 

alpha=0.05
n=length(pval)   #14 427
n_na=length(pval_na)   #8484

aux=pval_sorted<=alpha*(1:n)/n
aux_na=pval_na_sorted<=alpha*(1:n_na)/n_na
khat=max(which(aux))
khat_na=max(which(aux_na))

#Region de rejet : 

R=which(pval<=alpha*khat/n)
print(length(R))  #7050


R_na=which(pval_na<=alpha*khat_na/n_na)
print(length(R_na))  #4007


maxFP(pval,alpha)   #6306
maxFP(pval_na,alpha) #3831


load('final_tree.RData')
C=final_tree[[2]]
leaf_list=final_tree[[1]]

Z=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha)
V.star(1:n,C,Z,leaf_list)    #4407


load('final_tree_na.RData')
C_na=final_tree[[2]]
leaf_list_na=final_tree[[1]]

Z_na=zetas.tree(C_na,leaf_list_na,zeta.DKWM,pval_na,alpha)
V.star(1:n_na,C_na,Z_na,leaf_list_na)    #2703

