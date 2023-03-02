library(sanssouci)
library(cgwtools)
library(ggplot2)
library(reshape2)
library(doBy)
library(stringr)

bdd=read.csv("Fichier/peptides_YEASTUPS.txt",header=TRUE,sep="",blank.lines.skip=TRUE)

## Filtration ------------------------------------------------------------------

bdd[bdd==0] <- NaN
proteom_no_na = na.omit(bdd)

## Test ------------------------------------------------------------------------
obj = SansSouci(as.matrix(proteom_no_na[,-c(1,2)]),groups=c(rep(0,9),rep(1,9)),truth = as.numeric(str_detect(proteom_no_na[,'Leading_razor_protein'],'ups')))

cal = fit(obj,alpha=0.05,B=0,family='Simes')

pval = pValues(cal)

res =  data.frame(proteom_no_na[,c(1,2)],Pvalue = pval)


## Structure -------------------------------------------------------------------

Species = str_detect(proteom_no_na[,'Leading_razor_protein'],'ups')
Species[Species==TRUE] <- 'ups'
Species[Species==FALSE] <- 'levure'
res <- cbind(res,Species)         ## 77 peptides humaines seulement 

res_ordered = orderBy(~ Species + Leading_razor_protein,res)
proteom_no_na <- proteom_no_na[row.names(res_ordered),]

split_data=splitBy(formula = ~ Species + Leading_razor_protein,res_ordered)[]
#31 protéine humaines 

m=nrow(res_ordered)
n_leaf=length(split_data)

taille_leaf=rep(0,n_leaf) 
for (i in 1:n_leaf){
  taille_leaf[i]<-nrow(split_data[[i]])
}
taille_leaf_cum=cumsum(taille_leaf) 


leaf_list=list(1:(taille_leaf_cum[1]))
for (i in 2:n_leaf){
  j=i-1
  leaf_list=append(leaf_list,list((taille_leaf_cum[j]+1):(taille_leaf_cum[i])))
}
rm(taille_leaf,taille_leaf_cum,i,j)

aux=!str_detect(names(split_data),"ups")
n_levure=length(aux[aux])    
n_ups=length(aux[!aux])    
n_species_cum=c(n_levure,n_levure+n_ups)
rm(aux) 

C3=list(c(1,1))   
for (i in 2:n_leaf){
  C3<-append(C3,list(c(i,i)))
}

# C1 <- list(c(1,n_leaf))
# C2 <-list(c(1,n_levure) , c(n_levure+1,n_leaf))
C=list(list(c(1,n_leaf)),list(c(1,n_levure),c(n_levure+1,n_leaf)),C3)    
rm(C3)

final_tree=list(leaf_list,C)
rm(leaf_list,C,m,n_levure,n_species_cum,n_ups,n_leaf,i)


## Bornes ----------------------------------------------------------------------

pval = res_ordered[,'Pvalue']
C=final_tree[[2]]
leaf_list=final_tree[[1]]
m=nrow(proteom_no_na)
line_sorted_by_pval=order(pval)
# TP :

TP=cumsum(res_ordered[line_sorted_by_pval,]$Species=='ups')
df_bornes = data.frame(Index=1:m)

df_bornes <- cbind(df_bornes,TP)

alpha=0.05
K = 200 # Calculer tout les K lignes
ZL_DKWM=zetas.tree(C,leaf_list,zeta.DKWM,pval,alpha)
ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval,alpha)
#ZL_refined=zetas.tree.refined(C,leaf_list,zeta.DKWM,pval,alpha)

Vstar_DKWM=rep(0,m%/%K)
Vstar_HB=rep(0,m%/%K)
#Vstar_refined=rep(0,m%/%K)

for (i in 1:(m%/%K)){
  S=line_sorted_by_pval[1:(K*i)]
  Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
  Vstar_HB[i]<-V.star(S,C,ZL_HB,leaf_list)
  #Vstar_refined[i] <- V.star(S,C,ZL_refined,leaf_list)
  print(i)
}
# Vsimes : 
thr=alpha/m*(1:m)
Vsimes=sanssouci:::curveMaxFP(pval,thr)[1:(m%/%K)*K]

df <- cbind(df_bornes[1:(m%/%K)*K,],Vstar_DKWM,Vstar_HB,Vsimes)

resave(df,file='save/permut.RData')

#Plot : 
df_plot =  cbind(df[,1:2]
                 #,TP_Vsimes=df[,'Index']-df[,'Vsimes']
                 ,TP_Vstar_DKWM=df[,'Index']-df[,'Vstar_DKWM'])
#                 ,TP_Vstar_HB=df[,'Index']-df[,'Vstar_HB'])
#,TP_Vstar_refined=df[,'Index']-df[,'Vstar_refined'])

df_plot = melt(df_plot, id.vars = "Index")


ggplot(df_plot,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ggtitle('Lower Bound on True Positive in Yeast Data')


