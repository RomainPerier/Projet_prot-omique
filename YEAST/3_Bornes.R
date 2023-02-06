library(sanssouci)
library(doBy)

load('data.RData')

# On a fait correspondre les lignes de la base de données initiales à res_ordered
#> all(proteom[,'Sequence']==res_ordered[,'id'])
#[1] TRUE


pval = res_ordered[,'Pvalue']
proteom=cbind(proteom,pval)

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


#iter_pointeur = 0     
Vstar_DKWM=rep(0,m)


for (i in 1:m){
  S=line_sorted_by_pValue[1:i]
  Vstar_DKWM[i]<-V.star(S,C,ZL_DKWM,leaf_list)
  #iter_pointeur=iter_pointeur+1
  #print(iter_pointeur)
}

m=dim(res_ordered)[1]
plot(seq(m)-Vstar_DKWM,
     type='l',col='chocolate',lwd=3,
     main='V*',
     xlab='pval[1:n]',ylab='V.star')

# BH 


ZL_HB = zetas.tree(C,leaf_list,zeta.HB,pval,alpha)


#iter_pointeur = 0    
Vstar_HB=rep(0,m)


for (i in 1:m){
  S=line_sorted_by_pValue[1:i]
  Vstar_HB[i]<-V.star(S,C,ZL_HB,leaf_list)
  # iter_pointeur=iter_pointeur+1
  # print(iter_pointeur)
}

# Multiple zeta 

ZL_multiple = zeta.multiple(C,leaf_list,zeta.DKWM,zeta.HB,pval,0.05)


iter_pointeur = 0    
Vstar_mlt=rep(0,m)


for (i in 1:m){
  S=line_sorted_by_pValue[1:i]
  Vstar_mlt[i]<-V.star(S,C,ZL_multiple,leaf_list)
  iter_pointeur=iter_pointeur+1
  print(iter_pointeur)
}

plot(1:m-Vstar_mlt,
     type='l',col='seagreen',lwd=3,
     main='Vstar with multiple zeta',
     xlab='pval[1:n]',ylab='Vstar')


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
     main='Yeast Data',
     xlab='pval[1:n]')
lines(1:m-Vstar_DKWM,
     type='l',col='chocolate',lwd=3)
lines(1:m-Vstar_HB,
      type='l',col='chocolate4',lwd=3)
lines(1:m-Vsimes,
     type='l',col='darkorchid',lwd=3,)
legend(x="topright", legend=c("Oracle","k-Vstar with zeta.DKWM","k-Vstar with zeta.HB","k-Vsimes"), col=c("cadetblue",'chocolate',"chocolate4","darkorchid"), lty=c(1,1,1,1))


save(TP,Vstar_DKWM,Vstar_HB,Vstar_mlt,Vsimes,file='Borne.RData')
