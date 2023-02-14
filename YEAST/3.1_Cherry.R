library(sanssouci)
library(cherry)
library(doBy)

load('data.RData')
#Et aussi essayer une des méthodes de Goeman et Solari (nos "concurrents"), dans leur package à eux qui s'appelle cherry. La méthode "Hommel fast".

## True Positive ---------------------------------------------------------------

TP=cumsum(orderBy(~ Pvalue,res_ordered)$Species=='ups')

plot(TP,
     type='l',col='cadetblue',lwd=3,
     main='True Positif in Yeast Data',
     xlab='pval[1:n]',ylab='TP')

## Homel Fast ------------------------------------------------------------------

res_simes = hommelFast(res_ordered[,'Pvalue'])
res_hommel = hommelFast(res_ordered[,'Pvalue'],simes=FALSE)

m=dim(res_ordered)[1]
## Calcul Borne ----------------------------------------------------------------
# Curve Simes 


curveSimes(res_hommel)

curveSimes(res_simes)


# Curve Fisher


curveFisher(res_ordered[,'Pvalue'])
# sauvegarder dans Yeast car long à générer 