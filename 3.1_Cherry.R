library(sanssouci)
library(cherry)
library(doBy)
library(ggplot2)
load('save/data.RData')

## True Positive ---------------------------------------------------------------

TP=df_bornes[,'TP']

## Homel Fast ------------------------------------------------------------------

res_simes = hommelFast(res_ordered[,'Pvalue'])
res_hommel = hommelFast(res_ordered[,'Pvalue'],simes=FALSE)

m=nrow(res_ordered)
## Calcul Borne ----------------------------------------------------------------

curveSimes(res_hommel)

curveSimes(res_simes)

curveFisher(pval,order(pval))