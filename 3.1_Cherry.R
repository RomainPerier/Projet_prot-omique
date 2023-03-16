library(sanssouci)
library(cherry)
library(doBy)
library(ggplot2)
load('save/df_plot.RData')
load('save/res_ordered.RData')
load('save/final_tree.RData')
load('save/pval.RData')

## True Positive ---------------------------------------------------------------

TP=df_plot[df_plot$variable=='TP','value']

## Homel Fast ------------------------------------------------------------------

res_simes = hommelFast(res_ordered[,'Pvalue'])
res_hommel = hommelFast(res_ordered[,'Pvalue'],simes=FALSE)

m=nrow(res_ordered)
## Calcul Borne ----------------------------------------------------------------

#curveSimes(res_hommel)

curveSimes(res_simes,select=line_sorted_by_pval[1:100])

curveFisher(pval,order(pval))