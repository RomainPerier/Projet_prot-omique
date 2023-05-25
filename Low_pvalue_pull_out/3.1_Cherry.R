library(sanssouci)
library(cherry)
library(doBy)
library(ggplot2)
load('Low_pvalue_pull_out/save/df_plot.RData')
load('Low_pvalue_pull_out/save/res_ordered.RData')
load('Low_pvalue_pull_out/save/final_tree.RData')
load('Low_pvalue_pull_out/save/pval.RData')

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