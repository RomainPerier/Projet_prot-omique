library(sanssouci)
library(doBy)
library(stringr)

## Importation des données et Trie ---------------------------------------------

res=read.csv("ARATH/proteom_test_arath.txt",header=TRUE,sep="\t")[,c('id','Pvalue','Leading_razor_protein')]

Species=str_detect(res[,"Leading_razor_protein"],"ARATH") 
Species[Species==TRUE]<-"ARATH"
Species[Species==FALSE]<-"UPS1"
res=cbind(res,Species)
rm(Species)

res_ordered=orderBy(~ Species + Leading_razor_protein,res)

save(res_ordered,file='ARATH/save/res_ordered.RData')

## Groupement par Protéines ------------------------------------------------------------


split_data=splitBy(formula = ~ Species + Leading_razor_protein,res_ordered)[]
save(split_data, file = "ARATH/save/split_data.RData")

## Useful object --------------------------------------------------------------

pval = res_ordered[,'Pvalue']
line_sorted_by_pval = order(pval)
pval_sorted = sort(pval)
save(pval,line_sorted_by_pval,pval_sorted,file='ARATH/save/pval.RData')

