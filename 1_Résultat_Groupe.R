library(sanssouci)
library(doBy)
library(stringr)
library(ggplot2)

## Importation des données et Trie ---------------------------------------------
load('save/proteom.RData')
res=read.csv("Fichier/proteom_test.txt",header=TRUE,sep="\t")[,c('id','Pvalue','Leading_razor_protein')]

hist(res$Pvalue)

Species=!str_detect(res[,"Leading_razor_protein"],"ups")
Species[Species==TRUE]<-"levure"
Species[Species==FALSE]<-"ups"
res=cbind(res,Species)
rm(Species)

# 192 peptides de ups / 18 506 peptides de levure 

res_ordered=orderBy(~ Species + Leading_razor_protein,res)
proteom <- proteom[as.numeric(row.names(res_ordered)),]

save(proteom,file='save/proteom.RData')
save(res_ordered,file='save/res_ordered.RData')
rm(res)
## Groupement par Protéines ----------------------------------------------------


split_data=splitBy(formula = ~ Species + Leading_razor_protein,res_ordered)[]
# 2418 protéines, 42 humaines - 2376 levure

save(split_data, file = "save/split_data.RData")

## Useful object --------------------------------------------------------------

pval = res_ordered[,'Pvalue']
line_sorted_by_pval = order(pval)
pval_sorted = sort(pval)
save(pval,line_sorted_by_pval,pval_sorted,file='save/pval.RData')
