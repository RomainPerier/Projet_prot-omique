library(sanssouci)
library(doBy)
library(stringr)

## Importation des données et Trie ---------------------------------------------

res=read.csv("proteom_treated_test.txt",header=TRUE,sep="\t")[,c('id','Pvalue','Leading_razor_protein')]

Species=str_detect(res[,"Leading_razor_protein"],"ARATH") 
Species[Species==TRUE]<-"ARATH"
Species[Species==FALSE]<-"UPS1"
res=cbind(res,Species)
rm(Species)

res_ordered=orderBy(~ Species + Leading_razor_protein,res)

save(res_ordered,file='res_ordered.RData')

## Groupement par Protéines ------------------------------------------------------------


split_data=splitBy(formula = ~ Species + Leading_razor_protein,res_ordered)[]
save(split_data, file = "split_data.RData")


