library(sanssouci)
library(doBy)
library(stringr)

## Importation des données et Trie ---------------------------------------------
bdd=read.csv('proteom_treated.txt',header=TRUE,sep="\t")
res=read.csv("proteom_treated_test.txt",header=TRUE,sep="\t")[,c('id','Pvalue','Leading_razor_protein')]
#all(bdd[,c(1,2)]==res[,c(1,3)])
#[1] TRUE   -> pas d'erreur 



Species=!str_detect(res[,"Leading_razor_protein"],"ups") 
Species[Species==TRUE]<-"levure"
Species[Species==FALSE]<-"ups"
res=cbind(res,Species)
rm(Species)

# 192 peptides de ups / 18 506 peptides de levure 

res_ordered=orderBy(~ Species + Leading_razor_protein,res)
# res_ordered['i',]==res[i,]
save(res_ordered,file='res_ordered.RData')

## Groupement par Protéines ------------------------------------------------------------


split_data=splitBy(formula = ~ Species + Leading_razor_protein,res_ordered)[]
# 2418 protéines, 42 humaines - 2376 humaine 

save(split_data, file = "split_data.RData")


