library(sanssouci)
library(doBy)
library(stats)
library(stringr)
library(ggplot2)

## Importation des données et Trie ---------------------------------------------
load('save/proteom.RData')
res=read.csv("Fichier/proteom_treated_test.txt",header=TRUE,sep="\t")[,c('id','Pvalue','Leading_razor_protein')]

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

## Test sur les p-valeurs ------------------------------------------------------

pval_prostar = res_ordered$Pvalue

pval_sanssouci = as.numeric(pValues(fit(SansSouci(as.matrix(proteom[,4:21]),groups=c(rep(0,9),rep(1,9)),truth=as.numeric(proteom$Species=='ups')),B=0,alpha=0.05,family='Simes')))

df <- rbind(data.frame(Méthode="Prostar",Pvalue=pval_prostar),data.frame(Méthode="SansSouci",Pvalue=pval_sanssouci))

ggplot(data=df,aes(x=Pvalue,fill=Méthode))+geom_density(alpha=0.2)

plot(ecdf(pval_prostar))
lines(ecdf(pval_sanssouci),col='blue')
