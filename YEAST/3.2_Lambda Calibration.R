library(sanssouci)
library(ggplot2)
library(stringr)
library(cgwtools)
load('save/data.RData')

## Création de l'objet SansSouci -----------------------------------------------

Y <- as.matrix(proteom[,3:20]) #epr_all -> gene ligne - sample colonne donc idem 

groups <- c(rep(0,9),rep(1,9))

truth <- as.numeric(str_detect(proteom[,'Leading_razor_protein'],"ups"))

SS_obj <- SansSouci(Y,groups,truth)

resave(SS_obj,file='save/data.RData')

## Objet SansSouci  -----------------------------------------------------------

alpha = 0.05

m = nrow(proteom)
n = 18 #ncol

cat(" nHyp : ",nHyp(SS_obj))
cat(" nObs : ",nObs(SS_obj))
print(SS_obj)


## Calibration -----------------------------------------------------------------

cal <- fit(SS_obj,alpha=0.05,B=1000,family='Simes')
cal_Oracle <- fit(SS_obj,0.05,family="Oracle")
volcanoPlot(cal)

cat("label : ",label(cal))
pval_ss = pValues(cal)
cat("pValues : ",summary(pval_ss))
thr_ss = thresholds(cal)
cat("thresholds : ",summary(thr_ss))

resave(cal,cal_Oracle,file='save/data.RData')

## Pvaleur by sanssouci vs DAPAR -----------------------------------------------

pval_dapar = pval 

pval_ss = cal$output$p.value

df_pval = rbind(data.frame(pValeur=pval_dapar,Méthode='DAPAR'),data.frame(pValeur=pval_ss,Méthode='Sanssouci'))

ggplot(df_pval,aes(x=pValeur,fill=Méthode))+
  geom_density(alpha=0.2)+
  ggtitle("Densité des P-valeurs selon la méthode")

## Lambda Calibration ----------------------------------------------------------
# cal0 / cal_Oracle / cal

cal0 = fit(SS_obj,alpha=0.05,B=0,vamily='Simes')

confs <- list(Simes = predict(cal0, all = TRUE),
              "Simes+calibration" = predict(cal, all = TRUE),
              "Oracle" = predict(cal_Oracle, all = TRUE))

plotConfCurve(confs, xmax = 200)

## Comparaison avec les bornes précédentes -------------------------------------
# df_bornes 
aux = predict(cal,all=TRUE)


df_bornes <- cbind(df_bornes, Vsimes_cal = (1:m-data.frame(aux[aux$stat=='TP',]$bound)))
resave(df_bornes,file="save/data.RData")
df_plot <-  rbind(df_plot, data.frame(Index=1:m,variable='Vsimes_Cal',value=aux[aux$stat=='TP',]$bound))

resave(df_plot,file='save/data.RData')

ggplot(df_plot,aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ylim(c(0,200))+
  ggtitle('Lower Bound on True Positive in Yeast Data')