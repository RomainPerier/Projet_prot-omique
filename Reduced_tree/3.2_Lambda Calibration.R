library(sanssouci)
library(ggplot2)
library(stringr)
library(cgwtools)
load('Reduced_tree/save/df_plot.RData')
load('Reduced_tree/save/proteom.RData')
load('Reduced_tree/save/final_tree.RData')
load('Reduced_tree/save/pval.RData')

alpha = 0.01

m = nrow(proteom)
n = 18 #ncol

## Sans-souci
Y <- as.matrix(proteom[,3:20]) #epr_all -> gene ligne - sample colonne donc idem 

groups <- c(rep(0,9),rep(1,9))

truth <- as.numeric(str_detect(proteom[,'Leading_razor_protein'],"ups"))

SS_obj <- SansSouci(Y,groups,truth)

##save(SS_obj,file='ARATH/save/sanssouci_obj.RData')


## Calibration -----------------------------------------------------------------

cal <- fit(SS_obj,alpha=alpha,B=1000,family='Simes')
cal_Oracle <- fit(SS_obj,alpha,family="Oracle")
cal0 <- fit(SS_obj,alpha=alpha,B=0,family='Simes') #sans calibration

#pval_ss = pValues(cal)
#thr_ss = thresholds(cal)

save(cal,cal_Oracle,cal0,file='Reduced_tree/save/sanssouci_obj.RData')

confs <- list(Simes = predict(cal0, all = TRUE),
              "Simes+calibration" = predict(cal, all = TRUE),
              "Oracle" = predict(cal_Oracle, all = TRUE))

plotConfCurve(confs)

## Comparaison avec les df_plot précédentes -------------------------------------

aux = predict(cal,all=TRUE)



df_plot <-  rbind(df_plot, data.frame(Index=1:m,variable='Vsimes_Cal',value=aux[aux$stat=='TP',]$bound))

resave(df_plot,file='Reduced_tree/save/df_plot.RData')

ggplot(df_plot[df_plot$variable!='TP_Vstar_refined',],aes(x=Index,y=value,color=variable))+
  geom_line(lwd=1) +  
  ylim(c(0,300))+
  ggtitle('Lower Bound on True Positive in Yeast Data')

