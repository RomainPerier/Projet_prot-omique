#Dans le papier BNR il n'y a pas que la méthode Simes, il y a aussi une méthode dite de "lambda-calibration" qui est implémentée dans sanssouci et est censée s'adapter à la dépendance des données, vous pouvez l'essayer ?

library(sanssouci)
library(ggplot2)

jointFWERControl

bidule <- sanssouci:::jointFWERControl(pval,refFamily='kFWER',0.05) #?????????????
# jointFWERControl remplacer par CalibrateJER0 voir github 
# voir les objets sanssouci avec fit etc pour faire la lambda calibration # avec class comme en python mais version R, les méthodes c'est avec $

# res -> tableau de p-value 