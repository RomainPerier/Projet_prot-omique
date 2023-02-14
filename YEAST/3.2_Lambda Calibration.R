#Dans le papier BNR il n'y a pas que la méthode Simes, il y a aussi une méthode dite de "lambda-calibration" qui est implémentée dans sanssouci et est censée s'adapter à la dépendance des données, vous pouvez l'essayer ?

library(sanssouci)
library(ggplot2)

jointFWERControl

bidule <- sanssouci:::jointFWERControl(pval,refFamily='kFWER',0.05) #?????????????


# res -> tableau de p-value 