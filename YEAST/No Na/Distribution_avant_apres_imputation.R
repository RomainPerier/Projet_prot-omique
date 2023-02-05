library(ggplot2)

bdd_avant=read.csv("peptides_YEASTUPS.txt",header=TRUE,sep="",blank.lines.skip=TRUE)
bdd_apres=read.csv("proteom_treated.txt",header=TRUE,sep="",blank.lines.skip=TRUE)


## Low-Spike ----------------------------------------------------------------------------------
#on créée le vecteur contenant toutes les valeurs mesurées

tab_data_av=bdd_avant[,(3:11)]

data_av=c(tab_data_av[,1],
          tab_data_av[,2],
          tab_data_av[,3],
          tab_data_av[,4],
          tab_data_av[,5],
          tab_data_av[,6],
          tab_data_av[,7],
          tab_data_av[,8],
          tab_data_av[,9])

data_av=data_av[data_av!=0]

hist(data_av,breaks=10000,xlim=c(0,100000000))

tab_data_ap=bdd_apres[,(3:11)]

data_ap=c(tab_data_ap[,1],
          tab_data_ap[,2],
          tab_data_ap[,3],
          tab_data_ap[,4],
          tab_data_ap[,5],
          tab_data_ap[,6],
          tab_data_ap[,7],
          tab_data_ap[,8],
          tab_data_ap[,9])

data_ap=data_ap[data_ap!=0]

hist(data_ap,breaks=10000,xlim=c(0,100000000))
