no.pretraitement <- function(bdd){
  tab_data_av_low=bdd[,3:11]
  tab_data_av_high=bdd[,12:20]
  
  data_av_low = unlist(tab_data_av_low,use.names=FALSE)
  data_av_high = unlist(tab_data_av_high,use.names=FALSE)
  
  data_av=data.frame(Méthode='no_imputation',
                     rbind(data.frame(Condition="Low Spike",
                                    Value=data_av_low[data_av_low!=0])
                           ,data.frame(Condition="High Spike",
                                     Value=data_av_high[data_av_high!=0]))) 

  return(data_av)
}