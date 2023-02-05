function (C,leaf_list,method_node,method_leaf,pvalues,alpha){
  n_leaf = length(leaf_list) # = C[[length(C)]]
  H = length(C)
  K=0
  for (l in 1:H){ 
    K=K+length(C[[l]])
  }
  # Leaves : 
  zeta_leaves = rep(0,n_leaf)
  for (i in 1:n_leaf){
    leaf = leaf_list[[i]]
    zeta_leaves[i] <- method_leaf(pvalues[leaf],alpha/K)
  }
  #Nodes :
  zetas = list()
  for (l in H:1){ # On calcule pour chaque hauteur tout les zetas de la hauteur 
    Ch = C[[l]]
    zeta_nodes = rep(0,length(Ch))
    n = length(Ch)
    for (i in 1:n){# On calcule le zeta pour chaque noeud
      Chi = Ch[[i]]
      if (Chi[1]<Chi[2]){
        zeta_nodes[i] <- method_node(pvalues[unlist(leaf_list[Chi[1]:Chi[2]])],alpha/K)
      }
      else{
        zeta_nodes[i] <- zeta_leaves[Chi[1]]
      }
    }
    zetas[[l]] = zeta_nodes 
  }
  return(zetas)
}