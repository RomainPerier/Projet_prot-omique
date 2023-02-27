zetas.tree.refined <- function(C, leaf_list, method, pvalues, alpha) {
  H <- length(C)
  K <- sanssouci:::nb.elements(C)
  leaves <- length(leaf_list)
  zeta_leaves <- numeric(leaves)
  continue <- TRUE
  new_K <- K
  while (continue) {
    usage_K <- new_K
    new_K <- K
    CH <- C[[H]]
    for (i in 1:length(CH)) {
      CHi <- CH[[i]]
      if (CHi[1] == CHi[2]) {
        pvals <- pvalues[leaf_list[[CHi[1]]]]
        zeta_leaves[CHi[1]] <- method(pvals, alpha/usage_K)
        if (zeta_leaves[CHi[1]] == 0) {
          new_K <- new_K - 1
        }
      }
    }
    ZL <- list()
    for (h in H:1) {
      Ch <- C[[h]]
      len <- length(Ch)
      zeta_inter <- numeric(len)
      for (j in 1:len) {
        Chj <- Ch[[j]]
        if (Chj[1] < Chj[2]) {
          pvals <- pvalues[unlist(leaf_list[Chj[1]:Chj[2]])]
          zeta_inter[j] <- method(pvals, alpha/usage_K)
          if (zeta_inter[j] == 0) 
            new_K <- new_K - 1
        } else {
          zeta_inter[j] <- zeta_leaves[Chj[1]]
        }
      }
      ZL[[h]] <- zeta_inter
    }
    if (new_K == usage_K) {
      continue <- FALSE
    }
  }
  return(ZL)
}
