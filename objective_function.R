## Notation

# i is the intervened node
# j is any other node (j = 1, ...,q and j different from i)

# p_ij is the marginal posterior probability of edge inclusion (marginal PPI) for the undirected edge (i,j)

# Heterogeneity associated to edge (i,j)

H_j = function(p_ij){
  
  # i is fixed beforehand
  
  # i,j : nodes
  # p_ij, p_ji : marginal PPIs
  
  # 2 - 2*((p_ij/(p_ij + p_ji))^2 + (p_ji/(p_ij + p_ji))^2)
  # 1 - ((p_ij/(p_ij + p_ji))^2 + (p_ji/(p_ij + p_ji))^2)
  
  p_ij * p_ji
  
}

# Example

# H_j(1,0)       # homogeneity
# H_j(1,1)       # maximum heterogeneity
# H_j(0.5,0.5)   # maximum heterogeneity


# Total heterogeneity associated to intervention on node i

H_tot_j = function(i,P,lambda){
  
  # P : (q,q) matrix with marginal PPIs
  # i : the intervened node
  
  q = dim(P)[1]
  
  if(lambda == 0){
    
    h_i = (prod(P[i,-i]))^(1/(q-1))
    
  } else{
  
    h_i = (sum(P[i,-i]^lambda)/(q-1))^(1/lambda)
    
    }
  
  return(h_i)
  
}

# To find the optimal intervention target

opt.intervention = function(P,lambda){
  
  options("scipen" = 10)
  
  P[P == 0] = 0.00001
  
  H = sapply(1:dim(P)[1], H_tot_j, P, lambda = lambda)
  
  return(list(H = H, target = which.max(H)))
  
}
