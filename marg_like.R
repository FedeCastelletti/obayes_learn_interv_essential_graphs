chain_comp = function(ess){
  
  # find chain components of the EG ess
  
  ess = as(ess,"igraph")
  
  amat = as.matrix(get.adjacency(ess))  # if the argument is ess.obs (the essential graph, a graph object)
  wmat = matrix(as.integer((amat + t(amat)) > 1), nrow = nrow(amat))
  wg   = graph.adjacency(wmat, mode = "undirected")
  
  cc   = clusters(wg)
  neworder = order(cc$membership)
  a = matrix(0, nrow = length(cc$csize), ncol = length(cc$csize))
  b = cumsum(cc$csize)
  wmat = amat[neworder, neworder]
  
  for(i in 1: length(cc$csize)){
    for(j in 1: length(cc$csize)){
      if(j != i){
        a[i,j] = as.integer(sum(wmat[(max(b[i-1],0)+1):b[i],
                                      (max(b[j-1],0)+1):b[j]]) > 0)
      }
    }
  }
  
  rownames(a) = colnames(a) = as.character(1:length(b))
  
  chainorder = topOrder(a)
  vertorder  = c()
  chainsize  = c()
  
  for(k in 1:length(b)){
    vertorder = c(vertorder, which(cc$membership == chainorder[k]))
    chainsize = c(chainsize, cc$csize[chainorder[k]])
  }
  
  q = list(vert.order=vertorder,chain.size=chainsize)
  order = q$vert.order
  size  = c(1,q$chain.size)
  
  Tau = list(0)
  for(i in 2:(length(size))){
    Tau[[i-1]] = order[sum(size[1:(i-1)]):(sum(size[2:i]))]
  }
  
  Tau
  
}

Gamma = function(s, a){
  
  # Returns the logarithm of the multivariate Gamma function G(s,a) 
  
  (.25*s*(s-1))*log(pi) + sum(lgamma(a+.5*(1-1:s)))
  
}

Ehat = function(y, x){
  
  # Return the residuals of the least square estimate of a linear model
  
  beta = solve(t(x)%*%x)%*%t(x)%*%y
  y - x%*%beta
  
}

Int_G = function(G, i){
  
  # return the intervention graph of G for the target i
  
  A.int = G
  A.int[,i] = 0
  A.int
  
}

m_J = function(J, Y, tau_c, pa_tau, pa_tau_c, a_D, n_0){
  
  # compute the marginal likelihood m(Y_J) (selected columns of Y indexed by J)
  # J will be, from time to time, a clique or a separator
  
  if(identical(unlist(J), character(0))){
    m = 0
  }
  
  else{
    
    J   = as.numeric(J)
    J_c = length(J)
    Y_J = as.matrix(Y[,J])
    n   = nrow(Y_J)
    q   = ncol(Y_J)
    
    X_J = as.matrix(cbind(1,Y[,pa_tau]))
    
    m = (-.5*(n-n_0)*J_c)*log(pi) +
        Gamma(J_c,0.5*(a_D+n-pa_tau_c-1-tau_c+J_c)) - Gamma(J_c,0.5*(a_D+n_0-pa_tau_c-1-tau_c+J_c)) +
        .5*(J_c*(a_D+n_0-tau_c+J_c))*log(n_0/n) +
        (-.5*(n-n_0))*log(det(t(Ehat(Y_J,X_J))%*%Ehat(Y_J,X_J)))
    
    
  }
  
  m      
  
}

m_tau = function(tau, Y, G){
  
  # compute the marginal likelihood for a chain component tau given the data Y
  
  taug  = G[tau,tau]
  tau_c = length(tau)
  
  G = as (G,"graphNEL")     # simplify
  pa_tau = as.numeric(parents(tau,G))
  G = as (G,"matrix")
  
  pa_tau_c = length(pa_tau)
  
  n_0 = ncol(Y)+1
  a_D = tau_c-1  
  
  if(tau_c <= 2){
    
    m_J(tau,Y,tau_c,pa_tau,pa_tau_c,a_D,n_0)
    
  } else{
    
    C = mpd(taug)$cliques
    S = mpd(taug)$separators
    
    sum(unlist(lapply(C,m_J,Y,tau_c,pa_tau,pa_tau_c,a_D,n_0)))-
      sum(unlist(lapply(S,m_J,Y,tau_c,pa_tau,pa_tau_c,a_D,n_0)))
    
  }
  
}

marg_like = function(G, Y, D, I.cal){

  # compute the marginal likelihood of the I-EG G
  
  # G     : the adjacency matrix of the I-EG G
  # Y     : the complete (n,q) dataset (both observational and interventional data)
  # D     : (n,q) design matrix (associate each multivariate observationa to a target)
  # I.cal : family of intervention targets
  
  if(is.null(I.cal)){
    
    # if only observational data are available
    
    n = nrow(Y)
    A = G
    ess = as(G,"graphNEL")
    
    # need to specify first ess.obs (the essential graph)
    
    Tau = chain_comp(ess)
    
    return(sum(unlist(lapply(Tau, m_tau, Y = Y, G = G))))
    
  } else{
    
    # if (also) interventional data are available
    
    VI   = (1:q)[colSums(D) != 0]
    
    Tau  = chain_comp(G)
    TauI = setdiff(unlist(I.cal),integer(0))
    TauI = as.list(TauI)
    
    Tau0 = setdiff(Tau,VI)
    
    
    marg_I = function(i,Y,G){
      i = unlist(i)
      G = Int_G(G,i)
      m_tau(i,Y[which(D[,i]==1),],G)
    }
    
    marg_O = function(i,Y,G){
      i = unlist(i)
      m_tau(i,Y[which(D[,i]!=1),],G)
    }
    
    marg_Obs = function(i,Y,G){
      i = unlist(i)
      m_tau(i,Y,G)
    }
    
    sum(unlist(lapply(TauI,marg_I,Y,G)))+
      sum(unlist(lapply(TauI,marg_O,Y,G)))+
      sum(unlist(lapply(Tau0,marg_Obs,Y,G)))
    
  }

}
