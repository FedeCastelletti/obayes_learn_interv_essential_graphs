# Markov chain

# A_0     adjacency matrix of the initial CPDAG E_0 (for instance an empty matrix)
# alpha   accelerator parameter (fraction of operators to check)
# m       maximum number of edges (undirected or directed) implying a sparsity condition
# T       number of iterations

library(pcalg)
library(igraph)
library(graph)
library(gRbase)
library(mvtnorm)
library(ggm)

source("operators_i.r")
source("move_i_simple.r")
source("marg_like.r")

A_undirected = function(A){
  A_und = matrix(0, nrow(A), ncol(A))
  A_und[which(A == 1 & t(A) == 1,TRUE)] = 1
  A_und
}

n.edge = function(A){
  length(which(A[lower.tri(A)] == 1 | t(A)[lower.tri(A)] == 1))
}

post_i_ess_graphs = function(Y, I.cal, D, m, T, burn, verbose = FALSE){
  
  # Returns an MCMC chain over the I-EG space
  
  # INPUTS
  
  # Y       : the (n,q) data matrix
  # I.cal   : a (conservative) family of intervention targets
  # D       : (n,q) design matrix (associate each multivariate observationa to a target)
  # m       : the maximum number of edges (sparsity constraints on the I-EGs space)
  # T       : the number of MCMC iterations
  # burn    : burn-in period
  # verbose : TRUE/FALSE to see or not plot of posterior density of visited EGs
  # verbose:TRUE/FALSE to see or not plot of posterior density of visited EGs
  
  # OUTPUTS
  
  # A_chain : (q*q,T) matrix whose rows are the vectorized adjacency matrices of the graphs visited by the chain
  # PPI     : (q,q) matrix with marginal posterior probabilities of edge inclusion (frequency of visits associated to each edge)
  # PPI_und : (q,q) matrix with marginal posterior probabilities of inclusion of each undirected edge (frequency of visits associated to each undirected edge)
  
  
  q = ncol(Y)
  
  A_0 = matrix(0,q,q)
  rownames(A_0) = colnames(A_0) = 1:q
  
  A_chain = matrix(0,nrow=q*q,ncol=T+1)
  A_chain[,1] = A_0
  # matrix collecting the output of the MCMC (each column is the by-column-vectorized adjacency matrix of the accepted graph)
  
  # store space for logposterior
  logpost = rep(-Inf,T+1)
  
  A_old = A_0
  m_old = marg_like(A_old, Y, D, I.cal)
  G_und = G_undirected(A_old)
  
  colnames(A_0) = rownames(A_0) = 1:q
  
  PPI     = matrix(0, q, q)
  PPI_und = matrix(0, q, q)

  # MCMC iterations
  for(t in 2:T){
    
    # proposed move
    prop = move_i(A_old, m, G_und, q, I.cal = I.cal)
    A_new = prop[[1]]
    m_new = marg_like(A_new, Y, D, I.cal)
    
    # log prior ratio evaluation
    logprior.new = lgamma(n.edge(A_new)+1) + 
      lgamma(q*(q-1)/2-n.edge(A_new)+(2*q-2)/3-1)
    logprior.old = lgamma(n.edge(A_old)+1) + 
      lgamma(q*(q-1)/2-n.edge(A_old)+(2*q-2)/3-1)
    logprior = logprior.new - logprior.old
       
    # acceptance ratio
    ratio = min(0, m_new - m_old + logprior)
    
    if(log(runif(1,0,1)) < ratio){ # accept move
      
      A_old = A_new
      m_old = m_new
      G_und = G_undirected(A_old)
      logprior.old = logprior.new
      
      }
    
    if(t > burn){
      
      PPI     = PPI + A_old
      PPI_und = PPI_und + A_undirected(A_old)
      
    }
    
    # store chain value and posterior density
    A_chain[,t] = A_old
    logpost[t]  = m_old + logprior.old
    # print(ratio)
    
    # show plot of posterior if verbose
    if(verbose){
      if(t%%10 == 0 & t > burn) plot(logpost[burn:t], type = "l")
    }
 
  }
  
  return(list(A_chain = A_chain, PPI = PPI/(T-burn), PPI_und = PPI_und/(T-burn)))
  
}
