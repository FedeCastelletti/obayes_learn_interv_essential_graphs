gen.data = function(A.true, B, mu_0, Sigma, I = NULL, sigma_I = NULL, n){
  
  # Generate n q-variate observations from a Gaussian DAG model
  
  # INPUTS
  
  # A.true  : (q,q) adjacency matrix of the true DAG
  # B       : (q,q) matrix with regression coefficients in the observational data generating model Y = BY + eps
  # mu_0    : (q,1) mean vector of Y_1, ..., Y_q
  # Sigma   : (q,q) diagonal matrix with conditional variances of of Y_1, ..., Y_q
  # I       : intervention target, I in \{1, ..., q}
  # sigma_i : variance of the intervention density of y_I
  # n       : sample size
  
  # OUTPUTS
  
  # Y       : a (n,q) data of observations from Y | do(Y_I = y_i) 
  # D       : a (n,q) design matrix (associate each observation to the target node)
  
  
  A.int = A.true
  A.int[,I] = 0   # I remove all edges "pointing" to the intervened node (all j --> I)
  
  B_I = t(A.int*B) # remove the regression coefficients associated to  j -> I
  
  Sigma_I = Sigma
  Sigma_I[I,I] = sigma_I
  
  Id = diag(rep(1,q))
  
  S_I = solve(Id-B_I)%*%Sigma_I%*%t(solve(Id-B_I))
  
  Y = rmvnorm(n, mu_0, S_I)
  
  q = ncol(A.true)
  
  D = matrix(0, n, q)
  D[, I] = 1
  
  return(list(Y = Y, D = D))
  
}
