# Example

# Case 1: only observational data

library(pcalg)
library(mvtnorm)

q = 10

set.seed(123)
D = randomDAG(q, prob = 0.2)

A = as(D, "matrix")
A[A!=0] = 1

true.dag = as(A,"graphNEL")
true.dag.mat = as(true.dag,"matrix")

B = A*matrix(runif(q*q, 0.5, 1), q, q)

Sigma = diag(1, q, q)
mu_0  = rep(0, q)

# observational data

n.obs = n.int = 1000

source("gen_data.r")

data = gen.data(A.true = A, B = B, mu_0 = mu_0, Sigma = Sigma, n = n.obs)

Y.obs = data$Y
D.obs = data$D

# Apply OBES and find PPIs matrix

m = 2*q

T = 1000
burn = 200

source("posterior_i_ess.r")

t_0 = proc.time()
out.obs = post_i_ess_graphs(Y = Y.obs, I.cal = integer(), D = D.obs, m = m, T = T, burn = burn, verbose = TRUE)
t_1 = proc.time() - t_0

t_1

PPI     = out.obs$PPI
PPI_und = out.obs$PPI_und

med = round(PPI)
med = as(med, "graphNEL")
med = pdag2dag(med)
med = dag2essgraph(med$graph)

true.ess = dag2cpdag(true.dag)

shd(med, true.dag)
shd(med, true.ess)

plot(med)
plot(true.ess)

source("objective_function.r")

out.target = opt.intervention(PPI_und, lambda = 1)
I.opt = out.target$target


# Case 2: combine both observational and interventional data

data = gen.data(A.true = A, B = B, mu_0 = mu_0, Sigma = Sigma, I = I.opt, sigma_I = 0.1, n = n.int)

Y.int = data$Y
D.int = data$D

Y = rbind(Y.obs, Y.int)
D = rbind(D.obs, D.int)

I.cal = list(integer(), I.opt)

t_0 = proc.time()
out = post_i_ess_graphs(Y = Y, I.cal = I.cal, D = D, m = m, T = T, burn = burn, verbose = TRUE)
t_2 = proc.time() - t_0

PPI     = out$PPI
PPI_und = out$PPI_und

med = round(PPI)
med = as(med, "graphNEL")
med = pdag2dag(med)
med = dag2essgraph(med$graph, targets = I.cal)

true.ess = dag2essgraph(true.dag, targets = I.cal)

shd(med, true.dag)
shd(med, true.ess)

plot(med)
plot(true.ess)

out.target = opt.intervention(PPI_und, lambda = 1)
I.opt = out.target$target

