

#' Create one SBM (one community several group)

SBMatrix <- function(q = 4,N = 100,sparse = FALSE){
  require("Matrix")
  tmp <- runif(q)
  p = tmp/sum(tmp)
  nn = round(t(rmultinom(1, N, p)))
  n = nn/N
  
  p_ab <- Matrix(0,q,q)
  if(sparse == FALSE){
    p_ab[lower.tri(p_ab, diag = FALSE)] <- runif(q*(q-1)/2,min=0,max=1)
    p_ab = p_ab + t(p_ab) + diag(runif(q,min=0,max=1)) 
  }
  if(sparse == TRUE){
    p_ab[lower.tri(p_ab, diag = FALSE)] <- runif(q*(q-1)/2,min=0,max=1/N)
    p_ab = p_ab + t(p_ab) + diag(runif(q,min=0,max=1/N)) 
  }
  
  # construct Node matrix A 
  A = Matrix(0,N,N)
  # a is a vector for each block's number
  a = cbind(0,nn)
  for(k in q:1){
    a[k + 1] = sum(a[1:(k+1)])
  }
  
  for(i in 1:q){ # for different group
    for(j in 1:q){
      A[(a[i] + 1):a[i+1], (a[j] + 1):a[j+1]] <- Matrix(rbinom(nn[i]*nn[j],1,p_ab[i,j]),nrow= nn[i],ncol = nn[j])
    }
  }
  diag(A)<-0
  A[lower.tri(A, diag = FALSE)] <- 0
  A = A + t(A)
  data = list("Probability Matrix" = p_ab,"adjacency matrix" = A,"Group Probability" = n)
  return(data)
}

#' shuffle the square Matrix

shuffle <- function(A){
  N = dim(A)[1]
  shuffle <- sample(1: N)
  return(A[shuffle,shuffle])
}

#' generate SBM graph by igraph package
rNetwork <- function(n, pi, alpha = c(1)) {
  require("igraph")
  communities <- igraph::sample_sbm(n, pi, alpha)
}

#' Generate Probability Matrix

ProbMatrix <- function(b,N = 100){

  Q = sum(b)
  P = matrix(NA,Q,Q)
  a = c(0,b)
  for(k in length(b):1){
    a[k + 1] = sum(a[1:(k+1)])
  }
  
  for(i in 1:length(b)){
    M = matrix(runif(b[i]^2,min = 1/N),b[i],b[i])
    M = (M + t(M))/2
    P[(a[i]+1):(a[i+1]),(a[i]+1):(a[i+1])] <- M
  }
 
  P[is.na(P)] = runif(length(P[is.na(P)]),min= 0,max = 1/N)
  P = (P + t(P))/2
}

### clear boundary
ProbMatrix2 <- function(b,N = 100){
  
  Q = sum(b)
  P = matrix(NA,Q,Q)
  a = c(0,b)
  for(k in length(b):1){
    a[k + 1] = sum(a[1:(k+1)])
  }
  
  for(i in 1:length(b)){
    M = matrix(runif(b[i]^2,min = 0.8 ,max = 0.8),b[i],b[i])
    M = (M + t(M))/2
    P[(a[i]+1):(a[i+1]),(a[i]+1):(a[i+1])] <- M
  }
  
  P[is.na(P)] = runif(length(P[is.na(P)]),min= 0,max = 1/N)
  P = (P + t(P))/2
}



softthreshold <- function(x, lambda) {
  sapply(x, function(x){sign(x) * max(abs(x)-lambda, 0)})
}

hardthreshold <- function(x, lambda) {
  sapply(x, function(x){ifelse(x - lambda > 0, x,0)})
}



############## self written EM algorithm for SBM! #######################################

VEM_SBM <- function(G, Q, tau.init = NULL, nstart=5, eps=1e-5, maxIter=50){
  ## Initialization
  stopifnot(is.igraph(G))
  n <- gorder(G)
  J <- vector("numeric", maxIter)
  zero <- .Machine$double.eps
  A <-  as_adj(G)
  # m1 is connected edges, m0 is disconnected edges
  m1 <- t( which( (A==1) & (upper.tri(A, diag=FALSE)), arr.ind=TRUE) )
  m0 <- t( which( (A==0) & (upper.tri(A, diag=FALSE)), arr.ind=TRUE) )
  ### variational lower bound
  get_J <- function(theta, tau){
    ## fill up below
    # theta = list(); theta$pi = pi; theta$alpha = alpha/n
    pi =theta$pi; alpha = theta$alpha
    alpha[alpha < zero] <- zero
    pi[pi < zero] <- zero
    pi[pi > 1- zero] <- 1- zero
    tau[tau<zero] <- zero
    s1 = lapply(1:Q,function(x) lapply(1:Q, function(q,l=x){sum(tau[m1[1,],q] * tau[m1[2,],l] *log(pi[q,l]))}))
    s2 = lapply(1:Q,function(x) lapply(1:Q, function(q,l=x){sum(tau[m0[1,],q] * tau[m0[2,],l] *log(1-pi[q,l]))}))
    J = sum(tau * log(alpha)) - sum(tau * log(tau)) + sum(unlist(s1)) + sum(unlist(s2))
    ## fill up above
    J
  }
  ### M step: update theta (pi and alpha)
  M_step <- function(tau){
    ## fill up below
    tau[tau<zero] <- zero
    
    s1 = lapply(1:Q,function(x) lapply(1:x, function(q,l=x){sum(tau[m1[1,],q] * tau[m1[2,],l])}))
    s2 = lapply(1:Q,function(x) lapply(1:x, function(q,l=x){sum(tau[m0[1,],q] * tau[m0[2,],l])}))
    
    pi = matrix(0,Q,Q)
    pi[upper.tri(pi,diag = TRUE)] <- unlist(s1)/(unlist(s2) + unlist(s1))
    pi[lower.tri(pi)] <- t (pi[upper.tri(pi,diag = FALSE)])
    
    ## careless!! it is even not symmetric..
    #s1 = lapply(1:Q,function(x) lapply(1:Q, function(q,l=x){sum(tau[m1[1,],q] * tau[m1[2,],l])}))
    #s2 = lapply(1:Q,function(x) lapply(1:Q, function(q,l=x){sum(tau[m0[1,],q] * tau[m0[2,],l])}))
    #pi <- matrix(unlist(s1)/(unlist(s2) + unlist(s1)),Q,Q)
    
    alpha = colMeans(tau)
    ## fill up above
    alpha[alpha < zero] <- zero
    pi[pi < zero] <- zero
    pi[pi > 1- zero] <- 1- zero
    list(pi = pi, alpha = alpha)
  }
  ### E step: update the clustering parameters (tau)
  E_step <- function(theta, tau){
    ## fill this up
    pi =theta$pi; alpha = theta$alpha
    alpha[alpha < zero] <- zero
    pi[pi < zero] <- zero
    pi[pi > 1- zero] <- 1- zero
    tau[tau<zero] <- zero
    dbernoulli <- function(pi,X = A){pi^X*(1-pi)^(1-X)}
    
    s = exp(Reduce(cbind,lapply(1:Q,function(q) Reduce("+",lapply(1:Q, function(l) colSums(log(dbernoulli(pi[q,l])) * tau[,l]))) + log(alpha[q]) )))
    s / rowSums(s)
  }
  VEM <- function(tau) {
    cond <- FALSE; iter <- 0
    while(!cond){
      iter <- iter +1
      # M step
      theta <- M_step(tau)
      # E step
      tau <- E_step(theta, tau) 
      
      J[iter] <- get_J(theta, tau) # assess convergence
      if (iter > 1)
        cond <- (iter > maxIter) | (abs((J[iter] - J[iter-1])) < eps)
    }
    return(list(theta = theta, tau = tau, J = J[1:iter]))
  }

    if(is.null(tau.init) != TRUE){
      tau <- tau.init
      best <- VEM(tau)
    }
    else{
      out <- replicate(nstart, {
        cl0 <- sample(1:Q, n, replace = TRUE)
        tau <- matrix(0,n,Q); tau[cbind(1:n, cl0)] <- 1 
        VEM(tau)
      })
      best <- out[,which.max(sapply(out['J',], max))]
    }

  ## FILL THIS UPS FOR THE BEST MODEL
  # vBIC <-
  # vICL <-
  return(list(theta = best$theta,
              tau = best$tau, membership = apply(best$tau, 1, which.max),
              J = best$J))
}

###############################################

VEM_SBM2 <- function(G, Q, tau.init, theta.init, eps=1e-5, maxIter=50){
  ## Initialization
  stopifnot(is.igraph(G))
  n <- gorder(G)
  J <- vector("numeric", maxIter)
  zero <- .Machine$double.eps
  A <-  as_adj(G)
  # m1 is connected edges, m0 is disconnected edges
  m1 <- t( which( (A==1) & (upper.tri(A, diag=FALSE)), arr.ind=TRUE) )
  m0 <- t( which( (A==0) & (upper.tri(A, diag=FALSE)), arr.ind=TRUE) )
  ### variational lower bound
  get_J <- function(theta, tau){
    ## fill up below
    # theta = list(); theta$pi = pi; theta$alpha = alpha/n
    pi =theta$pi; alpha = theta$alpha
    alpha[alpha < zero] <- zero
    pi[pi < zero] <- zero
    pi[pi > 1- zero] <- 1- zero
    tau[tau<zero] <- zero
    s1 = lapply(1:Q,function(x) lapply(1:Q, function(q,l=x){sum(tau[m1[1,],q] * tau[m1[2,],l] *log(pi[q,l]))}))
    s2 = lapply(1:Q,function(x) lapply(1:Q, function(q,l=x){sum(tau[m0[1,],q] * tau[m0[2,],l] *log(1-pi[q,l]))}))
    J = sum(tau * log(alpha)) - sum(tau * log(tau)) + sum(unlist(s1)) + sum(unlist(s2))
    ## fill up above
    J
  }
  ### M step: update theta (pi and alpha)
  M_step <- function(tau){
    ## fill up below
    tau[tau<zero] <- zero
    
    s1 = lapply(1:Q,function(x) lapply(1:x, function(q,l=x){sum(tau[m1[1,],q] * tau[m1[2,],l])}))
    s2 = lapply(1:Q,function(x) lapply(1:x, function(q,l=x){sum(tau[m0[1,],q] * tau[m0[2,],l])}))
    
    pi = matrix(0,Q,Q)
    pi[upper.tri(pi,diag = TRUE)] <- unlist(s1)/(unlist(s2) + unlist(s1))
    pi[lower.tri(pi)] <- t (pi[upper.tri(pi,diag = FALSE)])
    
    ## careless!! it is even not symmetric..
    #s1 = lapply(1:Q,function(x) lapply(1:Q, function(q,l=x){sum(tau[m1[1,],q] * tau[m1[2,],l])}))
    #s2 = lapply(1:Q,function(x) lapply(1:Q, function(q,l=x){sum(tau[m0[1,],q] * tau[m0[2,],l])}))
    #pi <- matrix(unlist(s1)/(unlist(s2) + unlist(s1)),Q,Q)

    alpha = colMeans(tau)
    ## fill up above
    alpha[alpha < zero] <- zero
    pi[pi < zero] <- zero
    pi[pi > 1- zero] <- 1- zero
    list(pi = pi, alpha = alpha)
  }
  ### E step: update the clustering parameters (tau)
  E_step <- function(theta, tau){
    ## fill this up
    pi =theta$pi; alpha = theta$alpha
    alpha[alpha < zero] <- zero
    pi[pi < zero] <- zero
    pi[pi > 1- zero] <- 1- zero
    tau[tau<zero] <- zero
    dbernoulli <- function(pi,X = A){pi^X*(1-pi)^(1-X)}
    
    s = exp(Reduce(cbind,lapply(1:Q,function(q) Reduce("+",lapply(1:Q, function(l) colSums(log(dbernoulli(pi[q,l])) * tau[,l]))) + log(alpha[q]) )))
    s / rowSums(s)
  }
  VEM <- function(theta, tau) {
    cond <- FALSE; iter <- 0
    while(!cond){
      iter <- iter +1
      
      # E step
      tau <- E_step(theta, tau) 
      # M step
      theta <- M_step(tau)
      
      J[iter] <- get_J(theta, tau) # assess convergence
      if (iter > 1)
        cond <- (iter > maxIter) | (abs((J[iter] - J[iter-1])) < eps) 
    }
    return(list(theta = theta, tau = tau, J = J[1:iter]))
  }
  
  theta = theta.init
  tau = tau.init
  best <- VEM(theta, tau)

  
  ## FILL THIS UPS FOR THE BEST MODEL
  # vBIC <-
  # vICL <-
  return(list(theta = best$theta,
              tau = best$tau, membership = apply(best$tau, 1, which.max),
              J = best$J))
}



test_modul <- function(pi,n,alpha = c(n/8,n/8,n/4,n/2)){
  zero <- .Machine$double.eps
  library("igraph")
  library("Matrix")
  
  G = rNetwork(n,pi,alpha)
  A = as_adj(G)
  wc <- cluster_walktrap(G)
  memb = c(rep(1,n/4),rep(2,n/4),rep(3,n/2))
  memb1 = memb; memb1[sample(1:n/4,1)] = 2
  memb2 = memb; memb2[sample(1:n/4,1)] = 2; memb2[sample((n/4+1):n,1)] = 1
  memb3 = memb; memb3[sample(1:n/4,1)] = 4
  modul = c(modularity(G,memb),modularity(G,memb1),modularity(G,memb2),modularity(G,memb3))
  tmp = which.max(modul)
  a =  (tmp == 1)
  b = (modularity(G,membership(wc)) <= modularity(G,memb) + zero)
  data = list(a,b)
  return(data)
}


test_modul2 <- function(p1,p0,n,alpha =  c(n/2,n/2)){
  zero <- .Machine$double.eps
  library("igraph")
  library("Matrix")
  pi = Matrix(c(p1,p0,p0,p1),2,2) 
  
  G = rNetwork(n,pi,alpha)
  #A = as_adj(G)
  wc <- cluster_walktrap(G)
  memb = c(rep(1,n/2),rep(2,n/2))
  memb1 = memb; memb1[sample(1:n/2,1)] = 2
  memb2 = memb; memb2[sample(1:n/2,1)] = 2; memb2[sample((n/2+1):n,1)] = 1
  memb3 =memb; memb3[sample(1:n/2,1)] = 3
  
  modul = c(modularity(G,memb),modularity(G,memb1),modularity(G,memb2),modularity(G,memb3))
  tmp = which.max(modul)
  
  a =  (tmp == 1)
  b = (modularity(G,membership(wc)) <= modularity(G,memb) + zero)
  data = rbind(a,b)
  return(data)
}


#' evaluation of computation time 
#' @param n number of nodes
#' @param b b = c(2,3,4) means there are 3 groups, the first group has 2 communities etc.
#' @export t return computation time for our method(divide and conquer) and original method

mixer_time <- function(n, b){
  library("igraph")
  library("mixer")
  library("Matrix")
  library(nnet)
  library(plyr)
  source("SBMfunction.R")
  # assign block structure
  #b = c(2,2,2,2,2,2,2,2)
  #pi is the probability matrix
  
  pi = Matrix(ProbMatrix(b,N = 100)) 
  #pi <- Matrix(c(0.4,0.02,0.02,0.02,0.3,0.02,0.02,0.02,0.3),3,3)
  # n is the number of nodes
  #n = 240
  
  #alpha = rep(100,10) # length should be equal to dim(pi)[1]
  alpha = rep(n/sum(b),sum(b))
  
  stopifnot(length(alpha) == dim(pi)[1])
  
  G = rNetwork(n,pi,alpha)
  A = as_adj(G)
  #AA = shuffle(A)
  A_edge = as_data_frame(G, what="edges")
  # shuffle the edge
  AA_edge = A_edge[sample(1:dim(A_edge)[1]),]
  GG = graph_from_data_frame(AA_edge,directed = FALSE)
  AA = as_adj(GG)
  
  tmp = proc.time()
  mixer(matrix(AA,n,n),qmin= sum(b),qmax= sum(b))->xout # here we cannot use Matrix structure
  t1 = proc.time() - tmp
  m <- getModel(xout)
  #plot(xout)
  
  #m$q
  pi_1 = Matrix(m$Pis)
  #m$alphas
  
  round(pi_1,5)
  image(pi_1,main = "Estimated Probabilty Matrix")
  
  # divide by block
  
  tmp = proc.time()
  #wc <- cluster_walktrap(GG)
  fc <- cluster_fast_greedy(GG); wc = fc 
  AA_div=list()
  xout_div = list()
  m_div = list()
  pi2_div = list()
  for(i in 1:max(membership(wc))){
    AA_div[[i]] = AA[membership(wc) == i,membership(wc) == i]
    mixer(matrix(AA_div[[i]],sizes(wc)[i],sizes(wc)[i]),qmin=min(b),qmax=max(b))->xout_div[[i]]
    m_div[[i]] <- getModel(xout_div[[i]])
    #m_div[[i]] <- getModel(xout_div[[i]],q = sizes(wc)[[i]]/40)
    pi2_div[[i]] = m_div[[i]]$Pis
    #plot(xout_div[[i]])
  }
  
  pi_2 = Reduce(PutDiag,pi2_div)
  t2 = proc.time() - tmp
  t = c(t1[3],t2[3])
  names(t) <- c("whole","div")
  return(t)
}

