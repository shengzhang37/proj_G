#' Divide a graph into sub graph 
#' @param GG input graph
#' @param lower_size required lower bound of size
#' @param upper_size required upper bound of size (might not be attained because of control of Modularity)
#' @param lower_M required lower bound of Modularity for dividing sub-graph
#' @param lambda parameter for M_trend_peak function
#' @export mem assigned member for input graph
Graph_Div <- function(GG,lower_size = 50,upper_size = 2000,lower_M = 0.4,lambda = 0.01){
  GG = simplify(GG)
  AA = as_adj(GG)
  fc = cluster_fast_greedy(GG)
  m_trend = sapply(1:max(membership(fc)), function(i) {mem = cut_at(fc, no=i);modularity(GG,mem)})
  #plot(m_trend);lines(m_trend)
  n_g = M_trend_peak(m_trend,lambda)
  mem = cut_at(fc, no=n_g);table(mem);modularity(GG,mem)
  
  M_control = 1
  L_size_control = 1
  while( max(table(mem)) > upper_size  & M_control & L_size_control){
    adj = which(table(mem) > upper_size)
    n = length(adj)
    AA_div = lapply(adj, function(i) AA[mem == i,mem == i])
    GG_div = lapply(1:n, function(i) graph_from_adjacency_matrix(AA_div[[i]], mode = c("undirected")) )
    fc_div = lapply(1:n, function(i) cluster_fast_greedy(GG_div[[i]]))
    m_trend_div = lapply(1:n, function(j){
      sapply(1:max(membership(fc_div[[j]])), function(i) {mem = cut_at(fc_div[[j]], no=i);modularity(GG_div[[j]],mem)})
    } )
    n_g_div = lapply(m_trend_div, function(m_trend) M_trend_peak(m_trend,lambda = lambda))
    mem_div = lapply(1:n, function(i) cut_at(fc_div[[i]], no=n_g_div[[i]]))
    size_div = lapply(mem_div,table)
    mod_div = lapply(1:n,function(i) modularity(GG_div[[i]], mem_div[[i]]))
    
    M_control = length(which(mod_div > lower_M)) 
    L_size_control = sum( sapply(size_div, function(x) max(min(x) - lower_size, 0 )) )
    if(M_control){
      adj_bigM = adj[which(mod_div > lower_M)]
      mem_tmp <- mem
      for(i in 1:length(adj_bigM)){
        if(min(table( mem_div[[which(mod_div > lower_M)[i]]])) > lower_size){
          mem_tmp[mem_tmp == adj_bigM[i]] <- mem_div[[which(mod_div > lower_M)[i]]] +  adj_bigM[i] * 100}
      }
      ###### make discrete sequence into continous sequence ########
      df = count(mem_tmp)
      for(i in 1: dim(df)[1]){
        mem_tmp[mem_tmp == df[i,1]] = i
      }
      mem <- mem_tmp
    }
  }
  dat = list(mem = mem, fc = fc)
  return(dat)
}


#' Find the first Peak in M_trend
#' @param m_trend  Modularity trend as input
#' @param lambda if one cut has its neighbors within lambda, then we choose this cut
#' @export i as cut
M_trend_peak <- function(m_trend,lambda = 0.01){
  n = length(m_trend)
  i = 2
  diff = 1
  while(i<n & diff > lambda){
    diff =max(abs(m_trend[i+1] - m_trend[i]), abs(m_trend[i-1] - m_trend[i]))
    i = i+1 
  }
  return(i - 1)
}

#' Combine two block matrix into one
#' 
#' @param X,Y input Matrix
#' @export combined_matrix
#' @examples
#' A = Matrix(1,10,15)
#' B = Matrix(2,6,7)
#' PutDiag(A,B)
PutDiag <- function(X,Y){
  require("Matrix")
  n  = dim(X) + dim(Y)
  T = Matrix(0,n[1],n[2])
  T[1:dim(X)[1],1:dim(X)[2]]  = X
  T[(dim(X)[1]+1):n[1],(dim(X)[2]+1):n[2]] =Y
  return(T)
}
