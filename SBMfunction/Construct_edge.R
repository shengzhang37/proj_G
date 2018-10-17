Construct_General_P <- function(mem, P, A ,n = 100){
  n = dim(A)[1]
  #construct a special vector
  cata_prob = list()
  for(i in 1 : max(mem)){
    cata_prob[[i]] = rowSums(sapply(1:max(mem), function(j) (mem == j) * P[i,j] )) 
  }
  
  prob_row = unlist( lapply(1:n, function(i) cata_prob[[mem[i]]] ))
  edge_row = as.vector(t(A))
  data = list(prob_row = prob_row,edge_row = edge_row )
  return(data)
}


whole_edge = Construct_General_P(mem_whole,pi_1,AA)
div_edge = Construct_General_P(mem_div,pi_2,AA)

mem_cluster = membership(fc)
pi_3 = P_Matrix_fromclustering(mem_cluster,AA)
clu_edge = Construct_General_P(mem_cluster,pi_3,AA)


mem = mem 
P = pi_2
A = AA

edge_tmp = Construct_General_P(mem,P,A)

test <- data.frame(D = edge_tmp$edge_row,
                   A_div = edge_tmp$prob_row,  stringsAsFactors = FALSE)
basicplot <- ggplot(test, aes(d = D, m = A_div)) + geom_roc()
basicplot # shell中输入这个会画出图来

library(plotROC)
library(ggplot2)
library(pROC)

roc_obj2 <- roc(edge_tmp$edge_row, edge_tmp$prob_row)
auc(roc_obj2) # 返回结果
