P_Matrix_fromclustering <- function(mem,AA){
  
  P = Matrix(NA,max(mem),max(mem))
  
  for(i in 1:max(mem)){
    for(j in 1:max(mem)){
      P[i,j] = mean(AA[mem == i,mem == j])
    }
  }
  return(P)
}