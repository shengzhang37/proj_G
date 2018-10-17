
#' permutation the whole P as the same order as divided-and-combined P 
#' @param 
permutation_P <- function(P_whole, mem_whole, mem_div){
  n = length(mem_whole) 
  stopifnot(n == length(mem_div) )
  K = dim(P_whole)[1]
  mem_tmp = rep(0,n)
  mem_com = rbind(mem_div,mem_whole,mem_tmp)
  
  ## for the first time we match the majority, dont allow duplication
  
  count_div = count(mem_div)
  sort_div = sort.list(count_div[,2],decreasing = T)
  
  for(i in 1:K){
    df = count( mem_whole[mem_div == sort_div[i]])
    #df[which.max(df$freq),1]
    if(sum(mem_com[3,mem_whole == df[which.max(df$freq),1]]) == 0){
      mem_com[3,mem_whole == df[which.max(df$freq),1]] <- sort_div[i]}
    else if(!is.na( df$freq[2]) & sum(mem_com[3,mem_whole == df[which(df$freq == sort(df$freq,decreasing = T)[2])[1],1]]) == 0)
    {
      mem_com[3,mem_whole == df[which(df$freq == sort(df$freq,decreasing = T)[2])[1],1]] <- sort_div[i]}
    else if(!is.na( df$freq[3]) & sum(mem_com[3,mem_whole == df[which(df$freq == sort(df$freq,decreasing = T)[3])[1],1]]) == 0)
    {
      mem_com[3,mem_whole == df[which(df$freq == sort(df$freq,decreasing = T)[3])[1],1]] <- sort_div[i]}
  }
  
  ## for the second time, we assign the group to the nearest unassigned group
  mem_com2 = mem_com[,which(mem_com[3,]==0)]
  mem_div2 = mem_com2[1,]
  mem_whole2 = mem_com2[2,]
  count_div2 = count(mem_div2)
  sort_div2 = count_div2[sort.list(count_div2[,2],decreasing = T),1]
  
  mem_series = 1:K
  mem_series_left = mem_series[!mem_series %in% mem_com[3,]]
  for(i in 1:length(sort_div2)){
    df = count( mem_whole2[mem_div2 == sort_div2[i]])
    #df[which.max(df$freq),1]
    if(sum(mem_com2[3,mem_whole2 == df[which.max(df$freq),1]]) == 0){
      mem_com2[3,mem_whole2 == df[which.max(df$freq),1]] <- mem_series_left[which.min(abs(mem_series_left - sort_div2[i]))]
      mem_series_left = mem_series_left[-which.min(abs(mem_series_left - sort_div2[i])) ]
    }
  }
  mem_com[,which(mem_com[3,]==0)] = mem_com2
  ## group K-th is left, assign the rest one to it
  mem_com[3,mem_com[3,] == 0] = K
  
  ############## permutation ############################
  
  whole_after = sapply(1:K, function(ss)  mean(mem_com[3,mem_com[2,] == ss]))
  
  perm = rbind(1:K,whole_after)
  row.names(perm) <- c("whole_before","whole_after")
  
  
  pi_1_perm = P_whole
  
  for(i in perm[1,]){
    for(j in perm[1,]){
      pi_1_perm[perm[2,i],perm[2,j]] <- P_whole[i,j]
    }
  }
  
  return(pi_1_perm)
  
}



