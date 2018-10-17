require("igraph")
require("mixer")
require("Matrix")
library(nnet)
library(plyr)
library(plotROC)
library(ggplot2)

################################ Our Method #####################################
# config
Data_dir = "/Users/shengzhang/Desktop/third year I/SBM/Data/"
Function_dir = "/Users/shengzhang/Desktop/third year I/SBM/SBMfunction/"
Data_filename = "facebook_combined.txt"

# source useful function
source(paste0(Function_dir,"SBMfunction.R"))


# Read Data
AA_edge = read.table(paste0(Data_dir, Data_filename))
GG = igraph::graph_from_data_frame(AA_edge,directed = FALSE)
AA = as_adj(GG)


# use Graph_Div() to divide the graph, we can tune parameters
graph_div = Graph_Div(GG,upper_size = 800,lower_size = 30, lower_M = 0.03, lambda = 0.1)
mem = graph_div$mem
fc = graph_div$fc

# check the reult 
table(mem); max(table(mem));min(table(mem))
modularity(GG,mem) ## should be close(even better) to original fc 
modularity(GG,membership(fc)); table(cut_at(fc ,no =7))

# mixer_div
pi2_div = m_div = AA_div= xout_div = list()

for(i in 1:max(mem)){
  AA_div[[i]] = AA[mem == i,mem == i]
  mixer(matrix(AA_div[[i]],table(mem)[i],table(mem)[i]),qmin=2,qmax=30)->xout_div[[i]]
  m_div[[i]] <- getModel(xout_div[[i]])
  #m_div[[i]] <- getModel(xout_div[[i]],q = sizes(wc)[[i]]/40)
  pi2_div[[i]] = m_div[[i]]$Pis
  #plot(xout_div[[i]])
}

pi_2 = Reduce(PutDiag,pi2_div)

image(pi_2,main = "Estimated Probabilty Matrix(clustered)")

################################# Louvian ################################

pi_3 = P_Matrix_fromclustering(membership(fc),AA)
image(pi_3,main = "Estimated Probabilty Matrix(Louvian)")

################################# ROC ###################################






