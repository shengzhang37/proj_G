library(plotROC)
library(ggplot2)
library(pROC)

n = 400
stopifnot(n <= dim(AA)[1])
mem_div = mem

sub_index =  base::sample(1:(dim(AA)[1]-1),(n-1), replace = FALSE)
edge_row = lapply(sub_index, function(i) AA[i,-(1:i)])
prob2_row = lapply(sub_index, function(i){ j = (i + 1) : dim(AA)[1]
sapply(j, function(j)  pi_2[mem_div[i],mem_div[j]])
}
)

edge = unlist(edge_row) # 把matrix变成一个array
prob2 = unlist(prob2_row)

test <- data.frame(D = edge,
                   A_div =prob2,  stringsAsFactors = FALSE)
basicplot <- ggplot(test, aes(d = D, m = A_div)) + geom_roc()
basicplot # shell中输入这个会画出图来

roc_obj2 <- roc(edge, prob2)
auc(roc_obj2) # 返回结果


