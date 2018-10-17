D.ex = whole_edge$edge_row
M1 = whole_edge$prob_row
M2 = div_edge$prob_row
M3 = clu_edge$prob_row

test <- data.frame(D = D.ex, 
                   M1 = M1, M2 = M2, stringsAsFactors = FALSE)
basicplot <- ggplot(test, aes(d = D, m = M2)) + geom_roc()
basicplot

styledplot <- basicplot + style_roc()
styledplot

styledplot + geom_rocci()

head(test)

longtest <- melt_roc(test, "D", c("M1", "M2"))
head(longtest)

ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc() + style_roc()


roc_obj1 <- roc(D.ex, M1)
auc(roc_obj1) ## 0.993
roc_obj2 <- roc(D.ex, M2)
auc(roc_obj2) ## 0.975

roc_obj3 <- roc(D.ex, M3)
auc(roc_obj3)  ## 0.9314





