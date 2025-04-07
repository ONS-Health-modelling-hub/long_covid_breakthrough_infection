cov.dist.cat <- function(vars, dataset, exposure) {

out_list1 <- as.list(NULL)

for(i in 1:length(vars)) {
  
  dataset[[vars[i]]] <- as.factor(dataset[[vars[i]]])
  
  ds0 <- dataset[dataset[[exposure]]==0,]
  ds1 <- dataset[dataset[[exposure]]==1,]
  
  n_all <- table(dataset[[vars[i]]])
  n0 <- table(ds0[[vars[i]]])
  n1 <- table(ds1[[vars[i]]])
  
  p_all <- n_all / nrow(dataset)
  p0 <- n0 / nrow(ds0)
  p1 <- n1 / nrow(ds1)
  
  std_diff <- (p0 - p1) / sqrt((p0*(1-p0) + p1*(1-p1)) / 2)
  
  out_list1[[i]] <- as.data.frame(cbind(c(vars[i], rep("", nlevels(dataset[[vars[i]]])-1)),
                                        levels(dataset[[vars[i]]]),
                                        n_all, p_all, n0, p0, n1, p1, std_diff))
  
}

out_df1 <- out_list1[[1]]
if(length(vars)>1) {
  for(i in 2:length(vars)) {out_df1 <- rbind(out_df1, out_list1[[i]])}
}

rownames(out_df1) <- NULL
colnames(out_df1)[1:2] <- c("characteristic", "level")

return(out_df1)

}
