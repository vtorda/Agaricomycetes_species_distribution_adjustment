#!/usr/bin/Rscript
library(ape)
data <- read.table("kisb_species2.txt", header = T, row.names = 1, sep = "\t")
tree  <- read.tree("245_rooted2_partial_forTreePL_onlyone.tre")
data1 <- data[,-c(which(names(data) == "Dacrymyces"),which(names(data) == "Sistotrema"))] # exclude non-monophyletic groups
data2 <- data1
tree_specnumber <- length(tree$tip.label)
# one tailed biomial functions
binom_fun <- function(q, m, n, k) {phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE)}
binom_fun2 <- function(q, m, n, k) {phyper(q = q, m = m, n = n, k = k, lower.tail = TRUE)}
species_number <- data2[1,]
big_spec_number <- data2[2,]
ratio <- data2[2,] / 35471
samplesize <- tree_specnumber - sum(data2[1,]) + sum(data2[1,])
# crate a data frame with the starting species distribution and p-values
data2[3,] <- mapply(FUN = binom_fun, q = species_number, m = big_spec_number, n = 35471 - big_spec_number, MoreArgs = list(k = samplesize))
data2[3,] <- p.adjust(data2[3,], method = "BH")
data2[4,] <- mapply(FUN = binom_fun2, q = species_number, m = big_spec_number, n = 35471 - big_spec_number, MoreArgs = list(k = samplesize))
data2[4,] <- p.adjust(data2[4,], method = "BH")
data3<-data2

counter<-0 
while (any(data3[3,] < 0.05)) 
{
  resample <- 0 
  c1 <- rep(colnames(data3[3, data3[3,] < 0.05]), data3[1, data3[3,] < 0.05]) 
  i <- which(colnames(data3) == sample(c1,1)) 
  cat("firsti\t", names(data3)[i], "\t", data3[1,i], "\n")
  
  if (data3[1,i] > 1)
  {
    if (data3[1,i] / (tree_specnumber - sum(data2[1,]) + sum(data3[1,])) > data3[2, i] / 35471) 

    {
      data3[1,i] <- data3[1,i] - 1
      cat(names(data3)[i], "\n")
    }
    else 
    {
      cat("else", "\n")
      if (data3[1,i] < data2[1,i]) 
      {
        data3[1,i] <- data3[1,i] + 1
        cat("#####\nadd1\n#####\n") 
      }
    }
  }
  else
  {
   while (data3[1,i] <= 1)
    {
     cat("resample", "\t", resample, "\n")
      i <- which(colnames(data3) == sample(c1,1)) 
      resample <- resample + 1
      if (resample > tree_specnumber) break;
    } 
       if (data3[1,i] / (tree_specnumber - sum(data2[1,]) + sum(data3[1,])) > data3[2, i] / 35471) 

    {
      if (data3[1,i] > 1)
      {
        data3[1,i] <- data3[1,i] - 1
        cat("resamplenames", "\t", names(data3)[i], "\n")
      }
    }
    else 
    {
      if (data3[1,i] < data2[1,i]) 

      {
        data3[1,i] <- data3[1,i] + 1
        cat("#####\nadd1\n#####\n") 
      }
    }
  }
  data3[4,] <- mapply(FUN = binom_fun2, q = species_number, m = big_spec_number, n = 35471 - big_spec_number, 
                      MoreArgs = list(k = samplesize))
  data3[4,] <- p.adjust(data3[4,], method = "BH")
  
  if (any(data3[4,] < 0.05))
  {
    c2 <- which(data3[4,] < 0.5 & data3[1,] < data2[1,])
    if (length(c2) != 0)
    {
     data3[1,c2] <- data3[1,c2] + 1 
     cat("#####\nadd3\n#####\n")
    }
  }
  
  species_number <- data3[1,]
  samplesize <- tree_specnumber - sum(data2[1,]) + sum(data3[1,])
  data3[3,] <- mapply(FUN = binom_fun, q = species_number, m = big_spec_number, n = 35471 - big_spec_number, 
                      MoreArgs = list(k = samplesize))
  data3[3,] <- p.adjust(data3[3,], method = "BH")
  
  cat(length(which(data3[3,] > 0.05)),
      "\t", 
      length(which(data3[3,] < 0.05 & data3[1,i] / (tree_specnumber - sum(data2[1,]) + sum(data3[1,])) > data3[2, i] / 35471)),
      "\n")
  
  if (any(data3[3,] > 0.15))
  {
    numbers <- which(data3[3,] > 0.3 & data3[1,] < data2[1,])
    if (length(numbers > 0))
    {
      data3[1,numbers] <- data3[1,numbers] + 1
      cat("#####\nadd4\n#####\n")
    }
    
  }
  if (counter>0)
  {
    number2 <- which(data3[3,] < 0.05 & data3[1,] / (tree_specnumber-sum(data2[1,])+sum(data3[1,])) >= ratio & data3[1,] < data2[1,])
    if (length(number2) == 0) break;
    if (resample > tree_specnumber) break;
  }
  
  counter<-counter+1
}
save(data3, file="delfajok0208_3")
warnings()