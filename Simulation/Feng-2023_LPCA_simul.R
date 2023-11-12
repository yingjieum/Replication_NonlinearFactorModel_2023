#######################################################
#### This file replicates the simulation results ######
############## in Feng (2023)  ########################
## Optimal Est. of Large Dim. Nonlinear Factor Models #
########### Author: Yingjie Feng ######################
########## Last updated: 11/12/2023 ###################

rm(list = ls())
library(Rfast)
library(foreach)
library(doParallel)
library(doRNG)
library(irlba)
library(HDRFA)
library(RcppNumerical)
library(Hmisc)
source("Feng-2023_LPCA_simul_suppfuns.R")

# param
rep <- 2000
par <- read.csv("Feng-2023_LPCA_simul_models.csv", header = T, colClasses=c(rep("numeric", 6), "logical", "numeric"))

for (j in 1:9) {
n   <- par$n[j]
p   <- par$p[j]
model  <- funlist[[par$hdmodel[j]]]
nlam   <- par$nlam[j]
K      <- n^(2*1/(2*1+1))
Kseq   <- floor(K * c(0.5, 1, 1.5))
sigma  <- par$sigma[j]
binary <- par$binary[j]
prop.mat <- par$prop[j]

## simulation
# set # of cores
cl <- makeCluster(8)
registerDoParallel(cl)

output <- foreach (i = 1:rep, .options.RNG=1234, .packages=c('Rfast','irlba','RcppNumerical', 'HDRFA'),
                   .combine=rbind) %dorng% {
                     output <- sim(i, n, p, model, Kseq, nlam, sigma, binary, prop.mat)
                     output   # (5*rep) by (length(Kseq)+1) matrix
                   }

stopCluster(cl)

###################
write.table(output, paste("rawoutput_par", j, "txt", sep = "."), sep = ",", 
            row.names = F, col.names = F)
}


#################################
### Make tables in the paper
modlist <- list(c(1,2,3), c(4,5,6), c(7,8,9))
for (i in 1:3) {
  tab.out <- c()
  for (j in modlist[[i]]) {
    
    output <- as.matrix(read.table(paste("rawoutput_par", j, "txt", sep = "."), sep = ","))
    
    L  <- nrow(output)/rep
    aa <- cbind(output, rep(1:L, rep))
    
    ###############
    table <- rbind(
      colmeans(aa[aa[,ncol(aa)]==1,]),
      colmeans(aa[aa[,ncol(aa)]==3,]),
      colmeans(aa[aa[,ncol(aa)]==4,]),
      colmeans(aa[aa[,ncol(aa)]==5,])
    )
    
    tab.out <- rbind(tab.out, round(table[,c(1,2,3,6)], 3))
  }
  
  # Report main simulation results
  n.cgroup <- c(3, 1)
  colheads <- c("$K$=49", "$K$=99", "$K$=149", "")
  cgroup   <- c("LPCA", "GPCA")
  
  n.rgroup <- c(4, 4, 4)
  rgroup   <- c("Model 1", "Model 2", "Model 3")
  rowname  <- rep(c("MAE", "$q_\\alpha$=.1", "$q_\\alpha$=.5", "$q_\\alpha$=.9"), 3)
  latex(tab.out, file=paste("Table_MainPaper", i, ".txt", sep = ""), 
        append=FALSE, table.env=FALSE, center="none", title="",
        n.cgroup=n.cgroup, cgroup=cgroup, colheads=colheads,
        n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname
  )
}

