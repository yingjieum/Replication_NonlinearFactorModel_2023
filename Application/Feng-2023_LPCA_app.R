######################################################
#### This file replicates the empirical application ##
############## in Feng (2023)  #######################
########### Author: Yingjie Feng #####################
########## Last updated: 11/12/2023 ##################

rm(list = ls())
library(Rfast)
library(irlba)
library(RcppNumerical)
library(ggplot2)
library(CVXR)
source("Feng-2023_LPCA_app_suppfuns.R")

####################################################
#### Kansas example ################################
load("kansas.rda")

x <- as.data.frame(kansas[,c("abb","lngdpcapita","year_qtr")])
kansas.ts <- as.matrix(x[x[,1]=="KS",-1])
kansas.ts[,1] <- c(NA, diff(kansas.ts[,1])) * 100 
kansas.ts <- kansas.ts[-1,]
time <- kansas.ts[,2]

treated.ind <- kansas$treated[kansas$abb=="KS"]
T0 <- which.max(treated.ind) - 1

x[kansas$treated==1,2] <- NA
x <- reshape(x, idvar = "abb", timevar = "year_qtr", direction = "wide")
name <- x[,1]
id <- which(name=="KS")
x <- as.matrix(x[,-1])
x <- t(apply(x, 1, diff) * 100)

col.mean <- colMeans(x, na.rm=T)
x <- sweep(x, 2, col.mean)
n <- nrow(x); p <- ncol(x); K <- round(n^(2/3)); p1 <- 40; nlam <- 3
x[id, T0:p] <- 0

fullkmat  <- knn.index(A=x[,1:p1], Kseq=K)
mean.pred <- lpca(A=x[,(p1+1):p], index=fullkmat[,id], nlam=nlam, n=n, K=K, i=id)
mean.lpca <- mean.pred[1:(p-p1)] + col.mean[(p1+1):p]
gpca <- irlba(x, 3)
gpca.fit <- gpca$u[,1:2] %*% diag(gpca$d[1:2])%*% t(gpca$v[,1:2])
gpca.fit <- gpca.fit[id, (p1+1):p] + col.mean[(p1+1):p]
sc.fit   <- sc.pred(x, ind=id, post=(T0:p)); sc.fit <- sc.fit[(p1+1):p,]

count <- sum(mean.lpca[(T0-p1):(p-p1)] - kansas.ts[T0:p,1]>0)
# average post-treatment effect
ave.diff.lpca <- mean(kansas.ts[T0:p]) - mean(mean.lpca[(T0-p1):(p-p1)])
ave.diff.sc   <- mean(kansas.ts[T0:p]) - mean(sc.fit[(T0-p1):(p-p1)])

# Plotting
plot     <- ggplot()
plot.fit <- data.frame(y=mean.lpca, t=time[(p1+1):p], name="LPCA") 
plot.y   <- data.frame(y=kansas.ts[(p1+1):p,1], t=time[(p1+1):p], name="Kansas")
plot.sc <- data.frame(y=sc.fit, t=time[(p1+1):p], name="SC") 


plot <- plot + geom_line(data=plot.fit, aes(x=t, y=y, colour=name), linetype="dashed", size=.7)
plot <- plot + geom_line(data=plot.y,   aes(x=t, y=y, colour=name), linetype="solid", size=.7)
plot <- plot + geom_line(data=plot.sc,aes(x=t, y=y, colour=name), linetype="dotted", size=.7)
plot <- plot + scale_colour_manual(name="", values=c("black", "blue", "coral3"),
                                   guide=guide_legend(override.aes = list(
                                     linetype=c("solid", "dashed", "dotted"))))
plot <- plot + geom_vline(xintercept = 2012, linetype="dashed", color="black") + 
        xlab("time") + ylab("GDP per capita, growth rate") + 
        theme_bw() + 
        theme(legend.position = c(0.1,0.9),
              legend.background = element_rect(fill="transparent"), 
              panel.grid.minor = element_blank(),panel.grid.major = element_blank())
ggsave("kansas.growth.pdf", width = 5.5, height = 4)

######################
# recover the gdp path
gdp0 <- kansas$lngdpcapita[kansas$abb=="KS"][T0]
gdp.pre <- kansas$lngdpcapita[kansas$abb=="KS"][(p1+2):T0]
level.lpca <- c(rep(NA, T0-p1-2), gdp0, cumsum(mean.lpca[(T0-p1):(p-p1)])/100 + gdp0) 
level.gpca <- c(rep(NA, T0-p1-2), gdp0, cumsum(gpca.fit[(T0-p1):(p-p1)])/100 + gdp0)


x <- as.data.frame(kansas[,c("abb","lngdpcapita","year_qtr")])
x <- as.matrix(reshape(x, idvar = "abb", timevar = "year_qtr", direction = "wide")[,-1])
sc.fit.level <- sc.pred(x, ind=id, post=(T0+1):(p+1)); sc.fit.level <- sc.fit.level[(p1+2):(p+1),]

# Plotting
plot     <- ggplot()
plot.fit <- data.frame(y=level.lpca, t=time[(p1+1):p], name="LPCA") 
plot.y   <- data.frame(y=kansas$lngdpcapita[kansas$abb=="KS"][(p1+2):(p+1)], t=time[(p1+1):p], name="Kansas")
plot.sc <- data.frame(y=sc.fit.level, t=time[(p1+1):p], name="SC") 

plot <- plot + geom_line(data=plot.fit, aes(x=t, y=y, colour=name), linetype="dashed", size=.7)
plot <- plot + geom_line(data=plot.y,   aes(x=t, y=y, colour=name), linetype="solid", size=.7)
plot <- plot + geom_line(data=plot.sc,aes(x=t, y=y, colour=name), linetype="dotted", size=.7)
plot <- plot + scale_colour_manual(name="", values=c("black", "blue", "coral3"), 
                                   guide=guide_legend(override.aes = list(
                                     linetype=c("solid", "dashed", "dotted"))))
plot <- plot + geom_vline(xintercept = 2012, linetype="dashed", color="black") + 
        xlab("time") + ylab("GDP per capita") + 
        theme_bw() + 
        theme(legend.position = c(0.1,0.9),
              legend.background = element_rect(fill="transparent"), 
              panel.grid.minor = element_blank(),panel.grid.major = element_blank())
ggsave("kansas.level.pdf", width = 5.5, height = 4)

