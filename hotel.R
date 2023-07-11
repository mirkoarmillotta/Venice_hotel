
###############################################################################
## This file contains the code used in the paper "Unveiling Venice's hotel 
## competition networks from dynamic pricing digital markets" 
## by Mirko Armillotta, Konstantinos Fokianos and Andrea Guizzardi.
##
##
## Note: R version 4.2.1 used throughout.
##
###############################################################################


## Load packages

library(ks)
library(igraph)
library(pheatmap)
library(mvtsplot)
library(geodist)
library(crch)
library(car)
library(sf)
library(sp)
library(leaflet)

source("network estimation.R")      # Advance booking Network Autoregressive model (NAR) estimation (Algorithm 1)
source("adja_dist.R")               # Compute adjacency matrix from given probabilities
source("ols.nar.R")                 # Compute OLS of advance booking NAR model
source("genNAR.R")                  # Generate values from advance booking NAR model
source("data format.R")             # Format the data for the analysis


## Compute price distance

price_d <- aggregate( Price ~ id_hotel , Data, median)
price_dd <- data.matrix(price_d)
price_dd <- log(price_dd[,2])
price_dist <- abs(outer(price_dd, price_dd, '-'))


###################################################
## Format the data as follows:
## covariates in matrix z
## hotels per rows, time per column in matrix ho
## one ho matrix for each advance booking 
##################################################

ADV <- data_format(Data=Data)

View(ADV$ho[[1]])  # BAR for advance booking 0 for each hotel (by row) over time (by column)
View(ADV$geodist)  # geodesic distance between pairs of hotels
View(ADV$z)        # covariates

dim <- nrow(ADV$ho[[1]])         ## number of hotels
size <- ncol(ADV$ho[[1]])        ## temporal size
K <- length(ADV$ho)-1            ## maximal advance booking


###################### Descriptive statistics (Section 3) ##############################

# Table 1

BAR <- list()
summ <- matrix(0, K, 4)
for (k in 1:(K)){
  
  BAR[[k]] <- exp(ADV$ho[[k]])
  per5 <- quantile(BAR[[k]], 0.05)
  mm <- median(BAR[[k]])
  mea <- mean(BAR[[k]])
  per95 <- quantile(BAR[[k]], 0.95)
  summ[k,] <- c(per5, mm, mea, per95)
}
rownames(summ) <- seq(from=0, to=K-1, by=1)

View(summ)

# Figure 1
mvtsplot(t(ADV$ho[[1]]), sort=median)

# Figure A-1
mvtsplot(t(ADV$ho[[1]]), smooth.df = 15, sort = median)


####### Heatmap according to standard deviations (Figure 3)

sd_i <- matrix(0, dim, K)

for(k in 1:K){
  sd_i[,k] <- apply(ADV$ho[[k]],1, sd)
}


colnames(sd_i) <- seq(from=0, to=13, by=1)

outphe <- pheatmap(sd_i, scale="none", border_color=FALSE)
outphe$tree_row$order
clust_row <- cutree(outphe$tree_row, k=2)

outphe$tree_col$order
clust_col <- cutree(outphe$tree_col, k=3)

# Table 2-3

lattice <- ADV$geodist
lattice[lattice>200] <- 0
lattice[lattice!=0] <- 1

Desc <- cbind(ADV$z, rowSums(lattice))
summary(Desc)

Descr <- cbind(Desc, clust_row)
summary(Descr[Descr[,5]==1,c(1:4)])
summary(Descr[Descr[,5]==2,c(1:4)])

HO <- ADV$ho[[1]]
for(k in 2:15){
  HO <- cbind(HO, ADV$ho[[k]])
}
descr_ho <- cbind(exp(HO), clust_row)
descr_hoc1 <- descr_ho[descr_ho[,5161]==1,-5161]
descr_hoc2 <- descr_ho[descr_ho[,5161]==2,-5161]
summary(vec(descr_ho))
summary(vec(descr_hoc1))
summary(vec(descr_hoc2))

## Only from adv. book. 0 to 3

HO4 <- ADV$ho[[1]]
for(k in 2:4){
  HO4 <- cbind(HO4, ADV$ho[[k]])
}
descr_ho4 <- cbind(exp(HO4), clust_row)
descr_hoc14 <- descr_ho4[descr_ho4[,1377]==1,-1377]
descr_hoc24 <- descr_ho4[descr_ho4[,1377]==2,-1377]
summary(vec(descr_ho4))
summary(vec(descr_hoc14))
summary(vec(descr_hoc24))

##### plot cluster analysis (Figure 2)


plot(outphe$tree_row, labels=FALSE)
plot(outphe$tree_col)




#################### Compute competition probabilities (Section 4.2) ##############
#################### by Maximum Likelihood Estimation                ##################


d_g <- ADV$geodist[lower.tri(ADV$geodist, diag = FALSE)]  # geodesic distances
d_p <- price_dist[lower.tri(price_dist, diag = FALSE)]    # price distances


###### Truncated ligistic MLE for geodeisc distances

TLog = function(theta, y){
  m=exp(theta[1])
  s=exp(theta[2])
  
  l = dtlogis(x=y, location = m, scale = s, left = 0, right = 500, log = F)
  
  llik=-sum(log(l))
  llik
}

d_g[d_g==0] <- min(d_g[d_g!=0])/2 # avoid log(0) problem
d_g = d_g[d_g<=500]

x0 = c(mean(d_g), var(d_g))
x0 = log(x0)

TLog(theta=x0, y=d_g)

est_TL=optim(par=x0, fn=TLog, y=d_g, method="BFGS", hessian=T)
est_TL

tpar = exp(est_TL$par)
tpar                    # estimated parameters mu and lambda

plot(ecdf(d_g), xlim=c(0,500))
par(new=TRUE)

numtlr_d_g <- 1/(1+exp(-1*(d_g-tpar[1])/tpar[2])) - 1/(1+exp(-1*(0-tpar[1])/tpar[2]))
dentlr_d_g <- 1/(1+exp(-1*(500-tpar[1])/tpar[2])) - 1/(1+exp(-1*(0-tpar[1])/tpar[2]))
ptlr_d_g <- 1- numtlr_d_g/dentlr_d_g

plot(1-ptlr_d_g~d_g, col="orange", ylim=c(0,1), xlim=c(0,500))

ks.test(d_g, ptlogis, location=362.4541, scale = 124.4054, left = 0, right = 500)

####### Gamma MLE for price distances

d_p[d_p==0] <- min(d_p[d_p!=0])/2 # avoid log(0) problem

Gammal = function(theta, y){
  a=exp(theta[1])
  b=exp(theta[2])
  
  llik = sum(dgamma(x=y, shape = a, rate = b, log = T))
  -llik
}


x0 = c(1.22, 2.21)
x0 = log(x0)

Gammal(theta=x0, y=d_p)

est_G=optim(par=x0, fn=Gammal, y=d_p, method="BFGS", hessian=T)
est_G

exp(est_G$par)  # estimated parameters eta and nu 

a_shape_p <- 1.224101 
th_scale_p = 2.213598

plot(ecdf(d_p), xlim=c(0,2.5))
par(new=TRUE)
plot(pgamma(d_p, shape=a_shape_p, rate=th_scale_p)~d_p, col="pink", ylim=c(0,1), xlim=c(0,2.5))


ks.test(d_p, pgamma, shape=a_shape_p, rate=th_scale_p)

## Figure 4

plot(as.numeric(d_g), 1-ptlr_d_g, col="red", ylim=c(0,1), xlim=c(0,500),
     xlab="geodesic distances", ylab="F(x)",pch=20)
par(new=TRUE)
plot(ecdf(d_g), xlim=c(0,500),
     xlab="", ylab="", main="",
     verticals=T, do.points=F)
legend(0, 0.98, legend=c("c.d.f.", "e.c.d.f."),
       col=c("red", "black"), lty = 1, cex=0.5)
plot(pgamma(d_p, shape=a_shape_p, rate=th_scale_p)~d_p, col="red",
     ylim=c(0,1), xlim=c(0,2.5), xlab="price distances", ylab="", pch=20)
par(new=TRUE)
plot(ecdf(d_p), xlim=c(0,2.5), xlab="", ylab="", main="")

##### final probabilities (eq. (4))

p_clogp <- 1-pgamma(price_dist, shape=a_shape_p, rate=th_scale_p) 
diag(p_clogp) = 0

numtlr_d_g <- 1/(1+exp(-1*(ADV$geodist-tpar[1])/tpar[2])) - 1/(1+exp(-1*(0-tpar[1])/tpar[2]))
dentlr_d_g <- 1/(1+exp(-1*(500-tpar[1])/tpar[2])) - 1/(1+exp(-1*(0-tpar[1])/tpar[2]))
ptlr_d_g = numtlr_d_g/dentlr_d_g
ptlr_d_g[ptlr_d_g>1] = 1
p_cg <- 1 - ptlr_d_g
diag(p_cg) = 0

p_clog <- p_clogp * p_cg   # final probabilities
diag(p_clog) <- 0



################## Estimation results (Section 5) #######################

#### normalize regressors

lattice <- ADV$geodist
lattice[lattice>200] <- 0
lattice[lattice!=0] <- 1

zz <- ADV$z[,1:2]
n_room_s <- (ADV$z[,3]-mean(ADV$z[,3]))/sd(ADV$z[,3])
n_hotels_s <- rowSums(lattice)
n_hotels_s <- (n_hotels_s-mean(n_hotels_s))/sd(n_hotels_s)
zz <- cbind(zz,n_room_s,n_hotels_s)
View(zz)


#### Estimation (Algorithm 1)

set.seed(1234) #  for reproducibility

est0 <- NetEstgam(N=dim, TT=size, K=K, S=1000, ho=ADV$ho, z=zz, p_dist=p_clog)

## Table 4
View(round(est0$results[-K,], 4))

# check significance with Bonferroni correction

adim = 13*7
crit <- qnorm(0.025/adim, lower.tail = F)
abs(est0$results[-K,15:21]) > crit


# Figure 5

par(mfrow=c(1,1))
plot(est0$results[-K,2]~seq(from=0, to=12, by=1), xlab="Adv. booking, k", ylab="Estimated network effect")
abline(lm(est0$results[-K,2]~seq(from=0, to=12, by=1)), col="red")

## overall trend  effect
summary(lm(est0$results[-K,2]~seq(from=0, to=12, by=1)))



### ACF plot (Figure A-2)

err_p <- list()
sdd <- list()
for(k in 1:K){
  err_p[[k]] <- ADV$ho[[k]] - est0$hat_lambda[[k]]
  sdd[[k]] <- apply(err_p[[k]], 1, sd)
  err_p[[k]] <- err_p[[k]]/sdd[[k]]
}

par(mfrow=c(3,4))
par(mar=c(2,2,1.5,1.5))
acf(err_p[[1]][4,])
acf(err_p[[1]][8,])
acf(err_p[[1]][15,])
acf(err_p[[1]][25,])
acf(err_p[[1]][50,])
acf(err_p[[1]][89,])
acf(err_p[[1]][17,])
acf(err_p[[1]][85,])
acf(err_p[[1]][26,])
acf(err_p[[1]][82,])
acf(err_p[[1]][81,])
acf(err_p[[1]][94,])



#################### Hotels classification (Section 6) ################################

############################ K=13, whole booking window

## build frequency matrices row-to-col (RtC), col-to-row (CtR) and Mixed

tab1 <- list()
tab2 <- list()
tab3 <- list()
for(k in 1:K){
  tab1[[k]] <- matrix(0, dim, dim)
  tab2[[k]] <- matrix(0, dim, dim)
  tab3[[k]] <- matrix(0, dim, dim)
  for(i in 1:dim){
    for(j in 1:dim){
      if(est0$hat_A[[k]][i,j]==1 && est0$hat_A[[k]][j,i]==0) tab1[[k]][i,j] <- 1
      if(est0$hat_A[[k]][i,j]==0 && est0$hat_A[[k]][j,i]==0) tab2[[k]][i,j] <- 1
      if(est0$hat_A[[k]][i,j]==1 && est0$hat_A[[k]][j,i]==1) tab3[[k]][i,j] <- 1
    }
  }
}

tab1_sum <- Reduce('+', tab1) # row-to-col
tab2_sum <- Reduce('+', tab2) # no connections
tab3_sum <- Reduce('+', tab3) # mixed

# compute frequencies

f1 <- tab1_sum/K; f2 <- t(tab1_sum)/K; f3 <- tab3_sum/K; f4 <- tab2_sum/K

# compute mean and standard deviation of frequencies

n1 <- sum(f1!=0); n2 <- sum(f2!=0); n3 <- sum(f3!=0); n4 <- sum(f4!=0)
af1 <- sum(f1)/n1; af2 <- sum(f2)/n2; af3 <- sum(f3)/n3; af4 <- sum(f4)/n4
sdf1 <- sdf2 <- sqrt(sum((f1[f1!=0]-af1)^2)/n1); sdf3 <- sqrt(sum((f3[f3!=0]-af3)^2)/n3); sdf4 <- sqrt(sum((f4[f4!=0]-af4)^2)/n4)

# classify each pair of hotels 

bf1 <- bf2 <- max(0.5, af1+sdf1); bf3 <- max(0.5, af3+sdf3); bf4 <- max(0.5, af4+sdf4)
mf1 <- rowSums(f1>=bf1); mf2 <- rowSums(f2>=bf2); mf3 <- rowSums(f3>=bf3); mf4 <- rowSums(f4>=bf4)


# general classification of hotels 

# Follower = 1
# Leader = 2
# Mixed = 3
# Independent =4

FOLL <- ifelse((mf1>mf2), 1, 0) # classified as follower
LEAD <- ifelse((mf1<mf2), 2, 0) # classified as leader
IND <- ifelse((mf1+mf2+mf3==0), 4, 0) # classified as independent
MIX<- ifelse((LEAD + FOLL + IND==0),3,0) # classified as mixed

cl <- FOLL + LEAD + IND + MIX

## descriptive statistics for different profiles

# frequency of classification
sum(cl==1); sum(cl==2); sum(cl==3); sum(cl==4)

# add stars variable
Stars1 <- Data[,c(1,7)]
Stars <- Stars1[!duplicated(Stars1),]
Desc <- cbind(ADV$z, rowSums(lattice))
Descr_p <- cbind(Desc, Stars[,2])
colnames(Descr_p)[4] <- "Num. room 200 mt"
colnames(Descr_p)[5] <- "Stars"

# Descriptive statistics for each classified profile (Table 5)

De <- cbind(Descr_p, cl)
colMeans(De[De[,6]==1,c(1:5)])
colMeans(De[De[,6]==2,c(1:5)])
colMeans(De[De[,6]==3,c(1:5)])
colMeans(De[De[,6]==4,c(1:5)])

HOp <- ADV$ho[[1]]
for(k in 2:15){
  HOp <- cbind(HOp, ADV$ho[[k]])
}
descr_hop <- cbind(exp(HOp), cl)
descr_hoc1p <- descr_hop[descr_hop[,5161]==1,-5161]
descr_hoc2p <- descr_hop[descr_hop[,5161]==2,-5161]
descr_hoc3p <- descr_hop[descr_hop[,5161]==3,-5161]
descr_hoc4p <- descr_hop[descr_hop[,5161]==4,-5161]
mean(vec(descr_hop))
mean(vec(descr_hoc1p))
mean(vec(descr_hoc2p))
mean(vec(descr_hoc3p))
mean(vec(descr_hoc4p))



############################## K=4 Last minute

## build frequency matrices row-to-col (RtC), col-to-row (CtR) and Mixed

tab1.3 <- list()
tab2.3 <- list()
tab3.3 <- list()
for(k in 1:4){
  tab1.3[[k]] <- matrix(0, dim, dim)
  tab2.3[[k]] <- matrix(0, dim, dim)
  tab3.3[[k]] <- matrix(0, dim, dim)
  for(i in 1:dim){
    for(j in 1:dim){
      if(est0$hat_A[[k]][i,j]==1 && est0$hat_A[[k]][j,i]==0) tab1.3[[k]][i,j] <- 1
      if(est0$hat_A[[k]][i,j]==0 && est0$hat_A[[k]][j,i]==0) tab2.3[[k]][i,j] <- 1
      if(est0$hat_A[[k]][i,j]==1 && est0$hat_A[[k]][j,i]==1) tab3.3[[k]][i,j] <- 1
    }
  }
}

tab1_sum3 <- Reduce('+', tab1.3) # row-to-col
tab2_sum3 <- Reduce('+', tab2.3) # no connections
tab3_sum3 <- Reduce('+', tab3.3) # mixed


# compute frequencies

f13 <- tab1_sum3/4; f23 <- t(tab1_sum3)/4; f33 <- tab3_sum3/4; f43 <- tab2_sum3/4

# compute mean and standard deviation of frequencies

n13 <- sum(f13!=0); n23 <- sum(f23!=0); n33 <- sum(f33!=0); n43 <- sum(f43!=0)
af13 <- sum(f13)/n13; af23 <- sum(f23)/n23; af33 <- sum(f33)/n33; af43 <- sum(f43)/n43
sdf13 <- sdf23 <- sqrt(sum((f13[f13!=0]-af13)^2)/n13); sdf33 <- sqrt(sum((f33[f33!=0]-af33)^2)/n33); sdf43 <- sqrt(sum((f43[f43!=0]-af43)^2)/n43)

# classify each pair of hotels 

bf13 <- bf23 <- max(0.5, af13+sdf13); bf33 <- max(0.5, af33+sdf33); bf43 <- max(0.5, af43+sdf43)
mf13 <- rowSums(f13>=bf13); mf23 <- rowSums(f23>=bf23); mf33 <- rowSums(f33>=bf33); mf43 <- rowSums(f43>=bf43)


# general classification of hotels 

# Follower = 1
# Leader = 2
# Mixed = 3
# Independent =4

FOLL3 <- ifelse((mf13>mf23), 1, 0) # classified as follower
LEAD3 <- ifelse((mf13<mf23), 2, 0) # classified as leader
IND3 <- ifelse((mf13+mf23+mf33==0), 4, 0) # classified as independent
MIX3 <- ifelse((LEAD3 + FOLL3 + IND3==0),3,0) # classified as mixed

cl3 <- FOLL3 + LEAD3 + IND3 + MIX3

## Descriptive statistics for different profiles

# frequency of classification
sum(cl3==1); sum(cl3==2); sum(cl3==3); sum(cl3==4)


# Descriptive statistics for each classified profile (Table 6)

De3 <- cbind(Descr_p, cl3)
colMeans(De3[De3[,6]==1,c(1:5)])
colMeans(De3[De3[,6]==2,c(1:5)])
colMeans(De3[De3[,6]==3,c(1:5)])
colMeans(De3[De3[,6]==4,c(1:5)])


descr_hop3 <- cbind(exp(HOp), cl3)
descr_hoc1p3 <- descr_hop3[descr_hop3[,5161]==1,-5161]
descr_hoc2p3 <- descr_hop3[descr_hop3[,5161]==2,-5161]
descr_hoc3p3 <- descr_hop3[descr_hop3[,5161]==3,-5161]
descr_hoc4p3 <- descr_hop3[descr_hop3[,5161]==4,-5161]
mean(vec(descr_hop3))
mean(vec(descr_hoc1p3))
mean(vec(descr_hoc2p3))
mean(vec(descr_hoc3p3))
mean(vec(descr_hoc4p3))


##### Visualize hotel maps (Figure 6-7)

ADV2 <- data_format(Data=Data)
coord <- apply(ADV2$coord, 1:2, as.numeric)
lat <- coord[,1]
lon <- coord[,2]

meta <- data.frame("name"=seq(1:95), "lon"=lon, "lat"=lat)

Adg <- list()
for(k in 1:K){
  Adg[[k]]  <- graph.adjacency(est0$hat_A[[k]], mode="directed", weighted=NULL)
}

## Advance booking 0

Agg <- as_data_frame(Adg[[1]], what="both")
Adf <- data.frame("from"=Agg$edges[,1], "to"=Agg$edges[,2])

graphA <- graph.data.frame(Adf, directed=TRUE, vertices=meta)

ggt <- get.data.frame(graphA, "both")
vertt <- ggt$vertices

coordinates(vertt) <- ~lon+lat
edgest <- ggt$edges

edgest <- lapply(1:nrow(edgest), function(i) {
  as(rbind(vertt[vertt$name == edgest[i, "from"], ], 
           vertt[vertt$name == edgest[i, "to"], ]), 
     "SpatialLines")
})

for (i in seq_along(edgest)) {
  edgest[[i]] <- spChFIDs(edgest[[i]], as.character(i))
}

edgest <- do.call(rbind, edgest)

leaflet(vertt) %>% addTiles() %>% 
  addCircles(data = vertt, color="red", opacity=3,
             fillOpacity = 30) %>% addPolylines(data = edgest, weight=1)

## Advance booking 7

Agg <- as_data_frame(Adg[[8]], what="both")
Adf <- data.frame("from"=Agg$edges[,1], "to"=Agg$edges[,2])

graphA <- graph.data.frame(Adf, directed=TRUE, vertices=meta)

ggt <- get.data.frame(graphA, "both")
vertt <- ggt$vertices

coordinates(vertt) <- ~lon+lat
edgest <- ggt$edges

edgest <- lapply(1:nrow(edgest), function(i) {
  as(rbind(vertt[vertt$name == edgest[i, "from"], ], 
           vertt[vertt$name == edgest[i, "to"], ]), 
     "SpatialLines")
})

for (i in seq_along(edgest)) {
  edgest[[i]] <- spChFIDs(edgest[[i]], as.character(i))
}

edgest <- do.call(rbind, edgest)

leaflet(vertt) %>% addTiles() %>% 
  addCircles(data = vertt, color="red", opacity=3,
             fillOpacity = 30) %>% addPolylines(data = edgest, weight=1)



