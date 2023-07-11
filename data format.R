
###################################################
## Function which formats the data as follows:
## covariates in matrix z
## hotels per rows, time per column in matrix ho;
## one ho matrix for each advance booking 
##################################################

data_format <- function(Data){

Dat <- list()
ho <- list()

for (i in 1:15){
  Dat[[i]]<- Data[Data$adv.book == i-1,]
  ho[[i]] <- subset(Dat[[i]], select=-c(adv.book,Price,Stars))
  ho[[i]] <- tidyr::spread(ho[[i]], data_in, BAR, fill = 0)
  
  if(i==1){
    index <- as.matrix(ho[[i]][,1])                        # hotels id
    coord <- as.matrix(ho[[i]][,2:3])                      # coordinates
    dist <- as.matrix(dist(coord[,1:2], diag=T, upper=T))  # euclidean distance
    geodist <- geodist(coord[,1:2], measure = "haversine") # geodesic distance
    z <- as.matrix(ho[[i]][,4:6])                          # exogenous covariates
  }
  
  ho[[i]] <- as.matrix(ho[[i]][,-c(1:6)])
  
}
return(list(ho=ho, dist=dist, geodist=geodist, z=z, index=index, coord=coord))
}
