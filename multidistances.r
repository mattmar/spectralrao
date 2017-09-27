multieuclidean <- function(x) {
    tmp <- lapply(x, function(y) {
        (y[[1]]-y[[2]])^2
    })
    return(sqrt(Reduce(`+`,tmp)))
}

multimanhattan <- function(x) {
    tmp <- lapply(x, function(y) {
        abs(y[[1]]-y[[2]])
    })
    return(sqrt(Reduce(`+`,tmp)))
}

multicanberra <- function(x) {
    tmp <- lapply(x, function(y) {
        abs(y[[1]] - y[[2]]) / (abs(y[[1]]) + abs(y[[2]]))
    })
    return(Reduce(`+`,tmp))
}

multiminkowski <- function(x) {
  tmp <- lapply(x, function(y) {
    abs((y[[1]]-y[[2]])^lambda)
})
  return(Reduce(`+`,tmp)^(1/lambda))
}

multimahalanobis <- function(x){
  tmp<-matrix(NA,nrow=window^2,ncol=length(x))
  for(dimension in 1:length(x)){
    tmp[,dimension]<-matrix(x[[dimension]], ncol=1) 
    dimension=dimension+1
}
tmp <- tmp[!is.na(tmp[,1]),] 
if(rcond(cov(tmp)) <= 0.001){  
    tmp <- NA
    return(tmp)
} else {
    inverse<-solve(cov(tmp)) 
    tmp<-scale(tmp,center=T,scale=F)
    tmp<-as.numeric(t(tmp[1,])%*%inverse%*%tmp[1,])
    return(sqrt(tmp))
}
}