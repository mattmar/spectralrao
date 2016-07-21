######### SPECTRAL RAO #############################
## Code to calculate Rao's quadratic entropy on a
## numeric matrix, RasterLayer object (or lists)
## using a moving window. The function also calculates
## Shannon diversity index.
## Rao's Q Min = 0, if all pixel classes have
## distance 0. If the chosen distance ranges between
## 0 and 1, Rao's Max = 1-1/S (Simpson Diversity,
## where S is pixel classes).
## Latest update: 21th July
#####################################################
#Function
spectralrao<-function(matrix, distance_m="euclidean", p=NULL, window=9, mode="classic", shannon=FALSE, rescale=FALSE, na.tolerance=0.0, debugging=FALSE) {
#
#Load required packages
#
    require(raster)
#
#Initial checks
#
    if( !(is(matrix,"matrix") | is(matrix,"SpatialGridDataFrame") | is(matrix,"RasterLayer") | is(matrix,"list")) ) {
        stop("\nNot a valid input object.")
    }
#Change input matrix/ces names
    if( is(matrix,"SpatialGridDataFrame") ) {
        matrix <- raster(matrix)
    }
    if( is(matrix,"matrix") | is(matrix,"RasterLayer")) {
        rasterm<-matrix
    } else if( is(matrix,"list") ) {
        rasterm<-matrix[[1]]
    }
#Deal with matrix and RasterLayer in a different way
    if( is(matrix[[1]],"RasterLayer") ) {
        if( mode=="classic" ){
            rasterm<-round(as.matrix(rasterm),3)
            message("RasterLayer ok: \nRao and Shannon output matrices will be returned")
        }else if(mode=="multidimension" & shannon==FALSE){
            message(("RasterLayer ok: \nA raster object with multimension RaoQ will be returned"))
        }else if(mode=="multidimension" & shannon==TRUE){
            stop(("Matrix check failed: \nMultidimension and Shannon not compatible, set shannon=FALSE"))
        }
    }else if( is(matrix,"matrix") | is(matrix,"list") ) {
        if(mode=="classic"){
            message("Matrix check ok: \nRao and Shannon output matrices will be returned")
        }else if(mode=="multidimension" & shannon==FALSE){
            message(("Matrix check ok: \nA matrix with multimension RaoQ will be returned"))
        }else if(mode=="multidimension" & shannon==TRUE){
            stop("Matrix check failed: \nMultidimension and Shannon not compatible, set shannon=FALSE")
        }else{stop("Matrix check failed: \nNot a valid input, please provide a matrix, list or RasterLayer object")
    }
}
#
##Derive operational moving window
#
if( window%%2==1 ){
    w <- (window-1)/2
} else {stop("Moving window size must be odd")}
#
##Output matrices preparation
#
raoqe<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
shannond<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
#
#If classic RaoQ
#
if(mode=="classic"){
#Reshape values
    values<-as.numeric(as.factor(rasterm))
    rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
#Add fake columns and rows for moving window
    hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
    ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
    trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
#Derive distance matrix
    classes<-levels(as.factor(rasterm))
    d1<-dist(classes,method=distance_m)
#Loop over each pixel
    for (cl in (1+w):(dim(rasterm)[2]+w)) {
        for(rw in (1+w):(dim(rasterm)[1]+w)) {
            if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
                raoqe[rw-w,cl-w]<-NA
            } else {
                tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
                if( "NA's"%in%names(tw) ) {
                    tw<-tw[-length(tw)]
                }
                if(debugging) {
                    message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window)
                }
                tw_labels<-names(tw)
                tw_values<-as.vector(tw)
                if(length(tw_values) == 1) {
                    raoqe[rw-w,cl-w]<-NA
                } else {
                    p <- tw_values/sum(tw_values)
                    p1 <- diag(0,length(tw_values))
                    p1[upper.tri(p1)] <- c(combn(p,m=2,FUN=prod))
                    p1[lower.tri(p1)] <- c(combn(p,m=2,FUN=prod))
                    d2 <- unname(as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
                    raoqe[rw-w,cl-w]<-sum(p1*d2)
                }
            }
        }
    } # End classic RaoQ
#----------------------------------------------------#
} else if(mode=="multidimension"){
#
#If multimensional RaoQ
#
#Check if there are NAs in the matrices
    if ( is(rasterm,"RasterLayer") ){
        if(any(sapply(lapply(matrix, function(x) {as.matrix(x)}), is.na)==TRUE))
            message("\n Warning: One or more RasterLayers contain NA which will be threated as 0")
    } else if ( is(rasterm,"matrix") ){
        if(any(sapply(matrix, is.na)==TRUE) ) {
            message("\n Warning: One or more matrices contain NA which will be threated as 0")
        }
    }
#
##Check whether the distance is valid or not
#
    if( distance_m=="euclidean" | distance_m=="manhattan" | distance_m=="canberra" ) {
        #Define the distance functions
        #euclidean
        multieuclidean <- function(x) {
            tmp <- lapply(x, function(y) {
                (y[[1]]-y[[2]])^2
            })
            return(sqrt(Reduce(`+`,tmp)))
        }
        #manhattan
        multimanhattan <- function(x) {
            tmp <- lapply(x, function(y) {
                abs(y[[1]]-y[[2]])
            })
            return(sqrt(Reduce(`+`,tmp)))
        }
        #canberra
        multicanberra <- function(x) {
            tmp <- lapply(x, function(y) {
                abs(y[[1]] - y[[2]]) / (y[[1]] + y[[2]])
            })
            return(Reduce(`+`,tmp))
        }
#
##Decide what function to use
#
        if( distance_m=="euclidean") {
            distancef <- get("multieuclidean")
        } else if(distance_m=="manhattan") {
            distancef <- get("multimanhattan")
        } else if(distance_m=="canberra") {
            distancef <- get("multicanberra")
        }
    } else {
        stop("Distance function not defined for multidimensional Rao's Q; please chose among euclidean, manhattan or canberra")
    }
#
##Reshape values
#
    vls<-lapply(matrix, function(x) {as.matrix(x)})
#
##Rescale and add fake columns and rows for moving w
#
    hor<-matrix(NA,ncol=dim(vls[[1]])[2],nrow=w)
    ver<-matrix(NA,ncol=w,nrow=dim(vls[[1]])[1]+w*2)
    if(rescale) {
        trastersm<-lapply(vls, function(x) {
            t1 <- raster::scale(raster(cbind(ver,rbind(hor,x,hor),ver)))
            t2 <- as.matrix(t1)
            return(t2)
        })
    } else {
        trastersm<-lapply(vls, function(x) {
            cbind(ver,rbind(hor,x,hor),ver)
        })
    }
#
##Loop over all the pixels in the matrices
#
    if( (ncol(vls[[1]])*nrow(vls[[1]]))> 10000) {
        message("\n Warning: ",ncol(vls[[1]])*nrow(vls[[1]])*length(vls), " cells to be processed, may take some time... \n")
    }
    for (cl in (1+w):(dim(vls[[1]])[2]+w)) {
        for(rw in (1+w):(dim(vls[[1]])[1]+w)) {
            if( length(!which(!trastersm[[1]][c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
                raoqe[rw-w,cl-w] <- NA
            } else {

                tw<-lapply(trastersm, function(x) { x[(rw-w):(rw+w),(cl-w):(cl+w)]
            })
#
##Vectorize the matrices in the list and calculate
#the among matrix pairwase distances
                lv <- lapply(tw, function(x) {as.vector(t(x))})
                vcomb <- combn(length(lv[[1]]),2)
                vout <- c()
                for(p in 1:ncol(vcomb) ) {
                    lpair <- lapply(lv, function(chi) {
                        c(chi[vcomb[1,p]],chi[vcomb[2,p]])
                    })
                    vout[p] <- distancef(lpair)
                }
                raoqe[rw-w,cl-w] <- sum(rep(vout,2) * (1/(window)^4),na.rm=TRUE)
            }
        }
    }
} # end multimensional RaoQ
#----------------------------------------------------#
#
##ShannonD
#
if(shannon){
    #Reshape values
    values<-as.numeric(as.factor(rasterm))
    rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
#
##Add fake columns and rows for moving window
#
    hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
    ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
    trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
#
##Loop over all the pixels
#
    for (cl in (1+w):(dim(rasterm)[2]+w)) {
        for(rw in (1+w):(dim(rasterm)[1]+w)) {
            if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
                shannond[rw-w,cl-w]<-NA
            } else {
                tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]))
                if( "NA's"%in%names(tw) ) {
                    tw<-tw[-length(tw)]
                }
                tw_values<-as.vector(tw)
                p<-tw_values/sum(tw_values)
                p_log<-log(p)
                shannond[rw-w,cl-w]<-(-(sum(p*p_log)))
            }
        }
    } # End ShannonD
}
#----------------------------------------------------#
#
#Return the output
#
if( is(rasterm,"RasterLayer") ) {
    if( shannon) {
#Rasterize the matrices if matrix==raster
        rastertemp <- stack(raster(raoqe, template=matrix),raster(shannond, template=raster))
    } else if(shannon==FALSE){
        rastertemp <- raster(raoqe, template=rasterm)
    }
}
#
#Return multiple outputs
#
if( is(rasterm,"RasterLayer") ) {
    return(rastertemp)
} else if( !is(rasterm,"RasterLayer") & shannon ) {
    return(list(raoqe,shannond))
} else if( !is(rasterm,"RasterLayer") & shannon==FALSE ) {
    return(list(raoqe))
}
}
