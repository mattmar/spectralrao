######### SPECTRAL RAO #############################
## Code to calculate Rao's quadratic entropy on a
## numeric matrix, RasterLayer object (or lists)
## using a moving window. The function also calculates
## Shannon diversity index.
## Rao's Q Min = 0, if all pixel classes have
## distance 0. If the chosen distance ranges between
## 0 and 1, Rao's Max = 1-1/S (Simpson Diversity,
## where S is pixel classes).
## Last update: 29th February
#####################################################

#Function
spectralrao<-function(matrix,distance_m="euclidean",window=10,mode="classic",shannon=TRUE,debugging=F) {

#Load required packages
    require(raster)

#Initial checks
    if( !(is(matrix,"matrix") | is(matrix,"SpatialGridDataFrame") | is(matrix,"RasterLayer")) ) {stop("\nNot a valid input object.") }
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
            rasterm<-round(as.matrix(rasterm),4)
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

#Derive operational moving window
if( window%%2==1 ){
    w <- (window-1)/2
}else{stop("Moving window size must be odd")}

#Output matrices preparation
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
            if( !length(which(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2 )  {
                raoqe[rw-w,cl-w]<-NA
            }else{
                tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
                if("NA's" %in% names(tw)) {
                    tw<-tw[-length(tw)]
                }
                if(debugging){
                    message("Working on coords ",rw-w ,",",cl-w,". mw length: ",length(tw),", window size=",window)
                }
                tw_labels<-names(tw)
                tw_values<-as.vector(tw)
                if(length(tw_values) == 1) {
                    raoqe[rw-w,cl-w]<-NA
                }else{p<-tw_values/sum(tw_values)
                p1<-combn(p,m=2,FUN=prod)
                d2<-as.matrix(d1)
                d2[upper.tri(d1,diag=TRUE)]<-NA
                d3<-d2[as.numeric(tw_labels),as.numeric(tw_labels)]
                raoqe[rw-w,cl-w]<-sum(p1*d3[!(is.na(d3))])
            }
        }
    }
} # End classic RaoQ
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

#Reshape values
    vls<-lapply(matrix, function(x) {as.matrix(x)})
#Add fake columns and rows for moving w
    hor<-matrix(NA,ncol=dim(vls[[1]])[2],nrow=w)
    ver<-matrix(NA,ncol=w,nrow=dim(vls[[1]])[1]+w*2)
    trastersm<-lapply(vls, function(x) {cbind(ver,rbind(hor,x,hor),ver)})

# Loop over all the pixels in the matrices
    if( (ncol(vls[[1]])*nrow(vls[[1]]))> 10000) {
        message("\n Warning: ",ncol(vls[[1]])*nrow(vls[[1]])*length(vls), " cells to be processed, may take some time... \n")
    }

    for (cl in (1+w):(dim(vls[[1]])[2]+w)) {
        for(rw in (1+w):(dim(vls[[1]])[1]+w)) {
            tw<-lapply(trastersm, function(x) { x[(rw-w):(rw+w),(cl-w):(cl+w)]
        })
            distances<-lapply(tw, function(x) {
                out<-matrix(NA,nrow=length(x),ncol=length(x))
                for (i in 1:length(x)) {
                    out[,i] <- as.vector((x[i] - x)^2)
                    out[out %in% NA]<- 0
                }
                return(out)
            } )
            raoqe[rw-w,cl-w] <- sum(sqrt(Reduce('+',distances)) *(1/(window)^4))
        } # end multimensional RaoQ
    }
}
#
#ShannonD
#
if(shannon==TRUE){
    #Reshape values
    values<-as.numeric(as.factor(rasterm))
    rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])

#Add fake columns and rows for moving window
    hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
    ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
    trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)

#Loop over all the pixels
    for (cl in (1+w):(dim(rasterm)[2]+w)) {
        for(rw in (1+w):(dim(rasterm)[1]+w)) {
            if( !length(which(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2 )  { shannond[rw-w,cl-w]<-NA
        }else{
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
#
#Return the output
#
if( is(rasterm,"RasterLayer") ) {
    if( shannon==TRUE) {
#Rasterize the matrices if matrix==raster
        rastertemp <- stack(raster(raoqe, template=matrix),raster(shannond, template=raster))
    } else if(shannon==FALSE){
        rastertemp <- raster(raoqe, template=rasterm)
    }
}
#
#Return different outputs
#
if( is(rasterm,"RasterLayer") ) {
    return(rastertemp)
} else if( !is(rasterm,"RasterLayer") & shannon==TRUE ) {
    return(list(raoqe,shannond))
} else if( !is(rasterm,"RasterLayer") & shannon==FALSE ) {
    return(list(raoqe))
}
}
