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
spectralrao<-function(input, distance_m="euclidean", p=NULL, window=9, mode="classic", shannon=FALSE, rescale=FALSE, na.tolerance=0.0, simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE) {
#
#Load required packages
#
    require(raster)
#
#Initial checks
#
    if( !(is(input,"matrix") | is(input,"SpatialGridDataFrame") | is(input,"RasterLayer") | is(input,"list")) ) {
        stop("\nNot a valid input object.")
    }
#Change input input/ces names
    if( is(input,"SpatialGridDataFrame") ) {
        input <- raster(input)
    }
    if( is(input,"matrix") | is(input,"RasterLayer")) {
        rasterm<-input
    } else if( is(input,"list") ) {
        rasterm<-input[[1]]
    }
#Deal with matrices and RasterLayer in a different way
    if( is(input[[1]],"RasterLayer") ) {
        if( mode=="classic" ){
#If the data is float number transform it in integer
            isfloat<-FALSE
            if( !is.integer(getValues(rasterm)) ){
                isfloat<-TRUE
                mfactor<-100^simplify
                rasterm<-apply(as.matrix(rasterm), 1:2, function(y) round(y,simplify))
                rasterm<-apply(as.matrix(rasterm*mfactor), 1:2, function(x) as.integer(x))
                message("Converting in an integer matrix...")
            }else{
                rasterm<-as.matrix(rasterm)
            }
            message("RasterLayer ok: \nRao and Shannon output matrices will be returned")
        }else if(mode=="multidimension" & !shannon){
            message(("####\nRasterLayer ok: \nA raster object with multimension RaoQ will be returned\n####"))
        }else if(mode=="multidimension" & shannon){
            stop(("Matrix check failed: \nMultidimension and Shannon not compatible, set shannon=FALSE"))
        }
    }else if( is(input,"matrix") | is(input,"list") ) {
         if( mode=="classic" ){
#If the data is float number transform it in integer
            isfloat<-FALSE
            if( !is.integer(rasterm) ){
                isfloat<-TRUE
                mfactor<-100^simplify
                rasterm<-apply(as.matrix(rasterm), 1:2, function(y) round(y,simplify))
                rasterm<-apply(as.matrix(rasterm*mfactor), 1:2, function(x) as.integer(x))
                message("Converting in an integer matrix...")
            }else{
                rasterm<-as.matrix(rasterm)
            }
        }
        if(mode=="classic" & shannon){
            message("Matrix check ok: \nRao and Shannon output matrices will be returned")
        }else if(mode=="classic" & !shannon){
            message("Matrix check ok: \nRao output matrix will be returned")
        }else if(mode=="multidimension" & !shannon){
            message(("Matrix check ok: \nA matrix with multimension RaoQ will be returned"))
        }else if(mode=="multidimension" & shannon){
            stop("Matrix check failed: \nMultidimension and Shannon not compatible, set shannon=FALSE")
        }else{stop("Matrix check failed: \nNot a valid input | method | distance, please check all these options")
    }
}
if(nc.cores>1) {
    message("
    ###################### Starting parallel calculation ##########################")
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
if(shannon){
    shannond<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
}

if(mode=="classic") {
#
#If classic RaoQ - parallelized
#
    if(nc.cores>1) {
#
##required
#
        require(foreach)
        require(doSNOW)
        require(Rmpi)
#
##Reshape values
#
        values<-as.integer(as.factor(rasterm))
        if(debugging){cat("check_1")}
        rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
#
##Add fake columns and rows for moving window
#
        hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
        ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
        trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
        rm(rasterm_1,hor,ver)
#       
##Derive distance matrix
#
        classes<-levels(as.factor(rasterm))
        d1<-dist(classes,method=distance_m)
        rm(classes)
#       
##Create cluster object with given number of slaves
#
        plr<<-TRUE
        if(cluster.type=="SOCK" || cluster.type=="FORK") {
            cls <- parallel::makeCluster(nc.cores,type=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
        } else if(cluster.type=="MPI") {
            cls <- makeMPIcluster(nc.cores,outfile="",useXDR=FALSE,methods=FALSE,output="")
        }
        registerDoSNOW(cls)
        on.exit(stopCluster(cls)) # Close the clusters
        gc()
    #
    ##Start the parallelized loop over iter
    #
        raop <- foreach(cl=(1+w):(dim(rasterm)[2]+w)) %dopar% {
           if(debugging) {
            cat(cl)
        }
        vout<-c()
        for(rw in (1+w):(dim(rasterm)[1]+w)) {
            if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
                vv<-raoqe[rw-w,cl-w]
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
                if(length(tw_values) <= 2) {
                    vv<-raoqe[rw-w,cl-w]
                } else {
                    p <- tw_values/sum(tw_values)
                    p1 <- diag(0,length(tw_values))
                    p1[upper.tri(p1)] <- c(combn(p,m=2,FUN=prod))
                    p1[lower.tri(p1)] <- c(combn(p,m=2,FUN=prod))
                    d2 <- unname(as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
                    vv <- sum(p1*d2)
                }
            }
            vout<-append(vout,vv)
        }
        return(as.integer(vout*100^simplify)) #addedd bit of code
        gc()
    }
        raoqe<-do.call(cbind,raop) # End classic RaoQ - parallelized
#
##If classic RaoQ - sequential
#
    } else if(nc.cores==1) {
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
                if(length(tw_values) <= 2) {
                    raoqe[rw-w,cl-w]<-NA
                } else {
                    p <- tw_values/sum(tw_values)
                    p1 <- diag(0,length(tw_values))
                    p1[upper.tri(p1)] <- c(combn(p,m=2,FUN=prod))
                    p1[lower.tri(p1)] <- c(combn(p,m=2,FUN=prod))
                    d2 <- unname(as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
                     if(isfloat) {
                        raoqe[rw-w,cl-w]<-sum(p1*d2)/mfactor
                    } else {
                        raoqe[rw-w,cl-w]<-sum(p1*d2)
                    }
                }
            }
        } 
    } # End of for loop 
} # End classic RaoQ - sequential
#----------------------------------------------------#
} else if(mode=="multidimension"){
#
#If multimensional RaoQ
#
#Check if there are NAs in the matrices
    if ( is(rasterm,"RasterLayer") ){
        if(any(sapply(lapply(input, function(x) {as.matrix(x)}), is.na)==TRUE))
            message("\n Warning: One or more RasterLayers contain NA which will be threated as 0")
    } else if ( is(rasterm,"matrix") ){
        if(any(sapply(input, is.na)==TRUE) ) {
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
            return(Reduce(`+`,tmp))
        }
        #canberra
        multicanberra <- function(x) {
            tmp <- lapply(x, function(y) {
                abs(y[[1]] - y[[2]]) / (abs(y[[1]]) + abs(y[[2]]))
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
    vls<-lapply(input, function(x) {as.matrix(x)})
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
                tw[tw>1]<-1
                tw_values<-as.vector(tw)
                p<-tw_values/length(tw_values)
                p_log<-log(p)
                shannond[rw-w,cl-w]<-(-(sum(p*p_log)))
                #cat(paste(cl,rw,(-(sum(p*p_log)))),"\n")
            }
        }
    } # End ShannonD
}
#----------------------------------------------------#
#
#Return the output
#
# if( is(rasterm,"RasterLayer") ) {
#     if( shannon ) {
# #Rasterize the matrices if matrix==raster
#         rastertemp <- stack(raster(raoqe, template=input),raster(shannond, template=raster))
#     } else if(shannon==FALSE) {
#         rastertemp <- raster(raoqe, template=rasterm)
#     }
# }
#
#Return multiple outputs
#
if(debugging){
        message("check_2")
    }

if( shannon ) {
    outl<-list(raoqe,shannond)
    names(outl)<-c("Rao","Shannon")
    return(outl)
} else if( !shannon & mode=="classic") {
    if(isfloat & nc.cores>1) {
    raoqe<-do.call(cbind,raop)/mfactor/mfactor
    if(debugging){
        message("check_2.5")
    }
}
    outl<-list(raoqe)
    names(outl)<-c("Rao")
    return(outl)
} else if( !shannon & mode=="multidimension") {
    outl<-list(raoqe)
    names(outl)<-c("Multidimension_Rao")
    return(outl)
}

}