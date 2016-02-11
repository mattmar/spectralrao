#Function to calculate Rao's quadratic entropy on a numeric matrix o a RasterLayer object using a moving window. The function also calculate Shannon diversity index. Rao's Q Min = 0, if all pixel classes have distance 0. If the chosen distance ranges between 0 and 1, Rao's Max = 1-1/S (Simpson Diversity, where S is pixel classes).

#Function
spectralrao<-function(raster,distance_m="euclidean",window=3) {

#Change raster name
    rasterm<-raster

#Packages
    require(raster)

#Threating matrix and RasterLayer object in a different way
    if(is(rasterm,"RasterLayer")) {
        print("RasterLayer ok: Rao and Shannon output matrices will be returned")
        rasterm<-round(as.matrix(rasterm),2)
    }
    if(is(rasterm,"matrix")) {
        print("Matrix check ok: Rao and Shannon output matrices will be returned")
    }else{print("Not a valid input, please provide a matrix or a RasterLayer object")}

#Calculate true moving window
    oddn<-seq(1,window,2)
    oddn_pos<-which(oddn == 3)
    window=window-oddn_pos

#Matrix preparation
    raoqe<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
    shannond<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])

#Reshape values
    values<-as.numeric(as.factor(rasterm))
    rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])

#Add fake columns and rows for moving window
    hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=window)
    ver<-matrix(NA,ncol=window,nrow=dim(rasterm)[1]+window*2)
    trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)

#Derive distance matrix
    classes<-levels(as.factor(rasterm))
    d1<-dist(classes,method=distance_m)

#RaoQ
    for (cl in (1+window):(dim(rasterm)[2]+window)) {
        for(rw in (1+window):(dim(rasterm)[1]+window)) {
            twindow<-summary(as.factor(trasterm[c(rw-window):c(rw+window),c(cl-window):c(cl+window)]),maxsum=10000)
            if("NA's" %in% names(twindow)){
                twindow<-twindow[-length(twindow)]
            }
            twindow_labels<-names(twindow)
            twindow_values<-as.vector(twindow)
            if(length(twindow_values) == 1) {
                raoqe[rw-window,cl-window]<-0
            }else{p<-twindow_values/sum(twindow_values)
            p1<-combn(p,m=2,FUN=prod)
            d2<-as.matrix(d1)
            d2[upper.tri(d1,diag=TRUE)]<-NA
            d3<-d2[as.numeric(twindow_labels),as.numeric(twindow_labels)]
            raoqe[rw-window,cl-window]<-sum(p1*d3[!(is.na(d3))])
        }
    }
} # End Rao

#ShannonD
for (cl in (1+window):(dim(rasterm)[2]+window)) {
    for(rw in (1+window):(dim(rasterm)[1]+window)) {
        twindow<-summary(as.factor(trasterm[c(rw-window):c(rw+window),c(cl-window):c(cl+window)]))
        twindow_values<-as.vector(twindow)
        p<-twindow_values/sum(twindow_values)
        p_log<-log(p)
        shannond[rw-window,cl-window]<-(-(sum(p*p_log)))
    }
} # End Shannon

#Return the output
if(is(raster,"RasterLayer")){

# Rasterize the matrices
    rastertemp <- stack(raster(raoqe, template=raster),raster(shannond, template=raster))
}

if(is(raster,"RasterLayer")) {
    return(rastertemp)
}else{return(list(raoqe,shannond))}
}
