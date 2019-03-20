IndexComputation <- function(input, type, window=3, mode="single", BergerParker=FALSE, alpha=1, base=exp(1), rescale=FALSE, na.tolerance=0.0, simplify=3, nc.cores=1, cluster.type="MPI", debugging=FALSE, integer=FALSE, ...){
  #
  ## Load required packages
  #
  require(raster)
  require(svMisc)
  #
  ## Define function to check if a number is an integer
  #
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  Shannon <- FALSE
  Renyi <- FALSE
  Hill <- FALSE
  #
  ## Initial checks
  #
  if( !(is(input,"matrix") | is(input,"SpatialGridDataFrame") | is(input,"RasterLayer") | is(input,"list")) ) {
    stop("\nNot a valid input object.")
  }
  if( is(input,"SpatialGridDataFrame") ) {
    input <- raster(input) # Change input matrix/ces names
  }
  if( is(input,"matrix") | is(input,"RasterLayer")) {
    rasterm<-input
  } else if( is(input,"list") ) {
    rasterm<-input[[1]]
  }
  if(na.tolerance>1){
    stop("na.tolerance must be in the [0-1] interval. Exiting...")
  }
  if (base<.Machine$double.eps){
    stop("base value must be in the (0,+\u221E) interval. Exiting...")
  }
  if (!is.numeric(alpha)){
    stop("alpha must be a number or a vector. Exiting...")
  }
  if (mode=="single"){
    if (length(alpha)!=1){
        stop("In mode \"single\" alpha must be a single number. Exiting...")
    }
    if (alpha<0){
      stop("The alpha value must be a non-negative number. Exiting...")
    }
    if(abs(alpha-1)<.Machine$double.eps && !BergerParker){
      Shannon <- TRUE
    }
    if(alpha >=.Machine$integer.max){
      BergerParker <- TRUE
    }
    if(integer){
      alpha <- as.integer(alpha)
    }
  }
  else if(mode == "iterative"){
      if (length(alpha) != 2){
        stop("In mode \"iterative\" alpha must be a numeric vector containing the star and the stop value. Exiting...")
      }
    start <- as.integer(alpha[1])
     if (start<0){
       stop("The starting value must be a non-negative number. Exiting...")
     }
    stop <- as.integer(alpha[2])
    if (stop <= start){
      stop("Integer part of the starting value, alpha[1], must be strictly greater that the integer part of the stopping value, alpha[2]. Exiting...")
    }
    if (start <= 1 && 1 <= stop ){
        Shannon <- TRUE
      }
  }
  else if(mode == "sequential"){
    if ( length(alpha) < 2){
      stop("In mode \"sequential\" alpha must be a numeric vector containing at least two values. Exiting...")
    }
    if ( length(which(alpha < 0)) != 0 ){
      stop("The alpha values must be non-negative numbers. Exiting...")
    }
    if ( integer ){
      val <- unique( as.integer(alpha) )
    }
    else {
      val <- unique(alpha)
    }
  }
  else{
    stop("The choosen mode is not defined. Exiting...")
  }
  
  if(type=="Hill"){
    Hill <- TRUE
  }
  else if(type=="Rényi"){
    Renyi <- TRUE
  }
  else{
    stop("The choosen type is not defined. Exiting...")
  }
  # Deal with matrices and RasterLayer in a different way
  # If data are raster layers
  if( is(input[[1]],"RasterLayer") ) {
      isfloat <- FALSE # If data are float numbers, transform them in integer, this may allow for a shorter computation time on big datasets.
      if( !is.wholenumber(rasterm@data@min) | !is.wholenumber(rasterm@data@max) | is.infinite(rasterm@data@min) | !is.wholenumber(median(getValues(rasterm))) ){
        message("Converting input data in an integer matrix...")
        isfloat <- TRUE
        mfactor <- 100 ^ simplify
        rasterm <- getValues(rasterm) * mfactor
        rasterm <- as.integer(rasterm)
        rasterm <- matrix(rasterm, nrow(input), ncol(input), byrow = TRUE)
        gc()
      }
      else{
        rasterm <- matrix(getValues(rasterm), ncol = ncol(input), nrow = nrow(input), byrow=TRUE)
      }
    #Print user messages
    if (mode == "single"){
      if( BergerParker ){
        message("Matrix check OK: \nBerger-Parker output matrix will be returned")
      }else if( Shannon ){
        message("Matrix check OK: \nShannon output matrix will be returned")
      }else {
        message("Matrix check OK: \nRényi with parameter value=", alpha," output matrix will be returned")
      }
    }
    else if (mode == "iterative"){
        if( BergerParker && !Shannon){
          message("Matrix check OK: \nRényi output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix" )
        }else if( BergerParker && Shannon  ){
          message("Matrix check OK: \nRényi output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix (Shannon output matrix is included!)")
        }else if( !BergerParker && Shannon  ){
          message("Matrix check OK: \nRényi output matrix will be returned for parameters integer values in [",start,",",stop,"] (Shannon output matrix is included!)")
        }else {
          message("Matrix check OK: \nRényi output matrix will be returned for parameters integer values in [",start,",",stop,"]")
        }
      }
    # If data are a matrix or a list
  }else if( is(input,"matrix") | is(input,"list") ) {
      isfloat<-FALSE # If data are float numbers, transform them in integer
      if( !is.integer(rasterm) ){
        message("Converting input data in an integer matrix...")
        isfloat <- TRUE
        mfactor <- 100^simplify
        rasterm <- as.integer(rasterm*mfactor)
        rasterm <- matrix(rasterm,nrow(input),ncol(input),byrow=TRUE)
        gc()
      }else{
        rasterm<-as.matrix(rasterm)
      }
      #Print user messages
      if (mode == "single"){
        if( BergerParker ){
          message("Matrix check OK: \nBerger-Parker output matrix will be returned")
        }else if( Shannon ){
          message("Matrix check OK: \nShannon output matrix will be returned")
        }else {
          message("Matrix check OK: \nRényi with parameter value=", alpha," output matrix will be returned")
        }
      }
      else if (mode == "iterative"){
        if( BergerParker && !Shannon){
          message("Matrix check OK: \nRényi output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix" )
        }else if( BergerParker && Shannon  ){
          message("Matrix check OK: \nRényi output matrix will be returned for parameters integer values in [",start,",",stop,"] with Berger-Parker output matrix (Shannon output matrix is included!)")
        }else if( !BergerParker && Shannon  ){
          message("Matrix check OK: \nRényi output matrix will be returned for parameters integer values in [",start,",",stop,"] (Shannon output matrix is included!)")
        }else {
          message("Matrix check OK: \nRényi output matrix will be returned for parameters integer values in [",start,",",stop,"]")
        }
      }
  }
  #
  ## Derive operational moving window
  #
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of moving window must be an odd number. Exiting...")
  }
  
  if (nc.cores == 1){
    if(mode == "single") {
      if( Shannon ) {
        outS <- ShannonS(rasterm, w, base, debugging, type)
        message(("\nCalculation of Shannon's index is complete!\n"))
      } # End ShannonD
      else if( BergerParker ) {
        outS <- BergerParkerS(rasterm, w, base, debugging, type)
        message(("\nCalculation of Berger-Parker's index is complete!\n"))
      } # End BergerParker
      else{
        outS<- IndexS(rasterm, w, alpha, base, debugging, type)
        message(("\nCalculation of Rényi's index complete.\n"))
      }
    return (outS)
      }
    else if(mode == "iterative"){
      out <- list()
      for (ALPHA in start:stop){
        if(ALPHA == 1) {
          s <- "ShannonAlpha 1"
          out[[s]] <- ShannonS(rasterm, w, base, debugging, type)
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("RényiAlpha",as.character(ALPHA))
          out[[s]] <- IndexS(rasterm,w, ALPHA, base, debugging, type)
        }
      }
      message(("\nCalculation of Rényi's index complete.\n"))
        if (BergerParker){
          s<-"Berger-Parker"
          out[[s]] <- BergerParkerS(rasterm, w, base, debugging, type)
          message(("\nCalculation of Berger-Parker's index is also complete!\n"))
        }
      
      return(out)
    }
    else if(mode == "sequential"){
      out <- list()
      for (ALPHA in val){
        if(abs(ALPHA-1)<.Machine$double.eps) {
          s <- "ShannonAlpha 1"
          out[[s]] <- ShannonS(rasterm, w, base, debugging, type)
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("RényiAlpha",as.character(ALPHA))
          out[[s]] <- IndexS(rasterm,w, ALPHA, base, debugging, type)
        }
      }
      message(("\nCalculation of Rényi's index complete.\n"))
      if (BergerParker){
        s<-"Berger-Parker"
        out[[s]] <- BergerParkerS(rasterm, w, base, debugging, type)
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
      
      return(out)
    }
    
  }
  if (nc.cores>1){
    message("##################### Starting parallel calculation #######################")
    #
    ## Required packages for parallel calculation
    #
    require(foreach)
    require(doSNOW)
    require(parallel)
    if( cluster.type=="MPI" ){
      require(Rmpi)
    }
    #       
    ## Export variables in the global environment
    #
    if(isfloat) {
      sapply(c("mfactor"), function(x) {assign(x,get(x),envir= .GlobalEnv)})
    }
    #
    if(debugging){cat("#check: Rényi parallel function.")}
    plr<<-TRUE
    if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
      cls <- parallel::makeCluster(nc.cores,typedebugging=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
    } else if( cluster.type=="MPI" ) {
      cls <- makeMPIcluster(nc.cores,outfile="",useXDR=FALSE,methods=FALSE,output="")
    }
    registerDoSNOW(cls)
    clusterCall(cl=cls, function() library("parallel"))
     if(isfloat) {
       parallel::clusterExport(cl=cls, varlist=c("mfactor"))
     }
    on.exit(stopCluster(cls)) # Close the clusters on exit
    gc()

      if(mode == "single") {
      if( Shannon ) {
        outP <- ShannonP(rasterm, w, base, debugging, type)
      }
      else if (BergerParker){
        outP <- BergerParkerP(rasterm, w, base, debugging, type)
      }
      else{
        outP <- IndexP(rasterm, w, alpha, base, debugging, type)
      }
      return(do.call(cbind,outP))
      }
    else if(mode == "iterative"){
      outP <- list()
      for (ALPHA in start:stop){
        if(ALPHA == 1) {
          s <- "ShannonAlpha 1"
          out<- ShannonP(rasterm, w, base, debugging, type)
          outP[[s]] <- do.call(cbind,out)
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("RényiAlpha",as.character(ALPHA))
          out<- IndexP(rasterm, w, ALPHA, base, debugging, type)
          outP[[s]] <- do.call(cbind,out)
        }
      }
      message(("\nCalculation of Rényi's index complete.\n"))
      if (BergerParker){
        s <- "Berger-Parker"
        out <- BergerParkerP(rasterm, w, base, debugging, type)
        outP[[s]] <- do.call(cbind,out)
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
        return(outP)
    }
    else if( mode == "sequential"){
      outP <- list()
      for (ALPHA in val){
        if(abs(ALPHA-1)<.Machine$double.eps) {
          s <- "ShannonAlpha 1"
          out<- ShannonP(rasterm, w, base, debugging, type)
          outP[[s]] <- do.call(cbind,out)
          message(("\nCalculation of Shannon's index is also complete!\n"))
        } # End ShannonD
        else{
          s<-paste("RényiAlpha",as.character(ALPHA))
          out<- IndexP(rasterm, w, ALPHA, base, debugging, type)
          outP[[s]] <- do.call(cbind,out)
        }
      }
      message(("\nCalculation of Rényi's index complete.\n"))
      if (BergerParker){
        s <- "Berger-Parker"
        out <- BergerParkerP(rasterm, w, base, debugging, type)
        outP[[s]] <- do.call(cbind,out)
        message(("\nCalculation of Berger-Parker's index is also complete!\n"))
      }
      return(outP)
    }
  }
}


###nc=1

ShannonS <- function(rasterm, w, base, debugging, type){
  message("\nStarting Shannon-Wiener index calculation:\n")
  # Reshape values
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add "fake" columns and rows for moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  #
  ## Loop over all the pixels
  #
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      #if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
      #  out[rw-w,cl-w]<-NA
      #} else {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]))
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if(debugging) {
        message("Shannon-Wiener\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
      }
      tw_values<-as.vector(tw)
      p<-tw_values/sum(tw_values)
      if (type=="Rényi"){
        p_log<-log(p,base)
        out[rw-w,cl-w]<-(-(sum(p*p_log)))
        }
      else{
        p_log<-log(p)
        out[rw-w,cl-w]<-exp(-(sum(p*p_log)))
        }
      #}
    }   
    svMisc::progress(value=cl, max.value=(c((dim(rasterm)[2]+w)+(dim(rasterm)[1]+w))/2), progress.bar = FALSE)
  } 
  
  return(out)
}

BergerParkerS <- function(rasterm, w, base, debugging, type){
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  message("\nStarting Berger-Parker index calculation:\n")
  # Reshape values
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add "fake" columns and rows for moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  #
  ## Loop over all the pixels
  #
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      #if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
      #  out[rw-w,cl-w]<-NA
      #} else {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]))
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if(debugging) {
        message("Berger-Parker\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
      }
      tw_values<-as.vector(tw)
      p<-max(tw_values/sum(tw_values))
      if(type=="Rényi"){
        out[rw-w,cl-w]<-(log(1/p, base))
      }
      else{
        out[rw-w,cl-w]<-(1/p)
      }
      #}
    }   
    svMisc::progress(value=cl, max.value=(c((dim(rasterm)[2]+w)+(dim(rasterm)[1]+w))/2), progress.bar = FALSE)
  } 
  
  return(out)
}

IndexS <- function(rasterm, w, alpha, base, debugging, type){
  message("\nStarting Rényi index calculation with parameter value=",alpha,"\n")
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Reshape values
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Add fake columns and rows for moving window
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  # Loop over each pixel
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      #if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
      #  out[rw-w,cl-w]<-NA
      #} else {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if(debugging) {
        message("Rényi\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      #if clause to exclude windows with only 1 category
      #if(length(tw_values) < 2) {
      #  out[rw-w,cl-w]<-NA
      #} else {
      p <- tw_values/sum(tw_values)
      if (type=="Rényi"){
        out[rw-w,cl-w]<-1/(1-alpha) * drop(log(sum(p^alpha),base))
      }
      else {
        out[rw-w,cl-w]<-drop(sum(p^alpha))^(1/(1-alpha)) 
      }
      #}
      #} 
    } 
    svMisc::progress(value=cl, max.value=(c((dim(rasterm)[2]+w)+(dim(rasterm)[1]+w))/2), progress.bar = FALSE)
  } # End of for loop 
  return(out)
}

###nc>1

ShannonP<-function(rasterm, w, base, debugging, type){
  #
  ## Reshape values
  #
  values <- as.numeric( as.factor(rasterm) )
  rasterm_1 <- matrix(data = values, nrow = dim(rasterm)[1], ncol = dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor <- matrix(NA, ncol = dim(rasterm)[2], nrow = w)
  ver <- matrix(NA, ncol = w, nrow = dim(rasterm)[1]+ w * 2)
  trasterm <- cbind(ver, rbind(hor,rasterm_1,hor), ver)
  rm(hor, ver, rasterm_1, values); gc()
  #
  ## Progression bar
  #
  pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  ShannonOP <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.snow = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    ShannonOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      # if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
      #   vv<-NA
      #   return(vv)
      # } 
      #else {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if( debugging ) {
        message("Shannon - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      p <- tw_values/sum(tw_values)
      if (type=="Rényi"){
        vv <- (-(sum(p*log(p,base))))
      }
      else{
        vv <- exp(-(sum(p*log(p))))
      }
      return(vv)
      #}
    })
    return(ShannonOut)
  } # End Shannon - parallelized
  message(("\n\n Parallel calculation of Shannon's index complete.\n"))
  return(ShannonOP)
}

BergerParkerP<-function(rasterm, w, base, debugging, type){
  #
  ## Reshape values
  #
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  rm(hor,ver,rasterm_1,values); gc()
  #
  ## Progression bar
  #
  pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #
  ## Start the parallelized loop over iter
  #
  BergerParkerOP <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.snow = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    BergerParkerOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      # if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
      #   vv<-NA
      #   return(vv)
      # } 
      #else {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if( debugging ) {
        message("Berger-Parker - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", window)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      p <- max(tw_values/sum(tw_values))
      if (type=="Rényi"){
        vv <- (log(1/p, base))
      }
      else{
        vv <- 1/p
      }
      return(vv)
      #}
    })
    return(BergerParkerOut)
  } # End Berger-Parker - parallelized
  message(("\n\n Parallel calculation of Berger-Parker's index complete.\n"))
  return(BergerParkerOP)
}


IndexP<-function(rasterm, w, alpha, base, debugging, type){
  #
  ## Reshape values
  #
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  rm(hor,ver,rasterm_1,values); gc()
  #
  ## Progression bar
  #
  pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #
  ## Start the parallelized loop over iter
  #
  IndexOP <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.snow = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    IndexOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      # if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
      #   vv<-NA
      #   return(vv)
      # } 
      #else {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if( debugging ) {
        message("Rényi - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", window)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      p <- tw_values/sum(tw_values)
      if (type=="Rényi"){
        vv <- 1/(1-alpha) * drop(log(sum(p^alpha),base))
      }
      else{
        vv <- drop(sum(p^alpha))^(1/(1-alpha))
      }
      return(vv)
      #}
    })
    return(IndexOut)
  } # End classic Rényi - parallelized
  message(("\n\n Parallel calculation of Rényi's index complete.\n"))
  return(IndexOP)
}

m1<-matrix(0.5,ncol=3, nrow=3)
m2<-matrix(0.8,ncol=3, nrow=3)
m3<-matrix(0.8,ncol=3, nrow=3)
m4<-matrix(0.8,ncol=3, nrow=3)
M1<-cbind(rbind(m1,m2),rbind(m3,m4))
plot(raster(M1))
a<-IndexComputation (M1, type = "Rényi", mode = "sequential", nc.cores = 1, alpha=c(1,3), BergerParker = TRUE)
plot(raster(a$`ShannonAlpha 1`))
