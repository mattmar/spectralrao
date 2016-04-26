# spectralrao(...)
 Description: Applying Rao's index to remote sensing data

##Example:
###Random simulated spectral matrix
set.seed(26)
r1<-matrix(rpois(2500,lambda=5),nrow=50,ncol=50)

###Clustered simulated spectral matrix
xy1 <- matrix(rnorm(25, 0, .25), ncol=5,nrow=5)
xy2 <- matrix(rnorm(25, -.5, .1), ncol=5,nrow=5)
xy3 <- matrix(rnorm(25, .5, .1), ncol=5,nrow=5)
xy4 <- matrix(rnorm(25, .8, .25), ncol=5,nrow=5)
r2 <- cbind(rbind(xy1, xy2),rbind(xy3,xy4))

xy1 <- matrix(rnorm(20, 4.5, .01), ncol=5,nrow=4)
xy2 <- matrix(rnorm(30, 5.5, .1), ncol=5,nrow=6)
xy3 <- matrix(rnorm(40, 6.5, .1), ncol=5,nrow=8)
xy4 <- matrix(rnorm(10, 7.0, .05), ncol=5,nrow=2)
r4 <- t(cbind(rbind(xy1, xy2),rbind(xy3,xy4)))

###User simulated spectral matrix
r3<-matrix(data=c(-0.5,-0.5,-0.5,-0.5,0.8, -0.5,-0.5,-0.5,-0.5,0.8,0,0,0,0,0,0,0,0.8,0.8,0.8,0,0.5,0.8,0.8,0.5),nrow=5,ncol=5,byrow=F)


###Run the function on one dimension
raomatrix<-spectralrao(matrix=r2,distance_m="euclidean",window=3,shannon=TRUE)

###Comparison
par(mfrow=c(1,3))
raster::plot(raster(r2),main="landscape (NDVI)")
raster::plot(raster(raomatrix[[2]]),main="Shannon's H'")
raster::plot(raster(raomatrix[[1]]),main="Univariate Rao's Q")

###Raos vs Shannon
plot(raomatrix[[1]]~raomatrix[[2]],pch=16,col="grey",cex=2,xlab="ShannonD",ylab="RaoQ")

###Run the function as multidimensional RaoQ
raomatrix<-spectralrao(matrix=list(r2,r4),window=3,mode="multidimension",shannon=FALSE)

###Comparison
par(mfrow=c(1,3))
raster::plot(raster(r2),main="layer 1 (eg, NDVI)")
raster::plot(raster(r4),main="Layer 2 (eg, Soil pH)")
raster::plot(raster(raomatrix[[1]]),main="Multidimensional Rao's Q")

###Run the function as multidimensional RaoQ and RasterLayer as input
raomatrix<-spectralrao(matrix=list(ndvi,ndvi+7),window=3,mode="multidimension",shannon=FALSE)

###Plot the output
raster::plot(raster(raomatrix[[1]]))


###Import NDVI from GRASS and calculate shannon/rao (run only in a GRASS location with ndvi_2015_06_uk_5km as raster name)
library(rgrass7)
mm<-"ndvi_2015_06_uk_5km"
ndvi_2015_06_uk_5km<-readRAST(paste(mm,sep=""))
raomatrix<-spectralrao(get(mm),distance_m="euclidean",window=9,shannon=TRUE)

###Spatial matrices
raos <- raster(raomatrix[[1]],template=raster(mm))
shan <- raster(raomatrix[[2]],template=raster(mm))

###Plot the results
png("~/modis_ndvi_june_2015_2km.png",width = 300, height = 300,res=300,type = c("cairo"),units="mm",pointsize="16")
raster::plot(raster(mm))
dev.off()

png("~/shannon_9px_modis_ndvi_june_2015_2km.png",width = 300, height = 300,res=300,type = c("cairo"),units="mm",pointsize="16")
raster::plot(mask(shan,raster(mm)))
dev.off()

png("~/rao_9px_modis_ndvi_june_2015_2km.png",width = 300, height = 300,res=300,type = c("cairo"),units="mm",pointsize="16")
raster::plot(mask(raos,raster(mm)))
dev.off()

png("~/all_modis_ndvi_june_2015_2km.png",width = 900, height = 300,res=300,type = c("cairo"),units="mm",pointsize="35")
par(mfrow=c(1,3))
raster::plot(raster(mm))
raster::plot(shan)
raster::plot(mask(raos,raster(mm)))
dev.off()
