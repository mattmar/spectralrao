# spectralrao(...)
#Description: Applying Rao's index to remote sensing data'

##Example
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
raomatrix<-spectralrao(input=r1,distance_m="euclidean",window=3,shannon=TRUE,na.tolerance=1)

###Comparison
#Color palette
library('RColorBrewer')
clp<-brewer.pal(9,"RdYlGn")

png("~/spectralrao_monodimensional.png",width = 300, height = 105,res=300,type = c("cairo"),units="mm",pointsize="15")
par(mfrow=c(1,3),mar=c(5, 3, 4, 5),family="Arial")
raster::plot(raster(r2),main="Layer 1 (e.g., (NDVI)")
raster::plot(raster(raomatrix[[2]]),main="Shannon's H'",legend=F,col=clp)
raster::plot(raster(as.matrix(c(rep(2.1,10),rep(1.9,10)))), legend.only=TRUE,col=clp)
raster::plot(raster(raomatrix[[1]]),main="Univariate Rao's\n (Euclidean distance)",legend=T,col=clp)
dev.off()

###Check running time for parallelized and sequential functions on one dimension
r1<-matrix(rpois(25000,lambda=5),nrow=500,ncol=500)
system.time(raop<-spectralrao(input=r1,distance_m="euclidean",window=3,shannon=FALSE,na.tolerance=1, nc.cores=8)) #75.669
system.time(raos<-spectralrao(input=r1,distance_m="euclidean",window=3,shannon=FALSE,na.tolerance=1)) #89.064

###Comparison enhanced [standard legend]
#Color palette
library('RColorBrewer')
clp<-brewer.pal(9,"RdYlGn")

#Convert raster to dataframe
df <- melt(lapply(lapply(list(r2,raomatrix[[1]],raomatrix[[2]]),raster),as.data.frame,xy=T),id.vars=c("x", "y"), variable.name="index")

#Add correct label for facets
df$L1<-factor(df$L1,labels=c("NDVI","Rao","Shannon"),order=T)

#plot with ggplot2
ggplot(data=df, aes(x=x, y=y)) + 
geom_raster(aes(fill=value)) +
coord_equal() +
facet_wrap(~L1, nrow=1) + theme_bw() +
scale_fill_gradientn(colours=clp) +
geom_text(aes(label=round(df$value, digits = 2)),size=1.5) +
theme(legend.title=element_blank(),
	strip.text = element_text(family="Arial",size=12, colour = "black", angle = 0, face = "bold")) +
scale_x_continuous(expand = c(0,0)) +
scale_y_continuous(expand = c(0,0))

#Save with proper resolution
ggsave("~/raoshannon.png",dpi=300,height=3,width=9)

###Raos vs Shannon
plot(raomatrix[[1]]~raomatrix[[2]],pch=16,col="grey",cex=2,xlab="ShannonD",ylab="RaoQ")

###Run the function as multidimensional RaoQ
raomatrix<-spectralrao(input=list(r2,r4),window=3,mode="multidimension",distance_m="euclidean",na.tolerance=1,rescale=FALSE)

###Comparison
png("~/spectralrao_multidimensional.png",width = 300, height = 105,res=300,type = c("cairo"),units="mm",pointsize="15")
par(mfrow=c(1,3),mar=c(5, 3, 4, 5))
raster::plot(raster(r2),main="Layer 1 (e.g., NDVI)")
raster::plot(raster(r4),main="Layer 2 (e.g., Soil pH)")
raster::plot(raster(raomatrix[[1]]),main="Multidimensional Rao's Q\n(Euclidean distance)")
dev.off()

###Run the function as multidimensional RaoQ and RasterLayer as input
###Download NDVI from NASA at 0.1 degrees and derive Rao's and Shannon index; NASA likes to change urls, the one below may be outdated
cd ~
wget -O example_modis_ndvi_2015.tiff "http://neo.sci.gsfc.nasa.gov/servlet/RenderData?si=1690249&cs=rgb&format=TIFF&width=3600&height=1800"

wget "http://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/CNTR_2014_60M_SH.zip"
unzip CNTR_2014_60M_SH.zip

#R
#Load packages
library(rgdal)
library(raster)
library(rgeos)

#Import Europe boundary
world <- readOGR(dsn = "./CNTR_2014_60M_SH/Data/", layer = "CNTR_RG_60M_2014")

#Importa NDVI
ndvi2015 <- raster("~/example_modis_ndvi_2015.tiff", crs=CRS("+proj=longlat +datum=WGS84"))
extent(ndvi2015) <- extent(c(xmn=-180,xmx=180,ymn=-90,ymx=90))
ndvi2015_eu <- crop(ndvi2015,extent(-15,45,35,60))

# Change projection and mask
world <- spTransform(world,crs(ndvi2015_eu))
ndvi2015_eu_masked<-mask(ndvi2015_eu,world)

#Transform from 0-255 to -1-1
ndvi2015_final <- 2/255*ndvi2015_eu_masked-1 #Rescale in 0-1 interval

#Run raomatrix
raomatrix <- spectralrao(ndvi2015_final, distance_m="euclidean", window=9, shannon=TRUE, debugging=TRUE)

###Spatial matrices
raos <- raster(raomatrix[[1]],template=raster(ndvi2015_final))
shan <- raster(raomatrix[[2]],template=raster(ndvi2015_final))

###Plot the results
clrbrw<-colorRampPalette(c("#fc8d59","#ffffbf","#91cf60"))(10)

png("~/all_modis_ndvi_june_2015_2km_1.png",width = 149, height = 290,res=300,type = c("cairo"),units="mm",pointsize="15")
par(mfrow=c(3,1))

raster::plot(ndvi2015_final, main="MODIS NDVI June 2015",legend=FALSE)
raster::plot(ndvi2015_final, legend.only=TRUE, axis.args=list(at=seq(-1,0.99,length.out=5), labels=seq(-1,1,length.out=5)))

raster::plot(mask(shan,ndvi2015_final),main="Shannon's index",legend=FALSE)
raster::plot(shan, legend.only=TRUE, axis.args=list(at=seq(0,2.15,length.out=5), labels=round(seq(0,2,length.out=5),1)))

raster::plot(mask(raos,ndvi2015_final),main="Rao's Q"
,legend=FALSE)
raster::plot(raos, legend.only=TRUE, axis.args=list(at=seq(0.01,0.4,length.out=5), labels=round(seq(0,0.4,length.out=5),2)))

dev.off()


#Compare multiple distance, monodimensional rao
dst<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")

outl<-list()
for (d in 1:length(dst)) {
    outl[[d]] <- spectralrao(r2, distance_m=dst[d],window=3,shannon=FALSE)
    names(outl[[d]]) <- dst[d]
}

par(mfrow=c(2,3))
sapply(outl, function(x) {raster::plot(raster(x[[1]]),main=names(x))})


#Compare multiple distance, monodimensional rao
mdst<-c("euclidean","manhattan","canberra","minkowski","mahalanobis")

moutl<-list()
for (d in 1:length(mdst)) {
    moutl[[d]] <- spectralrao(list(r2,r4), distance_m=mdst[d],mode="multidimension",window=3, lambda=0.1,shannon=FALSE)
    names(moutl[[d]]) <- mdst[d]
}

par(mfrow=c(2,3))
sapply(moutl, function(x) {raster::plot(raster(x[[1]]),main=names(x))})


#Multidimensional with rasterlayer as input
#raomatrix<-spectralrao(matrix=list(ndvi,ndvi+7),window=3,mode="multidimension",shannon=FALSE)
###Plot the output
#raster::plot(raster(raomatrix[[1]]))