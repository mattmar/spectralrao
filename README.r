#spectralrao(...)
#Description: Applying Rao's index to remote sensing data

#For parallel computing in Ubuntu:
#you must have the following libraries installed before calling spectralrao()
apt-get update
apt-get upgrade
apt-get install mpi
apt-get install libopenmpi-dev
apt-get install r-cran-rmpi

#For parallel computing in MAC|Windows:
#Please check out this page http://www.stats.uwo.ca/faculty/yu/Rmpi/
#Rmpi spawn function is not supported in Windows
#https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/

##Minimum example to check that the function is working OK
###Result should be 0.222 in the central cell
a<-matrix(c(-0.5,1,0.5,1,0.5,1,1,1,0.5),ncol=3,nrow=3)
spectralrao(input=a,window=3,distance_m="euclidean",na.tolerance=0,shannon=FALSE)
spectralrao(input=a,window=3,distance_m="euclidean",na.tolerance=0,shannon=FALSE,nc.cores=2)

##Less minimum example
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
r3<-matrix(data=c(-0.5,-0.5,-0.5,-0.5,0.8,-0.5,-0.5,-0.5,-0.5,0.8,0,0,0,0,0,0,0,0.8,0.8,0.8,0,0.5,0.8,0.8,0.5),nrow=5,ncol=5,byrow=F)

###Run the function on one dimension
raomatrix<-spectralrao(input=r2,distance_m="euclidean",window=3)

###Plot results
library('RColorBrewer')
clp<-brewer.pal(9,"RdYlGn")

png("~/spectralrao_monodimensional.png",width = 300, height = 105,res=300,type = c("cairo"),units="mm",pointsize="15")
par(mfrow=c(1,3),mar=c(5, 3, 4, 5),family="Arial")
raster::plot(raster(r2),main="Layer 1 (e.g., (NDVI)")
raster::plot(raster(raomatrix[[2]]),main="Shannon's H'",legend=F,col=clp)
raster::plot(raster(raomatrix[[1]]),main="Univariate Rao's\n (Euclidean distance)",legend=T,col=clp)
dev.off()

###Check running time for parallelized and sequential functions on one dimension
r1<-matrix(rpois(25000,lambda=5),nrow=500,ncol=500)
system.time(raop<-spectralrao(input=r1,distance_m="euclidean",window=3,shannon=T,na.tolerance=1, nc.cores=2, cluster.type="SOCK")) #36.566
system.time(raos<-spectralrao(input=r1,distance_m="euclidean",window=3,shannon=FALSE,na.tolerance=1)) #146.668

###Plot results enhanced [ggplot+standard legend]
library('RColorBrewer')
clp<-brewer.pal(9,"RdYlGn")

###Convert raster to dataframe
df <- melt(lapply(lapply(list(r2,raomatrix[[1]],raomatrix[[2]]),raster),as.data.frame,xy=T),id.vars=c("x", "y"), variable.name="index")

###Add correct label for facets
df$L1<-factor(df$L1,labels=c("NDVI","Rao","Shannon"),order=T)

###Plot with ggplot2
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

###Save with proper resolution
ggsave("~/raoshannon.png",dpi=300,height=3,width=9)

###Rao vs Shannon
plot(raomatrix[[1]]~raomatrix[[2]],pch=16,col="grey",cex=2,xlab="ShannonD",ylab="RaoQ")

###Example for Rao with distance matrix defined by the user
###Distance matrix is a dataframe or matrix in the form: A | B | distance A-B
###All values in the matrix must have an associated distance.
set.seed(001)
n1=rep(c(1,3,2),each=70)
xmean=rnorm(210,0,0.4)
n1_clusterd=as.integer(n1+xmean)
n2<-matrix(n1_clusterd,ncol=15,nrow=14)
dist_user<-matrix(c(0,0,0,1,1,1,2,2,2,3,3,3,1,2,3,0,2,3,0,1,3,0,1,2,1,2,3,1,1,2,2,1,1,3,2,1),ncol=3)

raomatrix<-spectralrao(input=n2,window=3,distance_m=dist_user,na.tolerance=1)

par(mfrow=c(1,2),mar=c(5, 3, 4, 5),family="Arial")
raster::plot(raster(n2))
raster::plot(raster(raomatrix[[1]]))

###Example for Rao with distance function defined by the user
weir_distance <- function(x,y,...) {
	return(cos(y-x)-sin(y-x))
}

ds1<-spectralrao(input=r1,distance_m=weir_distance,window=3,shannon=FALSE,na.tolerance=1)
ds2<-spectralrao(input=r1,distance_m="euclidean",window=3,shannon=FALSE,na.tolerance=1)

###Comparison
par(mfrow=c(1,2))
raster::plot(raster(ds1[[1]]),main="Weir_distance")
raster::plot(raster(ds2[[1]]),main="Euclidean_distance")

###Run the function as multidimensional RaoQ
raomatrix<-spectralrao(input=list(r2,r4),window=3,mode="multidimension",distance_m="euclidean",na.tolerance=1,rescale=FALSE)

###Comparison
png("~/spectralrao_multidimensional.png",width = 300, height = 105,res=300,type = c("cairo"),units="mm",pointsize="15")
par(mfrow=c(1,3),mar=c(5, 3, 4, 5))
raster::plot(raster(n1),main="Layer 1 (e.g., NDVI)")
raster::plot(raster(r4),main="Layer 2 (e.g., Soil pH)")
raster::plot(raster(raomatrix[[1]]),main="Multidimensional Rao's Q\n(Euclidean distance)")
dev.off()

###Run the function as multidimensional RaoQ and RasterLayer as input
###Download NDVI from NASA at 0.1 degrees and derive Rao's and Shannon index; NASA likes to change urls, the one below may be outdated
cd ~
wget -O example_modis_ndvi_2019.tiff "https://neo.sci.gsfc.nasa.gov/servlet/RenderData?si=1772769&cs=rgb&format=FLOAT.TIFF&width=3600&height=1800"

wget "http://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/CNTR_2014_60M_SH.zip"
unzip CNTR_2014_60M_SH.zip

#R
#Load packages
library(rgdal)
library(raster)
library(rgeos)

#Import Europe boundary
world <- readOGR(dsn = "./CNTR_2014_60M_SH/Data/", layer = "CNTR_RG_60M_2014")

#Importa NDVI, set NA and crop on Europe
ndvi2019 <- raster("~/example_modis_ndvi_2019.tiff")
ndvi2019<-reclassify(ndvi2019, cbind(99999, NA), right=FALSE)
ndvi2019_eu <- crop(ndvi2019,extent(-15,45,35,60))

#Change projection and mask
world <- spTransform(world,crs(ndvi2019_eu))
ndvi2019_eu_masked<-mask(ndvi2019_eu,world)

#Run raomatrix
raomatrix <- spectralrao(ndvi2019_eu_masked, distance_m="euclidean", window=9, shannon=TRUE, nc.cores=8)

###Spatial matrices
raos <- raster(raomatrix[[1]],template=raster(ndvi2019_eu_masked))
shan <- raster(raomatrix[[2]],template=raster(ndvi2019_eu_masked))

###Plot results
clrbrw<-colorRampPalette(c("#fc8d59","#ffffbf","#91cf60"))(10)

png("~/all_modis_ndvi_june_2015_2km_1.png",width = 149, height = 290,res=300,type = c("cairo"),units="mm",pointsize="15")
par(mfrow=c(3,1))

raster::plot(ndvi2019_eu_masked, main="MODIS NDVI June 2015",legend=FALSE)
raster::plot(ndvi2019_eu_masked, legend.only=TRUE, axis.args=list(at=seq(-1,0.99,length.out=5), labels=seq(-1,1,length.out=5)))

raster::plot(mask(shan,ndvi2019_eu_masked),main="Shannon's index",legend=FALSE)
raster::plot(shan, legend.only=TRUE, axis.args=list(at=seq(0,2.15,length.out=5), labels=round(seq(0,2,length.out=5),1)))

raster::plot(mask(raos,ndvi2019_eu_masked),main="Rao's Q"
,legend=FALSE)
raster::plot(raos, legend.only=TRUE, axis.args=list(at=seq(0.01,0.4,length.out=5), labels=round(seq(0,0.4,length.out=5),2)))

dev.off()

###Compare multiple distance, monodimensional rao
dst<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")

outl<-list()
for (d in 1:length(dst)) {
    outl[[d]] <- spectralrao(r2, distance_m=dst[d],window=3,shannon=FALSE)
    names(outl[[d]]) <- dst[d]
}

par(mfrow=c(2,3))
sapply(outl, function(x) {raster::plot(raster(x[[1]]),main=names(x))})

###Compare multiple distance, monodimensional rao
mdst<-c("euclidean","manhattan","canberra","minkowski","mahalanobis")

moutl<-list()
for (d in 1:length(mdst)) {
    moutl[[d]] <- spectralrao(list(r2,r4), distance_m=mdst[d],mode="multidimension",window=3, lambda=0.1,shannon=FALSE)
    names(moutl[[d]]) <- mdst[d]
}

par(mfrow=c(2,3))
sapply(moutl, function(x) {raster::plot(raster(x[[1]]),main=names(x))})

#Check out how na.tolerance parameter works
set.seed(27)
rna<-matrix(c(rpois(2500,lambda=5),rep(NA,525)),nrow=55,ncol=55)

###Run the function with different na.tolerance
raona<-spectralrao(input=rna,distance_m="euclidean",window=9,shannon=FALSE,na.tolerance=0)
ranona<-spectralrao(input=rna,distance_m="euclidean",window=9,shannon=FALSE,na.tolerance=1)

###Comparison
par(mfrow=c(1,3),family="Arial")
raster::plot(raster(rna),main="Layer 1 (e.g., (NDVI)")
raster::plot(raster(raona[[1]]),main="Rao na.tolerance=0'",legend=F)
raster::plot(raster(ranona[[1]]),main="Rao na.tolerance=1'",legend=F)

#Multidimensional with rasterlayer as input, to be done
#raomatrix<-spectralrao(matrix=list(ndvi,ndvi+7),window=3,mode="multidimension",shannon=FALSE)
###Plot the output
#raster::plot(raster(raomatrix[[1]]))
