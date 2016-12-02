# Import NDVI from GRASS and calculate shannon/rao
library(rgrass7)
mm<-"ndvi_2015_06_uk_5km"
ndvi_2015_06_uk_5km<-readRAST(paste(mm,sep=""))
raomatrix<-spectralrao(get(mm),distance_m="euclidean",window=9,shannon=TRUE)

a <- raster(raomatrix[[2]],template=raster(mm))
b <- raster(raomatrix[[1]],template=raster(mm))

png("~/modis_ndvi_june_2015_2km.png",width = 300, height = 300,res=300,type = c("cairo"),units="mm",pointsize="16")
raster::plot(raster(mm))
dev.off()

png("~/shannon_10px_modis_ndvi_june_2015_2km.png",width = 300, height = 300,res=300,type = c("cairo"),units="mm",pointsize="16")
raster::plot(mask(a,raster(mm)))
dev.off()

png("~/rao_10px_modis_ndvi_june_2015_2km.png",width = 300, height = 300,res=300,type = c("cairo"),units="mm",pointsize="16")
raster::plot(mask(b,raster(mm)))
dev.off()

png("~/all_modis_ndvi_june_2015_2km.png",width = 900, height = 300,res=300,type = c("cairo"),units="mm",pointsize="35")
par(mfrow=c(1,3))
raster::plot(raster(mm))
raster::plot(a)
raster::plot(mask(b,raster(mm)))
dev.off()
