# Import NDVI from GRASS and calculate shannon/rao
library(rgrass7)
ndvi_2015_06_10km<-readRAST("ndvi_2015_06_10km")
raomatrix<-spectralrao(ndvi_2015_06_10km,distance_m="euclidean",window=10,shannon=TRUE)

a <- raster(raomatrix[[2]],template=raster(ndvi_2015_06_10km))
b <- raster(raomatrix[[1]],template=raster(ndvi_2015_06_10km))

png("~/modis_ndvi_june_2015_10km.png",width = 300, height = 300,res=300,type = c("cairo"),units="mm",pointsize="16")
raster::plot(raster(ndvi_2015_06_10km))
dev.off()

png("~/shannon_10px_modis_ndvi_june_2015_10km.png",width = 300, height = 300,res=300,type = c("cairo"),units="mm",pointsize="16")
raster::plot(mask(a,raster(ndvi_2015_06_10km)))
dev.off()

png("~/rao_10px_modis_ndvi_june_2015_10km.png",width = 300, height = 300,res=300,type = c("cairo"),units="mm",pointsize="16")
raster::plot(mask(b,raster(ndvi_2015_06_10km)))
dev.off()

png("~/all_modis_ndvi_june_2015_10km.png",width = 900, height = 300,res=300,type = c("cairo"),units="mm",pointsize="35")
par(mfrow=c(1,3))
raster::plot(raster(ndvi_2015_06_10km))
raster::plot(mask(a,raster(ndvi_2015_06_10km)))
raster::plot(mask(b,raster(ndvi_2015_06_10km)))
dev.off()
