library(raster)
ras = raster("/Users/mqwilber/Repos/rsf_swine/data/covariate_data/croplayer/tejon/tejon_fruit_and_nuts_2015.grd")
cras = crop(ras, extent(c(xmin=-118.80, xmax=-118.78, ymin=35, ymax=35.01)))

library(data.table)
pigdat = fread("/Users/mqwilber/Repos/rsf_swine/data/formatted/full_pig_data.csv")
tejon = pigdat[pigID == "tejonM302"]


pdf("tejon_ras.pdf", width=5, height=5)
plot(ras, col=terrain.colors(255)[c(255, 1)], xlab="Longitude", ylab="Latitude")
dev.off()

pdf("tejon_ras_w_pig.pdf", width=5, height=5)
plot(ras, col=terrain.colors(255)[c(255, 1)], xlab="Longitude", ylab="Latitude")
points(tejon$longitude, tejon$latitude, cex=0.09, pch=19)
dev.off()


source("/Users/mqwilber/Repos/rsf_swine/code/pigfxns.R")
oneinds = which(values(cras) == 1)
ones = sample(oneinds, 100, replace=F)
tras = cras
values(tras) = 0
values(tras)[ones] = 1

pdf("tejon_samp_ras.pdf", width=8, height=5)
plot(tras, col=terrain.colors(255)[c(255, 1)], xlab="Longitude", ylab="Latitude")
dev.off()


gxy = get_distance_to_gradient(cras)
long = rasterToPoints(cras)[, c('x')]
lat = rasterToPoints(cras)[, c('y')]
xoff = values(gxy$xgrad)
yoff = values(gxy$ygrad)

library(plotrix)
pdf("tejon_vect.pdf", width=8, height=5)
plot(cras, col=terrain.colors(255)[c(255, 1)], xlab="Longitude", ylab="Latitude")
vectorField(xoff, yoff, long, lat, headspan=0.05)
dev.off()


shp = shapefile("/Users/mqwilber/Repos/rsf_swine/data/covariate_data/croplayer/tejon/test.shp")

pdf("tejon_shapefile.pdf",  width=8, height=8)
plot(shp)
dev.off()

#crs(shp) = "+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

library(rgeos)
rsp = rasterToPoints(ras, spatial=TRUE)
short = rsp[sample(1:dim(rsp)[1], 1), ]

#crs(short) = "+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
#dists = gDistance(short, shp, byid=TRUE)
centers = gCentroid(shp[shp$DN == 1, ], byid=TRUE)

pdf("tejon_distance.pdf",  width=8, height=8)
plot(shp)
for(i in 1:nrow(centers@coords)){
	dat = rbind(centers[i, ]@coords, short@coords)
	#points(dat[, 1], dat[, 2], col="black", pch=19, cex=0.5)
	lines(dat[, 1], dat[, 2], col="red")
}
points(short@coords[, 1], short@coords[, 2], pch=19, cex=2, col="blue")
dev.off()
#points(gCentroid(shp[shp$DN == 1, ], byid=TRUE), col="blue", cex=0.3, pch=19)
