library(sf)
library(maps)

if(.Platform$OS.type=="windows") {
  maprootdir <- "~/Documents/usr/data/maps/"
} else {
  maprootdir <- "~/data/maps/"
}
#list.files(maprootdir)
mapdir <- paste0(maprootdir,"lds-new-zealand-5layers-SHP/")
mapdir <- paste0(mapdir,"nz-coastlines-and-islands-polygons-topo-1500k/")
coast <- read_sf(dsn=mapdir)
bbox.coast <- st_bbox(coast)
bbox.coast.polygon <- st_polygon(list(cbind(bbox.coast[c(1,3,3,1,1)],
                                            bbox.coast[c(2,2,4,4,2)])))

marlborough <- read.csv("data/marlborough.csv")
nmarlborough <- nrow(marlborough)

bbox.marlborough <- c(1650,5420,1710,5470)*1000  ## main one
xlim.marlborough <- bbox.marlborough[c(1,3)]
ylim.marlborough <- bbox.marlborough[c(2,4)]
bbox.marlborough.polygon <- st_polygon(list(cbind(bbox.marlborough[c(1,3,3,1,1)],
                                                  bbox.marlborough[c(2,2,4,4,2)])))
coast.marlborough <- st_crop(coast, bbox.marlborough.polygon)


fill.dmat <- function(dmat) {
  # fill the distance matrix by chaining points together
  # return the shortest distance between each pair of points
  nm <- nrow(dmat)
  for(i in 1:nm) {
    d1 <- dmat[i,]
    idx <- which(d1>0 & d1<Inf)

    for(j in idx) {
      dj <- dmat[j,]
      jdx <- which(dj>0 & dj<Inf)
      jdx1 <- jdx[jdx!=i]
      d1[jdx1] <- pmin(d1[jdx1], d1[j]+dj[jdx1])
    }
    dmat[i,] <- d1
    dmat[,i] <- d1
  }
  return(dmat)
}
