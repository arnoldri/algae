###################################################################
# SF functions

#' Extract XY coordinates from an object
#'
#' @export
ra_getxy <- function(x) unclass(x$geometry[[1]])

#' Convert global Longitude and Latitude to XY - defaults to NZ(XY)
#'
#' @export
ra_lonlat_to_xy <- function(lon, lat, crs=2193) {
  globalcrs <- 4326
  #marlcrs <- 2193
  gg <- st_as_sf(data.frame(lon=lon,lat=lat), coords=c("lon","lat"), crs=globalcrs, agr="constant")
  gg <- st_transform(gg, crs=crs)
  n <- nrow(gg)
  gg <- data.frame(lon=lon, lat=lat, t(array(unlist(gg),dim=c(2,n))))
  names(gg)[3:4] <- c("x","y")
  return(gg)
}

#' Convert XY to global Longitude and Latitude - defaults to NZ(XY)
#'
#' @export
ra_xy_to_lonlat <- function(x, y, crs=2193) {
  globalcrs <- 4326
  #marlcrs <- 2193
  gg <- st_as_sf(data.frame(x=x,y=y), coords=c("x","y"), crs=crs, agr="constant")
  gg <- st_transform(gg, crs=globalcrs)
  n <- nrow(gg)
  gg <- data.frame(t(array(unlist(gg),dim=c(2,n))), x=x, y=y)
  names(gg)[3:4] <- c("lon","lat")
  return(gg)
}


###################################################################
# Other Utilities

#' Make dates dd into the date of the preceding Monday
#'
#' @export
make.monday <- function(dates) {
  as.Date(-3+7*(as.integer(dates+3)%/%7), origin="1970-01-01")
}


#' Check for uniqueness
#'
#' @description Check for uniqueness
#'
#' @export
is.unique <- function(x) length(x)==length(unique(x))


#' Run the pacman 00LOCK file destroyer
#'
#' @description Stop 00LOCK files getting in the way of package builds
#'
#' @export
unlock <- function() {
  for(.lib in .libPaths()) {
    ff <- grep("00LOCK", list.files(.lib), value=TRUE)
    print(ff)
    if(length(ff)>0) unlink(paste0(.lib,"/",ff), recursive=TRUE, force=TRUE)
  }
  #pacman::p_unlock()
  invisible()
}
