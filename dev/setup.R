# More setting up

bbox.tasmar <- c(1570,5370,1725,5520)*1000
xlim.tasmar <- bbox.tasmar[c(1,3)]
ylim.tasmar <- bbox.tasmar[c(2,4)]

bbox.marlselect <- c(1645,5415,1720,5495)*1000
bbox.marlselect <- c(1645,5412,1725,5502)*1000
xlim.marlselect <- bbox.marlselect[c(1,3)]
ylim.marlselect <- bbox.marlselect[c(2,4)]
xlim.havelock <- c(1660,1675)*1000
ylim.havelock <- c(5425,5440)*1000
xlim.nydia <- c(1662,1677)*1000
ylim.nydia <- c(5437,5450)*1000
xlim.arapawa <- c(1695,1718)*1000
ylim.arapawa <- c(5432,5452)*1000
xlim.kp <- c(1662,1695)*1000
ylim.kp <- c(5427,5474)*1000


# Coastal polygons for mapping
bbox.marlselect.polygon <- st_polygon(list(cbind(bbox.marlselect[c(1,3,3,1,1)],
                                                 bbox.marlselect[c(2,2,4,4,2)])))
coast.marlselect <- st_crop(coast, bbox.marlselect.polygon)

bbox.tasmar.polygon <- st_polygon(list(cbind(bbox.tasmar[c(1,3,3,1,1)],
                                             bbox.tasmar[c(2,2,4,4,2)])))
coast.tasmar <- st_crop(coast, bbox.tasmar.polygon)

# polygons in model
# monitored polygons
habspoly_sf <- geojson_sf("data/HABS_exp_6_release_polygons.geojson")
habspoly_sf <- st_transform(habspoly_sf, crs=st_crs(coast))

# File 3
marlpoly_sf.latlon <- geojson_sf("data/Entire_marlb_polygons_volume_above_15m.geojson")
#st_crs(marlpoly_sf.latlon)
marlpoly_sf.latlon <- st_transform(marlpoly_sf.latlon, crs=st_crs(coast))
marlpoly_sf.latlon$volume <- marlpoly_sf.latlon$volume/1e9 # cubic km
marlpoly_sf.latlon$volume_15 <- marlpoly_sf.latlon$volume_15/1e9 # cubic km
# range(marlpoly_sf.latlon$CLUSTER_ID); nrow(marlpoly_sf.latlon) # 1 409

# use the CLUSTER_ID values from the XY file
#marlpoly_sf.latlon <- rename(marlpoly_sf.latlon, CLUSTER_ID.old=CLUSTER_ID)
#marlpoly_sf.latlon <- merge(marlpoly_sf.latlon, marlpoly_sf.xy[,c("id","CLUSTER_ID")] %>% st_drop_geometry())
#as.data.frame(marlpoly_sf.latlon)[marlpoly_sf.latlon$CLUSTER_ID!=marlpoly_sf.latlon$CLUSTER_ID.old,]
# alternative method without merging: reduce CLUSTER_ID values by 1 if above 400
marlpoly_sf.latlon <- rename(marlpoly_sf.latlon, CLUSTER_ID.old=CLUSTER_ID)
marlpoly_sf.latlon$CLUSTER_ID <- ifelse(marlpoly_sf.latlon$CLUSTER_ID.old>400,
                                        marlpoly_sf.latlon$CLUSTER_ID.old-1,marlpoly_sf.latlon$CLUSTER_ID)

marlpoly_sf <- marlpoly_sf.latlon
marlpoly_sf.trim <- ms_erase(marlpoly_sf, coast.marlselect)
marlpoly_sf.trim$area <- st_area(marlpoly_sf.trim)/1e6
#nrow(marlpoly_sf.trim) # 408

globalcrs <- 4326 # standard lat/lon global coordinates
marlcrs.latlon <- st_crs(marlpoly_sf.latlon)
marlcrs <- marlcrs.latlon

# Make adjacency matrix
adjlist.1 <- st_touches(marlpoly_sf.trim, marlpoly_sf.trim)
nn <- length(adjlist.1)
names(adjlist.1) <- marlpoly_sf.trim$CLUSTER_ID
adjlist <- lapply(adjlist.1, function (x) marlpoly_sf.trim$CLUSTER_ID[x])
adjlist <- adjlist[order(as.numeric(names(adjlist)))]
adjmat <- array(FALSE, dim=c(nn,nn))
dimnames(adjmat) <- list(names(adjlist),names(adjlist))
for(i in 1:nn) adjmat[i,adjmat[[i]]] <- TRUE
#all(adjmat==(t(adjmat))) # symmetric as it must be

# Connection matrices
load("data/conn2018.Rda") # creates conn
matdates <- as.Date(names(conn),format="%Y%m%d")

# Covariates for the same period
fname <- "data/Weekly_environmental_variables.xlsx"
covnames <- excel_sheets(fname)
covlist <- lapply(1:length(covnames), function(i) as.matrix(read_excel(path=fname,sheet=i,col_names=FALSE)))
names(covlist) <- covnames
for(vname in covnames) colnames(covlist[[vname]]) <- matdates

# Events
load(file="data/event2018.Rda") # creates object: event (indexed: CLUSTER_ID, week) # 6351 rows
load(file="data/polydat2018.Rda") # creates object: polydat_sf (indexed: CLUSTER_ID) # 408 rows
