require(RPostgreSQL)
require(raster)
require(rgeos)
require(rgdal)
require(maptools)
require(sp)
require(doMC)

#Establish connection with Boab
drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab
on.exit(dbDisconnect(con))

input <- raster("inst/extdata/Koala_HabSuit.tif")
input_res <- res(input)

output <- input
output[] <- NA

target <- "auscover.nsw_epsg3308_grid_bcul0"
target_res <- unname(unlist(dbGetQuery(con,paste0("SELECT ST_PixelHeight(rast), ST_PixelWidth(rast) FROM ",target," LIMIT 1;"))))

scale <- c(input_res[1]/target_res[1], input_res[2]/target_res[2])

#Get projection of target layer
proj.target.srid <- dbGetQuery(con,paste0("SELECT ST_SRID(rast) FROM ",target," WHERE rid=1;"))
proj.target <- dbGetQuery(con,paste0("SELECT proj4text FROM spatial_ref_sys WHERE srid=",proj.target.srid[[1]],";"))

#Create polygons overlay
overlay <- rasterToPolygons(input)

#Reproject to target layer
overlay_p <- spTransform(overlay, CRS(proj.target[[1]]))

#Load to boundary to server and assign projection
writeOGR(overlay_p, dsn="PG:dbname=qaeco_spatial user=qaeco password=Qpostgres15 host=boab.qaeco.com port=5432", layer="temp.overlay", driver="PostgreSQL", layer_options = "geometry_name=geom", overwrite_layer=TRUE)
dbGetQuery(con,paste0("SELECT UpdateGeometrySRID('temp','overlay','geom',",proj.target.srid[[1]],");"))

#Break raster into chunks for parallel processing
chunks <- split(seq_len(ncell(output)), ceiling(seq_along(seq_len(ncell(output)))/5000))

#Initialise parallel
registerDoMC(4)

#Perform spatial operation, extract value, and insert into output raster
output[] <- foreach(i = chunks, .combine = c) %dopar% {
  drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")
  unname(unlist(dbGetQuery(con,
    paste0("SELECT (ST_SummaryStats(ST_Union(ST_Clip(rast,geom)))).sum/",prod(scale)," 
    FROM auscover.nsw_epsg3308_grid_bcul0, temp.overlay 
    WHERE ST_Intersects(geom,rast) 
    AND ogc_fid BETWEEN ",min(i)," AND ",max(i)," 
    GROUP BY ogc_fid;")
    )))
}

plot(output)

writeRaster(output, filename="inst/extdata/Tree_Density.tif", format="GTiff", overwrite=TRUE)


##### Old code #####

# pb <- txtProgressBar(min = 0, max = ncell(input), style = 3)
# for (i in seq_len(ncell(input))) {
#   
#   #Create boundary from raster cell
#   boundary <- rasterToPolygons(rasterFromCells(input, i))
#   
#   #Reproject
#   boundary_p <- spTransform(boundary, CRS(proj.target[[1]]))
#   
#   #Load to boundary to server
#   writeOGR(boundary_p, dsn="PG:dbname=qaeco_spatial user=qaeco password=Qpostgres15 host=boab.qaeco.com port=5432", layer="temp.boundary", driver="PostgreSQL", layer_options = "geometry_name=geom", overwrite_layer=TRUE)
#   dbGetQuery(con,paste0("SELECT UpdateGeometrySRID('temp','boundary','geom',",proj.target.srid[[1]],");"))
# 
#   #Perform spatial operation, extract value, and insert into output raster
#   output[i] <- unname(unlist(dbGetQuery(con,paste0("SELECT (ST_SummaryStats(ST_Union(ST_Clip(raster.rast, polygon.geom), 'SUM'),1)).sum/",prod(scale)," AS rast FROM temp.boundary AS polygon, ",target," AS raster WHERE ST_Intersects(polygon.geom, raster.rast);"))))
#   
#   setTxtProgressBar(pb, i)
# }
# close(pb)
