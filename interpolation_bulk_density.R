## Packages 
install.packages(c("tidyverse","ggplot2","raster","sp","rgdal","sf","stars","gstat"))
my_packages <- c("tidyverse","ggplot2","raster","sp","rgdal","sf","stars","gstat")
lapply(my_packages, require, character.only = TRUE)
rm(my_packages)


## Making objects into spatial ----------------------------------------------------------------

# Loading objects
bulk <- read.csv("Input/BD_data_georef.csv")
bnd <- st_read("Input/vector/island_montreal.shp")

# Converting dataframes into spatial objects with 'sf' package
bulk_sf <- st_as_sf(bulk, coords = c("coord_X","coord_Y"), crs = 4326) # Making points vector into spatial object class in WGS84 coordinate system
bulk_reproj <- st_transform(bulk_sf, crs = 2950) # Transforming points vector coordinate system from WGS 84 into NAD83-MTM8
bulk_sp <- as_Spatial(bulk_reproj) # Converting sf spatial object into sp spatial object


## Interpolation: setup -----------------------------------------------------------

# Spatial grid according to extents of island object
extent(bnd)
x_range <- as.integer(c(267517,306653))
y_range <- as.integer(c(5029232,5062655))

# Spacing in spatial grid, 30m cells
grd <- expand.grid(x = seq(from = x_range[1], to = x_range[2], by = 30),
                   y = seq(from = y_range[1], to = y_range[2], by = 30))
rm(x_range); rm(y_range)

# Specifying coordinate system for grid
coordinates(grd) <- ~ x + y
proj4string(grd) <- CRS("+init=epsg:2950")
gridded(grd) <- TRUE
plot(grd, cex = 1.5)
points(bulk_sp, pch = 1, col = "red", cex = 1)

# Exploring variograms for kriging
variocloud <- variogram(BD_hybrid ~ 1, data = bulk_reproj, cloud = TRUE)
plot(variocloud) 
semivario <- variogram(BD_hybrid ~ 1, data = bulk_reproj)
plot(semivario)
semivario

# Acknowledging spatial correlation
hscat(BD_hybrid ~ 1, bulk_reproj, (0:9) * 100) # Scatterplot of lags (distances)

# Trying to fit best models onto variogram: manually
v <- variogram(BD_hybrid ~ 1, data = bulk_reproj)
vgm() # To check list of model functions
v_fit <- fit.variogram(v, vgm("Exp"))
v_fit
v_fit <- fit.variogram(v, vgm("Sph"))
v_fit
v_fit <- fit.variogram(v, vgm("Lin"))
v_fit
v_fit <- fit.variogram(v, vgm("Log"))
v_fit
# Range: value of last lag 
# Sill: point on y-axis where semivariogram starts to level off 
# Nugget: y-intercept 
# Partial sill: sill - nugget 

model_vario <- vgm(psill = 0.02119984,
                   model = "Lin",
                   nugget = 0.05460859,
                   range = 20943.01)
fit_vario <- fit.variogram(semivario, model_vario) # Fitting model to observed semi-variance
plot(semivario, fit_vario)


## Interpolation: Inverse Weight Distance ------------------------------------------------------

# Performing IDW
interpol_idw <- idw(formula = BD_hybrid ~ 1,
                    locations = bulk_sp,
                    newdata = grd)
output_idw <- as.data.frame(interpol_idw) # Saving output as dataframe for ggplot graphing
writeSpatialShape(interpol_idw, "Input/idw_ouput.shp") # Save interpolation result as vector shapefile

# Pre-graphing adjustments
class(interpol_idw) # Type spdf object
interpol_stars_idw <- st_as_stars(interpol_idw) # Converting spdf into as stars object 
bnd_idw <- st_transform(bnd, crs = st_crs(interpol_stars_idw)) # Change island polygon crs 
interpol_idw_crop <- interpol_stars_idw[bnd_idw] # Crop interpolation stars object to island polygon
output_idw <- as.data.frame(interpol_idw_crop) # Saving stars object as dataframe for ggplot
names(output_idw)[1:3] <- c("long","lat","bulk_density")

# Graphing IWD results
ggplot(data = output_idw, aes(x = long, y = lat)) +
  geom_tile(data = output_idw, aes(fill = bulk_density)) +
  scale_fill_continuous(low = "#D5CFD5", high = "#444244", na.value = "transparent") + 
  coord_equal() +
  theme_minimal() +
  labs(x = "\nLongitude",
       y = "Latitude\n",
       title = "Soil bulk density interpolated across \nthe island of Montreal using Inverse Distance Weighting",
       subtitle = "Projected in NAD83 MTM8") +
  theme(plot.title = element_text(color = "black", size = 16, face = "bold"),
        plot.title.position = "panel",
        axis.text.x = element_text(size = 10, face = "italic", angle = 45),
        axis.text.y = element_text(size = 10, face = "italic", angle = 45),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  guides(fill=guide_legend(title=expression(Bulk~density~(gcm^-3)))) 


## Interpolation: Ordinary Kriging --------------------------------------------------

# Creating 'gstat' formula object
ok <- gstat(formula = BD_hybrid ~ 1, model = model_vario, data = bulk_reproj)
# A formula object is created using the ~ operator, which separates names of dependent variables (to the left of the ~ symbol) and independent variables (to the right of the ~ symbol)
# Writing 1 to the right of the tilde means that there are no independent variables

# Performing interpolation
interpol_ok <- predict(ok, newdata = grd)
class(interpol_ok)
writeSpatialShape(interpol_ok, "Input/ok_ouput.shp") # Save interpolation result as vector shapefile

# Pre-graphing adjustments
class(interpol_ok) # Type spdf object
interpol_stars_ok <- st_as_stars(interpol_ok) # Converting spdf into as stars object 
bnd_ok <- st_transform(bnd, crs = st_crs(interpol_stars_ok)) # Change island polygon crs 
interpol_ok_crop <- interpol_stars_ok[bnd_ok] # Crop interpolation stars object to island polygon
output_ok <- as.data.frame(interpol_ok_crop) # Saving stars object as dataframe for ggplot
names(output_ok)[1:4] <- c("long","lat","bulk_density","variance")

# Graphing kriging results
ggplot(data = output_ok, aes(x = long, y = lat, fill = bulk_density)) +
  geom_tile(data = output_ok, aes(fill = bulk_density)) +
  scale_fill_continuous(low = "cyan", high = "orange", na.value = "transparent") + 
  coord_equal() +
  theme_minimal() +
  labs(x = "\nLongitude",
       y = "Latitude\n",
       title = "Soil bulk density interpolated across \nthe island of Montreal using Ordinary Kriging",
       subtitle = "Projected in NAD83 MTM8") +
  theme(plot.title = element_text(color = "black", size = 16, face = "bold"),
        plot.title.position = "panel",
        axis.text.x = element_text(size = 10, face = "italic", angle = 45),
        axis.text.y = element_text(size = 10, face = "italic", angle = 45),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  guides(fill=guide_legend(title=expression(Bulk~density~(gcm^-3))))

# Model evaluation
cv_ok <- gstat.cv(intpol) # Leave-One-Out cross validation, output type spdf
cv_ok <- st_as_sf(cv_ok) # Convert spdf into sf object
bubble(as(cv_ok[,"residual"], "Spatial")) # Bubble plot to examine residuals
sqrt(sum((cv_ok$var1.pred - cv_ok$observed)^2) / nrow(cv_ok))  # Calculating root mean square error (RMSE) = 0.2414772
# Interpreted as the standard deviation of the unexplained variance = variance of the residuals
# RMSE is in same units as response variable
# Lower RMSE indicates better fit


## Preparing Elevation Data ----------------------------------------------------------

# Definitions
# 1: DEM (Digital Elevation Model)
# Encapsulates both DSMs and DTMs, and can be generated from various methods
# 2: DSM (Digital Surface Model)
# An elevation model that captures both the environment’s natural and artificial features
# Includes the tops of buildings, trees, powerlines, and any other objects
# 3: DTM (Digital Terrain Model) 
# A bare-earth elevation model
# Does not contain any features above the bare-earth, even persistent ones
# Can be paired with DSMs to derive height information regarding objects on the surface

# Creating mosaic DTM, 1m resolution
dtm1 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_8_102.tif")
dtm2 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_8_103.tif")
dtm_mos1 <- st_mosaic(dtm1, dtm2)
rm(dtm1); rm(dtm2)
dtm3 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_9_103.tif")
dtm4 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_9_104.tif")
dtm_mos2 <- st_mosaic(dtm3, dtm4)
rm(dtm3); rm(dtm4)
dtm5 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_10_103.tif")
dtm6 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_10_104.tif")
dtm7 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_10_105.tif")
dtm_mos3 <- st_mosaic(dtm5, dtm6, dtm7)
rm(dtm5); rm(dtm6) ; rm(dtm7)
dtm_mos4 <- st_mosaic(dtm_mos1, dtm_mos2, dtm_mos3)
rm(dtm_mos1); rm(dtm_mos2); rm(dtm_mos3)
dtm8 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_11_103.tif")
dtm9 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_11_104.tif")
dtm_mos5 <- st_mosaic(dtm8, dtm9)
rm(dtm8); rm(dtm9)
dtm10 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_11_105.tif")
dtm11 <- read_stars("Input/raster/dtm/dtm_1m_utm18_e_11_106.tif")
dtm_mos6 <- st_mosaic(dtm10, dtm11)
rm(dtm10); rm(dtm11)
dtm_mos7 <- st_mosaic(dtm_mos5, dtm_mos6)
rm(dtm_mos5); rm(dtm_mos6)
dtm_mos <- st_mosaic(dtm_mos4, dtm_mos7)
rm(dtm_mos4); rm(dtm_mos7)
plot(dtm_mos)
write_stars(dtm_mos, "Input/raster/dtm/dtm_merged.tif")

# Creating mosaic DSM, 1m resolution
dsm1 <- read_stars("Input/raster/dsm/dsm_1m_utm18_e_8_102.tif")
dsm2 <- read_stars("Input/raster/dsm/dsm_1m_utm18_e_8_103.tif")
dsm_mos1 <- st_mosaic(dsm1, dsm2)
rm(dsm1); rm(dsm2)
dsm3 <- read_stars("Input/raster/dsm/dsm_1m_utm18_e_9_103.tif")
dsm4 <- read_stars("Input/raster/dsm/dsm_1m_utm18_e_9_104.tif")
dsm_mos2 <- st_mosaic(dsm3, dsm4)
rm(dsm3); rm(dsm4)
dsm5 <- read_stars("Input/raster/dsm/dsm_1m_utm18_e_10_103.tif")
dsm6 <- read_stars("Input/raster/dsm/dsm_1m_utm18_e_10_104.tif")
dsm7 <- read_stars("Input/raster/dsm/dsm_1m_utm18_e_11_103.tif")
dsm_mos3 <- st_mosaic(dsm5, dsm6, dsm7)
rm(dsm5); rm(dsm6) ; rm(dsm7)
dsm_mos4 <- st_mosaic(dsm_mos1, dsm_mos2, dsm_mos3)
rm(dsm_mos1); rm(dsm_mos2); rm(dsm_mos3)
dsm8 <- read_stars("Input/raster/dsm/dsm_1m_utm18_e_11_104.tif")
dsm9 <- read_stars("Input/raster/dsm/dsm_1m_utm18_e_11_105.tif")
dsm10 <- read_stars("Input/raster/dsm/dsm_1m_utm18_e_11_106.tif")
dsm_mos5 <- st_mosaic(dsm8, dsm9, dsm10)
rm(dsm8); rm(dsm9); rm(dsm10)
dsm_mos <- st_mosaic(dsm_mos4, dsm_mos5)
rm(dsm_mos4); rm(dsm_mos5)
plot(dsm_mos)
write_stars(dsm_mos, "Input/raster/dsm/dsm_merged.tif")


## Digital Terrain Model -------------------------------------------------------------

bulk_reproj <- st_transform(bulk_sf, crs = st_crs(bnd))
dtm <- read_stars("Input/raster/dtm/dtm_merged.tif")
borders <- st_read("Input/vector/island_montreal.shp")
st_bbox(dtm) # Extent of DTM
st_bbox(borders) # Extent of island
# Extents do not match, therefore will have to create new grid
grid <- st_as_sfc(st_bbox(borders)) # Setting new grid as extent of island shapefile

# Generating new grid 
grid_100 <- st_as_stars(grid, dx = 100, dy = 100) # New empty grid with cells 100m by 100m
dtm_avg_100 <- st_warp(src = dtm, grid_100, method = "average", use_gdal = TRUE) # Using average instead of nearest neighbour since appears smoother
# st_warp function: resample, transfers raster values from some original grid to a different grid
# Aligns several input rasters from different sources into same grid so they can be subject to spatial operators
# Also, reduces resolution of very detailed rasters to be more convenient to work with (in terms of memory use/processing time)
dtm_final_avg_100 <- dtm_avg_100[borders] # Cropping DEM raster to island polygon
names(dtm_final_avg_100) <- "terrain_100m"
plot(dtm_final_avg_100) 

# Joining new DtM grid to sf points spatial object
bulky <- st_join(bulk_reproj, st_as_sf(dtm_final_avg_100))
sum(is.na(bulky$terrain_100m)) # Only four data points with NA
# NAs are pixels whose centroids do no intersect with the polygon = outside
bulky <- bulky[!is.na(bulky$terrain_100m), ] # Removing NAs

# Graphing new DEM grid with points
plot(dtm_final_avg_100, breaks = "equal", col = terrain.colors(11), reset = FALSE)
plot(st_geometry(bulky), add = TRUE)


## Interpolation: Universal Kriging with DTM -------------------------------------------------------------

# Performing interpolation
uk <- gstat(formula = BD_hybrid ~ 1,
            model = model_vario,
            data = bulky)
interpol_uk <- predict(uk, dtm_final_avg_100)
write_stars(interpol_uk, dsn = "Input/uk_output.tif", layer = "var1.pred") # Saving raster as geotiff file

# Pre-graphing adjustments
class(interpol_uk) # Type stars object
bnd_uk <- st_transform(bnd, crs = st_crs(interpol_uk)) # Change island polygon crs 
interpol_uk_crop <- interpol_uk[bnd_uk] # Crop interpolation stars object to island polygon
output_uk <- as.data.frame(interpol_uk_crop) # Saving stars object as dataframe for ggplot
names(output_uk)[1:4] <- c("long","lat","bulk_density","variance")

# Graphing kriging results: ggplot
ggplot(data = output_uk, aes(x = long, y = lat)) +
  geom_tile(data = output_uk, aes(fill = bulk_density)) +
  scale_fill_continuous(low = "cyan", high = "orange", na.value = "transparent") +
  coord_equal() +
  theme_minimal() +
  labs(x = "\nLongitude",
       y = "Latitude\n",
       title = "Soil bulk density interpolated across \nthe island of Montreal using Universal Kriging",
       subtitle = "Projected in NAD83 MTM8") +
  theme(plot.title = element_text(color = "black", size = 16, face = "bold"),
        plot.title.position = "panel",
        axis.text.x = element_text(size = 10, face = "italic", angle = 45),
        axis.text.y = element_text(size = 10, face = "italic", angle = 45),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  guides(fill=guide_legend(title=expression(Bulk~density~(gcm^-3)))) 

# Model evaluation
cv_uk <- gstat.cv(uk) # Leave-One-Out cross validation, output type spdf
cv_uk <- st_as_sf(cv_uk) # Convert spdf into sf object
bubble(as(cv_uk[,"residual"], "Spatial")) # Bubble plot to examine residuals
sqrt(sum((cv_uk$var1.pred - cv_uk$observed)^2) / nrow(cv_uk)) # Calculating root mean square error (RMSE) = 0.2385589
# Interpreted as the standard deviation of the unexplained variance = variance of the residuals
# RMSE is in same units as response variable
# Lower RMSE indicates better fit 


## Digital Surface Model ---------------------------------------------------

bulk_reproj <- st_transform(bulk_sf, crs = st_crs(bnd))
dsm <- read_stars("Input/raster/dsm/dsm_merged.tif")
borders <- st_read("Input/vector/island_montreal.shp")
st_bbox(dsm) # Extent of DSM
st_bbox(borders) # Extent of island
# Extents do not match, therefore will have to create new grid
grid <- st_as_sfc(st_bbox(borders)) # Setting new grid as extent of island shapefile

# Generating new grid 
grid_100 <- st_as_stars(grid, dx = 100, dy = 100) # New empty grid with cells 100m by 100m
dsm_avg_100 <- st_warp(src = dsm, grid_100, method = "average", use_gdal = TRUE) # Using average instead of nearest neighbour since appears smoother
# st_warp function: resample, transfers raster values from some original grid to a different grid
# Aligns several input rasters from different sources into same grid so they can be subject to spatial operators
# Also, reduces resolution of very detailed rasters to be more convenient to work with (in terms of memory use/processing time)
dsm_final_avg_100 <- dsm_avg_100[borders] # Cropping DEM raster to island polygon
names(dsm_final_avg_100) <- "surface_100m"
plot(dsm_final_avg_100) 

# Joining new DSM grid to sf points spatial object
bulky <- st_join(bulk_reproj, st_as_sf(dsm_final_avg_100))
sum(is.na(bulky$surface_100m)) # Only four data points with NA
# NAs are pixels whose centroids do no intersect with the polygon = outside
bulky <- bulky[!is.na(bulky$surface_100m), ] # Removing NAs

# Graphing new DEM grid with points
plot(dsm_final_avg_100, breaks = "equal", col = terrain.colors(11), reset = FALSE)
plot(st_geometry(bulky), add = TRUE)


## Interpolation: Universal Kriging with DSM ---------------------------------------

# Performing interpolation
uk <- gstat(formula = BD_hybrid ~ 1,
            model = model_vario,
            data = bulky)
interpol_uk <- predict(uk, dsm_final_avg_100)
write_stars(interpol_uk, dsn = "Input/uk_output.tif", layer = "var1.pred") # Saving raster as geotiff file

# Pre-graphing adjustments
class(interpol_uk) # Type stars object
bnd_uk <- st_transform(bnd, crs = st_crs(interpol_uk)) # Change island polygon crs 
interpol_uk_crop <- interpol_uk[bnd_uk] # Crop interpolation stars object to island polygon
output_uk <- as.data.frame(interpol_uk_crop) # Saving stars object as dataframe for ggplot
names(output_uk)[1:4] <- c("long","lat","bulk_density","variance")

# Graphing kriging results: ggplot
ggplot(data = output_uk, aes(x = long, y = lat)) +
  geom_tile(data = output_uk, aes(fill = bulk_density)) +
  scale_fill_continuous(low = "cyan", high = "orange", na.value = "transparent") +
  coord_equal() +
  theme_minimal() +
  labs(x = "\nLongitude",
       y = "Latitude\n",
       title = "Soil bulk density interpolated across \nthe island of Montreal using Universal Kriging",
       subtitle = "Projected in NAD83 MTM8") +
  theme(plot.title = element_text(color = "black", size = 16, face = "bold"),
        plot.title.position = "panel",
        axis.text.x = element_text(size = 10, face = "italic", angle = 45),
        axis.text.y = element_text(size = 10, face = "italic", angle = 45),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  guides(fill=guide_legend(title=expression(Bulk~density~(gcm^-3)))) 

# Model evaluation
cv_uk <- gstat.cv(uk) # Leave-One-Out cross validation, output type spdf
cv_uk <- st_as_sf(cv_uk) # Convert spdf into sf object
bubble(as(cv_uk[,"residual"], "Spatial")) # Bubble plot to examine residuals
sqrt(sum((cv_uk$var1.pred - cv_uk$observed)^2) / nrow(cv_uk)) # Calculating root mean square error (RMSE) = 0.2385589
# Interpreted as the standard deviation of the unexplained variance = variance of the residuals
# RMSE is in same units as response variable
# Lower RMSE indicates better fit 

