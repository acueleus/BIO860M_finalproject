## Packages 
install.packages(c("tidyverse","ggplot2","raster","sp","rgdal","sf","stars","gstat"))
my_packages <- c("tidyverse","ggplot2","raster","sp","rgdal","sf","stars","gstat")
lapply(my_packages, require, character.only = TRUE)


## Making objects into spatial ----------------------------------------------------------------

# Loading objects
metal <- read.csv("Input/TE_data_georef.csv")
bnd <- st_read("Input/vector/island_montreal.shp")

# Converting dataframes into spatial objects with 'sf' package
metal_sf <- st_as_sf(metal, coords = c("coord_X","coord_Y"), crs = 4326) 
metal_reproj <- st_transform(metal_sf, crs = 2950)
metal_sp <- as_Spatial(metal_reproj) 

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
summary(grd) 
proj4string(grd) <- CRS("+init=epsg:2950")
summary(grd)
gridded(grd) <- TRUE
plot(grd, cex = 1.5)
points(bulk_sp, pch = 1, col = "red", cex = 1)

# Exploring variograms for kriging: lead
variocloud <- variogram(Pb_conc ~ 1, data = metal_reproj, cloud = TRUE)
plot(variocloud) 
semivario <- variogram(Pb_conc ~ 1, data = metal_reproj)
plot(semivario)
semivario

# Trying to fit best models onto variogram: lead
v <- variogram(Pb_conc ~ 1, data = metal_reproj)
vgm() # To check list of model functions
v_fit <- fit.variogram(v, vgm("Exp"))
v_fit
v_fit <- fit.variogram(v, vgm("Sph"))
v_fit
v_fit <- fit.variogram(v, vgm("Lin"))
v_fit
v_fit <- fit.variogram(v, vgm("Log"))
v_fit
model_vario <- vgm(psill = 11982.11,
                   model = "Log",
                   nugget = 10613.57,
                   range = 6353.064)
fit_vario <- fit.variogram(semivario, model_vario) # WARNING "singular model in variogram fit"
plot(semivario, fit_vario)

install.packages("automap")
library(automap)
model_vario_auto <- autofitVariogram(formula = Pb_conc ~ 1,
                                     input_data = metal_reproj,
                                     model = c("Exp","Log","Lin","Sph"),
                                     verbose = TRUE)
# Autofitted output:
# Selected:
#  model    psill    range
#1   Nug 6740.785    0.000
#2   Log 8371.972 4430.062
# Tested models, best first:
#  Tested.models kappa   SSerror
#1           Log     0 133038938
plot(semivario, model_vario_auto$var_model)


## Interpolation: Ordinary Kriging --------------------------------------------------

# Creating 'gstat' formula object: lead

ok <- gstat(formula = Pb_conc ~ 1, model = model_vario, data = metal_reproj) 

# Performing interpolation
interpol_ok <- predict(ok, newdata = grd)
# ERROR "In predict.gstat(ok, newdata = grd) : Covariance matrix singular at location [268987,5.02923e+06,0]: skipping"
# To address error message, would have to check the following:
# There is actually a spatial structure in dataset (via bubble plot)
# There are no duplicate locations
# The variogram model is not singular and has a good fit to the experimental variogram
# Try several values of range, sill, nugget and all the models in the gstat library
# The covariance matrix is positive definite and has positive eigen values = is it singular according to gstat and/or is.singular.matrix function
# There were enough pair of points to do the experimental variogram

# Pre-graphing adjustments
#class(interpol_ok) # Type spdf object
#interpol_stars_ok <- st_as_stars(interpol_ok) # Converting spdf into as stars object 
#bnd_ok <- st_transform(bnd, crs = st_crs(interpol_stars_ok)) # Change island polygon crs 
#interpol_ok_crop <- interpol_stars_ok[bnd_ok] # Crop interpolation stars object to island polygon
#output_ok <- as.data.frame(interpol_ok_crop) # Saving stars object as dataframe for ggplot
#names(output_ok)[1:3] <- c("long","lat","bulk_density")

# Graphing kriging results
#ggplot(data = output_ok, aes(x = long, y = lat, fill = bulk_density)) +
  #geom_tile(data = output_ok, aes(fill = bulk_density)) +
  #scale_fill_continuous(low = "green", high = "red", na.value = "transparent") + 
  #coord_equal() +
  #theme_minimal() +
  #labs(x = "\nLongitude",
       #y = "Latitude\n",
       #title = "Soil lead concentration interpolated across \nthe island of Montreal using Ordinary Kriging",
       #subtitle = "Projected in NAD83 MTM8") +
  #theme(plot.title = element_text(color = "black", size = 16, face = "bold"),
        #plot.title.position = "panel",
        #axis.text.x = element_text(size = 10, face = "italic", angle = 45),
        #axis.text.y = element_text(size = 10, face = "italic", angle = 45),
        #axis.title.x = element_text(size = 14),
        #axis.title.y = element_text(size = 14)) +
  #guides(fill=guide_legend(title=expression(Bulk~density~(gcm^-3))))

# Model evaluation
#cv_ok <- gstat.cv(intpol) # Leave-One-Out cross validation, output type spdf
#cv_ok <- st_as_sf(cv_ok) # Convert spdf into sf object
#bubble(as(cv_ok[,"residual"], "Spatial")) # Bubble plot to examine residuals
#sqrt(sum((cv_ok$var1.pred - cv_ok$observed)^2) / nrow(cv_ok))  # Calculating root mean square error (RMSE) = 0.2414772
# Interpreted as the standard deviation of the unexplained variance = variance of the residuals
# RMSE is in same units as response variable
# Lower RMSE indicates better fit

