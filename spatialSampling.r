## Simulate spatial designs for different polygons:

library(sf)
library(terra)
library(sdmTMB)
library(tidyverse)
library(ggeffects)

library(spbal) ## BAS
library(Spbsampling) ## Method used in Vilas et al. (I think for 'balanced sampling')
library(spsurvey) ## GRTS
library(BalancedSampling) ## LPM

## Use this to look up points on the grids:
lookuptable <- function(grid, xy = c("x", "y")){
  minxy <- c(min(grid[, xy[1]]), min(grid[, xy[2]]))
  ydiff <- grid[, xy[2]] |> diff() |> abs() |> round(7)
  yres <- min(ydiff[ydiff != 0])
  xdiff <- grid[, xy[2]] |> diff() |> abs() |> round(7)
  xres <- min(xdiff[xdiff != 0])
  minxy <- minxy - c(xres, yres)/2
  
  idx_state <- floor((grid[, xy[1]] - minxy[1]) / xres) + 1
  idy_state <- floor((grid[, xy[2]] - minxy[2]) / yres) + 1
  mat <- matrix(NA, nrow = nrow(grid), ncol = nrow(grid))
  mat[cbind(idx_state, idy_state)] <- 1:nrow(grid)
  return(list(lookup = mat, res = c(xres, yres), minxy = minxy))
}

## Matching sampled points to the grid.
matchpts <- function(xy, lookupinfo, dat){
  idx <- floor((xy[,1] - lookupinfo$minxy[1]) / lookupinfo$res[1]) + 1
  idy <- floor((xy[,2] - lookupinfo$minxy[2]) / lookupinfo$res[2]) + 1
  ids <- lookupinfo$lookup[cbind(idx,idy)]
  return(ids)
}

## Pacific Code from sdmTMB:
grid <- qcs_grid
grid$siteid <- 1:nrow(grid)
grid_sf <- st_as_sf(grid, coords = c("X", "Y"), crs = "EPSG:4269")
lookup <- lookuptable(grid, c("X", "Y"))

poly_sf <- grid %>% 
  st_as_sf(coords = c("X","Y")) %>% 
  st_buffer(dist = lookup$res[1]/2, endCapStyle = "SQUARE") %>% 
  st_union() %>%
  st_as_sf()

## Sample design
nsmp <- 50  # Number of samples
Ngrd <- nrow(grid) # Population size
iters <- 1000 # Simulations

## Create a surface to practice sampling:
# use a small number of knots for this example to make it fast:
mesh <- make_mesh(pcod, c("X", "Y"), n_knots = 60)

# fit a spatiotemporal model:
m <- sdmTMB(
 data = pcod,
 formula = density ~ 0 + as.factor(year) + s(depth, k = 5),
 time = "year", mesh = mesh, family = tweedie(link = "log"), spatial = "on"
)

# prepare a prediction grid:
nd <- replicate_df(grid, "year", unique(pcod$year))

# Note `return_tmb_object = TRUE` and the prediction grid:
predictions <- predict(m, newdata = nd, return_tmb_object = TRUE)

pred <- predictions$data %>% 
  filter(year == 2011) %>% 
  mutate(Estimate = exp(est))
  
pred %>% ggplot(aes(X, Y, fill = Estimate)) + 
  geom_raster() +
  scale_fill_viridis_c(trans = "sqrt") + 
  theme_bw() +
  coord_fixed()

grid <- grid %>% left_join(pred)

## Total Biomass is of interest:
## Note that sdmTMB has a different way of summing for log links etc via index.
total_biomass <- sum(grid$Estimate)

## Now try some different sampling schemes:
pincl <- nsmp/Ngrd
smp_ests <- data.frame()
for( i in 1:iters ){
  ## SRS
  pts.srs <- poly_sf %>% st_sample(nsmp)
  indx <- matchpts(pts.srs %>% st_coordinates(), lookup)
  smp_ests <- rbind(smp_ests, data.frame(method = "SRS", Estimate = sum(grid[indx, "Estimate"]/pincl)))
  ## Systematic
  pts.syst <- poly_sf %>% st_sample(nsmp, type = "regular")
  ## Just in case target sample size isn't reached.
  pincl_syst <- length(pts.syst)/Ngrd
  indx <- matchpts(pts.syst %>% st_coordinates(), lookup)
  smp_ests <- rbind(smp_ests, data.frame(method = "SYST", Estimate = sum(grid[indx,"Estimate"]/pincl_syst)))  
  ## BAS
  pts.bas <- BAS(poly_sf, n = nsmp)$sample
  indx <- matchpts(pts.bas %>% st_coordinates(), lookup)
  smp_ests <- rbind(smp_ests, data.frame(method = "BAS", Estimate = sum(grid[indx,"Estimate"]/pincl)))
  ## Could add GRTS or LPM here too. 
}

## Plot the estimate total biomass with uncertainty from design:
## Note that systematic might be the best but might be sampling slightly more
## than nsmp to get a full grid.
boxplot(data = smp_ests, Estimate~method)
abline(h = total_biomass, col = 'red')

## Show some stratification with depth:
grid <- grid %>% mutate(depth_grp = cut(depth, c(0, 100, 200, 300, 400, Inf)))
grid <- grid %>% group_by(depth_grp) %>% mutate(ngrp = n()) %>% ungroup()

poly_depth <- grid %>%
  st_as_sf(coords = c("X","Y")) %>% 
  st_buffer(dist = lookup$res[1]/2, endCapStyle = "SQUARE")  %>%
  group_by(depth_grp) %>% 
  summarize(geometry = st_union(geometry)) %>%
  st_as_sf()

## Equal Effort Stratification:
tab <- table(grid$depth_grp)
nlevels <- length(tab)
nstrat <- rep( floor(nsmp/nlevels), nlevels )
## Add in extra samples if nsmp %% nlevels != 0
nextra <- rmultinom(1, size = nsmp - sum(nstrat), prob = tab/sum(tab))
nstrat <- nstrat + nextra[,1]
names(nstrat) <- names(tab)

## Plot example:
pts.bas <- BAS(poly_depth, n = nstrat, stratum = "depth_grp")$sample
poly_depth %>% ggplot() + 
  geom_sf(aes(fill = depth_grp)) + 
  theme_bw()  + 
  geom_sf(data=pts.bas)

## Loop through drawing a BAS stratified sample:
for( i in 1:iters ){
  pts.bas <- BAS(poly_depth, n = nstrat, stratum = "depth_grp")$sample
  pts.bas$ndepth <- nstrat[pts.bas$depth_grp]
  indx <- matchpts(pts.bas %>% st_coordinates(), lookup)
  smp <- grid[indx, ]
  pincl <- pts.bas$ndepth/smp$ngrp
  smp_ests <- rbind(smp_ests, data.frame(method = "BAS-Depth", Estimate = sum(smp[,"Estimate"]/pincl)))
}

## Plot the estimate total biomass with uncertainty from design:
boxplot(data = smp_ests, Estimate~method)
abline(h = total_biomass, col = 'red')
## Yikes that made things worse!
