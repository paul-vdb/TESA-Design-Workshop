## Example for how to use BAS for a changing spatial extent:
## See: Using balanced acceptance sampling as a master sample for environmental surveys
## https://doi.org/10.1111/2041-210X.13003

library(sf)
library(terra)
library(sdmTMB)
library(tidyverse)
library(ggeffects)
library(spbal) ## BAS

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

## Define full spatial extent:
## Using pacific cod as an example:

grid <- qcs_grid
grid$siteid <- 1:nrow(grid)
grid_sf <- st_as_sf(grid, coords = c("X", "Y"), crs = "EPSG:4269")
lookup <- lookuptable(grid, c("X", "Y"))

poly_sf <- grid %>% 
  st_as_sf(coords = c("X","Y")) %>% 
  st_buffer(dist = lookup$res[1]/2, endCapStyle = "SQUARE") %>% 
  st_union() %>%
  st_as_sf()

## Full extent:
plot(poly_sf)

## Define the bounding box and random seed:
## box and seed are automatically returned when running BAS but here we want to be explicit.
bb_design <- BoundingBox(poly_sf, d = 2)
seed <- spbal:::findBASSeed(poly_sf, bb_design, n = 1)

## Now every sample is deterministic and the same, but you can ask for more or less whenever...
## Alternatively, you could ask for n=10000 and just use the order of the points and manually choose in and out.
smp1 <- BAS(poly_sf, 100, boundingbox = bb_design, seeds = seed)
smp2 <- BAS(poly_sf, 50, boundingbox = bb_design, seeds = seed)

poly_sf %>% ggplot() + 
  geom_sf() + 
  geom_sf(data = smp1$sample, col = 'red', shape = 3) + 
  geom_sf(data = smp2$sample, col = 'blue', shape = 2) +
  theme_bw()

## Now let's sample a smaller area, defined by shallow regions:
poly_sf_sub <- grid %>% 
  filter(depth < 200) %>%
  st_as_sf(coords = c("X","Y")) %>% 
  st_buffer(dist = lookup$res[1]/2, endCapStyle = "SQUARE") %>% 
  st_union() %>%
  st_as_sf()

poly_sf %>% ggplot() + 
  geom_sf() + 
  geom_sf(data = poly_sf_sub, fill = "blue") + 
  theme_bw()

## Sample designed for < 200 m at first:
smp_sub <- BAS(poly_sf_sub, n = 70, boundingbox = bb_design, seeds = seed)

poly_sf_sub %>% ggplot() + 
  geom_sf() + 
  theme_bw() + 
  geom_sf(data = smp_sub$sample)

## Now I will add the full extent. I can afford to do a total of 80 samples.
smp_full <- BAS(poly_sf, n = 80, boundingbox = bb_design, seeds = seed)

## These samples are no longer part of the survey.
samples_dropped <- smp_sub$sample %>% filter(!SiteID %in% smp_full$sample$SiteID)

## Red cross are the dropped points.
poly_sf %>% ggplot() + 
  geom_sf() + 
  theme_bw() + 
  geom_sf(data = smp_full$sample) + 
  geom_sf(data = smp_sub$sample, col = 'blue', shape = 17, size = 2) + 
  geom_sf(data = poly_sf_sub, fill = "blue", alpha = 0.1) +
  geom_sf(data = samples_dropped, shape = 3, col = 'red', size = 2) 
