## Comparing Spatially Balanced Sampling
## Impact of reduced sampling effort on inclusion probabilities.

## Paul van Dam-Bates
## RE: TESA Design Workshop

## Basic idea is to simulate a way the samples are generated and to check that you meet the
## target inclusion probabilities.
## Can put in more complicated scenarios for how samples are actually achieved and see that impact
## for how the samples meet the inclusion probs. E.g. for restricting nearby points in random sample.
## Strong differences in incl. probs -> biased Horvitz Thompson Estimator if not accounted for.

library(sf)
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

## Data from the Spbsamping package.
data("lucas_abruzzo", package = "Spbsampling")
grid <- lucas_abruzzo
grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = "EPSG:6875")
dis_sf <- st_distance(grid_sf)
lookup <- lookuptable(grid)

poly_sf <- grid %>% 
  st_as_sf(coords = c("x","y")) %>% 
  st_buffer(dist = lookup$res[1]/2, endCapStyle = "SQUARE") %>% 
  st_union() %>%
  st_as_sf()

## Setup sampling needed for pwd from Spbsampling.
con <- rep(0, nrow(dis_sf))
stand_dist <- stprod(mat = dis_sf, con = con)$mat

## Sample design
nsmp <- 50  # Number of samples
nmiss <- 10 # Reduced effort
Ngrd <- nrow(grid) # Population size
iters <- 1000 # Simulations

## Simulate inclusion probabilities for Spbsampling:
## Choose beta for strength of spatial spread.
smp_incl_bs <- numeric(Ngrd)
pb = txtProgressBar(min = 0, max = iters, initial = 0)
for( i in 1:iters ){
  # smp.bs <- pwd(dis = stand_dist, n = nsmp)$s[1,1:(nsmp-nmiss)] ## Missing last 10. Not random.
  smp.bs <- pwd(dis = stand_dist, n = nsmp, beta = 10)$s[1,sample(nsmp, nsmp-nmiss, replace = FALSE)] ## Missing last 10 at random?
  smp_incl_bs[smp.bs] <- smp_incl_bs[smp.bs] + 1
  setTxtProgressBar(pb,i)
}
close(pb)

## Simulate inclusion probabilities for Spbsampling:
## Choose beta for strength of spatial spread.
smp_incl_lpm <- numeric(Ngrd)
pb = txtProgressBar(min = 0, max = iters, initial = 0)
for( i in 1:iters ){
  smp.lpm <- lpm(prob = rep(nsmp/Ngrd, Ngrd), x=grid[,c('x','y')])
  smp.lpm <- smp.lpm[sample(nsmp, nsmp-nmiss, replace = FALSE)] ## Missing last 10 at random?
  smp_incl_lpm[smp.lpm] <- smp_incl_lpm[smp.lpm] + 1
  setTxtProgressBar(pb,i)
}
close(pb)

## Simulating inclusion probs from Simple Random Sample.
smp_incl_srs <- numeric(Ngrd)
pb = txtProgressBar(min = 0, max = iters, initial = 0)
for( i in 1:iters ){
  smp <- sample(Ngrd, nsmp)[1:(nsmp-nmiss)]
  smp_incl_srs[smp] <- smp_incl_srs[smp] + 1
  setTxtProgressBar(pb,i)
}
close(pb)

## Simulating inclusion probs from Balanced Acceptance Sampling.
## Sites were dropped in hierarchical order
smp_incl_bas <- numeric(Ngrd)
pb = txtProgressBar(min = 0, max = iters, initial = 0)
for( i in 1:iters ){
  smp <- BAS(poly_sf, n = nsmp)$sample
  ind <- smp %>% st_coordinates() %>% matchpts(lookup)
  ind <- ind[1:(nsmp-nmiss)]  # Removing sites in hierarchical order
  smp_incl_bas[ind] <- smp_incl_bas[ind] + 1
  setTxtProgressBar(pb,i)
}
close(pb)

## Simulating inclusion probs from GRTS
## Sites were dropped in hierarchical order
smp_incl_grts <- numeric(Ngrd)
pb = txtProgressBar(min = 0, max = iters, initial = 0)
for( i in 1:iters ){
  smp <- grts(sframe = grid_sf, n_base = nsmp)$sites_base$id[1:(nsmp-nmiss)] # Removing sites in hierarchical order
  smp_incl_grts[smp] <- smp_incl_grts[smp] + 1
  setTxtProgressBar(pb,i)
}
close(pb)

## Slowest method is PWD and it also is the least flexible that does not appear to meet
## the required incl probs.
inclprobs <- cbind(smp_incl_srs, smp_incl_bs, smp_incl_lpm, smp_incl_bas, smp_incl_grts)/iters
colnames(inclprobs) <- c("SRS", "PWD", "LPM", "BAS", "GRTS")

## Population should be close to equal inclusion probabilities centred around the
## target inclusion prob with some montecarlo error (increase iters to narrow that).
## Error is higher with PWD, not clear it actually reaches target inclusion probabilities.
## Also very slow to process all distances.
boxplot(inclprobs, ylab = "Inclusion Probability", xlab = "Method")
abline(h = (nsmp-nmiss)/Ngrd, col = 'red')


library(ggplot2)
## Baseline use SRS and compare to this:
ggplot(data = grid, aes(x = x, y = y)) + geom_tile(aes(fill = smp_incl_srs))

## Clear patterns of areas not sampled. Issue!!!!!
ggplot(data = grid, aes(x = x, y = y)) + geom_tile(aes(fill = smp_incl_bs))

## Nice random monte carlo error in incl. prob.
ggplot(data = grid, aes(x = x, y = y)) + geom_tile(aes(fill = smp_incl_lpm))

## Nice random monte carlo error in incl. prob.
ggplot(data = grid, aes(x = x, y = y)) + geom_tile(aes(fill = smp_incl_bas))

## Nice random monte carlo error in incl. prob.
ggplot(data = grid, aes(x = x, y = y)) + geom_tile(aes(fill = smp_incl_grts))
