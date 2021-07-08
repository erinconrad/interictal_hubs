library(spdep)
library(sp)
library(raster)
library(tidyverse)

dat <- read_csv('/Users/erinconrad/Desktop/research/interictal_hubs/results/for_r/r_HUP099.csv')

coords <- as.matrix(cbind(dat$x, dat$y, dat$z) )

ac <- autocov_dist(as.numeric(dat$abs), coords, nbs = 0.1, longlat = TRUE)

