# Simulation of Habitat Growth and Connectivity Metrics Calculation (for each growing step) in a Landscape

# Load libraries
library(terra)
library(parallel)
library(ggplot2)
library(reshape2)
library(scales)
library(stats)
library(landscapemetrics)
library(igraph)

# Load custom growing functions
sapply(list.files("./growing_functions", full.names = TRUE), source, .GlobalEnv)

# Load landscape raster
landscape <- rast("./TestMatrix.tif")

# Define fixed parameters for growth
params <- list(
  swfCover = 0.30,  # Target SWF coverage
  iterations = 500,  # Maximum number of iterations
  swfCat = 1,       # Cell value representing SWF
  agriCat = 4,      # Cell value representing crop fields
  boundaryCat = 3,  # Cell value representing crop boundaries
  ExpPriority = "mixed", # preferred directional priority
  ExpDirection = "mixed", # preferred direction 
  NNeighbors = 1, # Number of agri neighbours a SWF cells must have to expand
  maxDistance = 2, # Max distance (in pixel) from each SWF cell to look for agri neighbours
  queensCase = TRUE,
  np = detectCores() - 1,  # Parallel processes
  deBug = FALSE,
  gd = 3,  # Maximum search radius for growing SWF
  qn = 1,  # Number of agri cells to expand for each SWF cell
  rq = 0.50  # Proportion of SWF cells used for expansion at each iteration
)

# Define varying parameters for growth
weights <- c(0.05, 0.5, 0.99)  # Weights of boundary respect to interior (probability in a binomial random draw)
params.mt <- expand.grid(W = weights)

# Prepare the initial landscape matrix
landscapeMatrix <- t(as.matrix(landscape, wide = TRUE))

# Run growth simulation for each scenario
results.lt <- mclapply(seq_len(nrow(params.mt)), function(idx) {
  w <- params.mt[idx, "W"]
  message("Boundary Weights: ", w, "\nagri Weights: ", 1 - w, "\n")
  swf.molder(
  	Hmatrix = landscapeMatrix, 
		swfCover = params$swfCover, 
		swfCat = params$swfCat, 
		agriCat = params$agriCat, 
		boundaryCat = params$boundaryCat, 
		agriW=1-w, # This is the only varying parameter 
		boundaryW=w, # This is the only varying parameter
		Q=params$qn, 
		ExpPriority = params$ExpPriority, 
		ExpDirection = params$ExpDirection, 
		reduceQTo=params$rq, 
		iterations = params$iterations, 
		NNeighbors = params$NNeighbors, 
		maxDistance = params$maxDistance, 
		queensCase = params$queensCase, 
		maxGDistance= params$gd, 
		np = params$np)
}, mc.cores = 2, mc.preschedule = FALSE)

# Visualize a selected scenario
scenario <- 3  # Example with W=0.99
nrows <- 8
ncols <- ceiling(length(results.lt[[scenario]]) / nrows)

par(mfrow = c(nrows, ncols))
lapply(results.lt[[scenario]], function(i) {
  plot(rast(i))
})
dev.off()

# Calculate connectivity metrics (ECA) for each simulated matrix
source("./connectivity_functions/connectivityMetricsV2.0.R")
results.r <- mclapply(results.lt, connMetrics, simplify = TRUE, rasterRef = landscape, directions = 4, 
                      minArea = 0.0026, d = 10, iLandscape = FALSE, np = 4, mc.cores = 2)

# Process and visualize ECA results
resultsECA.df <- do.call(rbind.data.frame,
	lapply( 1:length(results.r), function(x) {
		y <- do.call(rbind.data.frame, results.r[[x]])
		y$scenario <- x
		y$iCov <- y$coverage[1]
		return(y)
	}
	)
	)
ggplot(resultsECA.df, aes(x = coverage, y = eca, group = scenario)) + 
  geom_line(aes(color=as.factor(scenario)))