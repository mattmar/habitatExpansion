#' Calculate Equivalent Connected Area (ECA) for a list of habitat rasters
#'
#' This function calculates the Equivalent Connected Area (ECA) based on a negative exponential kernel of dispersal probabilities for a given set of habitat rasters.
#' 
#' @param rasterList A list of rasters representing habitat distributions to be analyzed.
#' @param rasterRef A reference raster representing the initial state of habitat distribution.
#' @param directions The number of directions to be used for patch connectivity. Default is 8.
#' @param minArea The minimum area threshold for patches to be considered in the analysis. Default is 0.
#' @param alpha The alpha parameter for the negative exponential kernel. Default is 1/2.
#' @param b The b parameter for the negative exponential kernel. Default is 0.5.
#' @param d The maximum dispersal distance for the negative exponential kernel. Default is 1000.
#' @param iLandscape Logical indicating whether to calculate additional landscape metrics. Default is FALSE.
#' @param np The number of cores to use for parallel processing. Default is 1.
#' @param simplify Logical indicating whether to use a simplified approach for calculating ECA. Default is TRUE.
#' 
#' @return A list containing ECA, Hanski connectivity metric, habitat coverage, and other optional landscape metrics.
#' 
#' @examples
#' results <- connMetrics(rasterList=myRasters, rasterRef=myReferenceRaster)

connMetrics <- function(rasterList, simplify=TRUE, rasterRef, directions = 8, minArea = 0, alpha = 1/2, b = 0.5, d = 1000, iLandscape=FALSE, np=1) {
  
  # Derive the 0 p of dispersal
  lambda = -log(0.001) / d
  
  # Add the initial raster to the list of simulate matrix
  rasterList <- c(list(t(as.matrix(rasterRef,wide=TRUE))), rasterList)
  
  # Simple check if majority is NAs
  NAs <- length(which(is.na(rasterList[[1]])))/(1e3*1e3)

  # If check passed then proceed with ECA
  if( NAs<0.5 ) {
    results <- mclapply(rasterList, function(x) {

          # Reclassify matrix in habitat==1 and non-habitat==NA
          y <- terra::classify(
            x=terra::rast(t(x), crs = crs(rasterRef), ext = ext(rasterRef)), 
            rcl=matrix(c(2,1,3,NA,4,NA,5,NA,6,NA),ncol=2,byrow=T)
            )

          # Identify patches
          patches <- as.polygons(terra::patches(y, directions=8, zeroAsNA=TRUE, allowGaps=FALSE))

          # Add column for area
          patches$area <- expanse(patches, unit="m")

          # Filter patches by minimum area (minArea is in hectares)
          patchesSub <- subset(patches,  patches$area >= minArea * 10000)
          
          # Derive the number of patches
          patchN <- nrow(patchesSub)

          # Calculate distances between patches
          distM <- as.matrix(terra::distance(patchesSub, unit = "m"))
          # Calculate adjacency matrix, which the distance matrix transformed in a probability given the kernel
          adjM <- exp(-lambda*distM)

          # Make the graph from the adjacency matrix
          g <- graph_from_adjacency_matrix(adjM, mode = "undirected", weighted = TRUE)

          # Calculate ECA
          hanski <- NA
          ecaS <- calcECA(patchesSub$area, calcMaxPP(g, simplify=simplify)) / 10000

          # Calculate habitat coverage as a proportion of total landscape
          coverage <- sum(patches$area) / (ncell(x)*25)

          # Calculate SAI, clumpy and ENN
          landInd <- list(NA, NA, NA)
          if(iLandscape) landInd <- calculate_lsm(z, what= c("lsm_l_enn_mn", "lsm_l_shape_mn","lsm_c_clumpy", progress=T))

          # Return as a list
          list(hanski = hanski, eca = ecaS, coverage = coverage, patchN=patchN, na=NAs, enn=landInd[[1]], clumpy=landInd[[2]], sai=landInd[[3]])
          }, mc.cores = np)
    } else {list(hanski = NA, eca = NA, coverage = NA, na=NA, enn = NA, clumpy= NA, sai= NA)
    }
  }

#' Calculate the ECA metric
#'
#' @param areas A numeric vector of patch areas.
#' @param max_product_probs A matrix of maximum product probabilities between patches.
#' @return The ECA value.

  calcECA <- function(areas, max_product_probs) {
    eca_sum <- 0
    for (i in 1:length(areas)) {
      for (j in 1:length(areas)) {
        eca_sum <- eca_sum + areas[i] * areas[j] * max_product_probs[i, j]
      }
    }
    sqrt(eca_sum)
  }

#' Calculate the maximum product probabilities between habitat patches
#'
#' @param graph An igraph object representing habitat patches and their connectivity.
#' @param simplify Logical indicating whether to simplify the path calculations. Default is FALSE.
#' @param cutoff The cutoff parameter for path calculations. Used only if `simplify` is FALSE.
#' @return A matrix of maximum product probabilities.
  calcMaxPP <- function(graph, simplify = FALSE, cutoff=-1) {
    vcount <- vcount(graph)
    max_product_probs <- matrix(0, nrow = vcount, ncol = vcount, dimnames = list(V(graph)$name, V(graph)$name))
    components <- components(graph)

    for (comp in unique(components$membership)) {
      nodes <- which(components$membership == comp)
      for (i in nodes) {
        for (j in nodes) {
          if (i < j) {
            tryCatch({
              if (simplify) {
            # Here ECA can be "simplified" assuming the the shortest path is an approximation of the maximum product of node p among all possible paths between i and j (that otherwise would take forever to derive with all_simple_paths(link=-1)). This part uses the Dijkstra's algorithm. Where the "shortest" path is not necessarily the one with the shortest physical distance, but rather the path with the highest cumulative weight, which represents the most likely route for species dispersal.
            paths <- all_shortest_paths(graph, from = i, to = j, weights = -log(E(graph)$weight))
            if(!is.list(paths)) { paths_list <- list(paths$res) } else{ paths_list <- paths$res }
            } else {
            # Here complete ECA
            paths_list <- all_simple_paths(graph, from = i, to = j, cutoff=cutoff)
          }
          
          if (length(paths_list) > 0) {
            # This derives the probability of each path given the weights that are defined by the negative exponential distance kernel
            path_probs <- sapply(paths_list, function(path) {
              if (length(path) > 1) {
                edge_seq <- path_to_edge_seq(graph, path)
                prod(E(graph)$weight[edge_seq])
                } else {
                1  # Path to self or direct connection
              }
              })
            max_product_probs[i, j] <- max(path_probs, na.rm = TRUE)
            max_product_probs[j, i] <- max_product_probs[i, j]  # Fill symmetric entry
            } else {
              message("No valid paths found from ", i, " to ", j)
            }
            }, error = function(e) {
              message("Error in path calculation from ", i, " to ", j, ": ", e$message)
              })
            } else if (i == j) {
        max_product_probs[i, j] <- 1  # Self-connection
      }
    }
  }
}

return(max_product_probs)
}

# Helper function to convert a path to a sequence of edge IDs
path_to_edge_seq <- function(graph, path) {
  sapply(1:(length(path) - 1), function(k) {
    get.edge.ids(graph, c(path[k], path[k + 1]))
    })
}
