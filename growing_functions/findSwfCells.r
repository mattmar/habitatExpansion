# Helper functions
#' @export
findSwfCells <- function(matrix, swfCat, agriCat, boundaryCat, agriW, boundaryW, NNeighbors, queensCase, maxDistance = 1 , maxGDistance, Q = 1, reduceQTo, ExpPriority, ExpDirection, np = 1) {

    nrows <- nrow(matrix)
    ncols <- ncol(matrix)
    directions <- if (queensCase) { # Generate directions for neighbor checking
        expand.grid(row = -maxDistance:maxDistance, col = -maxDistance:maxDistance)
        } else {
            rbind.data.frame(cbind(-1:1, 0), cbind(0, -1:1))
        }

        directions <- directions[!apply(directions, 1, function(x) all(x == 0)), ]
        names(directions) <-c("row","col")
        
        indices <- expand.grid(row = 1:nrows, col = 1:ncols)
        indices <- indices[which(as.vector(matrix) == swfCat),]

        result <- parallel::mclapply(1:nrow(indices), function(i) {

            neighborsDF <- findAndSelectNeighbors(indices[i, ], matrix, swfCat, agriCat, boundaryCat, agriW, boundaryW, ExpPriority, ExpDirection, NNeighbors, directions, Q, maxGDistance)
            if( !is.null(dim(neighborsDF)) ) neighborsDF$index <- i
            return(neighborsDF)
            }, mc.cores=np, mc.preschedule=TRUE)

        allocatedCells <- result[which(sapply(result, function(x) !is.null(dim(x))))]

    if( length(allocatedCells)>0 && reduceQTo!=0 ) { # Subsample allocated cells
        sampledAllCells <- sample(1:length(allocatedCells),ceiling(length(allocatedCells)*reduceQTo))
        if(length(sampledAllCells)>0) {
            allocatedCells <- allocatedCells[sampledAllCells]
        }
    }

    allocatedCells <- do.call(rbind, allocatedCells)

    return(allocatedCells[,1:2])
}

findAndSelectNeighbors <- function(index, matrix, swfCat, agriCat, boundaryCat, agriW, boundaryW, ExpPriority, ExpDirection, NNeighbors, directions, Q, maxGDistance) { # Function to find and select closest agriCat neighbors
    r <- index$row
    c <- index$col
    
    if ( countAgriNeighbors(r, c, matrix, agriCat, boundaryCat, directions) >= NNeighbors ) {
        AgriCells <- findClosestAgriCells(matrix, c(r, c), agriCat, boundaryCat, agriW, boundaryW, Q, maxGDistance, ExpPriority, ExpDirection)
        return(AgriCells)
    }
    return(list())
}

# Helper function to count the number of neighbour boundaryCat for each SWF cat
countAgriNeighbors <- function(r, c, matrix, agriCat, boundaryCat, directions) {
    sum(sapply(1:nrow(directions), function(dir) {
        i <- r + directions[dir, "row"]
        j <- c + directions[dir, "col"]
        i >= 1 && i <= nrow(matrix) && j >= 1 && j <= ncol(matrix) && matrix[i, j] %in% c(boundaryCat)
        }))
}
