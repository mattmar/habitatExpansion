# # Find the Q closest neighbourns
# test <- matrix(c(2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1), nrow=4)
# gravity.pos <- findSwfCells(test,swfCat=2,agriCat=1,AgriNeighbors=4,TRUE,maxDistance=1)
# agri.pos <- as.data.frame(which(test==1,arr.ind=T))
# do.call(rbind.data.frame, ClosestDiagonalAgriCell(test, agri.pos, gravity.pos, 1, 2, 2, 2, np=1))

# ClosestDiagonalAgriCell(test, agri.pos, gravity.pos, 1, 2, 2, np=1)
#' @export
findClosestAgriCells <- function(matrix, targetPos, agriCat, boundaryCat, agriW, boundaryW, Q, maxGDistance, ExpPriority, ExpDirection) {
    # targetPos =c(81,8)
    nrows <- nrow(matrix)
    ncols <- ncol(matrix)
    r <- targetPos[1]
    c <- targetPos[2]
    maxRadius <- ifelse(is.numeric(maxGDistance), maxGDistance, max(nrows, ncols))

    # Precompute the logical matrix for agriCat cells
    agriCatMatrix <- matrix == agriCat | matrix == boundaryCat

    # Initialize variables
    selectedNeighbors <- data.frame(row = integer(), col = integer())
    radius <- 1

    while(nrow(selectedNeighbors) < Q && radius <= maxRadius) {
        rowRange <- max(1, r-radius):min(nrows, r+radius)
        colRange <- max(1, c-radius):min(ncols, c+radius)
        
        # Create a submatrix for the current radius
        subMatrix <- agriCatMatrix[rowRange, colRange]
        
        # Find row and column indices of agriCat cells in the submatrix
        trueElements <- which(subMatrix, arr.ind = TRUE)
        agriRows <- rowRange[trueElements[, "row"]]
        agriCols <- colRange[trueElements[, "col"]]

        # Update selectedNeighbors
        newNeighbors <- data.frame(row = agriRows, col = agriCols, stringsAsFactors = FALSE)

        radius <- radius + 1
    }

    # Order and select closest neighbors
    orderedNeighbors <- orderCellsByDistance(newNeighbors, r, c, ExpPriority, ExpDirection)
    
    if( !is.na(boundaryW) ) {
        prob.lulc <- ifelse(matrix[as.matrix(orderedNeighbors)]==agriCat, agriW, boundaryW)!=0

        if( any(prob.lulc) ) {
        # N_agri <- length(which(matrix[as.matrix(orderedNeighbors)] == agriCat)) +0.1
        # N_boundary <- length(which(matrix[as.matrix(orderedNeighbors)] == boundaryCat)) +0.1

        # # Calculate observed proportions
        # P_agri <- N_agri / (N_agri + N_boundary)
        # P_boundary <- N_boundary / (N_agri + N_boundary)

        # # Adjust weights based on the inverse of observed proportions
        # W_prime_agri <- agriW * (1 / P_agri)
        # W_prime_boundary <- boundaryW * (1 / P_boundary)

        # # Normalize adjusted weights to ensure they sum to 1 (or maintain their original sum)
        # sum_W_prime <- W_prime_agri + W_prime_boundary
        # W_agri <- W_prime_agri / sum_W_prime
        # W_boundary <- W_prime_boundary / sum_W_prime

            chosenNeighbors <- orderedNeighbors[sample.int(
                    n=nrow(orderedNeighbors),
                    size=min(Q,length(which(prob.lulc))),
                    prob=ifelse(matrix[as.matrix(orderedNeighbors)]==agriCat, agriW, boundaryW)),]
            } else { 
                chosenNeighbors <- c()
                } 
            } else {
                chosenNeighbors <- orderedNeighbors[1:min(nrow(orderedNeighbors), Q), ]
            }
            return(chosenNeighbors)
        }
