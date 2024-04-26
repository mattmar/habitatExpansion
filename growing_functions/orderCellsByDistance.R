orderCellsByDistance <- function(cellList, targetRow, targetCol, ExpPriority, ExpDirection) {
    if (nrow(cellList) == 0) {
        return(cellList)
    }

    # Calculate Manhattan distance for each cell from the target
    cellList$manhattanDist <- abs(cellList$row - targetRow) + abs(cellList$col - targetCol)
    
    # Identify diagonal, vertical, and horizontal cells
    cellList$isDiagonal <- abs(cellList$row - targetRow) == abs(cellList$col - targetCol)
    cellList$isVertical <- cellList$row != targetRow & cellList$col == targetCol
    cellList$isHorizontal <- cellList$row == targetRow & cellList$col != targetCol

    # Assign a random tiebreaker score
    cellList$tiebreaker <- runif(nrow(cellList))
    if (ExpPriority == "mixed") {
        if (ExpDirection == "orthogonal") {
            cellList$tiebreaker[cellList$isDiagonal] <- 0
        } else if (ExpDirection == "diagonal") {
            cellList$tiebreaker[!cellList$isDiagonal] <- 0
        }
    } else if (ExpPriority == "vertical") {
        cellList$tiebreaker[!cellList$isVertical] <- 0
    } else if (ExpPriority == "horizontal") {
        cellList$tiebreaker[!cellList$isHorizontal] <- 0
    }

    # Decide the ordering direction
    if (ExpDirection == "mixed") {
        ExpDirection <- sample(c("orthogonal", "diagonal", "vertical", "horizontal"), 1)
    }

    # Order by Manhattan distance, diagonal/vertical/horizontal, and tiebreaker
    if (ExpDirection == "diagonal") {
        cellList <- cellList[order(cellList$manhattanDist, -cellList$isDiagonal, cellList$tiebreaker), 1:2]
    } else if (ExpDirection == "orthogonal") {
        cellList <- cellList[order(cellList$manhattanDist, cellList$isDiagonal, cellList$tiebreaker), 1:2]
    } else if (ExpDirection == "vertical") {
        cellList <- cellList[order(cellList$manhattanDist, -cellList$isVertical, cellList$tiebreaker), 1:2]
    } else if (ExpDirection == "horizontal") {
        cellList <- cellList[order(cellList$manhattanDist, -cellList$isHorizontal, cellList$tiebreaker), 1:2]
    }

    return(cellList)
}
