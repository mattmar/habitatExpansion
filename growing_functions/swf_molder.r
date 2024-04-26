	#' SWF Molder Function
	#'
	#' This function modifies a habitat matrix to increase the cover of a specified category
	#' (swfCat) by clumping it within a defined kernel size, considering neighbor preferences
	#' and optionally reducing the selection by a given factor.
	#'
	#' @param Hmatrix A matrix representing the initial habitat state.
	#' @param swfCover The desired cover proportion for the swfCat category.
	#' @param swfCat The category value in the matrix to be increased.
	#' @param agriCat The category value in the matrix considered as non-habitat.
	#' @param Q The number of cells to be potentially modified in each iteration.
	#' @param ExpDir Function to determine the prioritization of diagonal neighbors.
	#' @param reduceQTo A factor by which to reduce the selection of cells during processing.
	#' @param iterations The number of iterations to run the modification process.
	#' @param NNeighbors The number of neighbor cells to consider for potential habitat clumping.
	#' @param maxDistance The maximum distance to look for neighbor cells.
	#' @param queensCase Logical; if TRUE, considers all 8 directions for neighbors; if FALSE, only orthogonal.
	#' @param np Number of cores for parallel processing.
	#' @param deBug Logical; if TRUE, prints debugging information during processing.
	#'
	#' @return A list of matrices representing the state of the habitat matrix after each iteration.
	#' @export
	swf.molder <- function(Hmatrix, swfCover=0.10, swfCat, agriCat, boundaryCat, agriW=NA, boundaryW=NA, Q, ExpPriority="mixed", ExpDirection="mixed", reduceQTo=0, iterations = 20, NNeighbors=0, maxDistance = 1, queensCase=FALSE, maxGDistance=1, np=1, deBug=FALSE) {

		matrices.list <- list()
		iteration = 0
		# SWFarea = length(which(Hmatrix%in%swfCat)) / length(which(Hmatrix%in%c(agriCat,swfCat)))
		SWFarea = length( which(Hmatrix%in%swfCat) ) / ( nrow(Hmatrix) * ncol(Hmatrix) )

		while(iteration < iterations && SWFarea<swfCover) {
			iteration=iteration+1
			
			agri.cells = as.data.frame(which(Hmatrix == agriCat | Hmatrix == boundaryCat, arr.ind = TRUE)) # Finds coords with agriCat and boundaryCat
			swf.cells <- as.data.frame(which(Hmatrix == swfCat, arr.ind = TRUE)) # Finds coords with swfCat

			if ( (nrow(swf.cells) != 0 ) ) { # Adds swf only if there is at least one swfCat pixel in the kernel
			# Chose a gravity cell and allocated NN all at once!

			allocated.cells <- findSwfCells(Hmatrix, swfCat, agriCat, boundaryCat, agriW, boundaryW, NNeighbors, queensCase, maxDistance, maxGDistance, Q, reduceQTo, ExpPriority, ExpDirection, np)

			if ( !is.null(allocated.cells) && nrow(allocated.cells) >= 1 ) {
				Hmatrix[as.matrix(allocated.cells)] <- swfCat
				} else {
					iteration=iteration+1;
				} 
				
				} else{
					iteration=iteration+1;
				}

				if( iteration>1 && identical(Hmatrix, matrices.list[[length(matrices.list)]]) ) {
					break
					} else {
						matrices.list[[iteration]] <- Hmatrix
						SWFarea = length( which(Hmatrix%in%swfCat[1]) ) / ( nrow(Hmatrix) * ncol(Hmatrix) )
						message(paste("Iteration: ", iteration, "; SWF cover:", round(SWFarea,2)))
					}

				}

				return(matrices.list)
			}