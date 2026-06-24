#' Result object from slow example 
#' 
#' This object is the result from a slow function call in the examples of 
#' \code{\link{multiResGrid}}, and is used for the remainder of the example 
#' section.
#' 
#' 
#' 
#' 
#' @format A data frame with 285 rows and 7 variables
#' \itemize{
#'   \item count The number of records in the grid cell
#'   \item countw The weighted number of farms in the grid cell
#'   \item UAA The total utilized agricultural area per grid cell
#'   \item weight_UAA The sum of weights used to calculate the UAA
#'   \item res Resolution of grid cells
#'   \item ID Grid cell ID
#'   \item geometry The coordinates of the grid cells 
#' }
#' @docType data
#' @keywords datasets
#' @name himg5
#' @usage data(himg5)
#' @format A data frame with 285 rows and 7 variables
#' 
#' @details
#' 
#' This is the result of a very slow function call in the example of 
#' \code{\link{multiResGrid}}:
#' 
#' himg5 = multiResGrid(fsl,  vars = c("UAA"), weights = "EXT_MODULE", ifg = fsg, 
#'                       strat = "STRA_ID_CORE", checkReliability = TRUE, 
#'                       reliabilitySplit = TRUE, rounding = FALSE, pseudoreg = "REGIONS")
#'                       save(himg5, file = "d:/git/jskoien/MRG/data/himg5.rda")
#' 
#' 
#' 
NULL
