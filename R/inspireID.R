#' Function to add inspireIDs to a multi-resolution grid, particular for 
#' European data sets
#' 
#' @eval MRGparam("himg")
#' @eval MRGparam("borders2")
#' @eval MRGparam("cntrCol")
#' 
#' @details
#' The function will attempt to add an ID column following the principles of 
#' INSPIRE (Infrastructure for Spatial Information in Europe) if the 
#' multi-resolution grid has the inspire coordinate reference system (CRS=3035). 
#' See also https://epsg.io/3035 and https://inspire.ec.europa.eu/crs/3035.
#' 
#' The function will fail if the grid has any other CRS.
#' 
#' The function assumes that the grid is a correct multi-resolution grid
#' without overlapping grid cells. 
#' If the column res is missing, the function will create this as the 
#' square root of the area of each grid cell. This will not work if the 
#' grid has already been clipped with coastal/country/region borders.
#' The function tests if the number of different grid cell sizes is reasonable 
#' (maximum 10) and will issue a warning if it is higher. 
#'  
#' The two first letters in the ID is usually the country code. This can be added 
#' in two different ways. Either by specifying the column with country code
#' (argument: \code{cntrCol}) or by passing a polygon with country code information
#'  (argument: \code{borders}). It is assumed that this polygon has a column
#'  \code{CNTR_CODE} with the country codes, as typical in country-objects downloaded
#'  from GISCO. 
#'  
#'  However, it should be noted that grid cells on the borders will be associated
#'  to one of the countries - which in some cases could be unwanted. The 
#'  method uses the \code{\link[sf]{st_nearest_feature}} functionality for 
#'  associating a grid cell with a country. The example shows how this information
#'  can be added later, for users who want more control over the country association.
#'  
#'  If neither borders nor country column is submitted, the function will use
#'  NA instead of country code in the ID.
#' 
#' 
#' @examples
#' \donttest{
#' library(sf)
#' library(dplyr)
#' # These are SYNTHETIC agricultural FSS data 
#' data(ifs_dk) # Census data
#' 
#' # Create spatial data
#' ifg = fssgeo(ifs_dk, locAdj = "LL")
#'
#' ress = c(1,5,10,20,40, 80, 160)*1000
#' # Gridding Utilized agricultural area (UAA)
#' ifl = gridData(ifg, "UAA",res = ress)
#' 
#' # Create a multi-resolution grid only with farm number as confidentiality rule, then plot results
#' himg = multiResGrid(ifl, checkReliability = FALSE, suppresslim = 0)
#' 
#' himg = inspireID(himg) 
#' 
#' # It is easy to modify the country information afterwards
#' if (require(giscoR)) {
#'   borders = gisco_get_nuts(nuts_level = 0, epsg = 3035)
#'   himg1 = inspireID(himg, borders)
#'   himg2 = st_join(himg, borders %>% select(CNTR_CODE), join = st_nearest_feature)
#'   himg2 = inspireID(himg2, cntrCol = "CNTR_CODE")
#' # The border issues cause some grid cell to be classified as German, although 
#' # all data is from Denmark 
#'   table(substr(himg1$ID, 1, 2))
#' }
#' }
#' 
#' @export
inspireID = function(himg, borders, cntrCol) {
  if (!st_crs(himg)$epsg == 3035) stop("Inspire ID can only be added 
                                               to grids with EPSG:3035")

  #' @importFrom units set_units
  if (!"area" %in% names(himg)) areas = set_units(st_area(himg), NULL) else areas = himg$area
  if (!"res" %in% names(himg)) ress = sqrt(areas) else ress = himg$res
  if (!missing(cntrCol)) {
    cntr = himg[[cntrCol]]
  } else if (!missing(borders)) {
  #' @importFrom sf st_nearest_feature
    cntr = st_join(himg, borders, join = st_nearest_feature)$CNTR_CODE
  } else cntr = NA
  if (length(unique(ress)) > 10) warning(paste("There are ", length(unique(ress)), 
            "unique grid cell sizes, are you sure this is a regular multi-resolution grid?"))
  X = Y = L2 = ii = NULL # Avoid R CMD check Notes
  lccheck = st_coordinates(himg) %>% as.data.frame()  %>% 
             group_by(L2) %>%
             mutate(ii = 1) %>% summarise(nn = sum(ii))
  if (!all(lccheck$nn == 5)) stop("Not all elements are polygons with 4 corners")
  lls = cbind(st_coordinates(himg)[seq(1,5*dim(himg)[1], 5), ], res = ress) %>%
       as.data.frame() %>% 
       mutate(inspireID = paste0(cntr, "_CRS3035RES", res, "MN", X, "E", Y))
  himg$ID = lls$inspireID
  himg
}

