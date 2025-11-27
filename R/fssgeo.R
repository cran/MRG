#' Function that creates an sf-object from IFS data
#'
#' @eval MRGparam("ifs")
#' @eval MRGparam("crsOut")
#' @eval MRGparam("geocol")
#' @eval MRGparam("locAdj")
#'
#' @details The geo-location in the FSS file has a particular format. For 2020, it includes country, coordinate reference system (CRS), resolution
#' (precision of coordinates) and coordinates
#' in one attribute ("GEO_LCT"). For past years, the FSS data structure differs and it includes three separate columns, like latitudes, longitudes and coordinate reference system.
#' This function splits the attribute in its individual parts, and creates an
#' sf-object with the correct coordinates and CRS.
#'
#' @returns An \code{\link[sf]{sf}}-object with the locations of the survey or census data
#'
#' @examples
#' data(ifs_dk)
#' ifg = fssgeo(ifs_dk)
#'
#'
#' @export
fssgeo = function(ifs, crsOut = NA, geocol = "GEO_LCT", locAdj = FALSE) {
  if (length(geocol) == 1 & geocol %in% names(ifs)) {
    geo = ifs[[geocol]]
    if (length(grep("^[A-Z][A-Z]", geo)) == 0) stop("The column ", geocol, "does not include inspire format coordinates")
    gspl = strsplit(geo, "_")
    countries = unlist(lapply(gspl, "[[", 1))
    country = unique(countries)
    gsplr = unlist(lapply(gspl, "[[", 2))
    gsplr = strsplit(gsplr, "RES")
    crsi = unlist(lapply(gsplr, "[[", 1))
    crsi = as.numeric(gsub("CRS", "", crsi))
    gsplr2 = unlist(lapply(gsplr, "[[", 2))
    gspl2 = strsplit(gsplr2, "MN")
    gsplr3 = unlist(lapply(gspl2, "[[", 2))
    clist = as.numeric(unlist(strsplit(gsplr3, "E")))
    coor = matrix(clist, ncol = 2, byrow = TRUE)
    ifs$xx = coor[,2]
    ifs$yy = coor[,1]
    ifs$crsi = crsi
    ifs$country = countries
  } else if (length(geocol) == 3) {
    ifs$xx = ifs[[geocol[1]]]
    ifs$yy = ifs[[geocol[1]]]
    ifs$crsi = ifs[[geocol[3]]]
  } else if(sum(c("A_1_1_NUMBER", "A_1_2_NUMBER", "A_1_3_CRD_REF") %in% names(ifs)) == 3) {
    ifs$yy=ifs$A_1_1_NUMBER
    ifs$xx=ifs$A_1_2_NUMBER
    ifs$crsi = ifs$A_1_3_CRD_REF
    ifs$crsi = ifelse(ifs$crsi %in% c("ETRS89","4"), 4258, 3035)
  } else {
    stop("geolocations are not given as an inspire format string or as three column names with x-coordinates, y-coordinates and crs")
  }
  sf_crs<-NULL
  crss = unique(ifs$crsi)
  for (icrs in crss){
    #' @importFrom dplyr filter bind_rows
    df <- ifs %>% filter(crsi == icrs)
    stf <- st_as_sf(df, coords = c("xx", "yy"), crs = st_crs(icrs))
    if (!is.na(crsOut)) stf = stf %>% st_transform(crs = crsOut)
    sf_crs <-  bind_rows(sf_crs, stf)
  }
  #' @importFrom sf st_crs
  if (!isFALSE(locAdj)) sf_crs = locAdjFun(sf_crs, locAdj)
  sf_crs
}
