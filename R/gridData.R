#' Function that converts point data to gridded data  (polygon values) or a list of gridded data
#' 
#'
#' @eval MRGparam("ifg")
#' @eval MRGparam("vars")
#' @eval MRGparam("weights")
#' @eval MRGparam("res")
#' @eval MRGparam("nclus")
#' @eval MRGparam("confrules")
#' @eval MRGparam("crsOut")
#' @eval MRGparam("verbose")
#' @eval MRGparam("locAdj", extra = "Please use with care in this function. It 
#'              will make it possible to produce the grid,but notice 
#'              that the coordinates of \\code{ifg}
#'              will be left untouched, which can cause problems if it is used in other functions.")
#' @eval MRGparam("centre")
#' 
#' 
#' @details This will create hierarchical grids of the selected variable(s), at the requested resolution(s),
#'          and using the requested function. In reality, the function will usually be sum,
#'          mean or max3, where the last one gives the average of the three highest numbers in the grid
#'          cell.
#'          
#'          Additionally, the function will always return the extrapolated number of farms per grid unit.
#'          The result will either be a set of sf-polygons (default) or a stars object.
#'
#' @returns A hierarchical list of gridded data, in the different resolutions requested.
#' Each grid also includes the count of records used for the gridding, and the
#' sum of the weights.
#'
#'
#' @details
#' Some of the sub functions are not sufficiently protected against too large memory usage, and 
#' in some cases this can lead to an aborted Rsession. If this happens, it is possible to use 
#' \code{\link{gridData2}}, which is somewhat slower, but better protected against memory issues.
#' 
#'
#'
#' @examples
#' \donttest{
#' library(sf)
#' library(dplyr)
#' if (!require(ggplot2)) print("Plotting of results will not work without installation of ggplot2")
#' if (!require(viridis)) print("Some of the plots will not work without installation of ggplot2")
#' if (require(giscoR)) {
#'   useBorder = TRUE 
#' } else {
#'   useBorder = FALSE
#'   print("You need to install giscoR for plotting borders and clipping the gridded maps")
#' }
#'
#' # These are SYNTHETIC agricultural FSS data 
#' data(ifs_dk) # Census data
#' # Create spatial data
#' ifg = fssgeo(ifs_dk, locAdj = "LL")
#' 
#' if (useBorder) {
#' # Read country borders, only used for plotting
#'   borders = gisco_get_nuts(nuts_level = 0)
#'   dkb = borders[borders$CNTR_CODE == "DK",] %>% st_transform(crs = 3035)
#' }
#' 
#' ress = c(1,5,10,20,40,80)*1000
#' ifl = gridData(ifg, vars = c("UAA", "UAAXK0000_ORG"), weights = "EXT_CORE", 
#'                res = ress)
#' ifl2 = gridData(ifg, vars = c("UAA", "UAAXK0000_ORG"), weights = "EXT_CORE", 
#'                res = ress, nclus = 2)
#' all.equal(ifl, ifl2)
#' if (require(ggplot2)) {
#' ifall = do.call("rbind", ifl)
#' g1 = ggplot() + geom_sf(data = ifall, aes(fill = count, color = count)) +
#'  scale_fill_viridis( name = "number of \n holdings", trans = "log10") +
#'  scale_color_viridis( name = "number of \n holdings", trans = "log10") +
#'  coord_sf(crs = 3035) +
#'  theme_bw() +
#'  ggtitle("Number of holdings for different resolutions") +
#'  facet_wrap(vars(res))
#'  if (useBorder) g1 = g1 + geom_sf(data = dkb, fill = NA, colour='black', lwd = 1)
#'  g1
#' }
#'#'
#'
#' MRGcluster(action = "stop")
#'
#'}
#'
#'
#' @export
gridData <-function (ifg, res = 1000, vars = NULL, weights = NULL,  
                   nclus = 1, confrules = "individual", crsOut = NA, verbose = FALSE, 
                   locAdj = FALSE, centre = FALSE) {
  
  checkVars(vars)
  if (verbose) print(paste("Function environment",  rlang::env_name(environment(fun = gridData))))
  if (centre & length(res) > 1) warning("centre should not be TRUE if res has more than one resolution, as the grids will not be overlapping")
  if ( !length(weights) %in% c(0,1,length(vars))) stop(paste("The length of weight should be 0,1 or equal to the length of vars"))
  if (!inherits(ifg, "sf"))  {
    tfss = system.time(ifg <- fssgeo(ifg, crsOut = crsOut, locAdj = locAdj))[3]
    if (verbose) print(paste("Spent ", tfss, " seconds on fssgeo"))
  } 
  if (!is.na(crsOut)) {
    if (st_crs(crsOut) != st_crs(ifg)) ifg = st_transform(ifg, crsOut)
  } else crsOut = st_crs(ifg)$epsg
  ifg$count = 1
  #' @importFrom terra rast rasterize ext xFromCol yFromRow res values
  #' @importFrom utils object.size
  if (length(res) > 1) {
    #' @importFrom dplyr lag
    rrat = res/lag(res, 1)
    if (nclus > 1) {
      cl = MRGcluster(nclus = nclus, action = "start")
      clusterEvalQ(cl, c(require(terra), require(sf), require(stars)))
      clusterExport(cl, varlist = c("res", "ifg", "weights", "vars", "addweights"), envir=environment())
      ret = parLapply(cl, res, fun = gridData, ifg = ifg, weights = weights, vars = vars, crsOut = crsOut, verbose = verbose)
    } else ret = lapply(res, FUN = gridData, ifg = ifg, weights = weights, vars = vars, crsOut = crsOut, verbose = verbose)
    return(ret)
  }

  rext = st_bbox(ifg)
  rext[c("xmin", "ymin")] = floor(rext[c("xmin", "ymin")]/res)*res - ifelse(centre, 0.5*res, 0)
  rext[c("xmax", "ymax")] = ceiling(rext[c("xmax", "ymax")]/res)*res + ifelse(centre, 0.5*res, 0)
  if (!is.na(crsOut) && length(grep("EPSG", crsOut)) == 0) rcrs = paste0("EPSG:", crsOut) else rcrs = crsOut
  tr1 = system.time(r0 <- rast(ext = ext(rext), resolution = res, crs = rcrs))[3]
  if (verbose) print(paste("Spent ", tr1, " seconds on creating base raster"))
  if (verbose) print(paste("Resolution, projections ifg, r0 and crsOut:", res, st_crs(ifg)$epsg, st_crs(r0)$epsg, crsOut, rcrs))  
    coors = st_coordinates(ifg)
    rx = xFromCol(r0)-res(r0)[1]/2
    ry = yFromRow(r0)-res(r0)[2]/2
    cx = unique(coors[,1])
    cy = unique(coors[,2])
    xdiffs = unlist(lapply(rx, FUN = function(rxx) min(abs(rxx-cx))))
    ydiffs = unlist(lapply(ry, FUN = function(ryy) min(abs(ryy-cy))))
    if (min(xdiffs)/res(r0)[1] < 1e-9 | min(ydiffs)/res(r0)[2] < 1e-9) {
        warning("One or more points are practically on the border between grid cells,
                                        it is advisible to shift the coordinates a bit (for 
                   example run \\code{ifg <- st_jitter(ifg)} beforehand, 
                or look at the \\code{locAdj} argument of \\code{\\link{fssgeo}}") 
    }

  if (verbose) print(paste("before rasterize - dim(ifg): ", dim(ifg)[1], "res:", res))
  if (verbose) print(paste(names(ifg)))

  tr2 = system.time(dnum <- rasterize(ifg, field = "count", r0, fun="sum"))[3]
  if (verbose) print(paste("Spent ", tr2, " seconds on rasterizing ifg"))
  if (verbose) print("rasterize(ifg, \"count\") succeeded")
  names(dnum) = "count"
  if (!is.null(vars)) {
    ifg = addweights(ifg, vars, weights)
    ifg$countw = ifg$count * st_drop_geometry(ifg[, paste0("weight_", vars[1])][[1]])
    if (verbose) print("rasterizing weighted count")
    if (verbose) print(paste("before rasterize countw - dim(ifg): ", dim(ifg)[1]))
    if (verbose) print(paste(names(ifg)))
    tr3 = system.time(dnumw <- rasterize(ifg, field = "countw", r0, fun = "sum"))[3]
    if (verbose) print(paste("succeeded rasterizing weighted count in ", tr3, "seconds"))
    names(dnumw) = "countw"
    dind = dweight = list()
    for (iw in 1:length(vars)) {
      if (confrules == "individual") ifg[st_drop_geometry(ifg[, vars[iw]]) == 0, paste0("weight_", vars[iw])]  = 0
      ifg[, paste(vars[iw],"_w", iw, sep="")] = 
              st_drop_geometry(ifg)[, vars[iw]] * 
                   st_drop_geometry(ifg[,paste0("weight_", vars[iw])])
      if (verbose) print(paste("rasterizing variable ", vars[iw]))
      tr4 = system.time(dind[[iw]] <- rasterize(ifg, field = paste0(vars[iw], "_w", iw), r0, fun = "sum")) [3]
      if (verbose) print(paste("succeeded rasterizing variable ", vars[iw], " in ", tr4, "seconds"))
      
      if (verbose) print(paste("rasterizing weight ", iw))
      tr5 = system.time(dweight[[iw]] <- rasterize(ifg, field = paste0("weight_", vars[iw]), r0, fun = "sum"))
      if (verbose) print(paste("succeeded rasterizing weight", iw, " in ", tr5, "seconds"))
      names(dweight[[iw]]) = paste0("weight_", vars[iw])
  #    names(dind[[iw]]) = c(vars[iw], paste0(vars[iw], "_w", iw))
      names(dind[[iw]]) = vars[iw]
    }
    #' @importFrom stars st_as_stars
    if (verbose) print("creating rast-object ")
    tr6 = system.time(ss <- rast(list(dnum, dnumw, rast(dind), rast(dweight))))[3]
    if (verbose) print(paste("succeeded creating rast-object in ", tr6, " seconds"))
    #' @importFrom terra values
    terra::values(ss) = terra::values(ss) # Forcing raster data into memory
    if (verbose) print(paste("succeeded forcing raster data into memory, with size in MB:", object.size(ss)/1e6))
    ts1 = system.time(ifsret <- st_as_stars(ss))
    if (verbose) print(paste("succeeded creating stars-object, with size in MB:", object.size(ifsret)/1e6 , " in ", ts1, "seconds"))
  }  else {
    dnumw = dnum
    names(dnumw) = "countw"
    if (verbose) print("rasterizing count ")
    tr11 = system.time(ss <- rast(list(dnum, dnumw)))
    if (verbose) print(paste("succeeded rasterizing count in ", tr11, "seconds"))
    terra::values(ss) = terra::values(ss)
    if (verbose) print(paste("succeeded forcing raster data into memory, with size in MB:", object.size(ss)/1e6))
    ts11 = system.time(ifsret <- st_as_stars(ss))
    if (verbose) print(paste("succeeded creating stars-object, with size in MB:", object.size(ifsret)/1e6 ), "in ", ts11, " seconds")
  }
    #' @importFrom sf st_as_sf
    ts2 = system.time(ifsret2 <- st_as_sf(ifsret))
    if (verbose) print(paste("succeeded creating sf-object, with size in MB:", object.size(ifsret2)/1e6, " in ", ts2, " seconds" ))
    rm(ifsret)
    ifsret2$res = res
    ifsret2$ID = 1:dim(ifsret2)[1]
    #' @importFrom sf st_crs st_transform
    if (!is.na(crsOut) & !is.na(st_crs(ifsret2)))  ifsret2 = st_transform(ifsret2, st_crs(crsOut))
  ifsret2
}

addweights = function(ifg, vars, weights) {
if (missing(weights) || is.null(weights) || (length(weights) == 1 && weights == 1)) {
  for (iw in 1:length(vars)) if (!paste0("weight_", vars[iw]) %in% names(ifg)) ifg[,paste0("weight_", vars[iw])] = 1 
} else if (length(weights) == 1) { 
  for (iw in 1:length(vars)) if (!paste0("weight_", vars[iw]) %in% names(ifg)) ifg[,paste0("weight_", vars[iw])] = as.numeric(data.frame(ifg)[, weights])
} else {
  for (iw in 1:length(weights)) if (!paste0("weight_", vars[iw]) %in% names(ifg)) ifg[,paste0("weight_", vars[iw])] = 
      as.numeric(data.frame(ifg)[, weights[iw]])
}
ifg
}

