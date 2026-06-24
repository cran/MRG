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
#' @eval MRGparam("tile_nx")
#' @eval MRGparam("tile_ny")
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
#'          The function is better protected  against memory issues than \code{\link{gridData}}, but somewhat slower.
#'
#' @returns A hierarchical list of gridded data, in the different resolutions requested.
#' Each grid also includes the count of records used for the gridding, and the
#' sum of the weights.
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
#' 
#' ress = c(1,5,10,20,40,80)*1000
#' system.time(ifl <- gridData(ifg, vars = c("UAA", "UAAXK0000_ORG"), weights = "EXT_CORE", 
#'                res = ress))
#' system.time(ifl2 <- gridData(ifg, vars = c("UAA", "UAAXK0000_ORG"), weights = "EXT_CORE", 
#'                res = ress, nclus = 2))
#' system.time(ifl11 <- gridData2(ifg, vars = c("UAA", "UAAXK0000_ORG"), weights = "EXT_CORE", 
#'                res = ress))
#' system.time(ifl22 <- gridData2(ifg, vars = c("UAA", "UAAXK0000_ORG"), weights = "EXT_CORE", 
#'                res = ress, nclus = 2))
#' system.time(ifl13 <- gridData2(ifg, vars = c("UAA", "UAAXK0000_ORG"), weights = "EXT_CORE", 
#'                res = ress, tile_nx = 1, tile_ny = 1))
#' system.time(ifl23 <- gridData2(ifg, vars = c("UAA", "UAAXK0000_ORG"), weights = "EXT_CORE", 
#'                res = ress, nclus = 2, tile_nx = 1, tile_ny = 1))
#' system.time(ifl14 <- gridData2(ifg, vars = c("UAA", "UAAXK0000_ORG"), weights = "EXT_CORE", 
#'                res = ress))
#' system.time(ifl24 <- gridData2(ifg, vars = c("UAA", "UAAXK0000_ORG"), weights = "EXT_CORE", 
#'                res = ress, nclus = 2))
#' all.equal(ifl, ifl2)
#' 
#' st_bbox(ifl[[1]])
#' st_bbox(ifl2[[1]])
#' st_bbox(ifl11[[1]])
#' st_bbox(ifl22[[1]])
#' st_bbox(ifl13[[1]])
#' st_bbox(ifl23[[1]])
#' st_bbox(ifl14[[1]])
#' st_bbox(ifl24[[1]])
#' 
#' summary(ifl[[1]])
#' summary(ifl2[[1]])
#' summary(ifl11[[1]])
#' summary(ifl22[[1]])
#' summary(ifl13[[1]])
#' summary(ifl23[[1]])

#' all.equal(as.data.frame(summary(ifl[[1]])), as.data.frame(summary(ifl2[[1]])))
#' all.equal(as.data.frame(summary(ifl[[1]])), as.data.frame(summary(ifl11[[1]])))
#' all.equal(as.data.frame(summary(ifl[[1]])), as.data.frame(summary(ifl22[[1]])))
#' all.equal(as.data.frame(summary(ifl[[1]])), as.data.frame(summary(ifl13[[1]])))
#' all.equal(as.data.frame(summary(ifl[[1]])), as.data.frame(summary(ifl23[[1]])))
#' all.equal(as.data.frame(summary(ifl[[1]])), as.data.frame(summary(ifl14[[1]])))
#' all.equal(as.data.frame(summary(ifl[[1]])), as.data.frame(summary(ifl24[[1]])))

#' 
#'#'
#'
#' MRGcluster(action = "stop")
#'
#'}
#'



#' @export
gridData2 <- function(ifg, res = 1000, vars = NULL, weights = NULL,
                     nclus = 1, confrules = "individual", crsOut = NA, verbose = FALSE,
                     locAdj = FALSE, centre = FALSE,
                     tile_nx = NULL, tile_ny = NULL) {
  checkVars(vars)
#' @importFrom rlang env_name
  if (verbose) print(paste("Function environment", env_name(environment(fun = gridData2))))
  if (centre & length(res) > 1) warning("centre should not be TRUE if res has more than one resolution, as the grids will not be overlapping")
  if (!length(weights) %in% c(0,1,length(vars))) stop("The length of weight should be 0,1 or equal to the length of vars")
  
  if (!inherits(ifg, "sf"))  {
    tfss = system.time(ifg <- fssgeo(ifg, crsOut = crsOut, locAdj = locAdj))[3]
    if (verbose) print(paste("Spent", tfss, "seconds on fssgeo"))
  }
  
  if (!is.na(crsOut)) {
    if (st_crs(crsOut) != st_crs(ifg)) ifg = st_transform(ifg, crsOut)
  } else crsOut = st_crs(ifg)$epsg
  
  ifg$count = 1
  
  if (length(res) > 1) {
#' @importFrom dplyr lag
    rrat = res / lag(res, 1)
    if (nclus > 1) {
      cl = MRGcluster(nclus = nclus, action = "start")
      clusterEvalQ(cl, c(require(terra), require(sf)))
      clusterExport(
        cl,
        varlist = c("res", "ifg", "weights", "vars", "addweights",
                    "spatRasterToSfTiles"),
        envir = environment()
      )
      ret = parLapply(
        cl, res, fun = gridData2, ifg = ifg, weights = weights, vars = vars,
        crsOut = crsOut, verbose = verbose, tile_nx = tile_nx, tile_ny = tile_ny
      )
    } else {
      ret = lapply(
        res, FUN = gridData2, ifg = ifg, weights = weights, vars = vars,
        crsOut = crsOut, verbose = verbose, tile_nx = tile_nx, tile_ny = tile_ny
      )
    }
    return(ret)
  }
  
  rext = st_bbox(ifg)
  rext[c("xmin", "ymin")] = floor(rext[c("xmin", "ymin")] / res) * res - ifelse(centre, 0.5 * res, 0)
  rext[c("xmax", "ymax")] = ceiling(rext[c("xmax", "ymax")] / res) * res + ifelse(centre, 0.5 * res, 0)
  
  if (!is.na(crsOut) && length(grep("EPSG", crsOut)) == 0) {
    rcrs = paste0("EPSG:", crsOut)
  } else {
    rcrs = crsOut
  }
#' @importFrom terra rast xFromCol yFromRow res rasterize values
  tr1 = system.time(r0 <- rast(ext = ext(rext), resolution = res, crs = rcrs))[3]
  if (verbose) print(paste("Spent", tr1, "seconds on creating base raster"))
  if (verbose) print(paste("Resolution, projections ifg, r0 and crsOut:",
                           res, st_crs(ifg)$epsg, st_crs(r0)$epsg, crsOut, rcrs))
  
  coors = st_coordinates(ifg)
  rx = xFromCol(r0) - res(r0)[1] / 2
  ry = yFromRow(r0) - res(r0)[2] / 2
  cx = unique(coors[,1])
  cy = unique(coors[,2])
  xdiffs = unlist(lapply(rx, FUN = function(rxx) min(abs(rxx - cx))))
  ydiffs = unlist(lapply(ry, FUN = function(ryy) min(abs(ryy - cy))))
  if (min(xdiffs) / res(r0)[1] < 1e-9 | min(ydiffs) / res(r0)[2] < 1e-9) {
    warning("One or more points are practically on the border between grid cells,
            it is advisible to shift the coordinates a bit (for example run ifg <- st_jitter(ifg)
            beforehand, or look at the locAdj argument of fssgeo")
  }
  
  if (verbose) print(paste("before rasterize - dim(ifg):", dim(ifg)[1], "res:", res))
  if (verbose) print(names(ifg))
  
  tr2 = system.time(dnum <- rasterize(ifg, field = "count", r0, fun = "sum"))[3]
  if (verbose) print(paste("Spent", tr2, "seconds on rasterizing ifg"))
  names(dnum) = "count"
  
  if (!is.null(vars)) {
    ifg = addweights(ifg, vars, weights)
    ifg$countw = ifg$count * st_drop_geometry(ifg[, paste0("weight_", vars[1])][[1]])
    
    tr3 = system.time(dnumw <- rasterize(ifg, field = "countw", r0, fun = "sum"))[3]
    if (verbose) print(paste("Succeeded rasterizing weighted count in", tr3, "seconds"))
    names(dnumw) = "countw"
    
    dind = dweight = list()
    for (iw in seq_along(vars)) {
      if (confrules == "individual") {
        ifg[st_drop_geometry(ifg[, vars[iw]]) == 0, paste0("weight_", vars[iw])] = 0
      }
      
      ifg[, paste(vars[iw], "_w", iw, sep = "")] =
        st_drop_geometry(ifg)[, vars[iw]] *
        st_drop_geometry(ifg[, paste0("weight_", vars[iw])])
      
      tr4 = system.time(
        dind[[iw]] <- rasterize(ifg, field = paste0(vars[iw], "_w", iw), r0, fun = "sum")
      )[3]
      if (verbose) print(paste("Succeeded rasterizing variable", vars[iw], "in", tr4, "seconds"))
      
      tr5 = system.time(
        dweight[[iw]] <- rasterize(ifg, field = paste0("weight_", vars[iw]), r0, fun = "sum")
      )[3]
      if (verbose) print(paste("Succeeded rasterizing weight", iw, "in", tr5, "seconds"))
      
      names(dweight[[iw]]) = paste0("weight_", vars[iw])
      names(dind[[iw]]) = vars[iw]
    }
    
    tr6 = system.time(ss <- rast(list(dnum, dnumw, rast(dind), rast(dweight))))[3]
    if (verbose) print(paste("Succeeded creating rast-object in", tr6, "seconds"))
  } else {
    dnumw = dnum
    names(dnumw) = "countw"
    tr11 = system.time(ss <- rast(list(dnum, dnumw)))[3]
    if (verbose) print(paste("Succeeded creating rast-object in", tr11, "seconds"))
  }
  
  values(ss) = values(ss)
  if (verbose) print(paste("Succeeded forcing raster data into memory, with size in MB:", object.size(ss) / 1e6))
  
  tile_choice <- choose_tiles(ss[[1]], nclus = nclus, tile_nx = tile_nx, tile_ny = tile_ny, verbose = verbose)
  tile_nx <- tile_choice$nx
  tile_ny <- tile_choice$ny
  ts2 = system.time(
    ifsret2 <- spatRasterToSfTiles(ss, nx = tile_nx, ny = tile_ny, na.rm = TRUE, verbose = verbose)
  )[3]
  if (verbose) print(paste("Succeeded creating sf-object, with size in MB:", object.size(ifsret2) / 1e6, "in", ts2, "seconds"))
  
  ifsret2$res = res
  ifsret2$ID = seq_len(nrow(ifsret2))
  
  if (!is.na(crsOut) & !is.na(st_crs(ifsret2))) {
    ifsret2 = st_transform(ifsret2, st_crs(crsOut))
  }
  
  ifsret2
}


addweights = function(ifg, vars, weights) {
  if (missing(weights) || is.null(weights) || (length(weights) == 1 && weights == 1)) {
    for (iw in 1:length(vars)) {
      if (!paste0("weight_", vars[iw]) %in% names(ifg)) {
        ifg[, paste0("weight_", vars[iw])] = 1
      }
    }
  } else if (length(weights) == 1) {
    for (iw in 1:length(vars)) {
      if (!paste0("weight_", vars[iw]) %in% names(ifg)) {
        ifg[, paste0("weight_", vars[iw])] = as.numeric(data.frame(ifg)[, weights])
      }
    }
  } else {
    for (iw in 1:length(weights)) {
      if (!paste0("weight_", vars[iw]) %in% names(ifg)) {
        ifg[, paste0("weight_", vars[iw])] = as.numeric(data.frame(ifg)[, weights[iw]])
      }
    }
  }
  ifg
}


choose_tiles <- function(r, nclus = 1, tile_nx = NULL, tile_ny = NULL, verbose = FALSE) {
  # Respect user-supplied values
  if (!is.null(tile_nx) && !is.null(tile_ny)) {
    return(list(nx = tile_nx, ny = tile_ny))
  }
  #' @importFrom terra ncell
  nc <- ncell(r)
  
  # Base choice from raster size
  if (nc <= 5e5) {
    nx <- ny <- 1
  } else if (nc <= 2e6) {
    nx <- ny <- 2
  } else if (nc <= 8e6) {
    nx <- ny <- 4
  } else if (nc <= 2e7) {
    nx <- ny <- 6
  } else {
    nx <- ny <- 8
  }
  
  # If parallelization is already used over resolutions, avoid too much tiling overhead
  if (nclus > 1) {
    nx <- min(nx, 4)
    ny <- min(ny, 4)
  }
  
  if (verbose) {
    message(sprintf("Automatic tiling: %s x %s for %s cells", nx, ny, format(nc, big.mark = ",")))
  }
  
  list(nx = nx, ny = ny)
}




spatRasterToSfTiles <- function(ss, nx = 4, ny = 4, na.rm = TRUE, verbose = FALSE) {
  stopifnot(inherits(ss, "SpatRaster"))
#' @importFrom terra nlyr ext crop as.polygons extract values "values<-"
  stopifnot(nlyr(ss) >= 1)
  
  r0 <- ss[[1]]
  e <- ext(r0)
  
  xbreaks <- seq(e$xmin, e$xmax, length.out = nx + 1)
  ybreaks <- seq(e$ymin, e$ymax, length.out = ny + 1)
  
  out_list <- vector("list", nx * ny)
  k <- 1L
  
  for (ix in seq_len(nx)) {
    for (iy in seq_len(ny)) {
      if (verbose) {
        message(sprintf("Processing tile (%s,%s) of (%s,%s)", ix, iy, nx, ny))
      }
      
      tile_ext <- ext(
        xbreaks[ix], xbreaks[ix + 1],
        ybreaks[iy], ybreaks[iy + 1]
      )
      
      ss_tile <- crop(ss, tile_ext, snap = "near")
      if (is.null(ss_tile) || ncell(ss_tile) == 0) {
        out_list[[k]] <- NULL
        k <- k + 1L
        next
      }
      
      # Keep only cells where at least one layer has a value
      vals_mat <- values(ss_tile, mat = TRUE)
      if (is.null(vals_mat) || nrow(vals_mat) == 0) {
        out_list[[k]] <- NULL
        k <- k + 1L
        next
      }
      
      keep <- rowSums(!is.na(vals_mat)) > 0
      if (!any(keep)) {
        out_list[[k]] <- NULL
        k <- k + 1L
        next
      }
      
      # Build a mask raster with 1 for cells to keep, NA otherwise
      mask_r <- ss_tile[[1]]
      values(mask_r) <- ifelse(keep, 1, NA)
      
      # Polygonize one cell at a time
      p <- as.polygons(
        mask_r,
        aggregate = FALSE,
        values = FALSE,
        na.rm = na.rm
      )
      
      if (is.null(p) || nrow(p) == 0) {
        out_list[[k]] <- NULL
        k <- k + 1L
        next
      }
      
      # Attach all layer values to the polygons
      vals <- extract(ss_tile, p, ID = FALSE)
      
      # Safety check
      keep2 <- rowSums(!is.na(vals)) > 0
      if (!any(keep2)) {
        out_list[[k]] <- NULL
        k <- k + 1L
        next
      }
      
      p <- p[keep2, ]
      vals <- vals[keep2, , drop = FALSE]
      
      p_sf <- st_as_sf(p)
      p_sf <- cbind(p_sf, vals)
      
      out_list[[k]] <- p_sf
      
      rm(ss_tile, vals_mat, keep, mask_r, p, vals, keep2, p_sf)
      gc()
      
      k <- k + 1L
    }
  }
  
  out_list <- Filter(Negate(is.null), out_list)
  
  if (length(out_list) == 0) {
    return(st_sf(geometry = sf::st_sfc(), crs = sf::st_crs(terra::crs(ss))))
  }
  
  do.call(rbind, out_list)
}
