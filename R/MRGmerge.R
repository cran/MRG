#' Merge two or more multi-resolution grids to a common resolution 
#' 
#' 
#' @eval MRGparam("himg1")
#' @eval MRGparam("himg2")
#' @eval MRGparam("vars1")
#' @eval MRGparam("vars2")
#' @eval MRGparam("postProcess")
#' @eval MRGparam("aggr")
#' @eval MRGparam("na.rm")
#' @eval MRGparam("ellipsisMerge")
#' 
#' @details
#' This function can merge different multi-resolution grids to a common 
#' resolution, i.e., it will select the grid cells with the lowest resolution,
#' as these are the ones defining the restrictions.  
#' 
#' The function will merge the variable names in \code{vars1, vars2, ...}
#' if they exist. If they are missing, the function will look for variable
#' names in the attributes of the grids (\code{attr(himg, "vars")}). These
#' are added by \code{\link{multiResGrid}}, but will often disappear if the
#' grid has been manipulated, or has been exported to another format for transmission.
#' 
#' If the variables are not given as \code{vars} or attributes, the function
#' will try to guess them from the column names. Typical column names used by
#' MRG (mostly temporary variables such as \code{small}, \code{confidential} etc)
#' will be ignored. If variable names partly coincide with any of these names,
#' or with \code{count}, \code{res}, \code{geometry}, it is necessary to specify vars.
#' 
#' The multi-resolution grids must be passed as named parameters if more than two 
#' are given. 
#' 
#' Common variable names in different grids should be avoided.
#' 
#' The default of the function is to treat NA-values as zeroes when merging 
#' (through \code{na.rm} in sums).
#' It will therefore not be possible to separate restricted grid cells
#' from grid cells with zero observations after merging, except for the ones that
#' have been left as they were. The alternative would
#' be a much higher number of NA-values in the merged grids.
#' 
#' The resulting grid will most likely not have exactly the same values as a 
#' multi-resolution grid produced 
#' directly from the microdata. If the input-grids have been post-processed
#' (the normal situation when not having access to the microdata), the 
#' grid cell values have usually been rounded, and some might have been 
#' suppressed. As these rounded and potentially suppressed values are summed,
#' their values are likely to deviate from those that are computed directly
#' from the microdata through a joint gridding process.
#' 
#' The argument \code{aggr} will decide on the direction of aggregation. 
#' If \code{aggr == "merge"}, The values in high resolution grid cells will 
#' be aggregated to match those of lower resolution grid cells in the 
#' second grid. If \code{aggr == "disaggr"}, the values of the lower resolution
#' grid cells will be redistributed equally among higher resolution grid cells,
#' according to their area. 
#' Note that this will most likely result in grid cell values that are apparently
#' confidential (for example having less than 10 individuals). These are still
#' not confidential values, but are average values from a larger area. 
#' This will in most cases be fine if the data is used for analyses, 
#' but publication of such values should be done with care.
#' 
#' Also note that if more than 2 MRG-grids are merged at the same time, 
#' then the redistribution will occur more than once. If the resolution of some grid cells
#' becomes higher for each redistribution, with some of the high resolution grid cells missing,
#' then the average values might differ for different high resolution grid cells coming from
#' the same low value grid cell. See the plotted examples of h2 and h22.
#' 
#' @returns
#' The function produces a new multiresolution grid, which is a
#' \code{\link[sf]{sf}}-object with polygons.
#' 
#' @examples
#' \donttest{
#' library(sf)
#' library(dplyr)
#' library(ggplot2)
#' library(viridis)
#' 
#' # These are SYNTHETIC agricultural FSS data 
#' data(ifs_dk) # Census data
#' ifs_weight = ifs_dk %>% dplyr::filter(Sample == 1) # Extract weighted subsample
#' 
#' # Create spatial data
#' ifg = fssgeo(ifs_dk, locAdj = "LL")
#' fsg = fssgeo(ifs_weight, locAdj = "LL")
#' 
#' # We use the numeric part of the farmtype to create a third variable. This 
#' # is done for the an example, the value does not have any meaning when treated 
#' # like this
#' ifg$ft = as.numeric(substr(ifg$FARMTYPE, 3, 4))^2
#' 
#' ress = c(1,5,10,20,40, 80, 160)*1000
#' # Create regular grid of the variables
#' ifl = gridData(ifg, vars = c("UAA", "UAAXK0000_ORG", "ft"), res = ress)
#'
#' # Create the different multi-resolution grids
#' himg1 = multiResGrid(ifl, vars = "UAA", ifg = ifg, postProcess = FALSE)
#' himg2 = multiResGrid(ifl, vars = "UAAXK0000_ORG", ifg = ifg, postProcess = FALSE)
#' himg3 = multiResGrid(ifl, vars = "ft", ifg = ifg, postProcess = FALSE)
#' 
#' # The grids have different number of polygons
#' dim(himg1)
#' dim(himg2)
#' dim(himg3)
#' 
#' hh1 = MRGmerge(himg1, himg2, himg3 = himg3)
#' dim(hh1)
#' # Postprocessing can also be done on the merged object
#' hh11 = MRGmerge(himg1, himg2, himg3 = himg3, postProcess = TRUE, rounding = -1)
#' dim(hh11)
#' summary(hh1$UAA-hh11$UAA)
#' 
#' # Here the merging will instead redistribute average values to 
#' # the higher resolution grid cells, and also seeing the effect
#' # of merging a third layer
#' hh2 = MRGmerge(himg1, himg2, aggr = "disaggr")
#' hh22 = MRGmerge(himg1, himg2, himg3 = himg3, aggr = "disaggr")
#' himg2$orgShare = himg2$UAAXK0000_ORG/himg2$res^2 * 10000
#' hh2$orgShare = hh2$UAAXK0000_ORG/hh2$res^2 * 10000
#' hh22$orgShare = hh22$UAAXK0000_ORG/hh22$res^2 * 10000
#' # Plot the organic share (organic area relative to grid cell area) for
#' # the original MRG grid for organic area, and after merging with the higher
#' # resolution maps.
#' p1 = ggplot(himg2) + geom_sf(aes(fill = orgShare)) + ggtitle("original") +
#'       scale_fill_viridis()
#' p2 = ggplot(hh2) + geom_sf(aes(fill = orgShare)) + ggtitle("merged two")+
#'       scale_fill_viridis() 
#' p3 = ggplot(hh22) + geom_sf(aes(fill = orgShare)) + ggtitle("merged three")+
#'       scale_fill_viridis() 
#' if (require(patchwork)) p1 + p2 + p3 + plot_spacer() + plot_layout(guides = 'collect')
#' 
#' # If two data sets share the same variable, one of them has to be renamed.
#' # (A comparison of the two can act as a indication of possible errors 
#' # introduced through the post-processing)
#' 
#' himg21 = multiResGrid(ifl, vars = c("UAA", "UAAXK0000_ORG"), ifg = ifg, postProcess = FALSE)
#' hh3 = try(MRGmerge(himg1, himg21, himg3 = himg3))
#' himg21 = himg21 %>% rename(UAA2 = UAA, weight_UAA2 = weight_UAA) 
#' hh3 = MRGmerge(himg1, himg21, himg3 = himg3)
#' 
#' 
#' summary(hh3[, c("UAA", "UAA2")])
#' 
#' himg4 = multiResGrid(ifl, vars = c("UAA", "ft", "UAAXK0000_ORG"), ifg = ifg, postProcess = FALSE)
#' summary(hh1[, c("UAA", "UAAXK0000_ORG", "ft")])
#' summary(himg4[, c("UAA", "UAAXK0000_ORG", "ft")])
#' }
#'            
#'            
#' @export
MRGmerge = function(himg1, himg2, vars1, vars2, na.rm = TRUE, postProcess = FALSE, aggr = "merge", ...) {
  # To avoid R CMD check notes for missing global variables
  if (!missing(vars1) && inherits(vars1, "data.frame")) stop("vars1 is a data.frame. Did you want
                                          to pass a third MRG-grid? Then it must be named, see the help file.")
  countw = ID = ID2 = area1 = area2 = NULL
  dots = list(...)
  #  Separate dots in himgs and vars
  hmgs = dots[grep("himg", names(dots))]
  vvs = dots[grep("vars", names(dots))]
  if ((inherits(himg1, "data.frame") | inherits(himg1, "sf")) & !missing(himg2)) {
    himgs = list(himg1, himg2)
    if (length(hmgs) > 0) himgs = c(himgs, hmgs)
  } else himgs = himg1
  if (length(himgs) <=1) stop("not enough grids to combine")  
  if (!missing(vars1) && is.list(vars1)) {
    vars = vars1
  }  else {
    if (!missing(vars1)) vars = list(vars1) else vars = list(NULL)
    if (!missing(vars2)) vars[[2]] = vars2 else vars = c(vars, list(NULL))
    if (length(vvs) > 0) vars = c(vars, vvs) else if (length(himgs) > 2) vars = c(vars, list(rep(NULL, length(himgs)-2))) 
  }
  
  h1 = himgs[[1]]
  if (is.null(vars[[1]])) {
    vars1 = attr(h1, "vars")
    vars1 = vars1[!vars1 %in% c("ID", "res", "area", attr(h1, "sf_column"))] 
  } else vars1 = vars[[1]]
  
  if (is.null(vars1)) vars1 = getVars(h1)  
  #' @importFrom dplyr rename
  vars1 = c("count1", "countw1", vars1, names(h1)[grep("weight_", names(h1))])
  if (!"count" %in% names(h1)) h1$count = NA
  if (!"countw" %in% names(h1)) h1$countw = NA
  h1 = h1 %>% rename(count1 = count, countw1 = countw)
  if (!"ID" %in% names(h1)) {
    h1 = h1 %>% mutate(ID = 1:dim(h1)[1])
  } else if (length(unique(h1$ID)) < length(h1$ID)) {
    h1 = h1 %>% mutate(ID = 1:dim(h1)[1])
    cat("Object with repeated IDs, it was fixed here, but this could indicate overlapping grid cells, please check \n")
    
  }
  sfcol = attr(h1, "sf_column")
  for (il in 2:length(himgs)){
    #' @importFrom units set_units
    h1 = h1 %>% mutate(area1 = set_units(st_area(h1), NULL))
    h2 = himgs[[il]]
    sfcol2 = attr(h2, "sf_column")
    #' @importFrom sf st_geometry "st_geometry<-"
    if (sfcol != sfcol2) st_geometry(h2) = sfcol
    if (!"ID" %in% names(h2)) {
      h2 = h2 %>% mutate(ID = 1:dim(h2)[1])
    } else if (length(unique(h2$ID)) < length(h2$ID)) {
      h2 = h2 %>% mutate(ID = 1:dim(h2)[1])
      cat("Object with repeated IDs, it was fixed here, but this could indicate overlapping grid cells, please check \n")
    }
    if (!"count" %in% names(h2)) h2$count = NA
    if (!"countw" %in% names(h2)) h2$countw = NA
    if (is.null(vars[[il]])) {
      vars2 = attr(h2, "vars") 
      vars2 = vars2[!vars2 %in% c("ID", "res", "area", attr(h2, "sf_column"))] 
    } else vars2 = vars[[il]]
    if (is.null(vars2)) vars2 = getVars(h2)  
    h2 = h2 %>% mutate(area2 = set_units(st_area(h2), NULL))    
    h2 = h2 %>% rename(!!paste0("count", il) := count, !!paste0("countw", il) := countw, ID2 = ID)
    vars2 = c(paste0("count", il), paste0("countw", il), vars2, names(h2)[grep("weight_", names(h2))])
    allvars = c(vars1, vars2)
    if (any(duplicated(allvars))) stop(cat("There are variables with same name, please rename:", allvars[duplicated(allvars)], "\n"))
    #' @importFrom sf st_intersection
    hm = st_intersection(h1, h2)
    hm$newArea = st_area(hm)
    units(hm$newArea) = NULL
    arange = diff(range(hm$newArea))
    hm = hm[hm$newArea > arange/1e6,]
    
    h1tab = aggregate(data.frame(Freq = rep(1, length(hm$ID)), cArea = hm$newArea), by = list(ID = hm$ID), sum)
    h2tab = aggregate(data.frame(Freq = rep(1, length(hm$ID)), cArea = hm$newArea), by = list(ID2 = hm$ID2), sum)
    h1u = h1tab$ID[h1tab$Freq == 1]
    h2u = h2tab$ID2[h2tab$Freq == 1]
    unx = which(hm$ID %in% h1u & hm$ID2 %in% h2u)   
    une = which(!(hm$ID %in% h1u) & !(hm$ID2 %in% h2u)) 
    unm1 = which(!h1$ID %in% hm$ID)
    unm2 = which(!h2$ID %in% hm$ID2)
    if (length(une) > 0) stop("Overlap error - could it be that at least one of the multiresolution grids has overlapping grid cells?")
    
    
    h1l = which(hm$area1 > hm$newArea)
    h2l = which(hm$area2 > hm$area1)

    h11 = hm[unx,]
    if (length(unm1) > 0) {
      h11 = bind_rows(h11, h1[unm1,])
      h1ia = which(hm$ID == h1$ID[unm1])
    } else h1ia = NULL
    if (length(unm2) > 0) {
      h11 = bind_rows(h11, h2[unm2,])
      h2ia = which(hm$ID2 == h2$ID2[unm2])
    }
    h1i = unique(c(which(hm$ID %in% h1tab$ID[h1tab$Freq > 1]), h1l))
    h2i = unique(c(which(hm$ID2 %in% h2tab$ID2[h2tab$Freq > 1]), h2l))
    
    if (aggr == "merge") {
      if (length(h1i) > 0) {
        h1a = hm[h1i,] %>% arrange(ID)
        h1aggr = aggregate(h1a[,vars2], by=list(ID = h1a$ID), FUN = sum, na.rm = na.rm )
        h1aggr$area1 = set_units(st_area(h1aggr), NULL)
        h1ids = unique(h1a$ID)
        h1b = h1[h1$ID %in% h1ids,] %>% arrange(ID)
        if (!all.equal(h1b$ID, h1aggr$ID)) stop("mismatch in aggregated IDs - h1a")
        h1aggr = cbind(h1b[, vars1], st_drop_geometry(h1aggr))
      } else h1aggr = NULL
      if (length(h2i) > 0) {
        h2a = hm[h2i,] %>% arrange(ID2)
        h2aggr = aggregate(h2a[,vars1], by=list(ID2 = h2a$ID2), FUN = sum, na.rm = na.rm ) %>% arrange(ID2)
        h2aggr$area2 = set_units(st_area(h2aggr), NULL)
        h2ids = unique(h2a$ID2)
        h2b = h2[h2$ID2 %in% h2ids,] %>% arrange(ID2)
        if (!all.equal(h2b$ID2, h2aggr$ID2)) stop("mismatch in aggregated IDs - h2b")
        h2aggr = cbind(h2b[, vars2], st_drop_geometry(h2aggr))
      } else h2aggr = NULL
      h1 = bind_rows(h11[, c(vars1, vars2, "area1", "area2")], h1aggr[,-which(names(h1aggr) == "ID")], h2aggr[,-which(names(h2aggr) == "ID2")])
    } else if (aggr == "disaggr") {
      #' @importFrom dplyr join_by 
      hm = hm %>% left_join(., h1tab, by = join_by(ID)) %>% left_join(., h2tab, by = join_by(ID2),
                                                                      suffix = c(".1", ".2"))
      h1a = h2a = NULL
      if (length(h1i) > 0) {
        h1a = hm[h1i,] %>% arrange(ID)
        #' @importFrom dplyr across
        h1a = h1a %>% mutate(across(!!vars1, ~ .x*cArea.2/cArea.1))
      }
      if (length(h2i) > 0) {
        h2a = hm[h2i,] %>% arrange(ID2)
        h2a = h2a %>% mutate(across(!!vars2, ~ .x*cArea.1/cArea.2))
      }
      h1 = rbind(hm[unx,], h1a, h2a)
    }
    vars1 = c(vars1, vars2)
    h1$ID = 1:dim(h1)[1]
    h1 = h1 %>% select(-matches("Freq|cArea|res.1|newArea|ID2"))
    h1$res = sqrt(st_area(h1)) %>% units::set_units(., NULL)
    if ("area2" %in% names(h1)) h1 = h1 %>% select(-area2)
  }
  h1 = h1 %>% select(-area1)
  vars = vars1[-grep("count|weight_|area1", vars1)]
  attr(h1, "vars") = vars
  if(postProcess) h1 = mergePP(h1, vars = vars1, ...)
  for (ii in 1:length(himgs)) {
    if (sum(h1[[paste0("count", ii)]], na.rm = TRUE) == 0) h1[[paste0("count", ii)]] = NULL
    if (sum(h1[[paste0("countw", ii)]], na.rm = TRUE) == 0) h1[[paste0("countw", ii)]] = NULL
  }
  h1
}


# To be able to extract only the relevant parameters from dots (ellipsis)
mergePP = function(himg, vars, remCols = TRUE, rounding = -1, ...){
  
  mc <- match.call()
  mc[[1]] <- as.name("MRGpostProcess")
  mc = mc[names(mc) %in% c("", "himg", "rounding", "remCols", "vars")]
  eval(mc, parent.frame())
  
}

if (FALSE) {
  
  sum(himg1$UAA)
  sum(himg2$UAAXK0000_ORG, na.rm = TRUE)
  sum(himg3$ft)
  colSums(st_drop_geometry(hh1)[, which(names(hh1) %in% c("UAA", "UAAXK0000_ORG", "ft"))], na.rm = TRUE)
  
  sum(st_area(himg1)/1e6)
  sum(st_area(himg2)/1e6)
  sum(st_area(himg3)/1e6)
  sum(st_area(hh1)/1e6)
  
}


getVars = function(h1, incCount = FALSE) {
  hnameso = names(h1)
  hnames = tolower(hnameso)
  if (incCount) {
    rids = grep("weight|geometry|res|small|reliability|idcount|idfail|vres|idRem|confidential|ufun|dom|freq|id", hnames)
  } else {
    rids = grep("count|weight|geometry|res|small|reliability|idcount|idfail|vres|idRem|confidential|ufun|dom|freq|id", hnames)
  }
  nonum = which(!unlist(lapply(h1, is.numeric) ))
  rids = unique(c(rids, nonum))
  if (length(rids) > 0) hnameso[-rids] else hnameso
}
