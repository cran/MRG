#' Merge two or more multi-resolution grids to a common resolution 
#' 
#' 
#' @eval MRGparam("himg1")
#' @eval MRGparam("himg2")
#' @eval MRGparam("vars1")
#' @eval MRGparam("vars2")
#' @eval MRGparam("postProcess")
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
#' @returns
#' The function produces a new multiresolution grid, which is a
#' \code{\link[sf]{sf}}-object with polygons.
#' 
#' @examples
#' \donttest{
#' library(sf)
#' library(dplyr)
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
MRGmerge = function(himg1, himg2, vars1, vars2, na.rm = TRUE, postProcess = FALSE, ...) {
  # To avoid R CMD check notes for missing global variables
  countw = ID = ID2 = NULL
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
    h2 = himgs[[il]]
    sfcol2 = attr(h2, "sf_column")
    #' @importFrom sf st_geometry "st_geometry<-"
    if (sfcol != sfcol2) st_geometry(fam) = sfcol
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
    h1tab = aggregate(rep(1, length(hm$ID)), by = list(ID = hm$ID), sum)
    h2tab = aggregate(rep(1, length(hm$ID)), by = list(ID2 = hm$ID2), sum)
    h1u = h1tab$ID[h1tab$x == 1]
    h2u = h2tab$ID2[h2tab$x == 1]
    unx = which(hm$ID %in% h1u & hm$ID2 %in% h2u)   
    une = which(!(hm$ID %in% h1u) & !(hm$ID2 %in% h2u)) 
    if (length(une) > 0) stop("Overlap error - could it be that at least one of the multiresolution grids has overlapping grid cells?")
    h11 = hm[unx,]
    h1i = which(hm$ID %in% h1tab$ID[h1tab$x > 1])
    h2i = which(hm$ID2 %in% h2tab$ID2[h2tab$x > 1])
    if (length(h1i) > 0) {
      h1a = hm[h1i,] %>% arrange(ID)
      h1aggr = aggregate(h1a[,vars2], by=list(ID = h1a$ID), FUN = sum, na.rm = na.rm )
      h1ids = h1a$ID[!duplicated(h1a$ID)]
      h1b = h1[h1$ID %in% h1ids,] %>% arrange(ID)
      if (!all.equal(h1b$ID, h1aggr$ID)) stop("mismatch in aggregated IDs - h1a")
      h1aggr = cbind(h1b[, vars1], st_drop_geometry(h1aggr))
    } else h1aggr = NULL
    if (length(h2i) > 0) {
      h2a = hm[h2i,] %>% arrange(ID2)
      h2aggr = aggregate(h2a[,vars1], by=list(ID2 = h2a$ID2), FUN = sum, na.rm = na.rm ) %>% arrange(ID2)
      h2ids = h2a$ID2[!duplicated(h2a$ID2)]
      h2b = h2[h2$ID2 %in% h2ids,] %>% arrange(ID2)
      if (!all.equal(h2b$ID2, h2aggr$ID2)) stop("mismatch in aggregated IDs - h2b")
      h2aggr = cbind(h2b[, vars2], st_drop_geometry(h2aggr))
    } else h2aggr = NULL
    h1 = rbind(h11[, c(vars1, vars2)], h1aggr[,-which(names(h1aggr) == "ID")], h2aggr[,-which(names(h2aggr) == "ID2")])
    vars1 = c(vars1, vars2)
    h1$ID = 1:dim(h1)[1]
  }
  h1$res = sqrt(st_area(h1))
  units(h1$res) = NULL
  vars = vars1[-grep("count|weight_", vars1)]
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
