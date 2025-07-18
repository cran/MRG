
MRGparam = function(par) {

  ret = switch(par,
               MRGobject = "@param MRGobject An object including all
                 the necessary variables and parameters, from a call to \\code{\\link{createMRGobject}}",
               gdl = "@param gdl A list of gridded data with different resolutions (from a call to \\code{\\link{gridData}}",
               MRGinp = "@param MRGinp Either an MRGobject (from a call to \\code{\\link{createMRGobject}}) or
                           a list of gridded data with different resolutions (from a call to \\code{\\link{gridData}} or
                           a gridded sf-object (typically from an earlier call to \\code{multiResGrid})",
              himg = "@param himg The grid resulting from a call to multiResGrid",
              himg1 = "@param himg1 Either a grid resulting from a call to multiResGrid, or a list of such grids",
              himg2 = "@param himg2 A grid resulting from a call to multiResGrid",
              ifs = "@param ifs A data.frame or tibble with the locations and the data of the survey or census data" ,
               crsOut = "@param crsOut The coordinate reference system (crs) to be used ",
               ifg = "@param ifg Either a data.frame or tibble or sf-object with the locations and the data of the survey or census data,
                          or a list of such objects." ,
               ress = "@param ress A vector with the different resolutions", 
              coords = "@param coords Names of the numeric columns holding coordinates",
              res = "@param res A resolution or a vector with the different resolutions", 
              lnames = "@param lnames Names for the different surveys or censuses if ifg is a list. 
                         Typically it could be survey years",
               geovar = "@param geovar Name of geodata variable in the objects. Must me the same for all of the surveys/censuses, if 
                              the data sets are not submitted as sf-objects",
               vars = "@param vars Variable(s) of interest that should be aggregated (necessary when ifg is
                        used for individual farm specific anonymization rules)",
              vars1 = "@param vars1 Variable(s) of interest that should be merged from the first grid, or a list of variables,
                                one for each grid in the list \\code{himg1}",
              vars2 = "@param vars2 Variable(s) of interest that should be merged from the second grid",
              weights = "@param weights Extrapolation factor(s) (weights) wi of unit i in the sample of units nc 
                      falling into
                               a specific cell c. Weights are used for disclosure control measures. 
                               A weight of 1 will be used if missing.
                             If only one weight is given, it will be used for all variables. If the length is more than one,
                             the length has to be equal to the number of variables. If the same weight is used for several variables,
                             it must be repeated in the weights-vector",
               mincount = "@param mincount The minimum number of farms for a grid cell (threshold rule)",
              countFeatureOrTotal = "@param countFeatureOrTotal Should the frequency limit be applied on records with a positive
                               value for a certain feature, or on all records, independent of value of feature", 
              minpos = "@param minpos Minimum number of positive values for a variable", 
               nlarge = "@param nlarge Parameter to be used if the nlarge(st) farms should count for maximum plim percent of
                        the total value for the variable in the grid cell (see details of \\code{\\link{gridData}})",
               plim = "@param plim See nlarge",
               verbose = "@param verbose Indicates if some extra output should be printed. Usually TRUE/FALSE, but can also have 
                         a value of 2 for \\code{\\link{multiResGrid}} for even more output.",
               nclus = "@param nclus Number of clusters to use for parallel processing. No parallelization is used
                           for \\code{nclus = 1}.",
               clusType = "@param clusType The type of cluster; see \\code{\\link[parallel]{makeCluster}} for more details.
                        The default of makeCluster is used if type is missing or NA",
               domEstat = "@param domEstat Should the dominance rule be applied as in the IFS handbook (TRUE), where 
                         the weights are rounded before finding the  first nlarge contributors, or should 
                        it be the first nlarge contributors*weight, where also fractions are considered (FALSE)?", 
               consistencyCheck = "@param consistencyCheck logical; whether consistency between the gridded values and
                         the similar values from ifg should be checked. The gridded value is derived
                         from rasterize and the second one from st_join. The two methods can in some 
                         cases treat border cases between grid cells differently.",       
               outfile = "@param  outfile File to direct the output in case of parallel processing, 
                         see \\code{\\link[parallel]{makeCluster}} for more details.",
               splitlim = "@param splitlim For large dataset - split the data set in batches of more or less splitlim size",
               checkDominance = "@param checkDominance Logical - should the dominance rule be applied?",
              checkPpercent = "@param checkPpercent Logical - should the p-percent rule be applied?",
              pPercent = "@param pPercent Which limit to use for the p-Percent rule?",
              checkReliability = "@param checkReliability Logical - should the prediction variance be checked, and used for the aggregation?
                         This considerably increases computation time",
               pseudoreg = "@param pseudoreg A column with regions to be used to define pseudostrata if checkReliability is TRUE.
                            This is used for the cases when one or more strata only has a single record (and the weight is 
                            different from one). This makes variance calculation impossible, so such strata are 
                            merged into a pseudostrata. If pseudoreg is given (for example a column with the country name,
                            or NUTS2 region), the pseudostrata will be created separately for each pseudoreg region.", 
               userfun = "@param userfun This gives the possibility to add a user defined function with additional confidentiality rules which 
                             the grid cell has to pass", 
              fargs = "@param fargs The name of the necessary variables of userfun", 
              strat = "@param strat Column name defining the strata for stratified sampling, used if checkReliability is TRUE",
               confrules = "@param confrules Should the frequency rule (number of holdings) refer to the number of holdings with 
                           a value of the individual vars above zero (\"individual\") or the total number of holdings in 
                           the data set (\"total\")? ",
               suppresslim = "@param suppresslim Parameter that can be used to avoid that almost empty grid cells are merged with cells 
                            with considerably higher number of observations. The value is a minimum share of the total
                            potential new cell for a grid cell to be aggregated. See below for more details.",  
               sumsmall = "@param sumsmall Logical; should the suppresslimSum value be applied on the sum of
                            small grid cells within the lower resolution grid cell?
                            Note that different combinations of suppreslim and suppreslimSum values 
                            might not give completely intuitive results.For instance, if both are equal, then
                            a higher value can lead to more grid cells being left unaggregated for smaller grid sizes, leading
                            to aggregation for a large grid cell",
               suppresslimSum = "@param suppresslimSum Parameter similar to suppreslim, but affecting the total
                              of grid cells to be suppressed",
               plotIntermediate = "@param plotIntermediate Logical or number - make a simple plot showing which grid cells have already 
                           passed the frequency rule. plotintermediate = TRUE, the function will wait 5 seconds after plotting 
                           before continuing, otherwise it will wait plotintermediate seconds.",
               addIntermediate = "@param addIntermediate Logical; will add a list of all intermediate himgs
                           and lohs (overlay of himg and the lower resolution grid) as an attribute to
                           the object to be returned",
               reliabilitySplit = "@param reliabilitySplit Logical or number - parameter to be used in
                            calculation of the reliability (if checkReliability = TRUE). It can either give the 
                            number of groups, or if TRUE, it will create groups of approdcimately 50,000 records per group. 
                            If FALSE, the data set will not be split, independent on the size.", 
               locAdj = "@param locAdj parameter to adjust the coordinates if they are exactly on the borders between grid cells. The values
                          can either be FALSE, or \"jitter\" (adding a small random value to the coordinates, essentially spreading
                          them randomly around the real location), \"UR\", \"UL\", \"LR\" or \"LL\", to describe which corner of the grid 
                          cell the location belong (upper right, upper left, lower right or lower left).",
              remZeroes = "@param remZeroes Set to TRUE if the gridding should only be done on a 
                            subset of the data set, for which the value(s) of the variable(s) of interest
                            is larger than zero. If there is more than one variable of interest, the subset
                            will include all records where at least one of the variables are above zero.",
              postProcess = "@param postProcess Logical; should the postprocessing be done as part
                             of creation of the multiresolution grid (TRUE), or be done in a separate 
                             step afterwards (FALSE). The second option is useful when wanting
                             to check the confidential grid cells of the final map",
              remCols = "@param remCols Logical; Should intermediate columns be removed? Can be set
                          to FALSE for further analyses. Temporary columns will not be removed if their names 
                          partly match the variable names of \\code{vars}",
              rounding = "@param rounding either logical (FALSE) or an integer indicating the number 
                      of decimal places 
                       to be used. Negative values are allowed (such as the default
                       value rounding to the closest 10). See also the details
                       for \\code{digits} in \\code{\\link{round}}.",
              na.rm = "@param na.rm Should NA values be removed when summing values (essentially
                          treating them equal to zero) ",
              action = "@param action How to treat the values of overlapping grid cells. Possible values are:
                             \\describe{
                             \\item{none}{return an \\code{sf} data.frame just with the overlapping grid cells}
                               \\item{sum}{sum the values of the overlapping grid cells - NAs are ignored unless both cells
                               are NA or one is NA and one is 0}
                               \\item{sumna}{sum the values of the overlapping grid cells - sum is NA if any of them is NA}
                               \\item{avg}{avg of the grid cells - NAs are ignored unless both cells
                               are NA or one is NA and one is 0}
                               \\item{avgna}{avg of the grid cells - avg is NA if any of them is NA}
                               \\item{replace}{replace the problematic grid cells with grid cells from \\code{himg2}}
                             } ",
              ellipsis = "@param ... Possible arguments to underlying functions",
              ellipsisMerge = "@param ... Additional grids (himg3, himg4, ...) and variables (vars3, vars4, ...) to be merged. 
                       Additional grids and variables must be named.",
              fill = "@param fill Which column to use for fill color. The column name should be unquoted, as in ggplot-calls",
              color = "@param color Which column to use for border color. The column name should be unquoted, as in ggplot-calls.
                        The color is not used if \\code{lwd = 0} (although some devices, such as pdf, might print a line anway).",
              lwd = "@param lwd Line width for the grid cells. Default is zero, to remove or minimize (for pdf) the line width",
              name = "@param name Name to be used for color scale. The default is to use the name of the fill/color column. 
                             \\code{name = NULL} will give no name. ",
              borders = "@param borders A polygon object with borders than can be drawn on top of the multi-resolution grid. The 
                              object will also be used to clip the grid if \\code{clip = TRUE}.",
              xlim = "@param xlim The limits for the x-axis. The default is to use the bounding box of the grid.",
              ylim = "@param ylim The limits for the y-axis. The default is to use the bounding box of the grid.",
              title = "@param title The title of the plot",
              crs = "@param crs The coordinate reference system (CRS) into which all data should be projected before plotting. 
                          If not specified, will use the CRS defined in the first sf layer of the plot. ",
              clip = "@param clip Logical; should the grid be clipped to the borders object (if exsisting)?",
              transform = "@param transform Possible transformation of the color scale, typical values can be  \\code{\"log\"},
                          \\code{\"log10\"} or \\code{\"sqrt\"}, based on available transformations in the 
                          \\code{scales} package. See for example \\code{\\link[scales]{transform_log}} and other transformations for
                              more details.",
              show.legend = "@param show.legend Logical; should the legend be shown or not.",
              lwdb = "@param lwdb The line width for the border",
              option = "@param option The color map option to use, see \\code{\\link[viridis]{scale_color_viridis}} for more details",
              resc = "@param res The column name with grid cell sizes for the different grid cells",
              df = "@param df A data.frame or name of a csv file with multi-resolution data, 
                            only specifying the lower left corner of the grid cells",
              ellipsisc = "@param ... Additional parameters to read.csv if csv is a file name",
              var = "@param var Which variable to plot",
              linecolor = "@param linecolor Which column or color to use for lines between grid cells. The default is not to plot lines",
              coordscale = "@param coordscale Multiplication scaling factor for coordinates",              
              dsn = "@param dsn Source name to be used by \\code{\\link[sf]{st_write}}.  
                         Interpretation varies by driver: can be a filename, a folder, a database name, 
                         or a Database Connection. ",
              layer = "@param layer Layer name to be used by \\code{\\link[sf]{st_write}}. 
                          Varies by driver, may be a file name without extension; for database connection, it is the name of the table. 
                            If layer is missing, the basename of dsn is taken.",
              driver = "@param driver Character; name of driver to be used by \\code{\\link[sf]{st_write}}; 
                           if missing and dsn is not a Database Connection, 
                           a driver name is guessed from dsn.
                           st_drivers() returns the drivers that are available with their properties. 
                          links to full driver documentation are found at https://gdal.org/drivers/vector/index.html .",
              Estat = "@param Estat Indicate if Eurostat is the source of the data set. This is currently the default, 
                            but this might be changed in the future if other providers will follow the same conventions",
              cignore = "@param cignore Logical; Should the function ignore if parameters appear to be 
                           neither lat-lon or projected"
              )
ret
}
               
