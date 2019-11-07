# FUNCTIONS TO CREATE, INSPECT AND MANIPULATE PCAVIZ OBJECTS
# ======================================================================
# Here we define the main function, "pcaviz", for creating a pcaviz
# object from PCA results, and some functions to manipulate pcaviz
# objects.
#
#   function              what is does
#   --------------        ------------------------------------
#   pcaviz                create pcaviz object
#   summary.pcaviz        generate summary of pcaviz object
#   print.pcaviz          display summary of pcaviz object
#   print.summary.pcaviz  display summary of pcaviz object
#   subset.pcaviz         return subset of observations
#   pcaviz_rotate         counter-clockwise plane rotation of PCs
#   pcaviz_reflect        reflect selected PCs about origin
#   pcaviz_scale          scale selected PCs
#   pcaviz_translate      translate points in selected PCs
#   pcaviz_transform2d    apply combination of 2-d transformations
#

# This function creates a "pcaviz" object from PCA results, and other
# accompanying data ("dat"). Note that the standard deviations of the
# PCs are equal to the square roots of the eigenvalues.
pcaviz <- function (out.pca, x = NULL, sdev = NULL, var = NULL,
                    rotation = NULL, dat = NULL, pc.cols) {

  # Extract results from input "out.pca", if provided. Note that it is
  # recommended to use the "inherits" function to check the class of
  # an object instead of "is" because it is faster and doesn't require
  # the "methods" package.
  if (!missing(out.pca)) {
    if (!is.null(x) | !is.null(sdev) | !is.null(rotation))
      stop(paste("When \"out.pca\" is provided, inputs \"x\", \"sdev\" and",
                 "\"rotation\" should not be provided"))
    if (inherits(out.pca,"prcomp") | inherits(out.pca,"rpca")) {

      # Extract results from principal components analysis using
      # "prcomp" or "rpca".
      x        <- out.pca$x
      sdev     <- out.pca$sdev
      rotation <- out.pca$rotation
      if (inherits(out.pca,"rpca"))
        var <- out.pca$var
      else
        var <- sum(sdev^2)
    } else if (inherits(out.pca,"princomp")) {

      # Extract results from principal components analysis using
      # "princomp".
      x        <- out.pca$scores
      sdev     <- out.pca$sdev
      var      <- sum(sdev^2)
      rotation <- as(out.pca$loadings,"matrix")
    } else
      stop(gettextf(paste("Argument \"out.pca\" should be a prcomp,",
                          "princomp or rpca object; got an object of",
                          "class \"%s\""),class(out.pca)))
  }

  # Verify and process input "dat".
  if (!is.null(dat)) {
    if (!is.data.frame(dat))
      stop("Argument \"dat\" should be a data frame")
    if (is.null(names(dat)))
      stop("Column names should be supplied for argument \"dat\"")
  }
  
  # Verify and process input "pc.cols".
  if (missing(pc.cols)) {
    if (is.null(x) & !is.null(dat))
      pc.cols <- names(dat)[grepl("PC",names(dat),ignore.case = TRUE)]
    else
      pc.cols <- NULL
  } else {
    if (is.null(dat))
      stop("Argument \"pc.cols\" cannot be used without providing \"dat\"")
    if (!is.null(x))
      stop(paste("Argument \"pc.cols\" cannot be specified when",
                 "\"x\" or \"out.pca$x\" is also provided"))
    if (!is.character(pc.cols))
      stop("Argument \"pc.cols\" should be a character vector")
    if (!all(is.element(pc.cols,names(dat))))
      stop("Argument \"pc.cols\" should select valid columns of \"dat\"")
  }

  # Extract the PCs from data frame "dat".
  if (!is.null(pc.cols)) {
    nonpc.cols <- setdiff(names(dat),pc.cols)
    x          <- dat[pc.cols]
    dat        <- dat[nonpc.cols]
  }
  
  # Verify and process input "x". Note that at least 2 PCs must be
  # supplied.
  if (is.null(x))
    stop(paste("Rotated data should be provided by arguments \"out.pca\",",
               "\"x\" or \"dat\""))
  if (is.matrix(x))
    x <- as.data.frame(x)
  if (!is.data.frame(x))
    stop("Argument \"x\" should be a matrix or data frame")
  if (is.null(names(x)))
    stop("Column names should be supplied for argument \"x\"")
  if (!all(sapply(x,is.numeric)))
    stop("All columns of \"x\" should be numeric")
  if (!is.null(dat))
    if (nrow(x) != nrow(dat))
      stop("Rows of argument \"dat\" must match rows of \"x\"")
  if (ncol(x) < 2)
    stop("Objects of class \"pcaviz\" require at least 2 principal components")

  # Verify and process input "sdev".
  if (!is.null(sdev)) {
    if (!(is.numeric(sdev) & is.vector(sdev)))
      stop("Argument \"sdev\" should be a numeric vector")
    if (length(sdev) != ncol(x))
      stop("Argument \"sdev\" should contain an entry for each PC/eigenvector")
  } else
    sdev <- rep(NA,ncol(x))
  names(sdev) <- names(x)

  # Verify and process input "var".
  if (!is.null(var)) {
    if (!(is.numeric(var) & length(var) == 1))
      stop("Argument \"var\" should be a numeric scalar")
  } else
    var <- NA
  
  # Verify and process input "rotation". Note that the columns of the
  # loadings matrix should match up with the columns of the rotated
  # data matrix.
  if (!is.null(rotation)) {
    if (is.matrix(rotation))
      rotation <- as.data.frame(rotation)
    if (!is.data.frame(rotation))
      stop("Argument \"rotation\" should be a matrix or data frame")
    if (is.null(names(rotation)))
      stop("Column names should be supplied for argument \"rotation\"")
    if (!all(is.element(names(x),names(rotation))))
      stop(paste("All principal components should also be",
                 "columns of argument \"rotation\""))
  }

  # At this point, we should have the following data frames: "x",
  # which contains the PCs; "rotation", which contains the
  # eigenvectors (possibly NULL); and "dat", which contains
  # accompanying data (possibly NULL).
  #
  # Get the number of columns in the data table (n) and the number of
  # principal components (m).
  m <- ncol(x)
  if (is.null(dat))
    n <- 0
  else
    n <- ncol(dat)
  
  # Get the column types of the combined data frame.
  data.coltypes <- c(rep("other",n),rep("pc",m))
  if (!is.null(dat)) {
    data.coltypes[which(sapply(dat,is.factor))]  <- "categorical"
    data.coltypes[which(sapply(dat,is.numeric))] <- "continuous"
  }
  data.coltypes <- factor(data.coltypes,data.coltypes.levels())

  # Create the pcaviz object.
  if (is.null(dat))
    dat <- x
  else
    dat <- cbind(dat,x)
  out <- list(data            = dat,
              data.coltypes   = data.coltypes,
              sdev            = sdev,
              var             = var,
              rotation        = rotation,
              basis           = rbind(0,diag(m)),
              transformed.pcs = FALSE)
  names(out$data.coltypes) <- names(out$data)
  colnames(out$basis)      <- names(x)
  rownames(out$basis)      <- c("origin",names(x))
  class(out) <- "pcaviz"
  return(out)
}

# Generate a summary of the data contained in the PCAviz object.
summary.pcaviz <- function (object, ...) {
  n   <- ncol(object$data)
  out <- data.frame(variable = names(object$data),
                    type     = object$data.coltypes,
                    n        = sapply(object$data,function (x) sum(!is.na(x))),
                    stats    = "",
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  for (i in 1:n) {
    a <- object$data.coltypes[i]
    x <- object$data[[i]]
    if (a == "pc") {
      y <- sprintf("(%0.4g,%+0.3g,%+0.3g,%+0.3g)",
                   object$sdev[out$variable[i]],
                   min(x,na.rm = TRUE),median(x,na.rm = TRUE),
                   max(x,na.rm = TRUE))
    } else if (a == "continuous") {
      y <- sprintf("(%0.3g,%0.3g,%0.3g)",min(x,na.rm = TRUE),
                   median(x,na.rm = TRUE),max(x,na.rm = TRUE))
    } else if (a == "categorical") {
      y <- sprintf("%d levels, largest=%s (%d)",nlevels(x),
                   names(which.max(table(x))),max(table(x)))
    } else if (a == "other")
      y <- NA
    out[i,"stats"] <- y
  }
  class(out) <- c("summary.pcaviz","data.frame")
  attr(out,"transformed.pcs") <- object$transformed.pcs
  return(out)
}

# Display summary of pcaviz object.
print.summary.pcaviz <- function (x, n = 4, ...) {
  headers      <- list(categorical = "categorical variables",
                       continuous  = "continuous variables",
                       other       = "other variables",
                       pc          = "principal components (PCs)")
  descriptions <- list(categorical = NULL,
                       continuous  = "# statistics are (min,median,max)",
                       other       = NULL,
                       pc = paste0("# statistics are (s.d.,min,median,max)\n",
                                   "# s.d.=sqrt(eigenvalue)"))
  if (!(is.numeric(n) & length(n) == 1 & n > 0))
    stop("Argument \"n\" should be a positive scalar")
  if (attr(x,"transformed.pcs"))
    headers$pc <- paste("transformed",headers$pc)
  class(x) <- "data.frame"
  for (a in levels(x$type)) {
    header <- headers[[a]]  
    rows   <- which(x$type == a)
    if (length(rows) > n) {
      header <- sprintf("first %d (of %d) %s",n,length(rows),header)
      rows   <- rows[1:n]  
    }
    cat(header,": ",sep="")
    if (length(rows) == 0)
      cat("none\n")
    else {
      cat("\n")
      if (!is.null(descriptions[[a]]))
        cat(descriptions[[a]],"\n")
      print(x[rows,c("variable","n","stats")],right = FALSE,row.names = FALSE)
    }
  }
  return(invisible(x))
}

# Display summary of pcaviz object.
print.pcaviz <- function (x, n = 4, ...)
  print(summary(x,n,...))

# Return a new pcaviz object containing a subset of the
# observations. This function also illustrates how "tryCatch" can be
# used to give the user a more informative error message indicating
# that there was an error evaluating the "subset" expression (at the
# cost of making the code a little more difficult to follow).
subset.pcaviz <- function (x, subset, ...) {
  e <- substitute(subset)
  tryCatch(r <- eval(e,x$data,parent.frame()),
    error = function (x)
      stop(paste("Failure in evaluating expression \"subset\" on",
                 "x$data in subset.pcaviz:\n",x$message),
           call. = FALSE))
  r <- eval(e,x$data,parent.frame())
  if (!is.logical(r))
    stop("Argument \"subset\" should be a logical expression")
  x$data <- droplevels(x$data[r,])
  return(x)
}

# Apply counter-clockwise plane rotation to PCs.
pcaviz_rotate <- function (x, angle, dims, units = c("degrees","radians")) {
  units <- match.arg(units)

  # Check input "x".
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                        "object of class \"%s\""),class(x)))

  # Choose the two PCs to rotate.
  if (missing(dims))
    dims <- get.pc.cols(x)[1:2]
  if (!valid.pc.dims(x,dims))
    stop("Argument \"dims\" is not a valid selection of PCs in \"x\"")
  
  # Apply the counter-clockwise rotation to the selected PCs, and
  # return the pcaviz object containing the transformed PCs.
  cols              <- get.pc.cols(x)
  x$data[cols]      <- matrix_rotate(as.matrix(x$data[cols]),angle,dims,units)
  x$basis           <- matrix_rotate(x$basis,angle,dims,units)
  x$transformed.pcs <- TRUE
  return(x)
}

# Reflect PCs about the origin, individually along each selected
# dimension.
pcaviz_reflect <- function (x, dims) {
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                          "object of class \"%s\""),class(x)))
  if (!valid.pc.dims(x,dims))
    stop("Argument \"dims\" is not a valid selection of PCs in \"x\"")
  cols              <- get.pc.cols(x)
  x$data[cols]      <- matrix_reflect(as.matrix(x$data[cols]),dims)
  x$basis           <- matrix_reflect(x$basis,dims)
  x$transformed.pcs <- TRUE
  return(x)
}

# Scale selected PCs.
pcaviz_scale <- function (x, scale, dims) {
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                          "object of class \"%s\""),class(x)))
  if (!valid.pc.dims(x,dims))
    stop("Argument \"dims\" is not a valid selection of PCs in \"x\"")
  cols              <- get.pc.cols(x)
  x$data[cols]      <- matrix_scale(as.matrix(x$data[cols]),scale,dims)
  x$basis           <- matrix_scale(x$basis,scale,dims)
  x$transformed.pcs <- TRUE
  return(x)
}

# Translate selected PCs by a.
pcaviz_translate <- function (x, a, dims) {
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                          "object of class \"%s\""),class(x)))
  if (!valid.pc.dims(x,dims))
    stop("Argument \"dims\" is not a valid selection of PCs in \"x\"")
  cols              <- get.pc.cols(x)
  x$data[cols]      <- matrix_translate(as.matrix(x$data[cols]),a,dims)
  x$basis           <- matrix_translate(x$basis,a,dims)
  x$transformed.pcs <- TRUE
  return(x)
}

# This is a convenience function that applies a sequence of
# transformations to 2 PCs: (1) counter-clockwise rotation, (2)
# reflection, (3) scaling and (4) translation. Note that exactly two
# dimensions must be selected.
pcaviz_transform2d <- function (x, dims, angle = 0,reflect.x = FALSE,
                                reflect.y = FALSE,scale = c(1,1),a = c(0,0),
                                units = c("degrees","radians")) {
  units <- match.arg(units)

  # Check input "x".
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                          "object of class \"%s\""),class(x)))

  # Choose the two PCs to transform.
  if (missing(dims))
    dims <- get.pc.cols(x)[1:2]
  if (!valid.pc.dims(x,dims))
    stop("Argument \"dims\" is not a valid selection of PCs in \"x\"")

  # Apply the sequence of 2-d transformations, and return the pcaviz
  # object containing the transformed PCs.
  cols         <- get.pc.cols(x)
  x$data[cols] <- matrix_transform2d(as.matrix(x$data[cols]),dims,angle,
                                     reflect.x,reflect.y,scale,a,units)
  return(x)
}

# FUNCTIONS TO PLOT PCAVIZ OBJECTS
# ======================================================================
# We have made three design decisions in the implementation of these
# plotting interfaces: (1) as much as possible, we try to provide
# sensible defaults so that interesting, visually compelling plots can
# be generated with minimal user intervention; (2) the defaults should
# include as many plotting features as possible so that the user is
# more immediately aware of these features---of course, all default
# features can be easily turned off manually by setting one of the
# input arguments; (3) the plotting interfaces should be as flexible
# as possible so that they can accomodate manual tuning by the user to
# produce visually attractive plots for presentations, publications,
# webpages, etc.
#
#   function                  what is does
#   --------------            ----------------------------------------
#   pcaviz_colors_categorical colors for plotting categorical variables
#   pcaviz_colors_continuous  colors for plotting continuous variables
#   pcaviz_shapes             shapes for plotting categorical variables
#   theme_pcaviz              default theme for plotting pcaviz objects
#   pcaviz_reduce_whitespace  rescale PCs to reduce "whitespace" in a plot
#   pcaviz_abbreviate_var     abbreviate variable in pcaviz object
#   create_abbreviated_label  abbreviate factor or character vector
#   pcaviz_violin             visualize PC data using violin plots
#   pcaviz_ggplot             main plotting interface (S3 method "plot")
#   pcaviz_screeplot          Plots variances against PCs.
#   pcaviz_loadingsplot       Plots eigenvector ("loadings").
#

# This function specifies the default colors for plotting categorical
# variables in the PCAviz plotting interface. If more than 7 colors
# are requested, the colors are repeated. Since LCM(5,7) = 35,
# functions "pcaviz_colors_categorical" and "pcaviz_shapes" can be
# easily combined to create 35 color-shape combinations.
pcaviz_colors_categorical <- function (n = 7) {
  if (!(is.numeric(n) & length(n) == 1 & n > 0))
    stop("Argument \"n\" should be a positive scalar")
  colors <- c("#E69F00","#56B4E9","#009E73","#F0E442",
              "#0072B2","#D55E00","#CC79A7")
  return(rep(colors,length.out = n))
}

# This function specifies the default color scheme for plotting
# continuous variables in the PCAviz plotting interface.
pcaviz_colors_continuous <- function()
  c("#0D0887FF","#5402A3FF","#8B0AA5FF","#B93289FF",
    "#DB5C68FF","#F48849FF","#FEBC2AFF","#F0F921FF")

# This function specifies the default shapes for plotting categorical
# variables in the PCAviz plotting interface. If more than 5 shapes
# are requested, the shapes are repeated. Since LCM(5,7) = 35,
# functions "pcaviz_colors_categorical" and "pcaviz_shapes" can be
# easily combined to create 35 color-shape combinations.
pcaviz_shapes <- function (n = 5) {
  if (!(is.numeric(n) & length(n) == 1 & n > 0))
    stop("Argument \"n\" should be a positive scalar")
  shapes <- c(19,17,8,1,3)
  return(rep(shapes,length.out = n))
}

# Returns the default theme for plotting pcaviz objects. The name of
# the function follows the convention used in the "ggthemes" package.
theme_pcaviz <- function()
  theme_cowplot(font_size = 10) +
  theme(panel.border = element_blank(),
        plot.title   = element_text(face = "plain"))

# This function heuristically rescales the PC data in a pcaviz object
# to reduce the amount of "whitespace" (that is, unused space) in a
# plot. The first dimension is always unscaled, and the remaining
# dimensions are scaled relative to the first.
pcaviz_reduce_whitespace <- function (x, dims) {

  # Check input argument "x".
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                          "object of class \"%s\""),class(x)))

  # Choose the PC dimensions to rescale.
  if (missing(dims))
    dims <- get.pc.cols(x)
  if (!valid.pc.dims(x,dims))
    stop("Argument \"dims\" is not a valid selection of PCs in \"x\"")

  # Get the data columns that will be rescaled.
  data <- x$data[dims]

  # Determine the scaling factor for each of the selected colors.
  i     <- dims[1]
  s1    <- diff(range(data[i]))
  scale <- sapply(data,function (x) s1/diff(range(x)))
  
  # Rescale the selected columns.
  x$data[dims]      <- matrix_scale(as.matrix(data),scale,dims)
  x$basis           <- matrix_scale(x$basis,scale,dims)
  x$transformed.pcs <- TRUE
  return(x)
}

# Returns a vector of abbreviated strings, allowing for either a
# customized abbreviation (e.g., specific to country names) or a
# generic abbreviation.
create_abbreviated_label <-
  function (x,
            abbrv.generic = function(x) abbreviate(gsub("[^[:alnum:] ]","",x)),
            abbrv.custom = NULL) {

  # Check input argument "x", and get "y", the character vector to
  # convert.
  if (!(is.factor(x) | is.character(x)))
    stop("Argument \"x\" should be a character or factor")
  if (is.factor(x))
    y <- levels(x)
  else
    y <- x

  # Try to generate the custom abbreviations.
  n <- length(y)
  a <- rep(as.character(NA),n)
  if (!is.null(abbrv.custom)) {
    if (!is.function(abbrv.custom))
      stop("Argument \"abbrv.custom\" should be a function")
    a <- abbrv.custom(y)
    if (!is.character(a))
      stop("Argument \"abbrv.custom\" should return a character vector")
  }

  # Generate generic abbreviations, if needed.
  if (any(is.na(a))) {
    if (!is.function(abbrv.generic))
      stop("Argument \"abbrv.generic\" should be a function")
    a <- abbrv.generic(y)
    if (!is.character(a))
      stop("Argument \"abbrv.generic\" should return a character vector")
  }

  # If x is a factor, replace the levels with the abbreviations.
  # Otherwise, return the abbreviations.
  if (is.factor(x)) {
    if (length(unique(a)) < nlevels(x))
      stop(paste("create_abbreviated_label did not generate unique",
                 "abbreviations for factor levels; consider setting",
                 "abbrv.custom = NULL to avoid this"))
    levels(x) <- a
    return(x)
  } else
    return(a)
}

# Generates an abbreviated version of a specified variable in a pcaviz
# object using "create_abbreviated_label".
pcaviz_abbreviate_var <-
  function (x, col,
            abbrv.generic = function (x)
              abbreviate(abbreviate(gsub("[^[:alnum:] ]","",x)),
                         minlength = 2,named = FALSE,method = "both.sides"),
            abbrv.custom = function (x)
              countrycode(x,"country.name.en","iso2c",warn = FALSE)) {

  # Check input "x".
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                          "object of class \"%s\""),class(x)))

  # Check input "col". Note that additional checks are made in the
  # call to function "create_abbreviated_label" below.
  if (!valid.data.col(x,col))
    stop("Argument \"col\" is not a valid name of a data column in \"x\"") 

  # Generate the abbreviated variable.
  y       <- create_abbreviated_label(x$data[[col]],abbrv.generic,abbrv.custom)
  col.new <- paste0(col,".abbrv")
  if (is.element(col.new,names(x$data)))
    stop(gettextf(paste("Cannot create abbreviated column; column \"%s\"",
                        "already exists"),col.new))

  # Add the abbreviated variable to the data frame.
  return(pcaviz.add.col(x,y,col.new,x$data.coltypes[col]))
}

# This is the default function for generating the group summary. It
# takes as input a data frame, and outputs a data frame with the same
# number of columns, and only one row: if a column is numeric, it
# outputs the median; if the column is a factor, it outputs the level
# with the highest count; otherwise, it just outputs the first row. 
pcaviz_summary_default <-
  function (x, stat.numeric = function (x) median(x,na.rm = TRUE)) {

  # Initialize the output.
  out <- x[1,]

  # Repeat for each column.
  n <- ncol(x)
  for (i in 1:n) {
    y <- x[[i]]
    if (is.numeric(y))
      out[[i]] <- stat.numeric(y)
    else if (is.factor(y))
      out[[i]] <- names(which.max(table(y)))
  }
  return(out)
}

# Overlay a map of the World onto the current plot.
overlay_map_world <- function (g) {
  if (!inherits(g,"ggplot"))
    stop("Argument \"g\" should be a ggplot object")
  return(g + geom_path(data = map_data("world"),
                       aes_string(x = "long",y = "lat",group = "group"),
                       color = "gray",size = 0.5,inherit.aes = FALSE) +
         coord_cartesian(xlim = c(-180,190),ylim = c(-80,80)))
}

# Overlay a map of Europe onto the current plot.
overlay_map_europe <- function (g) {
  if (!inherits(g,"ggplot"))
    stop("Argument \"g\" should be a ggplot object")
  return(g +
         geom_path(data = map_data("world"),
                   aes_string(x = "long",y = "lat",group = "group"),
                   color = "gray",size = 0.5,inherit.aes = FALSE) +
         coord_cartesian(xlim = c(-12,43),ylim = c(30,70)))
}

# Overlay a map of the US, including US state boundaries.
overlay_map_usa <- function (g) {
  if (!inherits(g,"ggplot"))
    stop("Argument \"g\" should be a ggplot object")
  return(g +
         geom_path(data = map_data("state"),
                   aes_string(x = "long",y = "lat",group = "group"),
                   color = "gray",size = 0.5,inherit.aes = FALSE) +
         geom_path(data = map_data("world"),
                   aes_string(x = "long",y = "lat",group = "group"),
                   color = "gray",size = 0.5) +
         coord_cartesian(xlim = c(-130,-60),ylim = c(16,50)))
}

# Generate one or several "violin" plots for visualizing the
# relationship between PCs and a categorical variable.
pcaviz_violin <-
    function (x, data.col, pc.dims, colors, sorted = TRUE,
              rank.fun = median, horizontal = FALSE, theme = theme_pcaviz(),
              violin.params = list(trim = FALSE, show.legend = FALSE),
              plot.grid.params = list()) {
   
  # Check input "x".
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                          "object of class \"%s\""),class(x)))

  # Choose which categorical variable to plot against the PCs.
  if (missing(data.col)) {
    data.col <- get.cols.by.type(x,"categorical")
    if (length(data.col) == 0)
      stop(paste("At least one categorical variable is needed for violin",
                 "plot; no categorical variables found in \"x\""))
    data.col <- data.col[1]
  }
  if (!valid.data.col(x,data.col))
    stop("Argument \"col\" is not a valid name of a data column in \"x\"") 

  # Choose with PCs to plot against the categorical variable.
  if (missing(pc.dims)) {
    pc.dims <- get.pc.cols(x)
    pc.dims <- pc.dims[seq(1,min(4,length(pc.dims)))]
  }
  if (!valid.pc.dims(x,pc.dims))
    stop("Argument \"pc.dims\" is not a valid selection of PCs in \"x\"")

  # Choose which colors to use for filling in the "violin" densities.
  n <- nlevels(x$data[[data.col]])
  if (missing(colors))
    colors <- pcaviz_colors_categorical(n)
  if (is.null(colors))
    colors <- rep("black",n)
  if (length(colors) < n)
    stop("Argument \"colors\" must have one entry for each level of \"col\"")

  # Check inputs "sorted" and "horizontal".
  if (!is.TorF(sorted))
    stop("Argument \"sorted\" should be TRUE or FALSE")
  if (!is.TorF(horizontal))
    stop("Argument \"horizontal\" should be TRUE or FALSE")
  
  # Check input "rank.fun".
  if (!is.function(rank.fun))
    stop("Argument \"rank.fun\" should be a function")
  
  # Check input "theme"
  if (!(inherits(theme,"gg") & inherits(theme,"theme")))
    stop("Argument \"theme\" is not a valid ggplot2 theme object")

  # Check input "violin.params".
  if (!is.list(violin.params))
    stop("Argument \"violin.params\" should be a list")
  
  # Create the violin plot, or arrangement of violin plots.
  # Generate a list object containing the arguments to plot_grid.
  n                <- length(pc.dims)
  plot.grid.args   <- rep(list(NULL),n)
  plot.grid.labels <- rep("",n)
  for (i in 1:n) {
    plot.grid.args[[i]] <-
      pcaviz_violin_helper(x,data.col,pc.dims[i],colors,sorted,rank.fun,
                           horizontal,theme,violin.params)
    plot.grid.labels[[i]] <- LETTERS[i]
  }
  
  # Return a single plot, or a grid of plots.
  if (length(pc.dims) == 1)
    return(plot.grid.args[[1]])
  else
    return(do.call("plot_grid",c(plot.grid.args,
                                 list(labels = plot.grid.labels),
                                 plot.grid.params)))
}

# Function for generating one plot, or an arrangement of plots, from a
# single pcaviz object (e.g., multiple plots showing combinations of
# PCs in the horizontal and vertical axes).
pcaviz_ggplot <-
  function (x, coords,
            arrange.coords = c("all.pairs","each.vs.pc1","consecutive.pairs"),
            plotly = FALSE, draw.points = plotly, label, group,
            color, shape, colors, shapes, abbreviated.label, 
            group.summary.fun = pcaviz_summary_default,
            group.summary.labels = TRUE, draw.pc.axes, hide.xy.axes,
            include.with.pc.axes, draw.linear.fit, draw.confint,
            confint.level = 0.95, show.r.squared, preserve.scale,
            overlay = NULL,
            geom.point.params = if (plotly)
              list(size = 3,stroke = 0,na.rm = TRUE)
            else
              list(size = 2,stroke = 1,na.rm = TRUE),
            geom.text.params = list(size = 3,fontface = "plain",na.rm = TRUE,
              alpha = 1),
            geom.point.summary.params = list(shape = 19,stroke = 1,size = 6,
                show.legend = FALSE,alpha = 0.8),
            geom.text.summary.params = list(size = 3.25,fontface = "plain",
              color = "black",show.legend = FALSE,alpha = 0.8),
            geom.segment.pc.axes = list(color = "black",linetype = "solid",
              arrow = arrow(length = unit(5,"pt"),ends = "both",type = "open"),
              size = 0.3),
            geom.text.pc.axes = list(hjust = "left",size = 3),
            geom.abline.params.linearfit = list(color = "dimgray",
              linetype = "dashed"),
            geom.abline.params.confint = list(color = "dimgray",
              linetype = "dotted"),
            scale.pc.axes = 0.6, theme = theme_pcaviz(), show.legend = TRUE,
            plot.title, plot.grid.params = list(), tooltip = NULL,
            plotly.file = NULL, ...) {

  # Check input "x".
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                          "object of class \"%s\""),class(x)))

  # Check input argument "plotly".
  if (!is.TorF(plotly))
    stop("Argument \"plotly\" should be TRUE or FALSE")
      
  # Choose which variables to plot in the horizontal (x) and vertical
  # (y) axes. The default is to plot the first PC against the second PC.
  if (missing(coords))
    coords <- get.pc.cols(x)[1:2]
  if (!(valid.data.cols(x,coords) &
        is.data.type(x,coords,c("continuous","pc"))))
    stop(paste("Argument \"coords\" does not specify valid names of",
               "PC or data columns in \"x\""))
  if (length(coords) < 2)
    stop(paste("Argument \"coords\" must specify at least 2 PCs or",
               "data columns in \"x\""))
  if (plotly & length(coords) != 2)
    stop(paste("For plotly graphcs, argument \"coords\" must specify",
               "exactly 2 PCs or data columns in \"x\""))

  # Determine which variables (PCs) are plotted against each
  # other. (This setting really only matters when the number of
  # selected PCs is three or greater.
  arrange.coords <- match.arg(arrange.coords)      
  
  # Check input "draw.points".
  if (!is.TorF(draw.points))
    stop("Argument \"draw.points\" should be TRUE or FALSE")
  if (plotly & !draw.points)
    stop("Argument \"draw.points\" must be TRUE when plotly = TRUE")
  
  # Choose which variable to use for the sample labels. If draw.points
  # is set to TRUE, don't choose a variable unless supplied by the
  # user. If draw.points is set to FALSE, choose the first factor
  # ("categorical") column column, and if no factor column is
  # available, choose the first character column.
  if (missing(label))
    if (plotly)
        label <- NULL
    else {
      if (draw.points)
        label <- NULL
      else {
        i <- c(which(sapply(x$data,is.factor)),
               which(sapply(x$data,is.character)))
        i <- names(x$data)[i]
        if (length(i) == 0)
          label <- NULL
        else
          label <- i[1]
    }
  }
  if (!is.null(label)) {
    if (plotly)
      stop("label = TRUE is not functional for plotly graphics")
    if (!valid.data.col(x,label))
      stop("Argument \"label\" is not a valid name of a data column in \"x\"") 
  }
  
  # Decide whether to plot the full labels or the abbreviated
  # labels. Show the abbreviated labels instead of the original labels
  # if either: (1) abbreviated labels are provided in the pcaviz
  # object, or (2) the label is a categorical variable and the level
  # names are more than 4 characters long.
  alabel <- paste0(label,".abbrv")
  if (missing(abbreviated.label)) {
    abbreviated.label <- FALSE
    if (!is.null(label))  {
      if (is.element(alabel,names(x$data)))
        abbreviated.label <- TRUE
      else if (is.data.type(x,label,"categorical"))
        if (max(nchar(levels(x$data[[label]]))) > 4)
          abbreviated.label <- TRUE
    }
  }
  if (!is.TorF(abbreviated.label))
    stop("Argument \"abbreviated.label\" should be TRUE or FALSE")
  if (is.null(label) & abbreviated.label)
    stop("Cannot set \"abbreviated.label = TRUE\" when \"label\" is NULL")

  # If necessary, create a new data column containing an abbreviated
  # variable, and summarize the abbreviations used in output to the
  # user.
  if (abbreviated.label & !is.element(alabel,names(x$data)))
    x <- pcaviz_abbreviate_var(x,label)
  if (!abbreviated.label)
    alabel <- NULL
  else {
    cat("Abbreviations used in plot:\n")
    abbrv.tbl <- data.frame(levels(x$data[[label]]),levels(x$data[[alabel]]))
    names(abbrv.tbl) <- c(label,alabel)
    print(abbrv.tbl,row.names = FALSE,right = FALSE)
  }
  
  # Choose which categorical variable to use for compiling and
  # plotting summary statistics. The default is to use the same
  # variable used to plot the labels, if it is a categorical variable,
  # otherwise take the first categorical variable, if available.
  if (missing(group)) {
    if (plotly)
      group <- NULL   
    else if (!is.null(alabel) & is.data.type(x,alabel,"categorical"))
      group <- alabel
    else if (!is.null(label) & is.data.type(x,label,"categorical"))
      group <- label
    else {
      i <- get.cols.by.type(x,"categorical")
      if (length(i) > 0)
        group <- i[1]
      else
        group <- NULL
    }
  }
  if (!is.null(group))
    if (!(valid.data.col(x,group) & is.data.type(x,group,"categorical")))
      stop(paste("Argument \"group\" is not the name of a categorical",
                 "data column in \"x\""))

  # Choose which variable(s) to plot as different colors and/or
  # shapes.
  if (!missing(color) & missing(shape))
    shape <- NULL
  if (missing(color) & !missing(shape))
    color <- NULL
  if (missing(color) & missing(shape)) {
    if (!is.null(alabel) & is.data.type(x,alabel,"categorical")) {

      # When "alabel" is specified and it is a categorical variable,
      # vary color and shape according to the abbreviated label.
      color <- alabel
      shape <- alabel
    } else if (!is.null(label) & is.data.type(x,label,"categorical")) {

      # When "label" is specified and it is a categorical variable,
      # vary color and shape according to label.
      color <- label
      shape <- label
    } else if (!is.null(group)) {
    
      # When "group" is specified, vary color and shape according to the
      # group variable.
      color <- group
      shape <- group
    } else if (!draw.points) {

      # When only labels are drawn, set the color to the first categorical
      # variable, and if unavailable, the first continuous variable.
      i <- c(get.cols.by.type(x,"categorical"),
             get.cols.by.type(x,"continuous"))
      if (length(i) == 0)
        color <- NULL
      else
        color <- i[1]
      shape <- NULL
    } else if (length(get.cols.by.type(x,"categorical")) > 0) {
    
      # When points are drawn, and there is at least one categorical
      # variable, vary shape and color according to the first
      # categorical variable that is encountered.
      color <- get.cols.by.type(x,"categorical")[1]
      shape <- color
    } else if (length(get.cols.by.type(x,"continuous")) > 0) {

      # When points are drawn, and there is at least one continuous
      # variable, varycolor according to the first continuous variable
      # that is encountered.
      color <- get.cols.by.type(x,"continuous")[1]
      shape <- NULL
    }
  } 
  if (!is.null(color))
    if (!(valid.data.col(x,color) &
          is.data.type(x,group,c("categorical","continuous","pc"))))
      stop(paste("Argument \"color\" is not the name of a categorical",
                 "or continuous data column in \"x\""))
  if (!is.null(shape))
    if (!(valid.data.col(x,shape) & is.data.type(x,group,"categorical")))
      stop(paste("Argument \"shape\" is not the name of a categorical",
                 "data column in \"x\""))

  # At least one of points, labels and group summaries must be drawn,
  # otherwise there will be nothing in the plot to draw.
  if (!draw.points & is.null(label) & is.null(group))
    stop(paste("At least one of points, labels and group summaries",
               "must be chosen for plot"))

  # To avoid plotting inconsistencies with unused levels, drop all
  # unused levels in any factors included in the plots.
  cols         <- unique(c(coords,label,alabel,group,color,shape))
  x$data[cols] <- droplevels(x$data[cols])

  # Determine how to plot the colors and/or shapes if these have not
  # been selected by the user.
  if (missing(colors)) {
    if (is.null(color))
      colors <- NULL
    else if (is.data.type(x,color,"categorical"))
      colors <- pcaviz_colors_categorical(nlevels(x$data[[color]]))
    else if (is.data.type(x,color,c("continuous","pc")))
      colors <- pcaviz_colors_continuous()
    else
      colors <- NULL
  }
  if (missing(shapes)) {
    if (is.null(shape))
      shapes <- NULL
    else
      shapes <- pcaviz_shapes(nlevels(x$data[[shape]]))
  }

  # Create the colors and shapes Scale objects.
  if (!is.null(colors)) {
    if (is.data.type(x,color,"categorical")) {
      colors.scale <- scale_color_manual(values = colors)
      if (!is.null(alabel))
        if (color == label | color == alabel)

          # Treat the special case when the label and color are
          # determined by the same (categorical) variable, and the
          # labels are abbreviated; in this case we can add the
          # abbreviations to the legend.
          colors.scale <- scale_color_manual(values = colors,
                                             labels = sprintf("%s (%s)",
                                               levels(x$data[[label]]),
                                               levels(x$data[[alabel]])))
    } else 
      colors.scale <- scale_color_gradientn(colors = colors)
  } else
    colors.scale <- NULL
  if (!is.null(shapes))
    shapes.scale <- scale_shape_manual(values = shapes)
  else
    shapes.scale <- NULL

  # Check input argument "group.summary.fun".
  if (!is.function(group.summary.fun))
    stop("Argument \"group.summary.fun\" should be a function")
  
  # Check input argument "group.summary.labels".
  if (!is.TorF(group.summary.labels))
    stop("Argument \"group.summary.labels\" should be TRUE or FALSE")

  # Decide whether to draw the transformed PC axes.
  if (missing(draw.pc.axes))
    draw.pc.axes <- valid.pc.dims(x,coords) & x$transformed.pcs
  if (!is.TorF(draw.pc.axes))
    stop("Argument \"draw.pc.axes\" should be TRUE or FALSE")
  if (draw.pc.axes & !valid.pc.dims(x,coords))
    stop("Setting draw.pc.axes = TRUE requires PC columns for \"coords\"")
  
  # Decide whether to draw the horizontal (x) and vertical (y) axes.
  if (missing(hide.xy.axes))
    hide.xy.axes <- draw.pc.axes
  if (!is.TorF(hide.xy.axes))
    stop("Argument \"hide.xy.axes\" should be TRUE or FALSE")

  # Decide whether to add some additional information (e.g., to the PC
  # axes (e.g., proportion of variance explained).
  if (missing(include.with.pc.axes)) {
    if (length(coords) <= 2 &
        valid.pc.dims(x,coords) &
        !any(is.na(x$sdev))) {
      if (any(is.na(x$var))) {
        message("Variance explained will be added to the axis labels.")
        include.with.pc.axes <- "var"
      } else {
        message(paste("Proportion of variance explained (PVE) will be added",
                      "to the axis labels."))
        include.with.pc.axes <- "pve"
      }
    } else
      include.with.pc.axes <- "none"
  }
  if (!(is.character(include.with.pc.axes) &
        is.element(include.with.pc.axes,c("none","var","pve","eigenvalue"))))
    stop(paste("Argument \"include.with.pc.axes\" must be one of",
               "\"none\", \"var\", \"pve\" or \"eigenvalue\""))
  if (include.with.pc.axes != "none" & any(is.na(x$sdev)))
    stop(paste("include.with.pc.axes = \"",include.with.pc.axes,
               "\" is not a valid setting because standard devations ",
               "are not provided by pcaviz object \"x\"",sep = ""))
  if (include.with.pc.axes == "pve" & any(is.na(x$var)))
    stop(paste("include.with.pc.axes = \"pve\" is not a valid setting because",
               "total variance is not provided by pcaviz object \"x\""))
  
  # Decide whether to draw the linear fit and/or confidence interval.
  if (!missing(draw.linear.fit) & missing(draw.confint))
    draw.confint <- draw.linear.fit
  else if (missing(draw.linear.fit) & !missing(draw.confint))
    draw.linear.fit <- draw.confint
  else if (missing(draw.linear.fit) & missing(draw.confint)) {
    if (sum(x$data.coltypes[coords] == "pc") == 1) {
      draw.linear.fit <- TRUE
      draw.confint    <- TRUE
    }
    else {
      draw.linear.fit <- FALSE
      draw.confint    <- FALSE
    }
  }

  # Decide whether to show the estimate of the proportion of variance
  # in y explained x; see help(summary.lm) for more details about this.
  if (missing(show.r.squared))
    show.r.squared <- draw.linear.fit | draw.confint
  if (!is.TorF(show.r.squared))
    stop("Argument \"show.r.squared\" should be TRUE or FALSE")

  # Decide whether to preserve the scale of the x and y axes if this
  # option is not set by the user. The heuristic is that if the span
  # of the data samples each dimension is on roughly the same scale
  # (here we use a factor of 1.5), then it is reasonable to fix the
  # scale.
  if (missing(preserve.scale)) {
    spans <- sapply(x$data[coords],function (x) diff(range(x,na.rm = TRUE)))
    preserve.scale <-
      is.null(overlay) & (max(c(spans/spans[1],spans[1]/spans)) < 1.5)
  }
  if (!is.TorF(preserve.scale))
    stop("Argument \"preserve.scale\" should be TRUE or FALSE")

  # Check input argument "overlay".
  if (!is.null(overlay))
    if (!is.function(overlay))
      stop("Argument \"overlay\" should be a function")

  # Check arguments giving lists of parameters passed to geom_ objects.
  args.lst <- c(quote(geom.point.params),
                quote(geom.text.params),
                quote(geom.point.summary.params),
                quote(geom.text.summary.params),
                quote(geom.segment.pc.axes),
                quote(geom.text.pc.axes),
                quote(geom.abline.params.linearfit),
                quote(geom.abline.params.confint))
  if (!all(sapply(args.lst,function (x) is.list(eval(x)))))
    stop(paste(strwrap(paste("These arguments should be lists:",
                             paste(args.lst,collapse = ", "))),
               collapse = "\n"))
  
  # Check input argument "scale.pc.axes".
  if (!(is.numeric(scale.pc.axes) &
        length(scale.pc.axes) == 1 &
        scale.pc.axes > 0))
    stop("Argument \"scale.pc.axes\" should be a positive scalar.")
  
  # Determine whether to auto-generate the panel titles.
  if (missing(plot.title)) {
    generate.panel.titles <- TRUE
    plot.title            <- NULL
  }
  else
    generate.panel.titles <- FALSE
  if (!is.null(plot.title))
    if (!(is.character(plot.title) & length(plot.title) == 1))
      stop("Argument \"plot.title\" should be a character vector of length 1")

  # Get the number of co-ordinates to plot in the x and y axes.
  n <- length(coords)

  # Hide the legend if we are plotting multiple combinations of
  # co-ordinates.
  if (missing(show.legend))
    show.legend <- (n <= 2)
  if (!is.TorF(show.legend))
    stop("Argument \"show.legend\" should be TRUE or FALSE")
  if (plotly & !show.legend)
    stop("show.legend = FALSE is not implemented when plotly = TRUE")
  
  # Check argument "theme".
  if (!is.null(theme))
    if (!(inherits(theme,"theme") & inherits(theme,"gg")))
      stop("Argument \"theme\" must be a ggplot2 theme object;",
           "e.g., output of theme_minimal()")

  # Check input argument "plotly.file".
  if (!plotly & !is.null(plotly.file))
    stop("Argument \"plotly.file\" only applicable when plotly = TRUE")
  
  # Check input argument "tooltip".
  if (!plotly & !is.null(tooltip))
    stop("Argument \"tooltip\" only applicable when plotly = TRUE")
  if (!is.null(tooltip))
    if (!valid.data.cols(x,tooltip))
      stop(paste("Argument \"tooltip\" does not specify valid names of PCs",
                 "or data columns in \"x\"") )
  
  # Generate two list objects containing the arguments to plot_grid:
  # one for the ggplot objects, and another for the plot labels.
  plot.grid.args     <- rep(list(NULL),(n-1)^2)
  plot.grid.labels   <- rep("",(n-1)^2)
  letters.for.labels <- LETTERS

  # Repeat for all combinations of the co-ordinates.
  k <- 0
  for (i in 1:(n-1))
    for (j in 2:n) {
      k <- k + 1
      if (i < j) {
        if (generate.panel.titles)
          plot.title <- paste(coords[i],"vs.",coords[j])
        plot.grid.args[[k]] <-
          pcaviz_ggplot_helper(x,coords[c(i,j)],draw.points,label,group,
                               color,shape,colors.scale,shapes.scale,alabel,
                               group.summary.fun,group.summary.labels,
                               draw.pc.axes,hide.xy.axes,include.with.pc.axes,
                               draw.linear.fit,draw.confint,confint.level,
                               show.r.squared,preserve.scale,overlay,
                               geom.point.params,geom.text.params,
                               geom.point.summary.params,
                               geom.text.summary.params,
                               geom.segment.pc.axes,geom.text.pc.axes,
                               geom.abline.params.linearfit,
                               geom.abline.params.confint,scale.pc.axes,
                               theme,show.legend,plot.title)
        plot.grid.labels[[k]] <- letters.for.labels[1]
        letters.for.labels    <- letters.for.labels[-1]
      }
    }

  if (plotly) {

    # Return a plotly graph.
    out <- suppressMessages(ggplotly(plot.grid.args[[1]],
                                     tooltip = paste0("var.",tooltip)))

    # Save the plotly graph to an HTML file, if requested.
    if (!is.null(plotly.file))
      saveWidget(out,plotly.file,selfcontained = TRUE)
    return(out)
  } else {
  
    # Return a single (non-interactive) plot, or a grid of plots.
    if (n == 2)
      return(plot.grid.args[[1]])
    else
      return(do.call("plot_grid",c(plot.grid.args,
                                   list(labels = plot.grid.labels),
                                   plot.grid.params)))
  }
}

# Plots variance (or a related quantity) of each PC. The design of
# this function is based on function ggscreenplot from the rsvd
# package.
pcaviz_screeplot <-
    function (x, type = c("var","pve","eigenvalue"),
              geom.point.params = list(size = 2,stroke = 1,color = "black",
                                       na.rm = TRUE),
              geom.line.params = list(size = 1,color = "black",na.rm = TRUE),
              theme = theme_pcaviz(), plot.title = NULL, ...) {
        
  # Check input "x".
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                          "object of class \"%s\""),class(x)))
      
  # Determine which quantity to show in the vertical axis ("type").
  type <- match.arg(type)

  # Check input "theme"
  if (!(inherits(theme,"gg") & inherits(theme,"theme")))
    stop("Argument \"theme\" is not a valid ggplot2 theme object")

  # Check input "plot.title".
  if (!is.null(plot.title))
    if (!(is.character(plot.title) & length(plot.title) == 1))
      stop("Argument \"plot.title\" should be a character vector of length 1")
  
  # Get the quantity to show in the vertical axis.
  if (any(is.na(x$sdev)))
    stop(paste("Cannot create scree plot because standard devations",
               "are not provided by pcaviz object \"x\""))
  out     <- get.screeplot.quantity(x,type)
  y       <- out$y
  y.label <- out$y.label

  # Initialize the ggplot object.
  pc.labels <- get.pc.cols(x)
  n         <- length(pc.labels)
  plot.data <- data.frame(PC = 1:n,y = y)
  out       <- ggplot(plot.data,aes_string(x = "PC",y = "y"),
                      environment = environment())

  # Draw the points and lines.
  out <- out +
         do.call("geom_point",geom.point.params) +
         do.call("geom_line",geom.line.params)

  # Set the horizontal and vertical labels.
  out <- out +
         scale_x_continuous(breaks = 1:n,labels = pc.labels) +
         xlab("") + ylab(y.label)

  # Set the limits of the vertical axis.
  if (type == "pve")
    out <- out + ylim(0:1)
  else
    out <- out + ylim(c(0,NA))
  
  # Set the theme.
  if (!is.null(theme))
    out <- out + theme

  # Add the plot title if one is provided.
  if (!is.null(plot.title))
    out <- out + labs(title = plot.title)
  
  # Return the ggplot object.
  return(out)
}

# Plots the absolute values of the loadings (eigenvectors), with
# several plotting options, including grouping variables by color.
pcaviz_loadingsplot <-
  function (x, pc.dim, color, colors, add.labels, min.rank = 0,
            gap = 0, geom.point.params = list(size = 1,na.rm = TRUE),
            theme = theme_pcaviz()) {

  # Check input "x".
  if (!inherits(x,"pcaviz"))
    stop(gettextf(paste("Argument \"x\" should be a pcaviz object; got an",
                          "object of class \"%s\""),class(x)))
  if (is.null(x$rotation))
    stop(paste("Cannot create loadings plot because rotation matrix",
               "is not provided by pcaviz object \"x\""))

  # Choose which PC to plot in the vertical axis. in the horizontal
  # (x) and vertical (y) axes. The default is to plot the first PC.
  if (missing(pc.dim))
    pc.dim <- get.pc.cols(x)[1]
  if (!(valid.pc.dims(x,pc.dim) & length(pc.dim) == 1))
    stop("Argument \"pc.dim\" is not a valid PC name in \"x\"")

  # Choose which column to show as different colors. The default is to
  # set the color to the first factor in the "rotation" table, when
  # one is available.
  if (missing(color)) {
    i <- which(sapply(x$rotation,is.factor))
    if (length(i) == 0)
      color <- NULL
    else {
      i     <- i[1]
      color <- names(x$rotation)[i]
    }
  }
  if (!is.null(color)) {
    if (!any(color == names(x$rotation)))
      stop(paste("Argument \"color\" should be a name of a factor column",
                 "in \"x$rotation\""))
    if (!is.factor(x$rotation[[color]]))
      stop(paste("Argument \"color\" should be a name of a factor column",
                 "in \"x$rotation\""))
  }
  
  # Determine how to plot the colors if these have not been selected
  # by the user.
  if (missing(colors)) {
    if (is.null(color))
      colors <- NULL
    else 
      colors <- pcaviz_colors_categorical(nlevels(x$rotation[[color]]))
  }

  # Create the colors Scale object.
  if (!is.null(colors))
    colors.scale <- scale_color_manual(values = colors)
  else
    colors.scale <- NULL
    
  # Determine whether to add labels to the points. The default is to
  # only add labels if there are at most 20 variables and the rows of
  # the x$rotation table have names.
  n <- nrow(x$rotation)
  if (missing(add.labels))
    add.labels <- (n <= 20 & !is.null(rownames(x$rotation)))
  if (!is.TorF(add.labels))
    stop("Argument \"add.labels\" should be TRUE or FALSE")

  # Check argument "min.rank".
  if (!(is.numeric(min.rank) & min.rank >=0 & min.rank <= 1))
    stop("Argument \"min.rank\" should be a number between 0 and 1")

  # Check argument "gap".
  if (gap > 0 & is.null(color))
    stop("Argument \"gap\" is only relevant when \"color\" is specified")
  
  # Check argument "geom.point.params".
  if (!is.list(geom.point.params))
    stop("Argument \"geom.pointparams\ should be a list")
  
  # Check argument "theme".
  if (!is.null(theme))
    if (!(inherits(theme,"theme") & inherits(theme,"gg")))
      stop("Argument \"theme\" must be a ggplot2 theme object;",
           "e.g., output of theme_minimal()")
  
  # Create the plot.
  return(pcaviz_loadingsplot_helper(x,pc.dim,color,colors.scale,add.labels,
                                    min.rank,gap,geom.point.params,theme))
}

# Define S3 methods for plotting pcaviz objects.
plot.pcaviz <- pcaviz_ggplot

# Define S3 method for pcaviz scree plots.
screeplot.pcaviz <- pcaviz_screeplot

# MATRIX TRANSFORMATION FUNCTIONS
# ======================================================================
# Here we define a suite of functions that apply basic transformations
# of points defined in a multi-dimensional Euclidean space. The points
# are stored as rows in a matrix. The functions defined here are:
#
#   function            transformation
#   --------------      ------------------------------------
#   matrix_rotate       counter-clockwise plane rotation
#   matrix_reflect      reflection along selected dim's about origin
#   matrix_scale        change of scale along selected dimensions
#   matrix_translate    translate points along selected dimensions
#   matrix_transform2d  apply sequence of 2-d transformations
#

# Internal function to convert degrees to radians. 
deg2rad <- function (x)
  x*pi/180

# Internal function to convert radians to degrees.
rad2deg <- function (x)
  x*180/pi

# Apply counter-clockwise plane rotation to row vectors in X; that is,
# return Y = X*R', where R is a plane rotation matrix. The rotation is
# applied to the first two dimensions by default.
matrix_rotate <-
  function (X, angle, dims = 1:2, units = c("degrees","radians")) {

  # Check input X.
  if (!is.matrix(X))
    stop("Argument \"X\" should be a matrix")
  n <- ncol(X)

  # Check input "angle".
  if (!(is.numeric(angle) & length(angle) == 1))
    stop("Argument \"angle\" should be a numeric scalar")

  # Process and check input "dims".
  if (length(dims) != 2)
    stop("Argument \"dims\" should be a vector of length 2")
  
  # Process and check input "units".
  units <- match.arg(units)
  if (units == "degrees")
    angle <- deg2rad(angle)

  # Construct the plane rotation matrix.
  a            <- angle
  R            <- diag(n)
  dimnames(R)  <- list(colnames(X),colnames(X))
  R[dims,dims] <- rbind(c(cos(a),-sin(a)),
                        c(sin(a),cos(a)))

  # Return the matrix of rotated points.
  return(X %*% t(R))
}

# Reflect points about the origin, individually along each selected
# dimension.
matrix_reflect <- function (X, dims) {

  # Check input X.
  if (!is.matrix(X))
    stop("Argument \"X\" should be a matrix")
  n <- ncol(X)

  # Reflect the points along the selected dimensions.
  X[,dims] <- (-X[,dims])
  return(X)
}

# Scale dimensions (columns) of matrix X; that is, return Y = X*S,
# where S is the scaling matrix. Note that it is also possible to set
# negative scale, which is equivalent to a reflection about the origin.
matrix_scale <- function (X, scale, dims) {

  # Check input X.
  if (!is.matrix(X))
    stop("Argument \"X\" should be a matrix")
  n <- ncol(X)

  # Check input "scale".
  if (!(is.numeric(scale) & length(scale) == length(dims)))
    stop(paste("Argument \"scale\" should be a numeric vector of the",
               "same length as \"dims\""))
  
  # Construct the scaling matrix.
  S           <- rep(1,n)
  names(S)    <- colnames(X)
  S[dims]     <- scale
  S           <- diag(S)
  dimnames(S) <- list(colnames(X),colnames(X))

  # Return the matrix of scaled points.
  return(X %*% S)
}

# Translate row vectors in X by a; that is return Y = X + A, where A
# is the translation matrix.
matrix_translate <- function (X, a, dims) {

  # Check input "X".
  if (!is.matrix(X))
    stop("Argument \"X\" should be a matrix")
  n <- ncol(X)

  # Check input "a".
  if (!(is.numeric(a) & length(a) == length(dims)))
    stop(paste("Argument \"a\" should be a numeric vector of the same",
               "length as \"dims\""))

  # Construct the translation matrix.
  A        <- rep(0,n)
  names(A) <- colnames(X)
  A[dims]  <- a
  A        <- matrix(A,nrow(X),ncol(X),byrow = TRUE)

  # Return the matrix of translated points.
  return(X + A)
}

# This is a convenience function that applies a sequence of 2-d
# transformations in the following order: (1) counter-clockwise
# rotation, (2) reflection, (3) scaling and (4) translation. Note that
# exactly two dimensions must be selected. Also note that the checks
# for some input arguments are skipped since it is assumed that they
# will be checked by matrix_rotate, etc.
matrix_transform2d <-
    function (X, dims = 1:2, angle = 0,reflect.x = FALSE, reflect.y = FALSE,
              scale = c(1,1),a = c(0,0), units = c("degrees","radians")) {
        
  # Process input "units".
  units <- match.arg(units)

  # Check inputs "reflect.x" and "reflect.y".
  if (!is.TorF(reflect.x))
    stop("Argument \"reflect.x\" should be TRUE or FALSE")
  if (!is.TorF(reflect.y))
    stop("Argument \"reflect.y\" should be TRUE or FALSE")
        
  # Apply: (1) counter-clockwise rotation, (2) reflection along first
  # dimension, (3) eflection along second dimension, (4) scaling, and
  # (5) translation.
  X <- matrix_rotate(X,angle,dims,units)
  if (reflect.x)
    X <- matrix_reflect(X,dims[1])
  if (reflect.y)
    X <- matrix_reflect(X,dims[2])
  X <- matrix_scale(X,scale,dims)
  X <- matrix_translate(X,a,dims) 
  return(X)
}

# ADDITIONAL PCAVIZ FUNCTIONS HIDDEN FROM USER
# ======================================================================
# Below are some additional functions used in the PCAviz package that
# are not directly accessible by the user. These functions typically
# forego input argument checks because we assume that the checks have
# already been taken care of by the function that calls them.

# Return TRUE if x is TRUE or FALSE.
is.TorF <- function (x)
  is.logical(x) & length(x) == 1

# Return the levels for the data.coltype field in the "pcaviz" class.
data.coltypes.levels <- function()
  c("pc","categorical","continuous","other")

# Return the names of the columns of a given type.
get.cols.by.type <- function (x, which.coltype)
  names(x$data.coltypes)[x$data.coltypes == which.coltype]

# Return the names of the PC columns.
get.pc.cols <- function (x)
  get.cols.by.type(x,"pc")

# Returns TRUE if input "cols" is a valid specification of data
# columns in pcaviz object x.
valid.data.cols <- function (x, cols) {
  if (length(cols) == 0)
    out <- FALSE
  else {
    out <- FALSE
    if (is.character(cols) & length(cols) > 0)
      if (all(is.element(cols,names(x$data))))
        out <- TRUE
  }
  return(out)
}

# Returns TRUE if input "col" is a valid specification of a data
# column in pcaviz object x.
valid.data.col <- function (x, col)
  valid.data.cols(x,col) & length(col) == 1


# Returns TRUE if input "dims" is a valid specification of PC dimensions
# in pcaviz object x.
valid.pc.dims <- function (x, dims) {
  out <- FALSE
  if (valid.data.cols(x,dims))
    if (all(x$data.coltypes[dims] == "pc"))
      out <- TRUE
  return(out)
}

# Returns TRUE if all columns ("cols") are of the specified data
# type(s), as specified by input "data.types".
is.data.type <- function (x, cols, data.types)
  all(is.element(x$data.coltypes[cols],data.types))

# Get the quantity to show in the vertical axis of the scree plot.
get.screeplot.quantity <- function (x, type) {
  if (type == "var") {
    y       <- x$sdev^2
    y.label <- "variance explained"
  } else if (type == "pve") {
    if (any(is.na(x$var)))
      stop(paste("Cannot plot proportion of variance explained because",
                 "total variance is not provided by pcaviz object \"x\""))
    y       <- x$sdev^2/x$var
    y.label <- "proportion of variance explained"
  } else if (type == "eigenvalue") {
    y       <- x$sdev^2
    y.label <- "eigenvalue"
  } else
    stop("Value for argument \"type\" is unsupported")
  return(list(y = y,y.label = y.label))
}

# Add a single column y to x$data.
pcaviz.add.col <- function (x, y, y.name, y.coltype) {

  # Add the abbreviated variable to the data frame.
  cols                   <- c(names(x$data),y.name)
  x$data                 <- cbind(x$data,data.frame(y))
  x$data.coltypes        <- factor(c(as.character(x$data.coltypes),
                                     as.character(y.coltype)),
                                   data.coltypes.levels())
  names(x$data)          <- cols
  names(x$data.coltypes) <- cols
  return(x)
}

# This function is used by "pcaviz_loadingsplot" to generate the
# loadings plots; this function should never be called directly by the
# user.
pcaviz_loadingsplot_helper <-
  function (x, pc.dim, color, colors.scale, add.labels, min.rank, gap,
            geom.point.params, ggtheme) {

  # Plot only the variable with rank greater than "min.rank", in which
  # the variables are ranked according to the absolute value of the PC
  # loading. 
  n          <- nrow(x$rotation)
  y          <- rank(abs(x$rotation[[pc.dim]]))/n
  rows       <- which(y >= min.rank)
  x$rotation <- x$rotation[rows,]

  # When gap > 0 and variables are plotted by "color", order the
  # variables by color.
  if (gap > 0 & !is.null(color)) {
    rows       <- order(x$rotation[[color]])
    x$rotation <- x$rotation[rows,]
  }
  
  # Determine the horizontal positions of the variables.
  n   <- nrow(x$rotation)
  pos <- 1:n
  if (!is.null(color)) {
    group <- x$rotation[[color]]
    if (nlevels(group) > 1) {
        r <- 0
        for (i in levels(group)) {
          j      <- which(group == i)
          m      <- length(j)
          pos[j] <- r + 1:m
          r      <- r + m + gap
        }
    } 
  }
      
  # Compile the plot data.
  plot.data           <- cbind(data.frame(x = pos),
                               x$rotation[c(pc.dim,color)])
  plot.data[[pc.dim]] <- abs(plot.data[[pc.dim]])
      
  # Initialize the ggplot object.
  out <- ggplot(plot.data,aes_string(x = "x",y = pc.dim,color = color),
                environment = environment())

  # Draw the points.
  out <- out + do.call("geom_point",geom.point.params)
    
  # Specify the colors.
  if (!is.null(colors.scale))
    out <- out + colors.scale

  # Label the points, if requested.
  if (add.labels) {
    out <- out + scale_x_continuous(breaks = 1:n,
                                    labels = rownames(x$rotation))
  } else
    out <- out + scale_x_continuous(breaks = NULL)

  # Adjust the vertical axis limits, but only if all the data are
  # being shown.
  if (min.rank == 0)
    out <- out + ylim(c(0,NA))
  
  # Specify the axis titles.
  out <- out + labs(x = "variable",
                    y = sprintf("abs(%s)",pc.dim))

  # Set the theme.
  if (!is.null(ggtheme))
    out <- out + ggtheme

  # Adjust the labels, if necessary.
  if (add.labels)
    out <- out + theme(axis.ticks.x = element_blank(),
                       axis.text.x  = element_text(angle = 45,hjust = 1))
  
  # Return the ggplot object.
  return(out)
}

# This function is used by "pcaviz_violin" to generate the violin
# plots; this function should never be called directly by the user.
pcaviz_violin_helper <- function (a, x, y, colors, sorted, rank.fun,
                                  horizontal, theme, violin.params) {

  # Get the columns of the data frame that are used to draw the plot.
  plot.data <- a$data[c(x,y)]

  # Sort the factor levels by the selected response variable, if this
  # is requested, so that the categories are drawn in increasing order
  # of the response variable. Ideally, we would add some exception
  # handling in the call to "tapply" (e.g., using "tryCatch") in order
  # to provide the user with a more informative message when an error
  # occurs. See "subset.pcaviz" for an example of this.
  if (sorted)
    plot.data[[x]] <-
      factor(as.character(plot.data[[x]]),
             names(sort(tapply(plot.data[[y]],plot.data[[x]],rank.fun))))

  # Initialize the ggplot object.
  out <- ggplot(data        = plot.data,
                mapping     = aes_string(x = x,y = y,fill = x),
                environment = environment())

  # Draw the violin densities. Ideally, we would add some exception
  # handling here (e.g., using "tryCatch") in order to provide the
  # user with a more informative message when an error occurs. See
  # "subset.pcaviz" for an example of this.
  out <- out + do.call("geom_violin",violin.params)
  
  # Specify the colors. Ideally, we would add some exception
  # handling here (e.g., using "tryCatch") in order to provide the
  # user with a more informative message when an error occurs. See
  # "subset.pcaviz" for an example of this.
  out <- out + scale_fill_manual(values = colors)
  
  # Set the theme. 
  if (!is.null(theme))
    out <- out + theme

  # Draw the violin densities vertically or horizontally.
  if (horizontal)
    out <- out + coord_flip()
  else
    out <- out +
      theme(axis.text.x  = element_text(angle = 45,hjust = 1),
            axis.ticks.x = element_line(size = 0))
  
  return(out)
}

# Return the two locations where the line segment defined by (x0, y0,
# x1, y1) intersects with the box defined by (xmin, ymin, xmax, ymax),
# or returns NULL if no intersections exist. Note this function
# assumes that there are either 0 intersections or 2 intersections
# (does not treat the situation in which there is exactly 1
# intersection). This function is used to draw the PC axes. 
# NOTE: This code is based on the Liang-Barsky algorithm.
collides <- function (x0, y0, x1, y1, xmin, ymin, xmax, ymax) {
  dx <- x1 - x0
  dy <- y1 - y0
  a  <- c(-dx,dx,-dy,dy)
  b  <- c(x0 - xmin,xmax - x0,y0 - ymin,ymax - y0)
  u1 <- (-Inf)
  u2 <- (+Inf)

  for (i in 1:4) {
    if (a[i] == 0 & b[i] < 0)
      return(NULL)
    r <- b[i]/a[i]
    if (a[i] < 0 & u1 < r)
      u1 <- r
    else if (a[i] > 0 & u2 > r)
      u2 <- r
  }
  
  if (u1 > u2 | u1 > 1 | u1 < 0)
    return(NULL)
  else
    return(list(x1 = x0 + u1*dx,
                y1 = y0 + u1*dy,
                x2 = x0 + u2*dx,
                y2 = y0 + u2*dy))
}

# Construct the data frame passed as the "data" argument to
# "geom_segment" to draw a segment on the line through points
# (x1,y1,x2,y2), but bounded by the extent of the data.
get_axis_data <- function (x, y, x1, y1, x2, y2, box.ratio = 0.6) {

  # Get the bounding box for the data points.
  xmin <- min(x,na.rm = TRUE)
  xmax <- max(x,na.rm = TRUE)
  ymin <- min(y,na.rm = TRUE)
  ymax <- max(y,na.rm = TRUE)
    
  # Create a box that delimits most of the data points (given by
  # vectors x and y), according to box.ratio.
  dx  <- xmax - xmin
  dy  <- ymax - ymin
  a   <- (1 - box.ratio)/2
  box <- list(x1 = xmin + a*dx,y1 = ymin + a*dy,
              x2 = xmax - a*dx,y2 = ymax - a*dy)

  # Determine the angle of the segment, in degrees.
  angle <- rad2deg(atan2(y2 - x1,x2 - x1))
  
  # Get the slope and direction of the axis segment.
  if (x1 == x2)
    pos.dir <- (y1 < y2)
  else {
    pos.dir <- (x1 < x2)
    m       <- (y2 - y1)/(x2 - x1)
    b       <- y1 - m*x1
  }

  # Create a new line segment for the axis with points falling outside
  # the box.
  if (x1 == x2)
    axis <- list(x1 = x1,y1 = ymin,x2 = x2,y2 = ymax)
  else {
    x1   <- 2*xmin
    x2   <- 2*xmax
    axis <- list(x1 = x1,y1 = m*x1 + b,x2 = x2,y2 = m*x2 + b)
  }

  # The axis segment is given by the two locations where the line
  # segment meets the box.
  out <- collides(axis$x1,axis$y1,axis$x2,axis$y2,
                  box$x1,box$y1,box$x2,box$y2)
  if (is.null(out))
    return(NULL)
  else {
    if (pos.dir)
      out <- with(out,list(x1 = x1,y1 = y1,x2 = x2,y2 = y2,angle = angle))
    else
      out <- with(out,list(x1 = x2,y1 = y2,x2 = x1,y2 = y1,angle = angle))
    return(as.data.frame(out))
  }
}

# This function is used by "pcaviz_ggplot" to generate plots from
# objects of class "pcaviz"; this function should never be called
# directly by the user.
pcaviz_ggplot_helper <-
  function (x, coords, draw.points, label, group, color, shape, colors.scale,
            shapes.scale, alabel, group.summary.fun, group.summary.labels,
            draw.pc.axes, hide.xy.axes, include.with.pc.axes, draw.linear.fit,
            draw.confint, confint.level, show.r.squared, preserve.scale,
            overlay, geom.point.params, geom.text.params,
            geom.point.summary.params, geom.text.summary.params,
            geom.segment.pc.axes, geom.text.pc.axes,
            geom.abline.params.linearfit, geom.abline.params.confint,
            scale.pc.axes, theme, show.legend, plot.title) {

  # Get the variables to plot on the horizontal (h) and vertical (v)
  # axes.
  h <- coords[1]
  v <- coords[2]

  # Set the default color for the points and labels, if necessary.
  if (is.null(color)) {
    default.color                      <- "royalblue"
    geom.text.params["color"]          <- default.color
    geom.point.params["color"]         <- default.color
    geom.point.summary.params["color"] <- default.color
  }
  
  # Determine which labels to plot.
  if (is.null(alabel))
    plot.label <- label
  else
    plot.label <- alabel
  
  # Fit a linear regression model y ~ x.
  fit.lm <- lm(x$data[coords[2:1]])
  
  # Initialize the ggplot object.
  aes.string.args        <- as.list(names(x$data))
  names(aes.string.args) <- paste0("var.",names(x$data))
  out <- ggplot(x$data,do.call("aes_string",aes.string.args),
                environment = environment())

  # Overlay points on additional graphic, if requested.
  if (!is.null(overlay))
    out <- overlay(out)
  
  # Draw the points, if requested.
  if (draw.points)
    out <- out + do.call("geom_point",
                         c(list(data = x$data,
                                aes_string(x = h,y = v,color = color,
                                           shape = shape),
                                inherit.aes = TRUE),geom.point.params))
  else
    out <- out + do.call("geom_point",
                         c(list(data = x$data,
                                aes_string(x = h,y = v,color = color),
                                shape = " ",inherit.aes = TRUE),
                           geom.point.params))
      
  # Draw the labels, if requested.
  if (!is.null(label)) 
    out <- out +
           do.call("geom_text",
                   c(list(data = x$data,
                          aes_string(x = h,y = v,label = plot.label,
                                     color = color),
                          show.legend = FALSE),geom.text.params))
    
  # Specify the colors.
  if (!is.null(colors.scale))
    out <- out + colors.scale

  # Specify the shapes.
  if (!is.null(shapes.scale))
    out <- out + shapes.scale

  # Add the group summary, if requested.
  if (!is.null(group)) {
      
    # Compile the summary data. Ideally, we would add some exception
    # handling in the call to "tapply" (e.g., using "tryCatch") in
    # order to provide the user with a more informative message when
    # an error occurs here. See "subset.pcaviz" for an example of
    # this.
    summary.data <- do.call("rbind",
                            by(x$data,x$data[[group]],group.summary.fun))
    summary.data <- cbind(data.frame(group.label = rownames(summary.data)),
                          summary.data)
    
    # Add the group summaries to the plot.
    out <- out + do.call("geom_point",
                         c(list(data = summary.data,
                                aes_string(x = h, y = v,color = color)),
                           geom.point.summary.params))

    # Add the group summary labels, if requested.
    if (group.summary.labels)
      out <- out + do.call("geom_text",
                           c(list(data = summary.data,
                                  aes_string(x = h,y = v,label = "group.label",
                                             color = color)),
                             geom.text.summary.params))
  }

  # Draw the PC axes, if requested.
  if (draw.pc.axes) {

    # Draw the two lines depicting the PC axes.
    pc.axis.h <-
      get_axis_data(x$data[[h]],x$data[[v]],x$basis["origin",h],
                    x$basis["origin",v],x$basis[h,h],x$basis[h,v],
                    scale.pc.axes)
    pc.axis.v <-
      get_axis_data(x$data[[h]],x$data[[v]],x$basis["origin",h],
                    x$basis["origin",v],x$basis[v,h],x$basis[v,v],
                    scale.pc.axes)
    if (is.null(pc.axis.h) | is.null(pc.axis.v))
      warning(paste("Unable to draw PC axes most likely because they do",
                    "not coincide with data points, or due to poor choice",
                    "of \"scale.pc.axes\""))
    else {
      pc.axes.data <-
        cbind(data.frame(label = paste0(" ",c(h,v)),stringsAsFactors = FALSE),
              rbind(pc.axis.h,pc.axis.v))

      # If requested, add the statistic (variance explained, PVE or
      # eigenvalue) to the axis labels.
      if (include.with.pc.axes != "none") {
        scree.data <- get.screeplot.quantity(x,include.with.pc.axes)
        if (include.with.pc.axes == "pve") {
          pc.axes.data[1,"label"] <-
            sprintf("%s (%0.2g%%)",pc.axes.data[1,"label"],100*scree.data$y[h])
          pc.axes.data[2,"label"] <- 
            sprintf("%s (%0.2g%%)",pc.axes.data[2,"label"],100*scree.data$y[v])
        } else {
          pc.axes.data[1,"label"] <-
            sprintf("%s (%0.3g)",pc.axes.data[1,"label"],scree.data$y[h])
          pc.axes.data[2,"label"] <- 
            sprintf("%s (%0.3g)",pc.axes.data[2,"label"],scree.data$y[v])
        }
      }
      
      # Draw the lines using geom_segment.
      out <- out + do.call("geom_segment",
                           c(list(data = pc.axes.data,
                                  mapping = aes_string(x = "x1",xend = "x2",
                                                       y = "y1",yend = "y2"),
                                  inherit.aes = FALSE),
                             geom.segment.pc.axes))

      # Label the PC axes. I rotate any labels that have an angle
      # between 90 and 270.
      pc.labels.data   <- pc.axes.data
      pc.labels.data$x <- pc.labels.data$x2
      pc.labels.data$y <- pc.labels.data$y2
      angle            <- pc.labels.data$angle %% 360
      rows             <- which(angle > 90 & angle < 270)
      pc.labels.data$angle[rows] <- pc.labels.data$angle[rows] + 180
      pc.labels.data$x[rows]     <- pc.labels.data$x1[rows]
      pc.labels.data$y[rows]     <- pc.labels.data$y1[rows]
      out <- out + do.call("geom_text",
                           c(list(data = pc.labels.data,
                                  mapping = aes_string(x = "x",y = "y",
                                                       label = "label",
                                                       angle = "angle"),
                                  inherit.aes = FALSE),
                             geom.text.pc.axes))
    }
  }

  # Draw the linear fit, if requested.
  if (draw.linear.fit) {
    out <- out + do.call("geom_abline",
                         c(list(intercept = coef(fit.lm)["(Intercept)"],
                                slope = coef(fit.lm)[h]),
                           geom.abline.params.linearfit))
  }

  # Draw the confidence intervals, if requested.
  if (draw.confint) {
    draw.line.for.confint <- function (p)
      do.call("geom_abline",
        c(list(intercept = qnorm(p,coef(fit.lm)["(Intercept)"],sigma(fit.lm)),
               slope = coef(fit.lm)[h]),
          geom.abline.params.confint))
    out <- out + draw.line.for.confint(0.5 - confint.level/2) +
                 draw.line.for.confint(0.5 + confint.level/2)
  }

  # Show the proportion of variance in y that is explained by x, if
  # requested.
  if (show.r.squared)
    plot.title <- sprintf("%s (r2 = %0.2g%%)",plot.title,
                          100*summary(fit.lm)$r.squared)

  # If requested, add the statistic (variance explained, PVE or
  # eigenvalue) to the axis labels.
  if (include.with.pc.axes != "none") {
    scree.data <- get.screeplot.quantity(x,include.with.pc.axes)
    if (include.with.pc.axes == "pve")
      out <- out + labs(x = sprintf("%s (%0.2g%%)",h,100*scree.data$y[h]),
                        y = sprintf("%s (%0.2g%%)",v,100*scree.data$y[v]))
    else
      out <- out + labs(x = sprintf("%s (%0.3g)",h,scree.data$y[h]),
                        y = sprintf("%s (%0.3g)",v,scree.data$y[v]))
  }
  
  # Preserve the scaling of the x and y axes, if requested.
  if (preserve.scale)
    out <- out + coord_fixed()
  
  # Add the plot title if one is provided.
  if (!is.null(plot.title))
    out <- out + labs(title = plot.title)

  # Show the legend, if requested.
  if (show.legend) {
    if (!draw.points)
      out <- out +
        guides(color = guide_legend(override.aes = list(shape = 20)))
  } else 
    out <- out + guides(color = "none",shape = "none")
      
  # Set the theme.
  if (!is.null(theme))
    out <- out + theme

  # Hide the horizontal (x) and vertical (y) axes, if requested.
  if (hide.xy.axes)
    out <- out +
      theme(axis.title = element_blank(),
            axis.line  = element_blank(),
            axis.ticks = element_blank(),
            axis.text  = element_blank())

  # Return the ggplot object.
  return(out)
}
