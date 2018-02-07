\name{pcaviz-plots}
\title{Plot pcaviz objects}
\alias{pcaviz-plots}
\alias{pcaviz-plot}
\alias{pcaviz-plotting}
\alias{pcaviz_colors_categorical}
\alias{pcaviz_colors_continuous}
\alias{pcaviz_shapes}
\alias{plot.pcaviz}
\alias{screeplot.pcaviz}
\alias{pcaviz_loadingsplot}
\alias{pcaviz_violin}
\alias{theme_pcaviz}
\alias{pcaviz_reduce_whitespace}
\alias{pcaviz_abbreviate_var}
\alias{pcaviz_summary_default}
\description{Interface for quickly creating visually compelling plots
  from Principal Component Analysis (PCA) and accompanying data. Default
  settings are designed to produce interesting plots with minimal user
  intervention; for example, running \code{plot(x)} on "pcaviz" object
  \code{x} with no additional arguments will produce a useful
  visualization of the PCA results in most cases (see also the examples
  accompanying this help page). The plotting interface is also 
  flexible and able accomodate fine-tuning for presentations,
  publications, webpages, etc.}
\usage{

  \method{plot}{pcaviz}(x, coords, plotly = FALSE, draw.points = plotly,
    label, group, color, shape, colors, shapes, abbreviated.label, 
    group.summary.fun = pcaviz_summary_default,
    group.summary.labels = TRUE, draw.pc.axes, hide.xy.axes,
    draw.linear.fit, draw.confint, confint.level = 0.95,
    show.r.squared, preserve.scale, overlay = NULL,
    geom.point.params = if (plotly)
      list(size = 3,stroke = 0,na.rm = TRUE)
    else
      list(size = 2,stroke = 1,na.rm = TRUE),
    geom.text.params = list(size = 3,fontface = "bold",na.rm = TRUE,
                            alpha = 1),
    geom.point.summary.params = list(shape = 19,stroke = 1,size = 10,
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
    scale.pc.axes = 0.8, theme = theme_pcaviz(), show.legend = TRUE,
    plot.title, plot.grid.params = list(), tooltip = NULL,
    plotly.file = NULL, \dots)

  \method{screeplot}{pcaviz}(x, type = c("var","pve","eigenvalue"),
    geom.point.params = list(size = 2,stroke = 1,color = "black",na.rm = TRUE),
    geom.line.params = list(size = 1,color = "black",na.rm = TRUE),
    theme = theme_pcaviz(), plot.title = NULL, \ldots)

  pcaviz_loadingsplot(x, pc.dim, color, colors, add.labels, min.rank = 0,
     gap = 0, geom.point.params = list(size = 1,na.rm = TRUE),
     theme = theme_pcaviz())

  pcaviz_violin(x, data.col, pc.dims, colors, sorted = TRUE,
    rank.fun = median, horizontal = FALSE, theme = theme_pcaviz(),
    violin.params = list(trim = FALSE, show.legend = FALSE),
    plot.grid.params = list())

  pcaviz_colors_categorical(n = 7)

  pcaviz_colors_continuous()

  pcaviz_shapes(n = 5)

  theme_pcaviz()

  pcaviz_reduce_whitespace(x, dims)

  pcaviz_abbreviate_var(x, col,
    abbrv.generic = function (x)
      abbreviate(abbreviate(gsub("[^[:alnum:] ]","",x)),
                  minlength = 2,named = FALSE,method = "both.sides"),
    abbrv.custom = function (x)
      countrycode(x,"country.name.en","iso2c",warn = FALSE))

  pcaviz_summary_default(x,
    stat.numeric = function (x) median(x,na.rm = TRUE))
}
\arguments{
  \item{x}{An object of class "pcaviz".}

  \item{coords}{The names of the PCs or data columns to plot in the
    horizontal and vertical axes. This may specify more than 2
    variables, in which case all combinations of the variables will be
    compared in different panels. For interactive graphics (\code{plotly
    = TRUE}), exactly two PCs or data columns must be chosen.}
  
  \item{col}{Name of PC or data column in pcaviz object \code{x}.}
  
  \item{data.col}{Name of categorical variable in \code{x} to plot
    against the PC(s) in \code{violin_plot}.}

  \item{dims}{Names of PCs in pcaviz object \code{x}.}

  \item{pc.dims}{Names of PCs in pcaviz object \code{x}.}

  \item{pc.dim}{Name of a PC in pcaviz object \code{x}.}
  
  \item{plotly}{If \code{TRUE}, an interactive graph will be generated
    using the \pkg{plotly} package. Note that some of the default plot
    settings are slightly different for interactive graphics.}

  \item{draw.points}{If \code{TRUE}, points are drawn by default instead
    of labels. (Both points and labels are drawn if labels is not
    \code{NULL}.)}

  \item{label}{The column to use for labeling the data points. The
    default is to label the data points according to the first
    categorical variable found, otherwise the first data column of type
    \code{character}.}

  \item{group}{Categorical variable (column of type \code{factor}) used
    to for compiling and plotting summary statistics.}

  \item{color}{PC or data column used to color the points or labels. For
    \code{pcaviz_loadingsplot}, this should be a column of the rotation
    matrix that is a factor (i.e., categorical). Otherwise, this may be
    a categorical or continuous variable.}

  \item{shape}{PC or data column used to determine shape of each point
    drawn. This option has no effect when \code{draw.points = FALSE}.}
  
  \item{colors}{Colors to use in the plots. This should be a vector of
    colors which is passed as input to
    \code{\link[ggplot2]{scale_fill_manual}},
    \code{\link[ggplot2]{scale_color_manual}} or
    \code{\link[ggplot2]{scale_color_gradientn}} depending on the
    variable used for \code{color}.}

  \item{shapes}{Shapes to use in the plots. This should be a should be a
    suitable vector specifying shapes which is passed as input to
    \code{\link[ggplot2]{scale_shape_manual}}.}

  \item{abbreviated.label}{Show the abbreviated labels instead of the
    original labels (as specified by argument \code{label}). This may
    either be provided by pcaviz object \code{x} in data column of the
    same name as \code{label}, but with an additional \code{".abbrv"}
    suffix, or if not provided, it will be automatically generated using
    function \code{pcaviz_abbreviate_var}. Abbreviated labels are
    drawn by default if the original label includes strings that are
    more than 4 characters long. When an abbreviated label is generated,
    for convenience the correspondence between full and abbreviated
    labels is printed.}

  \item{group.summary.fun}{Function used to generate summary
    statistics. It should be a function that takes a data frame as input
    and output as a data frame with the same number of columns, and only
    one row containing a summary of each of the table columns.}

  \item{group.summary.labels}{If \code{TRUE}, label the summary
    statistics by the \code{group} variable.}
  
  \item{draw.pc.axes}{If \code{TRUE}, draw the transformed PC axes.}

  \item{hide.xy.axes}{If \code{TRUE}, hide the axis lines and ticks.}

  \item{draw.linear.fit}{If \code{TRUE}, show the linear fit between the
    horizontal and vertical co-ordinates, as calculated by
    \code{\link[stats]{lm}}.}

  \item{draw.confint}{If \code{TRUE}, draw the confidence intervals of
    the linear fit based on the quantiles of the linear regression
    residual.}

  \item{confint.level}{A number between 0 and 1 specifying the quantile
    used to draw the confidence interval.}

  \item{show.r.squared}{If \code{TRUE}, append to the title the
    (unadjusted) proportion of variance in the vertical co-ordinate
    that is explained by the horizontal co-ordinate; see
    \code{\link[stats]{summary.lm}}.}

  \item{preserve.scale}{If \code{TRUE}, preserve the scaling of the
    axes using \code{\link[ggplot2]{coord_fixed}}.}

  \item{overlay}{A function that takes as input a ggplot object and
    returns a ggplot object, presumably a modified version of the input
    object. See Section "Overlays" below for more information about this
    argument.}

  \item{add.labels}{If \code{TRUE}, label variables along horizontal
    axis in \code{pcaviz_loadingsplot}. The default is to add labels
    only if a small number (<20) of variables are plotted}

  \item{min.rank}{A number between 0 and 1. When \code{min.rank > 0},
    plot loadings for variables only if the rank of the loading is
    greater than \code{min.rank} fraction of the variables. For example,
    if \code{min.rank = 0.9}, only the top 10\% of loadings will be
    shown in the loadings plot.}

  \item{gap}{Only relevant when \code{color} argument is specified.
    When \code{gap > 0}, the variables are grouped according to
    \code{color} (if they aren't already grouped), and an amount of
    space equal to \code{gap} is added between each of the groups.}

  \item{sorted}{If \code{TRUE}, the levels of the factor (categorical
    variable) in the violin plots are sorted in increasing order
    according to the selected PC(s). Note that if multiple violin plots
    are generated, the order of the factor levels may be different in
    each plot.}

  \item{rank.fun}{Function used to determine ranking for sorting
    factor levels in violin plots.}

  \item{horizontal}{If \code{TRUE}, draw the violin densities
    horizontally rather than vertically.}

  \item{theme}{A \pkg{ggplot2} theme.}

  \item{geom.point.params}{Additional parameters passed to
    \code{\link[ggplot2]{geom_point}} to customize plotting of points.}

  \item{geom.line.params}{Additional parameters passed to
    \code{\link[ggplot2]{geom_line}} to customize plotting of plots.}
  
  \item{geom.text.params}{Additional parameters passed to
    \code{\link[ggplot2]{geom_text}} to customize plotting of labels.}

  \item{geom.point.summary.params}{Additional parameters passed to
    \code{\link[ggplot2]{geom_point}} to customize plotting of points
    from group summary data.}

  \item{geom.text.summary.params}{Additional parameters passed to
    \code{\link[ggplot2]{geom_text}} to customize plotting of labels for
    group summaries.}

  \item{geom.segment.pc.axes}{Additional parameters passed to
    \code{\link[ggplot2]{geom_segment}} to customize drawing of PC axis
    lines.}

  \item{geom.text.pc.axes}{Additional parameters passed to
    \code{\link[ggplot2]{geom_text}} to customize text accompanying
    PC axis ines.}
  
  \item{geom.abline.params.linearfit}{Additional parameters passed to
    \code{\link[ggplot2]{geom_abline}} to customize drawing of linear
    fit.}

  \item{geom.abline.params.confint}{Additional parameters passed to
    \code{\link[ggplot2]{geom_abline}} to customize drawing of
    confidence intervals.}

  \item{violin.params}{Additional parameters passed to
    \code{\link[ggplot2]{geom_violin}} to customize plotting of violin
    densities.}

  \item{plot.grid.params}{Additional parameters passed to
    \code{\link[cowplot]{plot_grid}} to customize arrangement of plots.}
  
  \item{scale.pc.axes}{Number between 0 and 1 determining length of PC
    axis lines. Numbers close to 0 are generally not recommended as they
    may produce incorrectly displayed axis lines.}

  \item{show.legend}{If \code{FALSE}, hide the plot legend(s) that are
    normally generated by \code{ggplot2}. See also
    \code{\link[ggplot2]{guides}} for more fine-grained customization of
    plot legends.}

  \item{plot.title}{Plot title used, overriding the automatically
    generated plot title.}

  \item{tooltip}{Character vector specifying which PCs or data columns
    to display in the \pkg{plotly} "tooltip". For more details, see
    \code{\link[plotly]{ggplotly}}.}

  \item{plotly.file}{If not \code{NULL}, save plotly graph to this HTML
    file; see \code{\link[htmlwidgets]{saveWidget} for details.}}
  
  \item{n}{Number of colors or shapes to return.}

  \item{abbrv.generic}{Generic abbreviation for character strings; see
  \code{\link{create_abbreviated_label}} for details.}

  \item{abbrv.custom}{Custom abbreviation for character strings; see
    \code{\link{create_abbreviated_label}} for details. The default
    custom abbreviation is to convert country names in English to
    2-letter codes based on the based on the ISO 3166 standard; see
    \code{\link[countrycode]{countrycode}} and
    \code{\link[countrycode]{countrycode_data}} for more information;
    including other options for generating country abbreviations.}

  \item{stat.numeric}{Summary statistic for numeric vectors used by
    \code{pcaviz_summary_default}. It should be a function that accepts
    a numeric vector and returns a scalar.}

  \item{type}{Quantity to show in the vertical axis of the scree plot:
    \code{"var"}, variance explained; \code{"pve"}, proportion of
    variance explained; \code{"eigenvalue"}, which is equal to the
    variance explained.}
  
  \item{\dots}{Additional input arguments to generic plotting functions
    that are not used here.}
}
\details{
  Function \code{pcaviz_violin} generates one or several "violin" plots
  for visualizing the relationship between PCs and a categorical
  variable.
  
  Functions \code{pcaviz_colors_categorical} and \code{pcaviz_shapes}
  define default colors and shapes used for plotting categorical
  variables (i.e., factors) in the PCAviz plotting interface. These can
  be used as the \code{values} argument to
  \code{\link[ggplot2]{scale_color_manual}} and
  \code{\link[ggplot2]{scale_shape_manual}}. The colors are based on
  \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)} and the
  "Color blindness" paper by Wong (2008). If more than 7 colors are
  requested, the colors are repeated; if more than 5 shapes are
  requested, the shapes are repeated. Since LCM(5,7) = 35, functions
  \code{pcaviz_colors_categorical} and \code{pcaviz_shapes} can be
  easily combined to create as many as 35 unique color-shape
  combinations.

  Function \code{pcaviz_colors_continuous} specifies the default color
  scheme (turquoise--red) used for plotting continuous variables (i.e.,
  numeric data columns). This can be used as the "colors" argument for
  the \code{\link[ggplot2]{scale_color_gradientn}} function. These colors
  are obtained from \code{viridis::plasma(8)}.
  
  Function \code{theme_pcaviz} defines the default ggplot2 theme. It is
  based on the \pkg{cowplot} default theme.

  Function \code{pcaviz_reduce_whitespace} heuristically rescales the PC
  data to reduce the amount of "whitespace" (that is, unused space) in a
  PCAviz plot. This function often helps to make the plots more visually
  attractive. It also has the added benefit of minimizing the amount of
  automatic scaling done by \pkg{ggplot2}, so that, for example, rotated
  PC axes are less likely to be distorted. Dimension \code{dims[1]} is
  always unscaled, and the remaining dimensions are scaled relative to
  the first.

  Function \code{pcaviz_summary_default} is the default function for
  generating group summaries. It takes as input a data frame, and
  outputs a data frame with the same number of columns, and only one
  row: if a column is numeric, it outputs the median; if the column is
  a factor, it outputs the level with the highest count; otherwise, it
  just outputs the first row.
}
\value{
  Functions \code{plot.pcaviz} and \code{pcaviz_violin} return a
  \code{ggplot} object.
}
\section{Overlays}{
  For additional customization of PCAviz plots, we have included an
  "overlay" option. When not \code{NULL}, argument \code{overlay} should
  be a function that accepts a single input argument, a \code{ggplot}
  object. The return value is also a \code{ggplot} object, which is
  presumable a modified version of the input \code{ggplot} object. The
  overlay function is called before any graphics ("geom's") are drawn,
  so that any points, labels and so on are plotted on top of the overlay
  graphics.

  The overlay function can be very simple. For example, this overlay
  adds the x = y line to the plot: \code{overlay_abline <- function (p)
  p + geom_abline(slope = 1, intercept = 0, color = "orangered",
  linetype = "dotted")}.

  See \link{pcaviz-maps} for examples of overlays that draw world maps
  on the plot using the \pkg{maps} package.  
}
\note{
  Note that some features of \pkg{ggplot2} graphics are not fully
  implemented in \pkg{plotly}, so some of the plotting options are
  disabled when \code{plotly = TRUE}.

  All plotting variables must be specified as character strings. This is
  because the \pkg{ggplot2} aesthetic mappings are specified in PCAviz
  using \code{\link[ggplot2]{aes_string}}, not the more commonly used
  \code{\link[ggplot2]{aes}}. The latter is more flexible (e.g., can
  take expressions instead of column names) but is also more error-prone
  when used inside a function due to the fact that the expressions will
  be evaluated outside the function call.

  One important reason for using \pkg{ggplot2} graphics is that plots
  are easily extensible. However, due to ggplot's use of R's
  "lazy evaluation", the \code{inherit.aes} option should be set to
  \code{FALSE} when adding additional layers to the plot; e.g., using
  \code{\link[ggplot2]{geom_point}}. See the examples below and
  vignettes for demonstrations of how to correctly extend plots created
  using the PCAviz plotting functions.
}
\references{
  Wong., B. (2011). Points of view: color blindness. \emph{Nature
  Methods} \bold{8}, 441.

  Wong, B. and Krzywinski, M. (2013). Points of view: plotting symbols.
  \emph{Nature Methods} \bold{10}, 451.
}
\seealso{
  \code{\link{pcaviz}},
  \link{pcaviz-transformations},
  \link{pcaviz-maps},
  \code{\link{create_abbreviated_label}},
  \code{\link[stats]{lm}},
  \code{\link[graphics]{plot}},
  \code{\link[ggplot2]{ggplot}},
  \code{\link[ggplot2]{geom_point}},
  \code{\link[ggplot2]{geom_text}},
  \code{\link[ggplot2]{geom_violin}},
  \code{\link[ggplot2]{geom_abline}},
  \code{\link[ggplot2]{scale_shape_manual}},
  \code{\link[ggplot2]{scale_color_manual}},
  \code{\link[ggplot2]{scale_color_gradientn}},
  \code{\link[ggplot2]{theme}},
  \code{\link[plotly]{ggplotly}},
  \code{\link[htmlwidgets]{saveWidget}}
  \code{\link[cowplot]{theme_cowplot}},
  \code{\link[countrycode]{countrycode}},
  \code{\link[countrycode]{countrycode_data}}
}
\examples{

# Create pcaviz object from prcomp output.
data(iris)
x <- pcaviz(prcomp(iris[1:4]),dat = iris)

# Plot the samples projected onto the first two PCs, and label the
# samples by Species.
plot(x)

# Show Species as different colors and shapes.
plot(x,draw.points = TRUE)

# Show Petal Width as different colors, and Species as different shapes.
plot(x,draw.points = TRUE,color = "Petal.Width",shape = "Species")

# Plot the first PC vs. Petal Width.
plot(x,coords = c("PC1","Petal.Width"),draw.points = TRUE)

# Plot all combinations of PCs 1-4.
plot(x,coords = paste0("PC",1:4),draw.points = TRUE,group = NULL)

# Draw the rotated PC axes.
y <- pcaviz_reduce_whitespace(x)
y <- pcaviz_rotate(y,15)
plot(y,draw.points = TRUE)

# Show the proportion of variance explained by PCs 1-4.
screeplot(x)

# Plot Species against PCs 1-4.
pcaviz_violin(x)

\dontrun{

# One advantage of ggplot2 is that it easily allows for combining plots;
# here we give an example of overlaying additional data on top of the
# PCAviz scatterplot. Note that the inherit.aes = FALSE option is
# needed, otherwise it will give an error.
library(ggplot2)
d <- subset(x$data,Petal.Width > 2.25)
d <- with(d,data.frame(x = PC1,y = PC2))
plot(x,draw.points = TRUE,group = NULL) +
  geom_point(data = d,aes(x = x,y = y),
             shape = 1,size = 5.5,inherit.aes = FALSE)

# These the vignettes for extended demonstrations of the PCAviz
# plotting interface.
vignette("iris")
vignette("popres")
vignette("regmap")
}
}
