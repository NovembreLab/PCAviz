# PCAviz notes

Additional project notes go here.

## Project summary

A package that rewrites and expands upon the utility of a collection
of functions previously written for PCA plotting and visualizations
([prepSpatial.R](inst/original-files/code/prepSpatial.R) and
[PlotPCmaps.R](inst/original-files/code/PlotPCmaps.R)) is to be developed.
Additional plotting functions are to be developed (violin plots, scree
plots, PC loading plots). Additional utilities are to be added to
existing functions (mirroring of PCA plots, axis rotation for PCA
plots, hover-over descriptors and IDs for individual points, map inset
paired to PCA data).

## Data sets

Data sets we use for illustrating the PCAviz package:

+ Iris Data; type `data(iris)` in R.

+ POPRES data; type `data(popres)` in R.

+ RegMap data (see
  [the Nature Genetics paper](inst/original-files/etc/regmap.pdf)); type
  `data(regmap)` in R.

## Other useful resources

+ ["Genes mirror geography" paper](inst/original-files/etc/novembre2008.pdf)
  and [accompanying supplementary
materials](inst/original-files/etc/novembre2008-supplement.pdf).

+ ["Visualizing geography of genetic variants" paper](inst/original-files/etc/marcus2016.pdf)

+ [The rsvd package](http://github.com/Benli11/rSVD) for randomized
  singular value decomposition (and PCA) in R.

## Rough breakdown of project objectives

+ New world countries polygons dataset (from mapdata R package) to
replace outdated worldpolys.rda dataset.

+ Code from preSpatial.R to be modified for new dataset and made into
a processing function executed upon loading of the PCAviz library.

+ Functions to be updated to match new country polygons: `plotGeoPoints`
(preSpatial.R), `plotPCbiplot` (PlotPCmaps.R) and `plotPC` (PlotPCmaps.R).

+ plotPC function to be modified to enable mirroring of PCA plot axes
(PlotPCmaps.R).

+ Functions from PlotPCmaps.R to be modified to enable the viewing of
PCA point details upon hovering over the point, use of
[Shiny](http://shiny.rstudio.com/reference/shiny/latest/hoverOpts.html)
preferred.
  
+ New plotting functions to be developed/integrated:
[Violin plot](http://cran.r-project.org/package=vioplot),
[Scree plot](http://stat.ethz.ch/R-manual/R-devel/library/stats/html/screeplot.html) and PC loadings plot (Bar plot).

