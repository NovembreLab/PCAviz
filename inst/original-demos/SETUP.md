## Setup instructions

These are the steps I took to install R and the necessary packages in
order to be able to successfully run the R scripts in
[demo/R-2.5](demo/R-2.5).

These installation instructions worked on machines with Red Hat
Enterprise Linux 5.5 (kernel 2.6.18-194, gcc 4.1.2) and a version of
Scientific Linux 6 based on Red Hat Enterprise Linux 6.7
(kernel 2.6.32-573, gcc 4.4.7). **Note:** The latter is *midway1*.

First, install [R 2.5.1](http://stat.ethz.ch/pipermail/r-announce/2007/000465.html) (release date is June 26, 2007).

```bash
wget https://cran.r-project.org/src/base/R-2/R-2.5.1.tar.gz
tar -zvxf R-2.5.1.tar.gz
cd R-2.5.1
./configure --prefix=$INSTALL_ROOT/R-2.5.1
make 
make check
make install
```

Next, Install [Proj](http://proj4.org) 4.6.1.

```bash
wget https://github.com/OSGeo/proj.4/archive/proj_4_6_1.tar.gz 
tar -zxvf proj_4_6_1.tar.gz
cd proj.4-proj_4_6_1 
./configure --prefix=$INSTALL_ROOT/proj-4.6.1
make
make check
make install
```

Install [expat 2.0.1](https://sourceforge.net/projects/expat/files/expat/2.0.1).

```bash
wget --no-check-certificate https://superb-sea2.dl.sourceforge.net/project/expat/expat/2.0.1/expat-2.0.1.tar.gz
tar -zxvf expat-2.0.1.tar.gz
cd expat-2.0.1
./configure --prefix=$INSTALL_ROOT/expat-2.0.1
make all
make check
make install
```

Install [gdal](http://www.gdal.org) 1.5.4.

```bash
wget http://download.osgeo.org/gdal/gdal-1.5.4.tar.gz
tar -zxvf gdal-1.5.4.tar.gz
cd gdal-1.5.4
./configure --prefix=$INSTALL_ROOT/gdal-1.5.4 \
  --with-expat=$INSTALL_ROOT/expat-2.0.1
make
make install
```

Install R packages.

```bash
wget https://cran.r-project.org/src/contrib/Archive/RColorBrewer/RColorBrewer_1.0-2.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/sp/sp_0.9-28.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/maps/maps_2.0-40.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/rgdal/rgdal_0.5-33.tar.gz
R-2.5.1 CMD INSTALL RColorBrewer_1.0-2.tar.gz
R-2.5.1 CMD INSTALL sp_0.9-28.tar.gz
R-2.5.1 CMD INSTALL maps_2.0-40.tar.gz
R-2.5.1 CMD INSTALL \
  --configure-args='--with-proj-include=$INSTALL_ROOT/proj-4.6.1/include --with-proj-lib=$INSTALL_ROOT/proj-4.6.1/lib --with-gdal-config=$INSTALL_ROOT/gdal-1.5.4/bin/gdal-config' rgdal_0.5-33.tar.gz
```

Finally, before running R, make sure that your environment variables
are correctly set up to point to the installed libraries:

```bash
export LD_LIBRARY_PATH=$INSTALL_ROOT/proj-4.6.1/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$INSTALL_ROOT/gdal-1.5.4/lib:$LD_LIBRARY_PATH
```
