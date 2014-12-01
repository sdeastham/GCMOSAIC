#!/bin/csh -f
# compile.sh
set verbose

/bin/rm mosaic.x

cd  compile
/bin/rm  *.f90  Makefile*

cp -p  ../datamodules/*.f90  .
cp -p  ../main/*.f90  .
cp -p  ../gas/*.f90  .
cp -p  ../solver/*.f  .
cp -p  ../aerosol/*.f90  .
cp -p  ../cloud/*.f90  . 

#BSINGH mosaic_support.f90 has to go through C pre-processor (cpp)
/bin/rm module_mosaic_support.f90
cpp ../main/module_mosaic_support.f90 module_mosaic_support.f90

cp -p  ../Makefile  .

make  mosaic.x

/bin/mv  mosaic.x  ..

cd  ..

unset verbose
