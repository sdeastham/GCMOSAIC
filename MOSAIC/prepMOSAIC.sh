rm *.f90 *.o *.f

cp -p  BaseCode/datamodules/*.f90  .
cp -p  BaseCode/main/*.f90  .
cp -p  BaseCode/gas/*.f90  .
cp -p  BaseCode/solver/*.f  .
cp -p  BaseCode/aerosol/*.f90  .
cp -p  BaseCode/cloud/*.f90  . 

#BSINGH mosaic_support.f90 has to go through C pre-processor (cpp)
./cppSupport.sh

# Now ready to run Makefile
