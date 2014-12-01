#BSINGH mosaic_support.f90 has to go through C pre-processor (cpp)
rm module_mosaic_support.f90
cpp BaseCode/main/module_mosaic_support.f90 module_mosaic_support.f90

# Now ready to run Makefile
