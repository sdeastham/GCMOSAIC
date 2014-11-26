#!/bin/csh
source ~sing201/comp_scr/purge_all_load.csh
unsetenv LD_LIBRARY_PATH

if ($#argv == 0 ) then
    module load pgi/12.2
    set compiler = pgf90
    set flags    = "-g -C -byteswapio -Ktrap=fp -O0"
endif


#Load compiler
if ( $1 == 'pgf90' ) then
    module load pgi/12.2
    set compiler = $1
    set flags    = "-g -C -byteswapio -Ktrap=fp -O0"
else if ( $1 == 'lf95' ) then
    module load lahey/lf6481
    set compiler = $1
    set flags    = "--ap --chk a,e,s,u --dal  -O0 --trace --trap --pca -g"
else
    echo 'Please specify a valid compiler (pgf90 or lf95)'
    exit 1
endif


rm mosaic.x
./clean.sh
./compile.sh $compiler $flags |& tee outLog 
