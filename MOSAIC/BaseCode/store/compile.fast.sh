# compile.fast.sh
set verbose

cd  compile
/bin/rm  *.f90  Makefile*

cp -p  ../datamodules/*.f90  .
cp -p  ../main/*.f90  .
cp -p  ../gas/*.f90  .
cp -p  ../solver/*.f90  .
cp -p  ../aerosol/*.f90  .
cp -p  ../cloud/*.f90  . 

cp -p  ../Makefile.fast  .

make  -f Makefile.fast  mosaic.fast.x

/bin/mv  mosaic.fast.x  ..
cd  ..

unset verbose
