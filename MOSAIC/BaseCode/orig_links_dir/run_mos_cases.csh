#!/bin/csh

set cases = ( \
case2.inp \
case3.inp \
case5.inp \
case6.inp \
case9.inp \
case10.inp \
)

foreach case( $cases )
    echo processing ... $case
    mosaic.x $case
end
