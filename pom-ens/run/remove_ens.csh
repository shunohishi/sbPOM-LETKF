#!/bin/csh
#==============================================================
#   Remove ensemble member
#==============================================================

#--- Argument

set WORKDIR=$argv[1]; set REGION=$argv[2]; set MEMBER=$argv[3]

#--- Remove ensemble member

@ IMEM = 1

while($IMEM <= $MEMBER)

    set MEM=`printf "%03d" ${IMEM}`
    rm ${WORKDIR}/${MEM}/out/${REGION}.nc

    @ IMEM++

end
