#!/bin/csh
#--------------------------------------------------------------------------
# Submit Job cshfile with vcoordfile
#-------------------------------------------------------------------------

#---Argument---------------------------------------------------------------
set EPROC=$argv[1];set MEMBER=$argv[2]
set REGION=$argv[3];set EXE=$argv[4];set WORKDIR0=$argv[5]
#--------------------------------------------------------------------------

#---Submit Job------------------------------------------------------------
@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`

    cd ${WORKDIR0}/${MEM}
    echo "Submit Job: "${WORKDIR0}/${MEM}

    (mpiexec -n ${EPROC} -stdout stdout.${REGION} -stderr stderr.${REGION} --vcoordfile vcoordfile ./${EXE} && echo finished > ${WORKDIR0}/${MEM}/FINISHED) &
    cd ${WORKDIR0}

    @ IMEM++

end
#-------------------------------------------------------------------------

echo "Wait for Ensemble POM end"#-----------------------------------------
@ int = 10 #Interval [sec.]
@ isec = 0
while(${isec} >= 0)
    sleep ${int}s
    @ isec = ${isec} + ${int}
    @ sec = ${isec} % 60
    @ min = ${isec} / 60
    set FIN_NUM=`find ${WORKDIR0} -name FINISHED | wc -l`
    echo "sbPOM FINISHED [ ${FIN_NUM} / ${MEMBER} ], ${min}:${sec} passed"
    if(${FIN_NUM} == ${MEMBER})then
	echo "End ensemble sbPOM"
	break
    endif
end
#-------------------------------------------------------------------------
wait
exit

