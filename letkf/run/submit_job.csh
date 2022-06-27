#!/bin/csh
#--------------------------------------------------------------------------
# Submit Job cshfile with vcoordfile
#-------------------------------------------------------------------------

#---Argument---------------------------------------------------------------
set EPROC=$argv[1];set MEMBER=$argv[2]; set THREAD=$argv[3]
set REGION=$argv[4];set EXE=$argv[5];set WORKDIR=$argv[6]; set JOBID=$argv[7]
#--------------------------------------------------------------------------

#---Submit Job------------------------------------------------------------
@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`

    cd ${WORKDIR}/${MEM}
    echo "Submit Job: "${WORKDIR}/${MEM}
    echo ${JOBID} > JOBID     

    setenv OMP_NUM_THREADS ${THREAD}
    setenv PARALLEL ${THREAD}
    mpiexec -n ${EPROC} -stdout stdout.${REGION} -stderr stderr.${REGION} --vcoordfile vcoordfile ./${EXE} && echo finished > ${WORKDIR}/${MEM}/FINISHED_sbPOM &
    cd ${WORKDIR}

    @ IMEM++

end
#-------------------------------------------------------------------------

#echo "Wait for Ensemble POM end"#-----------------------------------------
@ int = 10 #Interval [sec.]
@ isec = 0
while(${isec} >= 0)
    sleep ${int}s
    @ isec = ${isec} + ${int}
    @ sec = ${isec} % 60
    @ min = ${isec} / 60
    set FIN_NUM=`find ${WORKDIR} -name FINISHED_sbPOM | wc -l`
    echo "sbPOM FINISHED [ ${FIN_NUM} / ${MEMBER} ], ${min}:${sec} passed"
    if(${FIN_NUM} == ${MEMBER})then
	echo "End ensemble sbPOM"
	break
    endif
end
-------------------------------------------------------------------------
wait
exit
