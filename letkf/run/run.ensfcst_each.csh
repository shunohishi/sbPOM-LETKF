#!/bin/csh
#====== ARGUMENT =========================#
set CDIR=${argv[1]};set WORKDIR=${argv[2]};set IMEM=${argv[3]}
set REGION=${argv[4]};set EPROC=${argv[5]};set ENODE=${argv[6]};set NPROC=${argv[7]}; set ETIME=${argv[8]};set PROG=${argv[9]}
set IT=${argv[10]}
#-----------------------------------------
set MEM=`printf "%03d" ${IMEM}`
set ETIME=`printf "%02d" ${ETIME}`
set date1=`perl juldays.prl ${IT}`
set yyyy1=`printf "%04d" ${date1[1]}`; set mm1=`printf "%02d" ${date1[2]}`; set dd1=`printf "%02d" ${date1[3]}`
#=========================================#

cd ${WORKDIR}
@ PROC = ${EPROC} * ${ENODE}
@ THREAD = ${NPROC} / ${EPROC}

    jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=${ENODE}
#JX --mpi proc=${EPROC}
#JX -L elapse=00:${ETIME}:00
#JX -N sbPOM_${REGION}_${yyyy1}${mm1}${dd1}_iau${MEM}
#JX -s

export PLE_MPI_STD_EMPTYFILE="off"
export OMP_NUM_THREAD=${THREAD}
export PARALLEL=${THREAD}

cd ${WORKDIR}/${MEM}
echo \${PJM_JOBID} > JOBID
mpiexec -n ${PROC} -stdout stdout.${REGION} -stderr stderr.${REGION} ./${PROG} && echo finished > FINISHED_sbPOM

EOF

exit 1

