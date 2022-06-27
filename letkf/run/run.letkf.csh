#!/bin/csh
#===== ARGUMENT ============================================
set switch_iau=${argv[1]}
set WORKDIR=${argv[2]}; set OUTPUT=${argv[3]}; set INFO=${argv[4]}
set MEMBER=${argv[5]}; set NPROC=${argv[6]}
set OBSDIR=${argv[7]}; set OBSFILE=${argv[8]}
set LETKFDIR=${argv[9]}; set PROG=${argv[10]}
set LPROC=${argv[11]}; set LNODE=${argv[12]}; set ETIME=${argv[13]}; set WTIME=${argv[14]}
set MODELDATADIR=${argv[15]}; set GRID=${argv[16]} 
set IT=${argv[17]}; set REGION=${argv[18]}

@ BT = ${IT} - 1
set date0=`perl juldays.prl ${BT}`
set yyyy0=`printf "%04d" ${date0[1]}`; set mm0=`printf "%02d" ${date0[2]}`; set dd0=`printf "%02d" ${date0[3]}`

set date1=`perl juldays.prl ${IT}`
set yyyy1=`printf "%04d" ${date1[1]}`; set mm1=`printf "%02d" ${date1[2]}`; set dd1=`printf "%02d" ${date1[3]}`

#=== TIMER =================================================

set ETIME=`printf "%02d" ${ETIME}`

#=== Setting ===============================================

echo "Link LETKF EXE file" #----------
rm -f ${WORKDIR}/${PROG}
ln -s ${LETKFDIR}/${PROG} ${WORKDIR}/${PROG}
#-------------------------------------

echo "Move Workdir"#---
cd ${WORKDIR}
#----------------------

echo "Link Observational data"#-------------------------
rm -f obs01.dat
if(-f ${OBSDIR}/${OBSFILE}${yyyy1}${mm1}${dd1}.nc)then
    ln -s ${OBSDIR}/${OBSFILE}${yyyy1}${mm1}${dd1}.nc obs01.dat
else
    echo "Not Find ${OBSDIR}/${OBSFILE}${yyyy1}${mm1}${dd1}.nc";exit 99
endif
#-------------------------------------------------------
	    
echo "Link Grid data"#----------------------------------
ln -s ${MODELDATADIR}/${GRID}.nc grd.nc
#-------------------------------------------------------

echo "Clear & restart --> gs01***.nc, anal***.nc" #-----
#gs01***.nc
@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`
    if(${switch_iau} == 0) set file="${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc"
    if(${switch_iau} == 1) set file="${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}_iau.nc"

    if(! -f ${file})then
	echo "***Error: Not Find ${file}"; exit 99
    else
	(rm -f gs01${MEM}.nc && \
	ln -s ${file} gs01${MEM}.nc &)
    endif
    @ IMEM++
end

#anal***.nc
@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`
    if(${switch_iau} == 0) set file="${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc"
    if(${switch_iau} == 1) set file="${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}_iau.nc"

    cp -v ${file} anal${MEM}.nc &

    @ IMEM++
end
#---------------------------------------------------------
	
echo "Copy anal001.nc --> gues_me,sp.nc anal_me,sp.nc"#---
if(${switch_iau} == 0) set file="${OUTPUT}/gues/001/restart.${yyyy0}${mm0}${dd0}.nc"
if(${switch_iau} == 1) set file="${OUTPUT}/gues/001/restart.${yyyy0}${mm0}${dd0}_iau.nc"
cp -v ${file} gues_me.nc &
cp -v ${file} gues_sp.nc &
cp -v ${file} anal_me.nc &
cp -v ${file} anal_sp.nc &
wait
#--------------------------------------------------------

#=== Submit LETKF Job ======================================
@ PROC = ${LNODE} * ${LPROC}
@ THREAD = ${NPROC} / ${LPROC}
echo "NODE: ${LNODE}"
echo "PROC: ${PROC}"
echo "THREAD: ${THREAD}"
echo "Submit LETKF Job"#--------------------------------
    jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=${LNODE}
#JX --mpi proc=${PROC}
#JX -L elapse=00:${ETIME}:00
#JX -N LETKF_${yyyy0}${mm0}${dd0}

export PLE_MPI_STD_EMPTYFILE="off"
export PARALLEL=${THREAD}
export OMP_NUM_THREADS=${THREAD}

echo \${PJM_JOBID} > JOBID
mpiexec -n ${PROC} -stdout stdout.letkf.${REGION} -stderr stderr.letkf.${REGION} ./${PROG} && echo finished > FINISHED_LETKF

EOF
#---------------------------------------------------------

#--- Wait for LETKF end ----------------------------------
@ isec = 0
@ imin = 0
@ jsec = 0
@ jmin = 0
@ iint = 10
@ ijob = 0
set JOBID=000000
while(${imin} >= 0)

    if(-f JOBID && ${ijob} == 0)then
	set JOBID=`head -n 1 JOBID`
	echo "Find JOBID: ${JOBID}"
	@ ijob++
    endif

    #TIMER
    sleep ${iint}s
    @ isec = ${isec} + ${iint}
    @ imin = ${isec} / 60
    @ imod = ${isec} % 60

    if(-f JOBID)then
	@ jsec = ${jsec} + ${iint}
	@ jmin = ${jsec} / 60
	@ jmod = ${jsec} % 60
    endif

    set FIN_NUM=`find ${WORKDIR} -name FINISHED_LETKF | wc -l`
    set STATS_NUM=`find ${WORKDIR} -name LETKF_${yyyy0}${mm0}${dd0}.${JOBID}.stats | wc -l`
    echo "LETKF FINISHED [${FIN_NUM}/1] STATS [${STATS_NUM}/1] ; ${imin}:${imod} passed"

    #BREAK
    if(${FIN_NUM} == 1)then
	sleep ${iint}s
	echo "LETKF FINISHED [${FIN_NUM}/1] STATS [${STATS_NUM}/1] ; ${imin}:${imod} ended"
	break
    endif

    #ERROR
    if(-f JOBID && ${ijob} != 0 && ! -f NOUT-00000)then
        set JOBID=`head -n 1 JOBID`
	jxdel ${JOBID}
	echo "***Error: Not working LETKF ${JOBID}"
	exit 88
    endif

    if(${jmin} == ${WTIME})then
	jxdel ${JOBID}
	exit 88
    endif

end
#----------------------------------------------------------

echo "Move LETKF output data" #----------------------------
if(-f gues_me.nc && -f gues_sp.nc && -f anal_me.nc && -f anal_sp.nc)then
    (mv gues_me.nc ${OUTPUT}/gues/mean/restart.${yyyy0}${mm0}${dd0}.nc &)
    (mv gues_sp.nc ${OUTPUT}/gues/sprd/restart.${yyyy0}${mm0}${dd0}.nc &)
    (mv anal_me.nc ${OUTPUT}/anal/mean/restart.${yyyy0}${mm0}${dd0}.nc &)
    (mv anal_sp.nc ${OUTPUT}/anal/sprd/restart.${yyyy0}${mm0}${dd0}.nc &)
else
    echo "***Error: Not Found gues/anal_me/sp.nc"; exit 99
endif

@ IMEM = 1
while(${IMEM} <= ${MEMBER})
    set MEM=`printf "%03d" ${IMEM}`
    if(-f anal${MEM}.nc)then
	(mv anal${MEM}.nc ${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc &)
    else
	echo "***anal${MEM}.nc not exist"; exit 99
    endif
    @ IMEM++
end
wait
#----------------------------------------------------------

#---Check restart file----------------------------------------------
@ ITRY = 0
@ NTRY = 3
while(${ITRY} <= ${NTRY})

    @ IMEM = 1
    @ FMEM = 0

    while(${IMEM} <= ${MEMBER})

        set MEM=`printf "%03d" ${IMEM}`
        set filename=${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc

        if(! -f ${filename})then
            echo "Wait to move ${filename}"
        else
            @ FMEM++
        endif

        @ IMEM++
    end

    if(${FMEM} == ${MEMBER})then
        break
    else
        sleep 60s
    endif

    if(${ITRY} == ${NTRY})then
        echo "***Error: Not found all restart file in  LETKF"
        exit 88
    endif

    @ ITRY++

end
wait

echo "Move & Remove Log" #---------------------------------
if(-f NOUT-00000)then
    (tail -n 17 NOUT-00000 && mv NOUT-00000 ${INFO}/letkf.${yyyy0}${mm0}${dd0}.log &)
else
    echo "***Error: Not found NOUT-00000"; exit 88
endif

if(-s stdout.letkf.${REGION})then
    (mv -f stdout.letkf.${REGION} ${INFO}/stdout.letkf.${REGION}.${yyyy0}${mm0}${dd0} &)
endif

if(-s stderr.letkf.${REGION})then
    (mv stderr.letkf.${REGION} ${INFO}/stderr.letkf.${REGION}.${yyyy0}${mm0}${dd0} &)
else
    (rm -f stderr.letkf.${REGION} &)
endif

if(-s LETKF_${yyyy0}${mm0}${dd0}.${JOBID}.stats)then
    (mv -f LETKF_${yyyy0}${mm0}${dd0}.${JOBID}.stats ${INFO}/ &)
    (rm -f LETKF_${yyyy0}${mm0}${dd0}.${JOBID}.out LETKF_${yyyy0}${mm0}${dd0}.${JOBID}.err &)
endif
wait
#----------------------------------------------------------

echo "Clean" #---------------------------------------------
@ IMEM = 1
while(${IMEM} <= ${MEMBER})
    set MEM=`printf "%03d" ${IMEM}`
    (rm -f gs01${MEM}.nc &)
    @ IMEM++
end

(rm -f ${PROG} obs01.dat grd.nc &)
(rm -f NOUT-* JOBID FINISHED_LETKF &)
#-----------------------------------------------------------
wait
exit 0

