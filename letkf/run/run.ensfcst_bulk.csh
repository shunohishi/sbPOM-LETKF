#!/bin/csh
#========= ARGUMENT ============================
set switch_iau=${argv[1]}
set CDIR=${argv[2]}; set WORKDIR=${argv[3]}; set OUTPUT=${argv[4]}; set INFO=${argv[5]}
set MEMBER=${argv[6]}; set NPROC=${argv[7]}
set MODELDIR=${argv[8]}; set MODELDATADIR=${argv[9]}; set TIDEDIR=${argv[10]}; set REGION=${argv[11]}
set EPROC=${argv[12]}; set ENODE=${argv[13]}; set ETIME=${argv[14]}; set WTIME=${argv[15]}
set PROG=${argv[16]}; set BUDGET=${argv[17]}
set TS_NUDGE=${argv[18]}; set TI_NUDGE=${argv[19]}; set SS_NUDGE=${argv[20]}; set SI_NUDGE=${argv[21]}
set GRID=${argv[22]}; set TSCLIM=${argv[23]}; set IC=${argv[24]}; set LBC=${argv[25]}
set ATM=${argv[26]}; set FFLUX=${argv[27]}; set TSDATA=${argv[28]}
set DT=${argv[29]}; set ST=${argv[30]}; set IT=${argv[31]}

#========= TIMER ================================

@ ETIME2 = 2 * ${ETIME}
set ETIME=`printf "%02d" ${ETIME}`

#========= Date =================================

set sdate=`perl juldays.prl ${ST}`
set yyyys=`printf "%04d" ${sdate[1]}`; set mms=`printf "%02d" ${sdate[2]}`; set dds=`printf "%02d" ${sdate[3]}`

@ BT = ${IT} - 1
set date0=`perl juldays.prl ${BT}`
set yyyy0=`printf "%04d" ${date0[1]}`; set mm0=`printf "%02d" ${date0[2]}`; set dd0=`printf "%02d" ${date0[3]}`

set date1=`perl juldays.prl ${IT}`
set yyyy1=`printf "%04d" ${date1[1]}`; set mm1=`printf "%02d" ${date1[2]}`; set dd1=`printf "%02d" ${date1[3]}`

#=== MODEL Setting ==============================
@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`

    echo "Make in out DIR: $MEM" #---------------------
    foreach dir(in out)
	if( ! -d ${WORKDIR}/${MEM}/${dir})then
	    mkdir -p ${WORKDIR}/${MEM}/${dir}
	endif
    end
    #---------------------------------------------

    echo "Clean: ${MEM}" #--------------------------------
    (rm -f ${OUTPUT}/gues/${MEM}/restart.${yyyy1}${mm1}${dd1}.nc &)
    (rm -f stdout.${REGION} stderr.${REGION} &)
    (rm -f restart.nc &)
    (rm -f ${REGION}.grid.nc ${REGION}.tsclim.nc ${REGION}.ic.nc ${REGION}.tsdata.nc ${REGION}.lbc.nc ${REGION}.atm.nc &)
    (rm -f ${WORKDIR}/${MEM}/TIDE &)
    (rm -f ${WORKDIR}/${MEM}/${PROG} &)
    #----------------------------------------------

    @ IMEM++
end
wait

@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`

    echo "Check restart and iau file link: ${MEM}" #-----------------
    #Restart file

    if(${switch_iau} == 0) set ag="anal"
    if(${switch_iau} == 1) set ag="gues"

    if(! -f ${OUTPUT}/${ag}/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc)then
	echo "***Error: Not Find ${OUTPUT}/${ag}/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc";exit 99
    endif

    #IAU file
    if(! -f ${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc)then
	echo "***Error: Not Find ${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc"; exit 99
    endif

    @ IMEM++

end
wait

@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`

    echo "Make Netcdf file & TIDE dir & EXE file link: ${MEM}" #-----------------
    #Restart file

    if(${switch_iau} == 0) set ag="anal"
    if(${switch_iau} == 1) set ag="gues"

    (ln -s ${OUTPUT}/${ag}/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc ${WORKDIR}/${MEM}/in/restart.nc &)
    #IAU file
    (ln -s ${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc ${WORKDIR}/${MEM}/in/iau.nc &)
    (ln -s ${MODELDATADIR}/${GRID}.nc ${WORKDIR}/${MEM}/in/${REGION}.grid.nc &)
    (ln -s ${MODELDATADIR}/${TSCLIM}.nc ${WORKDIR}/${MEM}/in/${REGION}.tsclim.nc &)
    (ln -s ${MODELDATADIR}/${IC}.nc ${WORKDIR}/${MEM}/in/${REGION}.ic.nc &) 
    (ln -s ${MODELDATADIR}/${TSDATA}.${MEM}.nc ${WORKDIR}/${MEM}/in/${REGION}.tsdata.nc &)
    (ln -s ${MODELDATADIR}/${LBC}.${MEM}.nc ${WORKDIR}/${MEM}/in/${REGION}.lbc.nc &)
    (ln -s ${MODELDATADIR}/${ATM}.${MEM}.nc ${WORKDIR}/${MEM}/in/${REGION}.atm.nc &)
    (ln -s ${MODELDATADIR}/${FFLUX}.nc ${WORKDIR}/${MEM}/in/${REGION}.fflux.nc &)

    (ln -s ${TIDEDIR} ${WORKDIR}/${MEM}/TIDE &)

    (ln -s ${MODELDIR}/${PROG} ${WORKDIR}/${MEM}/${PROG} &)
    #----------------------------------------------

    cd ${WORKDIR}/${MEM}

    if(${switch_iau} == 0) set assim=3
    if(${switch_iau} == 1) set assim=1

    echo "Make NAME list: ${MEM}" #------------------------
cat <<EOF > pom.nml
&pom_nml
    title = "${REGION}"
    netcdf_file = "${REGION}"
    mode = 3
    assim = ${assim}
    nadv = 2
    nitera = 2
    sw = 1.0
    npg = 2
    dte = 12.0
    isplit = 25
    time_start = "${yyyys}-${mms}-${dds} 00:00:00 +00:00"
    nread_rst = 1
    read_rst_file = "restart.nc"
    read_iau_file = "iau.nc"
    write_rst = ${DT}
    write_rst_file = "restart"
    write_ens = 999.
    write_ens_file = "ens"
    budget = ${BUDGET}
    days =  ${DT}.
    prtd1 = 1.0
    prtd2 = 1.0
    swtch = 9999.
    ts_nudge = ${TS_NUDGE}
    ti_nudge = ${TI_NUDGE}
    ss_nudge = ${SS_NUDGE}
    si_nudge = ${SI_NUDGE}
/
EOF
    #-----------------------------------------------

    @ IMEM++

end
wait

cd ${WORKDIR}

echo "Submit Job" #-------------------------
@ PROC = ${EPROC} * ${ENODE}
@ THREAD = ${NPROC} / ${EPROC}
jxsub --bulk --sparam "1-${MEMBER}"<<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=${ENODE}
#JX --mpi proc=${PROC}
#JX -L elapse=00:${ETIME}:00
#JX -N sbPOM_${REGION}_${yyyy1}${mm1}${dd1}_${MEM}

export PLE_MPI_STD_EMPTYFILE="off"
export OMP_NUM_THREAD=${THREAD}
export PARALLEL=${THREAD}

echo \${PJM_JOBID} > JOBID
cd ${WORKDIR}/\$(printf "%03d" "\${PJM_BULKNUM}")
mpiexec -n ${PROC} -stdout stdout.${REGION} -stderr stderr.${REGION} ./${PROG} && finished > FINISHED_sbPOM
cd ${WORKDIR}

EOF
#------------------------------------------------    

cd ${WORKDIR}
echo "Wait for sbPOM end"#---------------------------------------------
@ isec = 0
@ imin = 0
@ jsec = 0
@ jmin = 0
@ iint = 10
@ ijob = 0
@ ITRY = 0
@ NTRY = 3
set JOBID=000000
set LMEM=`printf "%03d" ${MEMBER}`
while(${imin} >= 0)

    sleep ${iint}s
    @ isec = ${isec} + ${iint}
    @ imin = ${isec} / 60
    @ imod = ${isec} % 60

    if(-f ${LMEM}/JOBID && ${ijob} == 0)then
        set JOBID=`head -n 1 ${LMEM}/JOBID`
        echo "Find JOBID: ${JOBID}"
        @ ijob++
    endif

    if(-f ${LMEM}/JOBID)then
        @ jsec = ${jsec} + ${iint}
        @ jmin = ${jsec} / 60
        @ jmod = ${jsec} % 60
    endif

    if(${jmin} >= ${ETIME2})then
        @ IMEM = 1
        while(${IMEM} <= ${MEMBER})
            set MEM=`printf "%03d" ${IMEM}`
            if(! -f ${MEM}/FINISHED_sbPOM || ! -f ${MEM}/out/restart.nc)then
                echo "Execute ${MEM} each ensemble fcst"
                cd ${CDIR}
                csh run.ensfcst_each.csh ${CDIR} ${WORKDIR} ${IMEM} ${REGION} ${EPROC} ${ENODE} ${NPROC} ${ETIME} ${PROG} ${IT}
                cd ${WORKDIR}
                @ jsec = 0
                @ jmin = 0
                set LMEM=`printf "%03d" ${IMEM}`
            endif
            @ IMEM++
        end
        @ ITRY++
    endif

    set FIN_NUM=`find ${WORKDIR} -name FINISHED_sbPOM | wc -l`
    set RESTART_NUM=`find ${WORKDIR}/*/out -name restart.nc | wc -l`
    echo "sbPOM_iau FINISHED JOB[ ${FIN_NUM} / ${MEMBER} ] RESTART FILE[ ${RESTART_NUM} / ${MEMBER}] ${imin}:${imod} passed"

    #set STATS_NUM=`find ${WORKDIR} -name "sbPOM_${REGION}_${yyyy1}${mm1}${dd1}_iau${JOBID}\[*\].stats" | wc -l`
    #echo "sbPOM_iau FINISHED JOB[ ${FIN_NUM} / ${MEMBER} ] RESTART FILE[ ${RESTART_NUM} / ${MEMBER}] STATS[ ${STATS_NUM} / ${MEMBER} ] ${imin}:${imod} passed"

    if(${FIN_NUM} == ${MEMBER} && ${RESTART_NUM} == ${MEMBER})then
        sleep ${iint}s
        break
    endif

    if(${imin} >= ${WTIME} || ${jmin} >= ${ETIME2} || ${ITRY} == ${NTRY})then
        jxdel ${JOBID}
        exit 88
    endif

end
#-----------------------------------------------------------------------

echo "Move & Remove Log" #-----------------------------------------
@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`

    #Restart file
    if(-f ${MEM}/out/restart.nc)then
	(mv ${MEM}/out/restart.nc ${OUTPUT}/gues/${MEM}/restart.${yyyy1}${mm1}${dd1}.nc &)
    else
	echo "***Error: Not Found ${MEM}/out/restart.nc"; exit 88
    endif

    #Output file
    if(-f ${MEM}/out/${REGION}.nc)then
	(mv ${MEM}/out/${REGION}.nc ${OUTPUT}/gues/${MEM}/${REGION}.${yyyy1}${mm1}${dd1}.nc &)
    else
	echo "***Error: Not Found ${MEM}/out/${REGION}.nc"; exit 88
    endif

    #Log
    if(-f ${MEM}/stdout.${REGION})then
	(mv -f ${MEM}/stdout.${REGION} ${INFO}/stdout.${REGION}.${yyyy1}${mm1}${dd1}.${MEM} &)
    else
	echo "Not Found ${MEM}/stdout.${REGION}"
    endif

    if(-s ${MEM}/stderr.${REGION})then
	(mv -f ${MEM}/stderr.${REGION} ${INFO}/stderr.${REGION}.${yyyy1}${mm1}${dd1}.${MEM} &)
    else
	(rm -f ${MEM}/stderr.${REGION} &)
    endif

    if(-s sbPOM_${REGION}_${yyyy1}${mm1}${dd1}${JOBID}\[${IMEM}\].stats)then
	(mv -f sbPOM_${REGION}_${yyyy1}${mm1}${dd1}${JOBID}\[${IMEM}\].stats ${INFO}/ &)
    else
	(rm -f sbPOM_${REGION}_${yyyy1}${mm1}${dd1}${JOBID}\[${IMEM}\].stats &)
    endif

    (rm -f sbPOM_${REGION}_${yyyy1}${mm1}${dd1}${JOBID}\[${IMEM}\].out &)
    (rm -f sbPOM_${REGION}_${yyyy1}${mm1}${dd1}${JOBID}\[${IMEM}\].err &)

    @ IMEM++

end
wait
#-------------------------------------------------------------------

#---Check restart file----------------------------------------------
@ ITRY = 0
@ NTRY = 3
while(${ITRY} <= ${NTRY})

    @ IMEM = 1
    @ FMEM = 0

    while(${IMEM} <= ${MEMBER})

        set MEM=`printf "%03d" ${IMEM}`
	set filename=${OUTPUT}/gues/${MEM}/restart.${yyyy1}${mm1}${dd1}.nc

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
        echo "***Error: Not found all restart file in sbPOM run"
        exit 88
    endif

    @ ITRY++

end
wait

echo "Daily & Ensemble mean"#---------------------------------------
cd ${CDIR}

(echo ${OUTPUT}/gues > ens_info.txt && \
 echo ${REGION} >> ens_info.txt && \
 echo ${yyyys} ${mms} ${dds} >> ens_info.txt && \
 echo ${yyyy1} ${mm1} ${dd1} >> ens_info.txt && \
 echo ${MEMBER} ${DT} >> ens_info.txt && \
 echo ${BUDGET} >> ens_info.txt &)
#---------------------------------------------------------------------
wait
exit 0
