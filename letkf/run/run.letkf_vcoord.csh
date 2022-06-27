#!/bin/csh
#===== ARGUMENT ============================================
set switch_iau=${argv[1]}
set WORKDIR=${argv[2]}; set OUTPUT=${argv[3]}; set INFO=${argv[4]}
set MEMBER=${argv[5]}; set EPROC=${argv[6]}
set OBSDIR=${argv[7]}; set OBSFILE=${argv[8]}
set LETKFDIR=${argv[9]}; set PROG=${argv[10]}  
set MODELDATADIR=${argv[11]}; set REGION=${argv[12]}
set GRID=${argv[13]}; set IT=${argv[14]}

echo "IT:"${IT}
@ BT = ${IT} - 1
set date0=`perl juldays.prl ${BT}`
set yyyy0=`printf "%04d" ${date0[1]}`;set mm0=`printf "%02d" ${date0[2]}`;set dd0=`printf "%02d" ${date0[3]}`

set date1=`perl juldays.prl ${IT}`
set yyyy1=`printf "%04d" ${date1[1]}`;set mm1=`printf "%02d" ${date1[2]}`;set dd1=`printf "%02d" ${date1[3]}`

#=== Setting ===============================================

echo "Link LETKF EXE file" #----------
(rm -f ${WORKDIR}/${PROG} && \
ln -s ${LETKFDIR}/${PROG} ${WORKDIR}/${PROG} &)
#-------------------------------------

echo "Move Workdir"#---
cd ${WORKDIR}
#----------------------

echo "Link Observational data"#-------------------------
rm -f obs01.dat
if(-f ${OBSDIR}/${OBSFILE}${yyyy1}${mm1}${dd1}.nc)then
    (ln -s ${OBSDIR}/${OBSFILE}${yyyy1}${mm1}${dd1}.nc obs01.dat &)
else
    echo "Not Find ${OBSDIR}/${OBSFILE}${yyyy1}${mm1}${dd1}.nc";exit 99
endif
#-------------------------------------------------------
	    
echo "Link Grid data"#----------------------------------
(ln -s ${MODELDATADIR}/${GRID}.nc grd.nc &)
#-------------------------------------------------------

echo "Clear & restart --> gs01***.nc, anal***.nc" #-----
#gs01***.nc
@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`
    if(${switch_iau} == 0) set file="${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc"
    if(${switch_iau} == 1) set file="${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}_iau.nc"

    if(! -f ${file})then
	echo "Not Find ${file}"
	exit 99
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
wait
#---------------------------------------------------------

echo "Copy anal001.nc --> gues_me,sp.nc anal_me,sp.nc"#---
cp -v anal001.nc gues_me.nc &
cp -v anal001.nc gues_sp.nc &
cp -v anal001.nc anal_me.nc &
cp -v anal001.nc anal_sp.nc &
wait
#---------------------------------------------------------

#=== Submit LETKF Job ======================================
echo "Submit LETKF Job"#--------------------------------
@ PROC = ${MEMBER} * ${EPROC}
echo ${PROC}
(mpiexec -n ${PROC} -stdout stdout.letkf.${REGION} -stderr stderr.letkf.${REGION} ./${PROG}; echo finished > FINISHED_LETKF) &
#---------------------------------------------------------

#--- Wait for LETKF end ----------------------------------
echo "Wait for LETKF end"#---------------------------------------------
@ imin = 0
@ itime = 0
while(${imin} >= 0)

    if(-f FINISHED_LETKF)then
	echo "Finished LETKF; ${imin} min passed"
        break
    else
        sleep 5s
        @ itime = ${itime} + 5
    endif

    if(${itime} % 60 == 0)then
        @ imin++
	echo "LETKF; ${imin} min passed"
    endif

end
#-----------------------------------------------------------------------

echo "Move LETKF output data" #----------------------------
(mv gues_me.nc ${OUTPUT}/gues/mean/restart.${yyyy0}${mm0}${dd0}.nc &)
(mv gues_sp.nc ${OUTPUT}/gues/sprd/restart.${yyyy0}${mm0}${dd0}.nc &)
(mv anal_me.nc ${OUTPUT}/anal/mean/restart.${yyyy0}${mm0}${dd0}.nc &)
(mv anal_sp.nc ${OUTPUT}/anal/sprd/restart.${yyyy0}${mm0}${dd0}.nc &)

@ IMEM = 1
while(${IMEM} <= ${MEMBER})
    set MEM=`printf "%03d" ${IMEM}`
    if(-e anal${MEM}.nc)then
	(mv anal${MEM}.nc ${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc &)
    else
	echo "***anal${MEM}.nc not exist";exit 99
    endif
    @ IMEM++
end
#----------------------------------------------------------
	
echo "Move & Remove Log" #---------------------------------
tail -n 17 NOUT-00000
(mv NOUT-00000 ${INFO}/letkf.${yyyy0}${mm0}${dd0}.log &)
(mv stdout.letkf.${REGION} ${INFO}/stdout.letkf.${REGION}.${yyyy0}${mm0}${dd0} &)
if(-s stderr.letkf.${REGION})then
    (mv stderr.letkf.${REGION} ${INFO}/stderr.letkf.${REGION}.${yyyy0}${mm0}${dd0} &)
else
    (rm -f stderr.letkf.${REGION} &)
endif
#----------------------------------------------------------

wait

echo "Clean" #---------------------------------------------
@ IMEM = 1
while(${IMEM} <= ${MEMBER})
    set MEM=`printf "%03d" ${IMEM}`
    (rm -f gs01${MEM}.nc &)
    @ IMEM++
end

(rm -f ${PROG} obs01.dat grd.nc &)
(rm -f NOUT-* FINISHED_LETKF &)
#-----------------------------------------------------------    

wait
exit 0

