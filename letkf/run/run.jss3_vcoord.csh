#!bin/csh
#======= DATE ====================#

set sdate=(2011 1 1)
set idate=(2015 7 7)
set edate=(2015 7 7)
#========== Switch(1:execute, other: skip)================#

set switch_prepare_ens=1 #Prepare ensemble member
set switch_iau=1         #Integration of sbPOM for IAU procedure
set switch_anal=0        #LETKF/Make analysis
set switch_gues=0        #sbPOM integration/Make guess 
set switch_rmens=0       #Remove ensemble member (except end day of each month)

#========= GENERAL ================#

set DIR0=/data/hp120279/u00399/TEST
set LDIR0=${DIR0}/letkf
set CDIR=`pwd`             #Current DIR
set WORKDIR=${LDIR0}/work  #Work DIR
set OUTPUT=${LDIR0}/output #Output DIR
set INFO=${LDIR0}/run/info #INFORMATION 
set MEMBER=100             #Ensemble Member
set EPROC=48               #Number of each Processor
set ENODE=1                #Number of each Node
set ETIME=10               #Elapse time [min]
set DT=1                   #Delta T [day]
set EMAIL="shun.ohishi@riken.jp"

@ NODE = ${ENODE} * ${MEMBER}
@ PROC = ${EPROC} * ${MEMBER}
echo "Total Node:"${NODE}
echo "Total Processor:"${PROC}

echo "Check OUTPUT Directory"
if(! -d ${OUTPUT})then
    mkdir -p ${OUTPUT}
endif

echo "Make Work directory"
rm -rf ${WORKDIR}; mkdir -p ${WORKDIR}

#========= OBSERVATION ==================#

set OBSDIR=${DIR0}/prep/obs #Observational DIR for Assimilation
set OBSFILE=obs             #Observation filename: ${OBSFILE}${yyyy}${mm}${dd}.nc

#========== LETKF ==========================#

set LETKFDIR=${LDIR0}/run #LETKF DIR
set LPROG=letkf.exe       #LETKF EXE File

echo "Make LETKF EXE & Log directory"
cd ${LETKFDIR}
make
cd ${CDIR}
if(! -f ${LETKFDIR}/${LPROG})then
    echo "***Error: Not found ${LPROG}"
    exit
endif
if(! -d ${INFO}) mkdir -p ${INFO}

#========== sbPOM ===============#

set MODELDIR=${DIR0}/pom-ens/run          #sbPOM DIR
set MODELDATADIR=${DIR0}/prep/in          #grid,tsclim,ic,lbc,atm,tsdata DIR
set MODELOUTPUTDIR=${DIR0}/pom-ens/output #sbPOM output DIR
set TIDEDIR=${DIR0}/pom-ens/TIDE          #Tide DIR
set REGION=test
set PPROG=pom.exe                         #EXE file
set NEST=0                                #One-way nest 1:on, 0:off
set BUDGET=1                              #1:execute, 0:skip
set TS_NUDGE=90.                          #SST nuding [day] (*0 --> not execute)
set TI_NUDGE=90.                          #T nuding [day]
set SS_NUDGE=90.                          #SSS nuding [day]
set SI_NUDGE=90.                          #S nuding [day]

echo "Make sbPOM EXE"
cd ${MODELDIR}
make
cd ${CDIR}

if(! -f ${MODELDIR}/${PPROG})then
    echo "***Error: Not Find $PPROG"; exit
endif

#========= Fortran =================#

set cfortran=frtpx    #Computation node
set lfortran=gfortran #Login node
set option="-I/vol0001/apps/oss/spack-latest/opt/spack/linux-rhel7-haswell/gcc-4.8.5/netcdf-fortran-4.5.2-j2y55tgvcbna4c33sokcgzhiwr5pb5vs/include -L/vol0001/apps/oss/spack-latest/opt/spack/linux-rhel7-haswell/gcc-4.8.5/netcdf-fortran-4.5.2-j2y55tgvcbna4c33sokcgzhiwr5pb5vs/lib -lnetcdff -L/vol0001/apps/oss/spack-latest/opt/spack/linux-rhel7-haswell/gcc-4.8.5/netcdf-c-4.7.3-inuxrtaypqvhkczzj4n7mhe3o5zmtwnh/lib -lnetcdf -lnetcdf -ldl -lm" #Option for login node compiler 

#========= Make output dir =========#

@ IMEM = 1
while(${IMEM} <= ${MEMBER})
    set MEM=`printf "%03d" ${IMEM}`
    if(! -d ${OUTPUT}/anal/${MEM}) mkdir -p ${OUTPUT}/anal/${MEM}
    @ IMEM++
end

if(! -d ${OUTPUT}/gues/mean) mkdir -p ${OUTPUT}/gues/mean 
if(! -d ${OUTPUT}/gues/sprd) mkdir -p ${OUTPUT}/gues/sprd
if(! -d ${OUTPUT}/anal/mean) mkdir -p ${OUTPUT}/anal/mean 
if(! -d ${OUTPUT}/anal/sprd) mkdir -p ${OUTPUT}/anal/sprd

#========= Set DATE ===============#

@ DT = ${DT}

#Integration julian day
@ ST = `perl caldays.prl ${sdate[1]} ${sdate[2]} ${sdate[3]}` 

#Assimilation start julian day
@ IT = `perl caldays.prl ${idate[1]} ${idate[2]} ${idate[3]}`
set yyyy1=`printf "%04d" ${idate[1]}`;set mm1=`printf "%02d" ${idate[2]}`;set dd1=`printf "%02d" ${idate[3]}`

#Assimilation end julian day
@ ET = `perl caldays.prl ${edate[1]} ${edate[2]} ${edate[3]}` 
set yyyy2=`printf "%04d" ${edate[1]}`;set mm2=`printf "%02d" ${edate[2]}`;set dd2=`printf "%02d" ${edate[3]}`
#_____________________________________________________________________________#

#=== Prepare ENSENBLE MEMBER ============================
if(${switch_prepare_ens} == 1)then
    echo "Prepare ensemble member"
    csh prepare_ens.csh ${idate} ${MEMBER} ${MODELOUTPUTDIR} ${OUTPUT}
    set iexit=$?
    if(${iexit} == 99) exit 99
endif

#=== Make vcoordfile ====================================
@ IMEM = 1
while(${IMEM} <= ${MEMBER})
    set MEM=`printf "%03d" ${IMEM}`
    if(! -d ${OUTPUT}/vcoord/${MEM}) mkdir -p ${OUTPUT}/vcoord/${MEM}
    @ IMEM++
end
csh make_vcoord.csh ${OUTPUT}/vcoord ${MEMBER} ${ENODE} ${EPROC} ${lfortran}
set iexit=$?
if(${iexit} == 99) exit 99

#=== Compile check elapse time ==========================
#${cfortran} -o check_elapsetime.out ens/check_elapsetime.f90

#=== Make mesp_ens.out =========================================
${lfortran} ens/julian.f90 ens/mesp_ens.f90 -o mesp_ens.out ${option}
(rm -f julian_day.mod &)
#_______________________________________________________________________#

@ itime = 0
while(${itime} >= 0)    

    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo ">>> Start (${yyyy1}.${mm1}.${dd1} - ${yyyy2}.${mm2}.${dd2})>>>"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

    #=== Submit Loop Job ===========================================
    rm -f  ${LETKFDIR}/FINISHED ${LETKFDIR}/ERROR

    pjsub <<EOF
#PJM -L "node=${NODE}"
#PJM -L "rscunit=rscunit_ft01"
#PJM -L "rscgrp=eap-small"
#PJM -L "elapse=00:10:00"
#PJM --mpi "shape=${MEMBER}"
#PJM --mpi "max-proc-per-node=${EPROC}"
#PJM --name "sbPOM-LETKF_${REGION}_${yyyy1}${mm1}${dd1}_${yyyy2}${mm2}${dd2}"
#PJM -s

csh run.loop_vcoord.csh ${switch_iau} ${switch_anal} ${switch_gues} ${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${MEMBER} ${EPROC} ${ENODE} ${ETIME} ${OBSDIR} ${OBSFILE} ${LETKFDIR} ${LPROG} ${MODELDIR} ${MODELDATADIR} ${TIDEDIR} ${REGION} ${PPROG} ${NEST} ${BUDGET} ${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} ${DT} ${ST} ${IT} ${ET}

if [ \$? -eq 99 ]; then
    echo "***Error: Abnormal termination"    
    echo error > ${LETKFDIR}/ERROR
fi

EOF

    #=== Wait for LETKF end =======================================
    @ imin = 0
    while(${imin} >= 0)

	if(-f ${LETKFDIR}/ERROR)then
	    echo "***Error: Abnormal termination"
	    echo "${LETKFDIR}/sbPOM-LETKF_${REGION}_${yyyy1}${mm1}${dd1}_${yyyy2}${mm2}${dd2}" | mail -s "Error" ${EMAIL}
	    exit 99
	endif

	set date=`perl juldays.prl ${IT}`
	set yyyy=`printf "%04d" ${date[1]}`;set mm=`printf "%02d" ${date[2]}`;set dd=`printf "%02d" ${date[3]}`
	
	#---Make mean/sprd file & Remove Previous Ensemble Member---------
	if(-f ${LETKFDIR}/FINISHED${yyyy}${mm}${dd} && ${switch_gues} == 1 \
	    && ${switch_rmens} == 1 && ${IT} != ${ET})then
	    (cp mesp_ens.out mesp_ens${yyyy}${mm}${dd}.out && \
	    ./mesp_ens${yyyy}${mm}${dd}.out && \
	    rm -f mesp_ens${yyyy}${mm}${dd}.out && \
	    csh remove_ens.csh ${OUTPUT} ${MEMBER} ${IT} ${DT} ${REGION} &)
	else if(-f ${LETKFDIR}/FINISHED${yyyy}${mm}${dd} && ${switch_gues} == 1)then
	    (cp mesp_ens.out mesp_ens${yyyy}${mm}${dd}.out && \
	    ./mesp_ens${yyyy}${mm}${dd}.out && \
	    rm -f mesp_ens${yyyy}${mm}${dd}.out &)
	endif #switch_gues

	#----------------------------    
    
	#---Manage Timer--------------------------------------------
	if(-f ${LETKFDIR}/FINISHED${yyyy}${mm}${dd})then
	    (rm -f ${LETKFDIR}/FINISHED${yyyy}${mm}${dd} &)
	    @ IT++
	endif

	
	if(-f ${LETKFDIR}/FINISHED)then
	    #---Break loop
	    (rm -f ${LETKFDIR}/FINISHED &)
	    (mv ${LETKFDIR}/sbPOM-LETKF_${REGION}_${yyyy1}${mm1}${dd1}_${yyyy2}${mm2}${dd2}.* ${INFO}/ &)
	    echo "=== Finished sbPOM-LETKF_${REGION}_${yyyy1}${mm1}${dd1}_${yyyy2}${mm2}${dd2} until ${yyyy}${mm}${dd}; ${imin} min ==="
	    break
	else
	    #---Continue loop & Count time 
	    sleep 60s
	    @ itime++
	    @ imin++
	    echo "All job; ${itime} min passed"
	    echo "sbPOM-LETKF_${REGION}_${yyyy1}${mm1}${dd1}_${yyyy2}${mm2}${dd2}; ${imin} min passed"
	endif

    end
    #===================================================================

    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo ">>> END (${yyyy1}.${mm1}.${dd1} - ${yyyy}.${mm}.${dd}) >>>"
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

    if(${IT} <= ${ET})then
	#---Modify idate
	set idate=`perl juldays.prl ${IT}`
	set yyyy1=`printf "%04d" ${idate[1]}`;set mm1=`printf "%02d" ${idate[2]}`;set dd1=`printf "%02d" ${idate[3]}`
    else
	#---Break loop
	echo "===Finished All job; ${itime} min ==="
	break
    endif

end

#=== Clean & Wait ===========================================================
#rm -f mesp_ens.out check_elapsetime.out
wait

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo ">>>>    All    JOB    END    >>>>"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
exit 0

