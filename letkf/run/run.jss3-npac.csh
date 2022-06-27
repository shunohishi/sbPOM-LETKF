#!bin/csh
#======= DATE ===========================================#

set sdate=(2011 1 1) #Integration start date
set idate=(2017 4 22) #Assimilation start date
set edate=(2018 12 31) #Assimilation end date

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo ">>> Start (${idate} - ${edate})>>>"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

#========== Switch(1:execute, other: skip)================#

set switch_prepare_ens=1 #Prepare ensemble member
set switch_iau=1         #Integration of sbPOM for IAU procedure
set switch_anal=1        #LETKF/Make analysis
set switch_gues=1        #sbPOM integration/Make guess 
set switch_rmens=1       #Remove ensemble member (except end day of each month)

#========= GENERAL ======================================#

set DIR0=/data/R/R2402/ohishi/NPAC
set LDIR0=${DIR0}/letkf-jra55do
set CDIR=`pwd`             #Current DIR
set WORKDIR=${LDIR0}/work  #Work DIR
set OUTPUT=${LDIR0}/output #Output DIR
set INFO=${LDIR0}/run/info #JSS2 INFORMATION 
set MEMBER=100             #Ensembel member
set NPROC=48               #Proceccor for each node
set DT=1                   #Delta T (unit: day)
set NTRY=3                 #Number of trial
set SJOB_TYPE=2            #1:Serial Job, #2:BULK Job,#3 VCOORD Job

#========= OBSERVATION ====================================#

set OBSDIR=${DIR0}/prep/obs #Observational DIR for Assimilation
set OBSFILE=obs-jra55do             #Observation filename: ${OBSFILE}${yyyy}${mm}${dd}.nc

#========== LETKF =========================================#

set LETKFDIR=${LDIR0}/run #LETKF DIR
set LPROG=letkf.exe       #LETKF EXE File
set LPROC=4               #Number of Processor for each thread
set LNODE=100              #Number of Node
set LETIME=30             #Elapse time [min.]
set LWTIME=360             #Wait time[min.]

#========== sbPOM =========================================#

set MODELDIR=${DIR0}/pom-ens-jra55do/run          #sbPOM DIR
set MODELDATADIR=${DIR0}/prep/in          #grid,tsclim,ic,lbc,atm,tsdata DIR
set MODELOUTPUTDIR=${DIR0}/pom-ens-jra55do/output #sbPOM output DIR
set TIDEDIR=${DIR0}/pom-ens-jra55do/TIDE          #Tide DIR
set REGION=npac
set PPROC=4                               #Number of Processor
set PNODE=1                               #Number of Node
set PETIME=10                              #Elapse time [min.]
set PWTIME=360                            #Wait time [min.]
set PPROG=pom.exe                         #EXE file
set NEST=0                                #One-way nest 1:on, 0:off
set BUDGET=1                              #1:execute, 0:skip
set TS_NUDGE=548.                          #SST nuding [day] (*0 --> not execute)
set TI_NUDGE=548.                          #T nuding [day]
set SS_NUDGE=90.                         #SSS nuding [day]
set SI_NUDGE=548.                          #S nuding [day]
set ATMNAME="jra55do"                   #jra55/jra55do

#========= Fortran =========================================#

set lfortran=ifort #Login node
set option="-assume byterecl -convert big_endian -mcmodel=medium -shared-intel -I/opt/JX/oss/x86_64/netcdf-fortran/4.5.3/include -I/opt/JX/oss/x86_64/netcdf/4.7.3/include -I/opt/JX/oss/x86_64/hdf5/1.12.0/include -L/opt/JX/oss/x86_64/netcdf-fortran/4.5.3/lib -lnetcdff -L/opt/JX/oss/x86_64/hdf5/1.12.0/lib -L/opt/JX/oss/x86_64/netcdf/4.7.3/lib -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -lm -lz -lcurl" #Option for login node compiler 

#========= Make dir =========================================#

###OUTPUT
echo "Check OUTPUT Directory"
if(! -d ${OUTPUT})then
    mkdir -p ${OUTPUT}
endif

@ IMEM = 1
while(${IMEM} <= ${MEMBER})
    set MEM=`printf "%03d" ${IMEM}`
    if(! -d ${OUTPUT}/anal/${MEM}) mkdir -p ${OUTPUT}/anal/${MEM}
    @ IMEM++
end

if(! -d ${OUTPUT}/gues/mean) mkdir -p ${OUTPUT}/gues/mean 
if(! -d ${OUTPUT}/gues/sprd) mkdir -p ${OUTPUT}/gues/sprd
if(! -d ${OUTPUT}/gues/eens) mkdir -p ${OUTPUT}/gues/eens
if(! -d ${OUTPUT}/anal/mean) mkdir -p ${OUTPUT}/anal/mean 
if(! -d ${OUTPUT}/anal/sprd) mkdir -p ${OUTPUT}/anal/sprd

###WORKDIR
echo "Make Work directory"
rm -rf ${WORKDIR} && mkdir -p ${WORKDIR}

###LETKF
echo "Make LETKF EXE & Log directory"
cd ${LETKFDIR}
make
cd ${CDIR}
if(! -f ${LETKFDIR}/${LPROG})then
    echo "***Error: Not found ${LPROG}"
    exit 99
endif
if(! -d ${INFO}) mkdir -p ${INFO}

###sbPOM
echo "Make sbPOM EXE"
cd ${MODELDIR}
make
cd ${CDIR}

if(! -f ${MODELDIR}/${PPROG})then
    echo "***Error: Not Find $PPROG"
    exit 99
endif

#========= Set DATE =====================================#
@ DT = ${DT}

#Integration julian day
@ ST = `perl caldays.prl ${sdate[1]} ${sdate[2]} ${sdate[3]}` 
set yyyys=`printf "%04d" ${sdate[1]}`;set mms=`printf "%02d" ${sdate[2]}`;set dds=`printf "%02d" ${sdate[3]}`

#Assimilation start julian day
@ IT = `perl caldays.prl ${idate[1]} ${idate[2]} ${idate[3]}` 

#Assimilation end julian day
@ ET = `perl caldays.prl ${edate[1]} ${edate[2]} ${edate[3]}` 


#=== Prepare ENSENBLE MEMBER ============================#
if(${switch_prepare_ens} == 1)then
    echo "Prepare ensemble member"
    csh prepare_ens.csh ${idate} ${MEMBER} ${MODELOUTPUTDIR} ${OUTPUT}
    set iexit=$?
    if(${iexit} == 99)then
	echo "***Error: Prepare ensemble member"
	exit 99
    endif
endif

#=== Make mesp_ens.out ==================================#
${lfortran} ens/julian.f90 ens/mesp_ens.f90 -o mesp_ens.out ${option}
if(! -f mesp_ens.out)then
    echo "***Error: Not found mesp_ent.out"
    exit 99
endif
(rm -f julian_day.mod &)

#=== Main Loop ==========================================#
while(${IT} <= ${ET})

    set date=`perl juldays.prl ${IT}`
    if(${date[2]} == 12)then
	@ yyyy2 = ${date[1]} + 1
	@ mm2 = 1
    else
	@ yyyy2 = ${date[1]}
	@ mm2 = ${date[2]} + 1
    endif
    set yyyy1=`printf "%04d" ${date[1]}`;set mm1=`printf "%02d" ${date[2]}`; set dd1=`printf "%02d" ${date[3]}`
    set yyyy2=`printf "%04d" ${yyyy2}`;set mm2=`printf "%02d" ${mm2}`

    #===MODEL Setting data===#

    set GRID=grid               #Model Grid data
    set TSCLIM=tsclim           #T/S Climatology from WOA
    set IC=ic.woa18.${mms}      #Initial Condition Netcdf file
    if($NEST == 0)then
	set LBC=lbc_mclim       #Lateral Boundary Condition Netcdf file
	set TSDATA=tsdata_mclim # Nudging Netcdf file
    else if($NEST == 1)then
	set LBC=lbc-assim.${yyyy1}${mm1}01.${yyyy2}${mm2}01
	set TSDATA=tsdata-assim.${yyyy1}${mm1}01.${yyyy2}${mm2}01
    endif
    set FFLUX=fflux.${yyyy1}${mm1}01.${yyyy2}${mm2}01 #Freshwater flux Netcdf file

    foreach file(${GRID}.nc ${TSCLIM}.nc ${IC}.nc ${FFLUX}.nc)
	if(-f ${MODELDATADIR}/${file})then
	    echo "Find ${file}"
	else
	    echo "***Error: Not find ${file}"
	    exit 99
	endif
    end

    #Atmospheric Condition Netcdf file
    set ATM=${ATMNAME}.${yyyy1}${mm1}01.${yyyy2}${mm2}01

    @ IMEM = 1
    while($IMEM <= $MEMBER)
	set MEM=`printf "%03d" ${IMEM}`
	foreach file(${LBC} ${TSDATA} ${ATM})
	    set file=${MODELDATADIR}/${file}.${MEM}.nc
	    if(-f ${file})then
		echo "Find ${file}"
	    else
		echo "***Error: Not find ${file}"
		exit 99    
	    endif
	end
	@ IMEM++
    end

    echo ">>>"
    echo ">>>BIGIN ${date} / ${edate}"
    echo ">>>"

    #=== ENSEMBLE FORECAST (run sbPOM) for IAU ===============
    if(${switch_iau} == 1)then
    
	echo " ====================================================="
	echo " === START ENSEMBLE PREDICTION FOR IAU at ${date}  ==="
	echo " === NUMBER of ENSEMBLE: ${MEMBER}                 ==="
	echo " ====================================================="
	 
	echo "Submit sbPOM Shell" #------------------------------------------
	if(${SJOB_TYPE} == 1)then
	    set CSHFILE=run.ensfcst_iau.csh
	else if(${SJOB_TYPE} == 2)then
	    set CSHFILE=run.ensfcst_iau_bulk.csh
	else if(${SJOB_TYPE} == 3)then
	    set CSHFILE=run.ensfcst_iau_vcoord.csh
	else
	    echo "***Error: csh file"
	    exit 99
	endif

	@ iexe = 1
	while(${iexe} > 0)
	    csh ${CSHFILE} ${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${MEMBER} ${NPROC} ${MODELDIR} ${MODELDATADIR} ${TIDEDIR} ${REGION} ${PPROC} ${PNODE} ${PETIME} ${PWTIME} ${PPROG} ${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} ${GRID} ${TSCLIM} ${IC} ${LBC} ${ATM} ${FFLUX} ${TSDATA} ${DT} ${ST} ${IT}
	    set iexit=$?
	    if(${iexit} == 88)then
                if(${iexe} == ${NTRY}) exit 99
		rm -rf ${WORKDIR}/*
		@ iexe++
	    else if(${iexit} == 99)then
		echo "***Error: sbPOM-IAU"
		exit 99
	    else
		break
	    endif
	end
	#------------------------------------------------------------------

	echo "Clean"#--------------
	rm -rf ${WORKDIR}/*
	#--------------------------

	echo "Move Current DIR"#---
	cd ${CDIR}
	#--------------------------

	echo " ==================================================="
	echo " === END ENSEMBLE PREDICTION FOR IAU at ${date}  ==="
	echo " === NUMBER of ENSEMBLE: ${MEMBER}               ==="
	echo " ==================================================="

    endif #switch_iau
    #=========================================================
    

    #=== LETKF ===============================================
    if(${switch_anal} == 1)then

	echo "========================================================="
	echo "=== START LETKF DATA ASSIMILATION at ${date}          ==="
	echo "========================================================="

	echo "Submit LETKF Shell" #----------------------------------
	set CSHFILE=run.letkf.csh

	@ iexe = 1
	while(${iexe} > 0)
	    csh ${CSHFILE} ${switch_iau} ${WORKDIR} ${OUTPUT} ${INFO} ${MEMBER} ${NPROC} ${OBSDIR} ${OBSFILE} ${LETKFDIR} ${LPROG} ${LPROC} ${LNODE} ${LETIME} ${LWTIME} ${MODELDATADIR} ${GRID} ${IT} ${REGION}
	    set iexit=$?
	    if(${iexit} == 88)then
                if(${iexe} == ${NTRY}) exit 99
		rm -rf ${WORKDIR}/*
		@ iexe++
	    else if(${iexit} == 99)then
		echo "***Error: LETKF"
		exit 99
	    else
		break
	    endif
	end
	#----------------------------------------------------------

	echo "Clean"#--------------
	rm -rf ${WORKDIR}/*
	#--------------------------

	echo "Move Current DIR"#---
	cd ${CDIR}
	#--------------------------

	echo "======================================================"
	echo "=== END LETKF DATA ASSIMILATION at ${date}        ===="
	echo "======================================================"

    endif #switch_anal
    #================================================================

    #=== ENSEMBLE FORECAST (run sbPOM) ==============================
    if(${switch_gues} == 1)then

	echo " ====================================================="
	echo " === START ENSEMBLE PREDICTION at ${date}          ==="
	echo " === NUMBER of ENSEMBLE: ${MEMBER}                 ==="
	echo " ====================================================="

	echo "Submit sbPOM Shell" #------------------------------------------
	if(${SJOB_TYPE} == 1)then
	    set CSHFILE=run.ensfcst.csh
	else if(${SJOB_TYPE} == 2)then
	    set CSHFILE=run.ensfcst_bulk.csh
	else if(${SJOB_TYPE} == 3)then
	    set CSHFILE=run.ensfcst_vcoord.csh
	else
	    echo "***Error: csh file"
	    exit 99
	endif

	@ iexe = 1
	while($iexe > 0)
	    csh ${CSHFILE} ${switch_iau} ${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${MEMBER} ${NPROC} ${MODELDIR} ${MODELDATADIR} ${TIDEDIR} ${REGION} ${PPROC} ${PNODE} ${PETIME} ${PWTIME} ${PPROG} ${BUDGET} ${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} ${GRID} ${TSCLIM} ${IC} ${LBC} ${ATM} ${FFLUX} ${TSDATA} ${DT} ${ST} ${IT}
	    set iexit=$?
	    if(${iexit} == 88)then
                if(${iexe} == ${NTRY}) exit 99
		rm -rf ${WORKDIR}/*
		@ iexe++
	    else if(${iexit} == 99)then
		echo "***Error: sbPOM"
		exit 99
	    else
		break
	    endif
	end

	echo "Clean"#--------------
	rm -rf ${WORKDIR}/*
	#--------------------------

	echo "Move Current DIR"#---
	cd ${CDIR}
	#--------------------------

	echo " ==================================================="
	echo " === END ENSEMBLE PREDICTION at ${date}          ==="
	echo " === NUMBER of ENSEMBLE: ${MEMBER}               ==="
	echo " ==================================================="
	
    endif #switch_gues
    #===================================================================

    #=== Ensemble mean/sprd & Remove each ensemble member #============
    cp mesp_ens.out mesp_ens${yyyy1}${mm1}${dd1}.out
    if(${switch_gues} == 1 && ${switch_rmens} == 1 && ${IT} != ${ET})then
	echo "Ensmble mean/sprd & Remove each Ensmble member"
	(./mesp_ens${yyyy1}${mm1}${dd1}.out && \
	 rm -f mesp_ens${yyyy1}${mm1}${dd1}.out && \
	 csh remove_ens.csh ${OUTPUT} ${MEMBER} ${IT} ${DT} ${REGION} &)
    else if(${switch_gues} == 1)then
	echo "Ensmble mean/sprd"
	(./mesp_ens${yyyy1}${mm1}${dd1}.out && \
	 rm -f mesp_ens${yyyy1}${mm1}${dd1}.out &)
    endif #switch_rmens

    echo "Move Current DIR"#---
    cd ${CDIR}
    #--------------------------
    #=====================================================================


    @ IT = ${IT} + ${DT}

end

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo ">>> END (${idate} - ${edate}) >>>"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
rm -f mesp_ens.out
wait
exit 0
