#!/bin/csh
#===Argument===================================================================
set switch_iau=${argv[1]};set switch_anal=${argv[2]};set switch_gues=${argv[3]}
set CDIR=${argv[4]};set WORKDIR=${argv[5]};set OUTPUT=${argv[6]};set INFO=${argv[7]}
set MEMBER=${argv[8]};set EPROC=${argv[9]};set ENODE=${argv[10]}; set ETIME=${argv[11]}
set OBSDIR=${argv[12]};set OBSFILE=${argv[13]}
set LETKFDIR=${argv[14]};set LPROG=${argv[15]}
set MODELDIR=${argv[16]};set MODELDATADIR=${argv[17]};set TIDEDIR=${argv[18]}
set REGION=${argv[19]};set PPROG=${argv[20]};set NEST=${argv[21]};set BUDGET=${argv[22]}
set TS_NUDGE=${argv[23]};set TI_NUDGE=${argv[24]};set SS_NUDGE=${argv[25]};set SI_NUDGE=${argv[26]}
set DT=${argv[27]};set ST=${argv[28]};set IT=${argv[29]};set ET=${argv[30]}

#Julian day --> date
set sdate=`perl juldays.prl ${ST}`;set edate=`perl juldays.prl ${ET}`
set yyyys=`printf "%04d" ${sdate[1]}`;set mms=`printf "%02d" ${sdate[2]}`;set dds=`printf "%02d" ${sdate[3]}`
#==============================================================================

while(${IT} <= ${ET})

    set date=`perl juldays.prl ${IT}`
    if(${date[2]} == 12)then
	@ yyyy2 = ${date[1]} + 1
	@ mm2 = 1
    else
	@ yyyy2 = ${date[1]}
	@ mm2 = ${date[2]} + 1
    endif
    set yyyy1=`printf "%04d" ${date[1]}`;set mm1=`printf "%02d" ${date[2]}`;set dd1=`printf "%02d" ${date[3]}`
    set yyyy2=`printf "%04d" ${yyyy2}`;set mm2=`printf "%02d" ${mm2}`

    #===MODEL Setting data===#
    set GRID=grid #Model Grid data
    set TSCLIM=tsclim #T/S Climatology from WOA
    set IC=ic.woa18.${mms} #Initial Condition Netcdf file
    if($NEST == 0)then
	set LBC=lbc_mclim #Lateral Boundary Condition Netcdf file
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
	    echo "***Error: Not find ${file}"; exit 99
	endif
    end

    #Atmospheric Condition Netcdf file
    set ATM=jra55.${yyyy1}${mm1}01.${yyyy2}${mm2}01

    @ IMEM = 1
    while(${IMEM} <= ${MEMBER})
	set MEM=`printf "%03d" ${IMEM}`
	foreach file(${LBC} ${TSDATA} ${ATM})
	    set file=${MODELDATADIR}/${file}.${MEM}.nc
	    if(-f ${file})then
		echo "Find ${file}"
	    else
		echo "***Error: Not find ${file}"; exit 99
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
	csh run.ensfcst_iau_vcoord.csh ${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${MEMBER} ${EPROC} ${MODELDIR} ${MODELDATADIR} ${TIDEDIR} ${REGION} ${PPROG} ${NEST} ${BUDGET} ${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} ${GRID} ${TSCLIM} ${IC} ${LBC} ${ATM} ${FFLUX} ${TSDATA} ${DT} ${ST} ${IT}
	set iexit=$?
	if(${iexit} == 99) exit 99

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
    #=============================================================

    #=== LETKF ===================================================
    if(${switch_anal} == 1)then

	echo "========================================================="
	echo "=== START LETKF DATA ASSIMILATION at ${date}          ==="
	echo "========================================================="

	echo "Submit LETKF Shell" #----------------------------------
	csh run.letkf_vcoord.csh ${switch_iau} ${WORKDIR} ${OUTPUT} ${INFO} ${MEMBER} ${EPROC} ${OBSDIR} ${OBSFILE} ${LETKFDIR} ${LPROG} ${MODELDATADIR} ${REGION} ${GRID} ${IT}
	set iexit=$?
	if(${iexit} == 99) exit 99
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
	csh run.ensfcst_vcoord.csh ${switch_iau} ${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${MEMBER} ${EPROC} ${MODELDIR} ${MODELDATADIR} ${TIDEDIR} ${REGION} ${PPROG} ${NEST} ${BUDGET} ${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} ${GRID} ${TSCLIM} ${IC} ${LBC} ${ATM} ${FFLUX} ${TSDATA} ${DT} ${ST} ${IT}
	set iexit=$?
	if(${iexit} == 99) exit 99

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

    #=== Write daily FINISHED log===========================================
    echo finished > ${LETKFDIR}/FINISHED${yyyy1}${mm1}${dd1}
    #=======================================================================

    #=== Check Elapse time==================================================
    #./check_elapsetime.out
    #set elapsetime=`gawk '{print $1}' check_elapsetime.dat`
    #(rm -f check_elapsetime.dat &)
    #if(${elapsetime} <= ${ETIME}) break
    #=======================================================================

@ IT = ${IT} + ${DT}

end

#=== Write FINISHED log =================================================
echo finished > ${LETKFDIR}/FINISHED
#========================================================================

wait
exit 0

