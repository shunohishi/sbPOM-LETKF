#!/bin/csh
#=== ARGUMENT =================================
set OUTPUT=${argv[1]}; set MEMBER=${argv[2]}
set IT=${argv[3]};set DT=${argv[4]};set REGION=${argv[5]}

#=== JULIAN DATE ==============================

@ BT = ${IT} - ${DT}
set date0=`perl juldays.prl ${BT}`
set yyyy0=`printf "%04d" ${date0[1]}`;set mm0=`printf "%02d" ${date0[2]}`;set dd0=`printf "%02d" ${date0[3]}`

set date1=`perl juldays.prl ${IT}`
set yyyy1=`printf "%04d" ${date1[1]}`;set mm1=`printf "%02d" ${date1[2]}`;set dd1=`printf "%02d" ${date1[3]}`

#== Check ENSEBMLE MEMBER =====================

@ IMEM = 1
@ CMEM = 0
while(${IMEM} <= ${MEMBER})
    set MEM=`printf "%03d" ${IMEM}`

    if(-f ${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc \
    && -f ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc \
    && -f ${OUTPUT}/gues/${MEM}/restart.${yyyy1}${mm1}${dd1}.nc)then
	@ CMEM++
    endif
    
    if(! -f ${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc)then
	echo "***Not FOUND ${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc"
    endif

    if(! -f ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc)then
	echo "***Not FOUND ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc"
    endif

    if(! -f ${OUTPUT}/gues/${MEM}/restart.${yyyy1}${mm1}${dd1}.nc)then
	echo "***Not FOUND ${OUTPUT}/gues/${MEM}/restart.${yyyy1}${mm1}${dd1}.nc"
    endif

    @ IMEM++
end

echo "Check ENSEMBLE MEMBER: ${CMEM}/${MEMBER}"

#=== Remove ENSEMBLE ==========================
if(${CMEM} != ${MEMBER})then
    echo "***Do not remove because of non-exist all files***"
else

    echo "Remove ensemble file"
    @ IMEM = 1
    while(${IMEM} <= ${MEMBER})
	set MEM=`printf "%03d" ${IMEM}`

	if(-f ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}_iau.nc)then
	    #echo "Remove gues/${MEM}/restart.${yyyy0}${mm0}${dd0}_iau.nc"
	    (rm ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}_iau.nc &)
	endif

	if(${date1[3]} == 1)then
            #echo "(Skip to remove because of the 1st day)"
        else if(${date1[3]} == 16)then
            #echo "(Skip to remove because of the 16th day)"
	else if(-f ${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc)then
	    #echo "Remove anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc"
	    (rm ${OUTPUT}/anal/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc &)
	endif

	if(${date1[3]} == 1)then
	    #echo "(Skip to remove because of the 1st day)"
	else if(${date1[3]} == 16)then
	    #echo "(Skip to remove because of the 16th day)"
	else if(-f ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc)then
	    #echo "Remove gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc"
	    (rm ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc &)
	endif

	if(-f ${OUTPUT}/gues/${MEM}/${REGION}.${yyyy0}${mm0}${dd0}.nc)then
	    #echo "Remove gues/${MEM}/${REGION}.${yyyy0}${mm0}${dd0}.nc"
	    (rm ${OUTPUT}/gues/${MEM}/${REGION}.${yyyy0}${mm0}${dd0}.nc &)
	endif

	@ IMEM++
    end

    if(-f ${OUTPUT}/gues/mean/restart.${yyyy0}${mm0}${dd0}.nc)then
	#echo "Remove gues/mean/restart.${yyyy0}${mm0}${dd0}.nc"
	(rm ${OUTPUT}/gues/mean/restart.${yyyy0}${mm0}${dd0}.nc &)
    endif

    if(-f ${OUTPUT}/gues/sprd/restart.${yyyy0}${mm0}${dd0}.nc)then
	#echo "Remove gues/sprd/restart.${yyyy0}${mm0}${dd0}.nc"
	(rm ${OUTPUT}/gues/sprd/restart.${yyyy0}${mm0}${dd0}.nc &)
    endif

    if(${date1[3]} == 1)then
	#echo "(Skip to remove because of the 1st day)"
    else if(${date1[3]} == 16)then
       	#echo "(Skip to remove because of the 16th day)"
    else if(-f ${OUTPUT}/anal/mean/restart.${yyyy0}${mm0}${dd0}.nc)then
	#echo "Remove anal/mean/restart.${yyyy0}${mm0}${dd0}.nc"
	(rm ${OUTPUT}/anal/mean/restart.${yyyy0}${mm0}${dd0}.nc &)
    endif
	
    if(${date1[3]} == 1)then
        #echo "(Skip to remove because of the 1st day)"
    else if(${date1[3]} == 16)then
        #echo "(Skip to remove because of the 16th day)"
    else if(-f ${OUTPUT}/anal/sprd/restart.${yyyy0}${mm0}${dd0}.nc)then
	#echo "Remove anal/sprd/restart.${yyyy0}${mm0}${dd0}.nc"
	(rm ${OUTPUT}/anal/sprd/restart.${yyyy0}${mm0}${dd0}.nc &)
    endif

endif

