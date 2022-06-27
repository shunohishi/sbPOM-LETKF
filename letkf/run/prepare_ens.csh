#!/bin/csh
#-----------------------------------------------------------------
# Prepare Ensemble member 
#-----------------------------------------------------------------

set date=(${argv[1]} ${argv[2]} ${argv[3]})
set MEMBER=${argv[4]}
set MODELOUTPUTDIR=${argv[5]}
set OUTPUT=${argv[6]}

#-----------------------------------------------------------------

set julian=`perl caldays.prl ${date[1]} ${date[2]} ${date[3]}`
@ ijul = ${julian} - 1
set date0=`perl juldays.prl ${ijul}`
set yyyy0=`printf "%04d" ${date0[1]}`; set mm0=`printf "%02d" ${date0[2]}`; set dd0=`printf "%02d" ${date0[3]}`

#Copy restart file
@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`

    if(! -d ${OUTPUT}/gues/${MEM}) mkdir -p ${OUTPUT}/gues/${MEM}
    if(-f ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc)then
	echo "Exist: ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc"
    else if(${IMEM} == ${MEMBER})then
	cp -v ${MODELOUTPUTDIR}/${yyyy0}${mm0}/${MEM}/out/ens.${dd0}.nc ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc &
	wait
    else
	cp -v ${MODELOUTPUTDIR}/${yyyy0}${mm0}/${MEM}/out/ens.${dd0}.nc ${OUTPUT}/gues/${MEM}/restart.${yyyy0}${mm0}${dd0}.nc &
    endif

@ IMEM++

end
wait
exit 0
