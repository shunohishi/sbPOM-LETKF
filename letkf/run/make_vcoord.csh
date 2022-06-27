#!/bin/csh
#--------------------------------------------------------------------------
# Make vcoordfile |
#------------------
# 2020.03 Created by S.Ohishi
#--------------------------------------------------------------------------
set WORKDIR=$argv[1];set MEMBER=$argv[2];set ENODE=$argv[3]; set EPROC=$argv[4]
set fortran="gfortran"
#--------------------------------------------------------------------------

echo "Make vcoordfile" #--------------------------------------------------
echo ${WORKDIR} > vcoordfile_info.txt
echo ${MEMBER} >> vcoordfile_info.txt
echo ${ENODE} >> vcoordfile_info.txt
echo ${EPROC} >> vcoordfile_info.txt
${fortran} ens/make_vcoordfile.f90 -o make_vcoordfile.out
./make_vcoordfile.out
rm -f make_vcoordfile.out
#--------------------------------------------------------------------------

echo "Check for Making vcoordfile" #------------------------------------------
@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`
    set vcoordfile=${WORKDIR}/${MEM}/vcoordfile
    set nline=`wc -l ${vcoordfile} | awk '{print $1}'`

    if(${EPROC} == ${nline})then
	echo "vcoordfile FINISHED [ ${IMEM} / ${MEMBER}]"
	@ IMEM++
    else
	echo "***Error @ ${MEM} vcoordfile";exit 99
    endif

end    
#--------------------------------------------------------------------------
exit 0
