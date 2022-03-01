#!/bin/csh
#--------------------------------------------------------------------------
# Make vcoordfile |
#------------------
# 2020.03 Created by S.Ohishi
#--------------------------------------------------------------------------
set WORKDIR0=$argv[1];set MEMBER=$argv[2];set ENODE=$argv[3]; set EPROC=$argv[4]
#--------------------------------------------------------------------------

#echo "Clean vcoordfile" #--------------------------------------------------
#@ IMEM = 1
#while(${IMEM} <= ${MEMBER})
#    set MEM=`printf "%03d" ${IMEM}`
#    set vcoordfile=${WORKDIR0}/${MEM}/vcoordfile
#    rm -f ${vcoordfile}
#    @ IMEM++
#end
#-------------------------------------------------------------------------

#echo "Make vcoordfile" #--------------------------------------------------
#@ IMEM = 1
#@ PN = $EPROC / $ENODE
#@ INUM = 0
#while(${IMEM} <= ${MEMBER})
#
#    set MEM=`printf "%03d" ${IMEM}`
#    set vcoordfile=${WORKDIR0}/${MEM}/vcoordfile
#
#    @ INODE = 1
#    while(${INODE} <= ${ENODE})
#        @ IPN = 1
#	while(${IPN} <= ${PN})
#	    echo  "(${INUM})" >> ${vcoordfile}
#	    @ IPN++
#	end
#	@ INODE++
#	@ INUM++
#    end
#    @ IMEM++
#end
echo "Make vcoordfile" #--------------------------------------------------
echo ${WORKDIR0} > vcoordfile_info.txt
echo ${MEMBER} >> vcoordfile_info.txt
echo ${ENODE} >> vcoordfile_info.txt
echo ${EPROC} >> vcoordfile_info.txt
gfortran ens/make_vcoordfile.f90 -o make_vcoordfile.out
./make_vcoordfile.out
rm -f make_vcoordfile.out
#--------------------------------------------------------------------------

echo "Check for Making vcoordfile" #------------------------------------------
@ IMEM = 1
while(${IMEM} <= ${MEMBER})

    set MEM=`printf "%03d" ${IMEM}`
    set vcoordfile=${WORKDIR0}/${MEM}/vcoordfile
    set nline=`wc -l ${vcoordfile} | awk '{print $1}'`

    if(${EPROC} == ${nline})then
	echo "vcoordfile FINISHED [ ${IMEM} / ${MEMBER}]"
	@ IMEM++
    else
	echo "***Error @ ${MEM} vcoordfile"
	exit 99
    endif

end    
#--------------------------------------------------------------------------
