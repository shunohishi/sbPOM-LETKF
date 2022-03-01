#=======================================================================
# LONG INTEGRATION VERSION
#=======================================================================

set REGION=test6
set EXE=pom.exe
set DIR0=/data/R/R2402/ohishi/TEST6/pom-ens
set TIDEDIR=${DIR0}/TIDE
set NCDIR=/data/R/R2402/ohishi/TEST6/prep/in
set CURDIR=`pwd`
set EPROC=32 #Each processor
set MEMBER=40 #Number of ensemble
set ENODE=1 #Each Node
set NEST=0 #One-way nest 1: On, 0: Off
set BUDGET=1 #Heat/Salinity budget convervation 1:On, 0:Off
set RM_ENS=1 #Remove ensemble member 1: On, 0: Off

set sdate=(2011 1) #start time
set idate=(2016 1) #initial time
set edate=(2016 1) #end time

#=========================================================================

echo "Compile"
#make clean
make
if(! -f ${EXE})exit

@ iyr = ${idate[1]}
@ imon = ${idate[2]}

set syear=`printf "%04d" ${sdate[1]}`;set smonth=`printf "%02d" ${sdate[2]}` 

while(${iyr} <= ${edate[1]})

if(${iyr} == ${edate[1]})then
    @ emon = ${edate[2]}
else
    @ emon = 12
endif

    while(${imon} <= ${emon})

    #-----------------------------------------
    # Setting year,month,day
    #-----------------------------------------

    if(${imon} == 12)then
	@ amon = 1
	@ ayr = ${iyr} + 1
    else
	@ amon = ${imon} + 1
	@ ayr = ${iyr}
    endif

    if(${imon} == 1)then
	@ bmon = 12
	@ byr = ${iyr} - 1
    else
	@ bmon = ${imon} - 1
	@ byr = ${iyr}
    endif

    set byear=`printf "%04d" ${byr}`;set bmonth=`printf "%02d" ${bmon}`
    set year=`printf "%04d" ${iyr}`;set month=`printf "%02d" ${imon}`
    set ayear=`printf "%04d" ${ayr}`;set amonth=`printf "%02d" ${amon}` 

    if(${imon} == 2 && ${iyr} % 4 == 0)then
	set nday=29.
    else if(${imon} == 2 && ${iyr} % 4 != 0)then
	set nday=28.
    else if(${imon} == 4 || ${imon} == 6 || ${imon} == 9 || ${imon} == 11)then
	set nday=30.
    else
	set nday=31.
    endif

    #-----------------------------------
    echo "Start ${year}${month}"
    #-----------------------------------

    echo "Make Workdir"

    set WORKDIR0=${DIR0}/output/${year}${month}
    set WORKDIRB=${DIR0}/output/${byear}${bmonth}
    rm -rf ${WORKDIR0}

    if(! -d ${WORKDIR0}/mean/out) mkdir -p ${WORKDIR0}/mean/out
    if(! -d ${WORKDIR0}/sprd/out) mkdir -p ${WORKDIR0}/sprd/out

	@ IMEM = 1
	while(${IMEM} <= ${MEMBER})

	    set MEM=`printf "%03d" ${IMEM}`

	    echo "Make in/out DIR: $MEM" #----------------
	    foreach dir(in out)
		if( ! -d ${WORKDIR0}/${MEM}/${dir})then
		    mkdir -p ${WORKDIR0}/${MEM}/${dir}
		endif
	    end
	    #----------------------------------------------

	    echo "Link Netcdf file" #----------------------
	    #rm -f ${WORKDIR0}/in/*.nc
	    #---grid.nc
	    (ln -s ${NCDIR}/grid.nc ${WORKDIR0}/${MEM}/in/${REGION}.grid.nc &)
	    
	    #---tsclim.nc
	    (ln -s ${NCDIR}/tsclim.nc ${WORKDIR0}/${MEM}/in/${REGION}.tsclim.nc &)
	    
	    #---ic.nc
	    (ln -s ${NCDIR}/ic.woa13.${smonth}.nc ${WORKDIR0}/${MEM}/in/${REGION}.ic.nc &)
	    #ln -s ${NCDIR}/ic.soda.${smonth}.nc ${WORKDIR0}/in/${REGION}.ic.nc
	    
	    #---lbc/tsdata.nc
	    if(${NEST} == 0)then
		(ln -s ${NCDIR}/tsdata_mclim.${MEM}.nc ${WORKDIR0}/${MEM}/in/${REGION}.tsdata.nc &)
		(ln -s ${NCDIR}/lbc_mclim.${MEM}.nc ${WORKDIR0}/${MEM}/in/${REGION}.lbc.nc &)
	    else if(${NEST} == 1)then
		(ln -s ${NCDIR}/lbc.${year}${month}01.${ayear}${amonth}01.nc ${WORKDIR0}/${MEM}/in/${REGION}.lbc.nc &)
		(ln -s ${NCDIR}/tsdata.${year}${month}01.${ayear}${amonth}01.nc ${WORKDIR0}/${MEM}/in/${REGION}.tsdata.nc &)
	    endif
    
	    #ensemble atm.nc
	    (ln -s ${NCDIR}/jra55.${year}${month}01.${ayear}${amonth}01.${MEM}.nc ${WORKDIR0}/${MEM}/in/${REGION}.atm.nc &)
	    
	    #fflux.nc
	    (ln -s ${NCDIR}/fflux.${year}${month}01.${ayear}${amonth}01.nc ${WORKDIR0}/${MEM}/in/${REGION}.fflux.nc &)

	    #tide dir
	    (ln -s ${TIDEDIR} ${WORKDIR0}/${MEM}/TIDE &)

	    #switch_rst
	    if(${iyr} == ${sdate[1]} && ${imon} == ${sdate[2]})then
		@ switch_rst = 0
	    else
		@ switch_rst = 1
		(ln -s ${WORKDIRB}/${MEM}/out/restart.nc ${WORKDIR0}/${MEM}/in/restart.nc &)
	    endif

	    #exe
	    if(${IMEM} == ${MEMBER})then
		ln -s ${DIR0}/run/${EXE} ${WORKDIR0}/${MEM}/${EXE}
	    else
		(ln -s ${DIR0}/run/${EXE} ${WORKDIR0}/${MEM}/${EXE} &)
	    endif

	    #------------------------------------------------- 
	    cd ${WORKDIR0}/${MEM}

            set nday=1.

	    echo "Make Name List: ${MEM}"#---------
cat <<EOF > pom.nml
&pom_nml
  title = "${REGION}"
  netcdf_file = "${REGION}"
  mode = 3
  assim = 0
  nadv = 2
  nitera = 2
  sw = 1.0
  npg = 2
  dte = 30.0
  isplit = 20
  time_start = "${syear}-${smonth}-01 00:00:00 +00:00"
  nread_rst = ${switch_rst}
  read_rst_file = "restart.nc"
  read_iau_file = "iau.nc"
  write_rst = ${nday}
  write_rst_file = "restart"
  write_ens = 99.
  write_ens_file = "ens"
  budget = ${BUDGET}
  days = ${nday}
  prtd1 = 1.0
  prtd2 = 1.0
  swtch = 9999.
  ts_nudge = 90.
  ti_nudge = 90.
  ss_nudge = 90.
  si_nudge = 90.
/
EOF
	    #----------------------------------------------

	    @ IMEM++

	end

echo "Make vcoordfile"#-------------------------
cd ${DIR0}/run/
csh make_vcoord.csh ${WORKDIR0} ${MEMBER} ${ENODE} ${EPROC}

cd ${WORKDIR0}
@ NODE = ${ENODE} * ${MEMBER}
@ PROC = ${EPROC} * ${MEMBER}
echo "Total Node:"${NODE}
echo "Total Processor:"${PROC}	
echo "Submit Job"#------------------------------
jxsub<<EOF 
#JX --bizcode R2402
#JX --usecode SP
#JX -L rscunit=MA
#JX -L node=${NODE}
#JX --mpi proc=${PROC}
#JX -L node-mem=29184Mi
#JX -L elapse=00:10:00
#JX -N sbPOM_${REGION}_${year}${month}
#JX -j
#JX -s

csh ${DIR0}/run/submit_job.csh ${EPROC} ${MEMBER} ${REGION} ${EXE} ${WORKDIR0}

EOF
#------------------------------------------------

        echo "Wait for Ensemble POM end"
        @ imin = 1
        while(${imin} > 0)
            set FIN_NUM=`find ${WORKDIR0} -name FINISHED | wc -l`
            echo "sbPOM FINISHED [ ${FIN_NUM} / ${MEMBER} ], ${imin} min passed"
            if(${FIN_NUM} == ${MEMBER})then
                @ imin = -999
            else
                sleep 60
                @ imin++
            endif
        end

	cd ${CURDIR}

	#Calculate Mean/Sprd
	echo ${WORKDIR0} > ens_info.txt
	echo ${REGION} >> ens_info.txt
	echo ${syear} ${smonth} 01 >> ens_info.txt
	echo ${year} ${month} ${nday} >> ens_info.txt
	echo ${MEMBER} >> ens_info.txt
	echo ${BUDGET} >> ens_info.txt 

	set option="-assume byterecl -convert big_endian -mcmodel=medium -shared-intel -I/opt/JX/oss/netcdf-fortran-4.2/x86_64/icc/include -L/opt/JX/oss/netcdf-fortran-4.2/x86_64/icc/lib -lnetcdff -I/opt/JX/oss/netcdf-4.2.1.1/x86_64/icc/include -L/opt/JX/oss/netcdf-4.2.1.1/x86_64/icc/lib -lnetcdf -L/opt/JX/oss/hdf5-1.8.9/x86_64/icc/lib -lhdf5_hl -lhdf5 -lz"
	ifort ens/mesp_ens.f90 -o mesp_ens_${year}${month}.out ${option}
	
	if(${RM_ENS} == 0)then
	    ./mesp_ens_${year}${month}.out && rm -f mesp_ens_${year}${month}.out &
	else if(${RM_ENS} == 1)then
	    ./mesp_ens_${year}${month}.out && rm -f mesp_ens_${year}${month}.out && csh remove_ens.csh ${WORKDIR0} ${REGION} ${MEMBER} &
	endif

	@ imon++

    end

@ iyr++
@ imon=1

end

#=========================================================================#
echo "#####ALL END######"
#=========================================================================#
