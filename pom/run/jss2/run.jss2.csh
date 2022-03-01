#=======================================================================
#LONG INTEGRATION VERSION
#=======================================================================

set REGION=test6
set EXE=pom.exe
set DIR0=/data/R/R2402/ohishi/TEST6/pom
set TIDEDIR=${DIR0}/TIDE
set NCDIR=/data/R/R2402/ohishi/TEST6/prep/in
set CURDIR=`pwd`
set PROC=32

set sdate=(2011 1) #start time
set idate=(2011 1) #initial time
set edate=(2015 12) #end time

@ iyr = ${idate[1]}
@ imon = ${idate[2]}

set syear=`printf "%04d" ${sdate[1]}`
set smonth=`printf "%02d" ${sdate[2]}` 

#=========================================================================

echo "Compile"
#make clean
make
if(! -f ${EXE})exit

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
    mkdir -p ${WORKDIR0}
    mkdir -p ${WORKDIR0}/in
    mkdir -p ${WORKDIR0}/out

    echo "Link Netcdf file"
    #rm -f ${WORKDIR0}/in/*.nc
    #---grid.nc
    ln -s ${NCDIR}/grid.nc ${WORKDIR0}/in/${REGION}.grid.nc

    #---tsclim.nc
    ln -s ${NCDIR}/tsclim.nc ${WORKDIR0}/in/${REGION}.tsclim.nc
    #---ic.nc
    ln -s ${NCDIR}/ic.woa13.${smonth}.nc ${WORKDIR0}/in/${REGION}.ic.nc
    #ln -s ${NCDIR}/ic.soda.${smonth}.nc ${WORKDIR0}/in/${REGION}.ic.nc
    #---tsdata.nc
    ln -s ${NCDIR}/tsdata_mclim.nc ${WORKDIR0}/in/${REGION}.tsdata.nc
    #ln -s ${NCDIR}/tsdata.${year}${month}01.${ayear}${amonth}01.nc ${WORKDIR0}/in/${REGION}.tsdata.nc
    #lbc.nc
    ln -s ${NCDIR}/lbc_mclim.nc ${WORKDIR0}/in/${REGION}.lbc.nc
    #ln -s ${NCDIR}/lbc.${year}${month}01.${ayear}${amonth}01.nc ${WORKDIR0}/in/${REGION}.lbc.nc
    #atm.nc
    ln -s ${NCDIR}/jra55.${year}${month}01.${ayear}${amonth}01.nc ${WORKDIR0}/in/${REGION}.atm.nc
    #fflux.nc
    ln -s ${NCDIR}/fflux.${year}${month}01.${ayear}${amonth}01.nc ${WORKDIR0}/in/${REGION}.fflux.nc

    ln -s ${TIDEDIR} ${WORKDIR0}/TIDE

    #switch_rst
    if(${iyr} == ${sdate[1]} && ${imon} == ${sdate[2]})then
    @ switch_rst = 0
    else
    @ switch_rst = 1
    ln -s ${WORKDIRB}/out/restart.nc ${WORKDIR0}/in/restart.nc
    endif

    echo "Make NML"
    cp ${EXE} ${WORKDIR0}/
    cd ${WORKDIR0}
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
  budget = 0
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

    #scs 25 min./month
    #wpac 70 min./month, con = 1 15min./2day
    echo "Submit Job"
jxsub <<EOF
#JX --bizcode R2402
#JX --usecode SP
#JX -L rscunit=MA
#JX -L node=1
#JX --mpi proc=${PROC}
#JX -L node-mem=29184Mi
#JX -L elapse=20:00
#JX -N sbPOM_${REGION}_${year}${month}
#JX -S

mpiexec -np ${PROC} -stdout stdout.${REGION} -stderr stderr.${REGION} ./${EXE}

echo finished > FINISHED

EOF

    echo "Wait for Job finish"
    @ imin = 1
    while(${imin} > 0)
	if(-f FINISHED)then
	sleep 1m
	echo "End ${year}${month}"
	@ imin = -999
	else
	sleep 1m
	echo "${imin} minitues have passed" 
	@ imin++
	endif
    end

    cd ${CURDIR}

    @ imon++

    end

@ iyr++
@ imon=1

end

#=========================================================================#
echo "#####ALL END######"
#=========================================================================#
