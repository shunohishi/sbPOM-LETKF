#!/bin/csh
set sdate=(2016 1)
set edate=(2016 1)

@ iyr = $sdate[1]
@ imon = $sdate[2]

while($iyr <= $edate[1])

    if($iyr == $edate[1])then
	@ emon = $edate[2]
    else
	@ emon = 12
    endif

    while($imon <= $emon)

    if($iyr == 2015 && $imon == 7)then
	csh run.jss2-anal.csh $iyr $imon 7 &
	csh run.jss2-anal.csh $iyr $imon 16 &
    else
	csh run.jss2-anal.csh $iyr $imon 1 &
	csh run.jss2-anal.csh $iyr $imon 16 &
    endif

    @ imon++

    end

    @ imon = 1
    @ iyr++

end
