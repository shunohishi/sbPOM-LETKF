#!/bin/csh

@ i =  426075
@ n =  426175

while($i <= $n)

if($i != 426108)then
jxdel $i &
endif

@ i++

end
