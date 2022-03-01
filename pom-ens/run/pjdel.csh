#!/bin/csh

@ i = 1707508
@ n = 1707607 

while($i <= $n)

pjdel $i &

@ i++

end
