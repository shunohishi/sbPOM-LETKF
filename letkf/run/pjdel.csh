#!/bin/csh

@ i = 4424392
@ n = 4424497

while($i <= $n)

pjdel $i &

@ i++

end
