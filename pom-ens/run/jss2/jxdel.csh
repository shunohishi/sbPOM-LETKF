#!/bin/csh

@ i = 53429975
@ n = 53430016

while($i <= $n)

jxdel $i &

@ i++

end
