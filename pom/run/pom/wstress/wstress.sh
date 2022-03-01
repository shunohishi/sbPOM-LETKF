#! /bin/sh


TXT=w5.txt
TXT2=w.txt
EPS=w15.eps
PNG=w15.png

gnuplot <<EOF
set term postscript eps color enhanced 32
set output '$EPS'
set xrange [0:20]
set yrange [0:2.5]
set grid
set nokey
set xtics 0,2
set xlabel 'U(m/s)'
set ytics 0,0.5
set ylabel 'CD_{10} x 10^{-3}'
plot \
'$TXT'  with lines lt 1 lw 20 lc rgb 'blue',\
'$TXT2' with lines lt 1 lw 20 lc rgb 'black'
EOF

convert $EPS $PNG
echo $PNG

