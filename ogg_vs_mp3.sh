#!/bin/bash

gnuplot -p << EOF
set title "Ogg Vorbis vs MP3 LAME"
set xlabel "File size in bytes"
set ylabel "RMSE (Root Mean Square Error)"
plot "./ogg-RD.sh.dat" with lines title "Ogg Vorbis", \
     "./mp3-RD.sh.dat" with lines title "MP3 LAME"
EOF
