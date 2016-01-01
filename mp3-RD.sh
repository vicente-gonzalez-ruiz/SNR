#!/bin/bash

sound_file=$1
ffmpeg -i $sound_file.wav -f s16le -y $sound_file.pcm 2> /dev/null

br=32
while [ $br -le 320 ]; do
  echo $br
  rm -f $sound_file.mp3
  lame -b $br $sound_file.wav -o $sound_file.mp3 2> /dev/null
  lame --decode $sound_file.mp3 $0.wav 2> /dev/null > /dev/null
  ffmpeg -i $0.wav -f s16le -y $0.pcm 2> /dev/null
  RMSE=`snr --block_size=4096 --FFT --file_A=$sound_file.pcm \
      --file_B=$0.pcm 2> /dev/null | grep RMSE | cut -d "=" -f 2`
  bitrate=`wc -c < $sound_file.mp3`
  echo $bitrate $RMSE >> $0.dat
  let br=br+32
done
