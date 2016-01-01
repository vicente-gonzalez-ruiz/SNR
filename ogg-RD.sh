#/bin/bash

sound_file=$1
ffmpeg -i $sound_file.wav -f s16le -y $sound_file.pcm

quality=-1
while [ $quality -le 10 ]; do
  echo $quality
  rm -f $sound_file.ogg
  oggenc -q $quality $sound_file.wav -o $sound_file.ogg 2> /dev/null
  oggdec $sound_file.ogg -o $0.wav 2> /dev/null > /dev/null
  ffmpeg -i $0.wav -f s16le -y $0.pcm 2> /dev/null
  RMSE=`snr --block_size=4096 --FFT --file_A=$sound_file.pcm \
--file_B=$0.pcm 2> /dev/null | grep RMSE | cut -d "=" -f 2`
  bitrate=`wc -c < $sound_file.ogg`
  echo $bitrate $RMSE >> $0.dat
  let quality=quality+1
done
