CC	= g++
FLAGS	= -pipe -O3 -static -D _FFT_
LIBS	= -lm -lfftw3

$(HOME)/bin/% :: %.c
	$(CC) $(FLAGS) $< -o $@ $(LIBS)

default:	$(HOME)/bin/snr

clean:
	rm -rf $(HOME)/bin/snr

