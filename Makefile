CC	= g++
FLAGS	= -pipe -O3 -static -D _FFT_
LIBS	= -lm -lfftw3

$(HOME)/bin/% :: %.c
	$(CC) $(FLAGS) $< -o $@ $(LIBS)

default:	$(HOME)/bin/snr	$(HOME)/bin/snr2D $(HOME)/bin/snr3D

clean:
	rm -rf $(HOME)/bin/snr $(HOME)/bin/snr2D $(HOME)/bin/snr3D

