objects= test.o simple_jepg.o

main: $(objects)
	g++ -o main $(objects)

test.o :simple_jepg.h
simple_jepg.o: simple_jepg.h JPEGSamples.h

.PHONY : clean

clean:
	rm main $(objects)
