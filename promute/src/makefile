objects = atom.o chain.o mutator.o protein.o utility.o

.PHONY: all clean

all: proMute proMuteBatch

proMute: proMute.o $(objects)
	g++ proMute.o $(objects) -o proMute

proMuteBatch: proMuteBatch.o $(objects)
	g++ proMuteBatch.o $(objects) -o proMuteBatch

%.o: %.cpp
	g++ -c $<

clean:
	-rm -f $(objects) proMute.o proMuteBatch.o proMute proMuteBatch
