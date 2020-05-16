objects = atom.o chain.o mutator.o protein.o utility.o

.PHONY: all clean

all: proMute makeScript

proMute: proMute.o $(objects)
	g++ proMute.o $(objects) -o proMute

makeScript: makeScript.o $(objects)
	g++ makeScript.o $(objects) -o makeScript

%.o: %.cpp
	g++ -c $<

clean:
	-rm -f $(objects) proMute.o makeScript.o proMute makeScript
