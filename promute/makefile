.PHONY: all compile userBuild clean

all: compile userBuild

compile:
	make -C src

userBuild:
	-rm -rf build
	mkdir build
	cd build && mkdir external
	cp -r external/em build/external
	cp -r external/scwrl build/external
	cp -r external/surfaceRacer build/external
	cp src/proMute build
	cp src/proMuteBatch build
	cp README build

clean:
	-rm -rf build
	make clean -C src
