#ifndef MUTATOR_H
#define MUTATOR_H

#include <unistd.h>
#include "protein.h"
#include "utility.h"
using namespace std;

#define SURF_RACER_FOLDER "./external/surfaceRacer/"
#define SURF_RACER_EXE "surfrace5_0_linux_64bit"

class Mutator {
	private:
		Mutator();

		constexpr static const double SURFACE_AREA_LIMIT = 70.0;

		static bool canMutateLocally(char originResidue, char mutant);
		static bool toRemove(Atom originAtom, char mutant);
		static bool toRemoveSerine(Atom originAtom);
		static bool largeToSmallMutation(char originResidue, char mutant);
		static bool shouldPerformEM(char originResidue, char mutant, int chainNum, int residue, string PDBFilename, bool alreadyRun, char performSR);
		static int sendToScwrl(Protein inputProtein, char chain, int residue, char mutant, string PDBFilename, string outputBase);
		static int mutateLocally(Protein inputProtein, char chain, int residue, char mutant, string outputBase);
		static float getSurfaceArea(int chainNum, int residue, string PDBFilename, bool alreadyRun);
		static Atom mutateAtom(Atom origAtom, char mutant);
		static Atom mutateAtomToSerine(Atom origAtom);
		static int performEnergyMinimization(string outputBase);

	public:

		static int performMutation(Protein inputProtein, char chain, int residue, char mutant, string PDBFilename, string outputBase, char performEM, char performSR);
};

#endif
