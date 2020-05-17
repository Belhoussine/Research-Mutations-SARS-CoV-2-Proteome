// proMute.cpp
// Katie Hursh
// Gideon Wolfe
// Last Modified: 17 May 2020

#include "mutator.h"

int main(int argc, char* argv[]) {
	if (argc < 5 || argc > 7) {
		cout << "Invalid number of arguments" << endl;
		cout << "Usage: " << argv[0] << " protein chain residueNum target [em/no]" << endl;
		return 1;
	} else if (strlen(argv[2]) > 1 || strlen(argv[4]) > 1) {
		cout << "Chain and Target should be 1 character" << endl;
		cout << "Usage: " << argv[0] << " protein chain residueNum target [em/no]" << endl;
		return 1;
	} else if (argc == 6 && (strlen(argv[5]) > 4 || (strcmp(argv[5],"no") && strcmp(argv[5],"em") && strcmp(argv[5],"srem")))) {
		cout << "Final argument should be 'em', 'srem' or 'no'" << endl;
		cout << "Usage: " << argv[0] << " protein chain residueNum target [em/no]" << endl;
		return 1;
	}

	// holders for the mutation information
	char chain = toUpper(argv[2][0]);
	char mutant = toUpper(argv[4][0]);
	int residue = atoi(argv[3]);
	char performEM = 'm';
	char performSR = 'm';

	// String holders for the files names and commands
	string proteinName = toUpper(string(argv[1]));
	string PDBFilename = proteinName + ".pdb";
	string outputBase = proteinName + "." + chain + string(argv[3]) + mutant;

	ifstream inputPDB;
	// Open the input and output files.  Quit program on failure
	if (!openPDBFile(PDBFilename, inputPDB)) {
		return -1;
	}

	cout << "Source File: " + PDBFilename << endl << endl;

	Protein inputProtein(inputPDB);
	cleanup(inputPDB);

	cout << "PDB File read, proceeding with mutation." << endl;

	if (argc >= 6 && strcmp(argv[5],"no") == 0) {
		performEM = 'n';
	} else if (argc == 6 && strcmp(argv[5],"em") == 0) {
		performEM = 'y';
	} else if (argc == 6 && strcmp(argv[5],"srem") == 0) {
		performEM = 'y';
		performSR = 'y';
	}

    return Mutator::performMutation(inputProtein, chain, residue, mutant, PDBFilename, outputBase, performEM, performSR);
}
