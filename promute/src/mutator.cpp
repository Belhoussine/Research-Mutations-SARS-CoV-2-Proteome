// mutator.cpp
// Katie Hursh
// Gideon Wolfe
// Last Modified: 17 May 2020
#include "mutator.h"

Mutator::Mutator() {}

bool Mutator::canMutateLocally(char originResidue, char mutant) {
  // FBJ: code being deprecated
  // local mutation no longer an option
  // forcing use of SCWRL in all instances for consistency
  /*if (originResidue == 'B' || originResidue == 'J' || originResidue == 'X' || originResidue == 'Z') {
    return false;
    } else if (mutant == 'G') {
    return true;
    } else if (mutant == 'A' && originResidue != 'G') {
    return true;
    } else if (mutant == 'S' && originResidue != 'G' && originResidue != 'A') {
    return true;
    } else {
    return false;
    }
  */
  return false;
}

bool Mutator::toRemove(Atom originAtom, char mutant) {
	string atomType = originAtom.getAtomName();

	// Any target other than G, A, and S are not handled
	if (mutant != 'G' && mutant != 'A' && mutant != 'S') {
		return false;
	}

	if (atomType == "H") {
		return true;
	} else if (atomType == "N" || atomType == "CA" || atomType == "C" || atomType == "O") {
		return false;
	} else if (atomType == "CB" && (mutant == 'A' || mutant == 'S')) {
		return false;
	}

	if (mutant == 'S') {
		return toRemoveSerine(originAtom);
	}

	return true;
}

bool Mutator::toRemoveSerine(Atom originAtom) {
	string atomType = originAtom.getAtomName();
	char aminoChar = originAtom.getAminoChar();

	switch (aminoChar) {
		case 'D':
			if (atomType == "OD2" || atomType == "OD1")
				return true;
			break;
		case 'E':
			if (atomType == "CD" || atomType == "OE2" || atomType == "OE1")
				return true;
			break;
		case 'F':
			if (atomType == "CD1" || atomType == "CE1" || atomType == "CZ" || atomType == "CE2"
				|| atomType == "CD2")
				return true;
			break;
		case 'H':
			if (atomType == "CD1" || atomType == "NE2" || atomType == "CE2" || atomType == "ND1"
				|| atomType == "CD2" || atomType == "CE1")
				return true;
			break;
		case 'I':
			if (atomType == "CD1" || atomType == "CG1")
				return true;
			break;
		case 'K':
			if (atomType == "CD" || atomType == "CE" || atomType == "NZ")
				return true;
			break;
		case 'L':
			if (atomType == "CD1" || atomType == "CD2")
				return true;
			break;
		case 'M':
			if (atomType == "SD" || atomType == "CE")
				return true;
			break;
		case 'N':
			if (atomType == "ND2" || atomType == "OD1")
				return true;
			break;
		case 'P':
			if (atomType == "CD")
				return true;
			break;
		case 'Q':
			if (atomType == "CD" || atomType == "OE1" || atomType == "NE2")
				return true;
			break;
		case 'R':
			if (atomType == "CD2" || atomType == "CD" || atomType == "NE" || atomType == "CZ"
				|| atomType == "NH2" || atomType == "NH1")
				return true;
			break;
		case 'T':
			if (atomType == "OG1")
				return true;
			break;
		case 'V':
			if (atomType == "CG1")
				return true;
			break;
		case 'W':
			if (atomType == "CD1" || atomType == "NE1" || atomType == "CE2" || atomType == "CZ2"
				|| atomType == "CH2" || atomType == "CZ3" || atomType == "CE3" || atomType == "CD2")
				return true;
			break;
		case 'Y':
			if (atomType == "CD2" || atomType == "CD1" || atomType == "CE2" || atomType == "CE1"
				|| atomType == "CZ" || atomType == "OH")
				return true;
			break;
	}

	// All other atoms not specifically mentioned
	return false;
}

bool Mutator::largeToSmallMutation(char originResidue, char mutant) {
	int originCount = Atom::getAtomCount(originResidue);
	int mutantCount = Atom::getAtomCount(mutant);

	if (originCount == -1 || mutantCount == -1) {
		return false;
	}
	return originCount > mutantCount;
}

bool Mutator::shouldPerformEM(char originResidue, char mutant, int chainNum, int residue, string PDBFilename, bool alreadyRun, char performSR) {
  cout << "Checking accessable surface area" << endl;
  if (getSurfaceArea(chainNum, residue, PDBFilename, alreadyRun) > SURFACE_AREA_LIMIT) {
    return false;
  } else {
    cout << "Running large to small mutation check" << endl;
    return !largeToSmallMutation(originResidue, mutant);
  }
}

int Mutator::sendToScwrl(Protein inputProtein, char chain, int residue, char mutant, string PDBFilename, string outputBase) {
	string FastaFilename = outputBase + ".fasta.txt";
	string FinalPDBFilename = outputBase + ".pdb";
	string scwrlCommand;

	ofstream outputFASTA;
	// Open the input and output files.  Quit program on failure
	if (!openFile(FastaFilename, outputFASTA)) {
		return -1;
	}

	inputProtein.createFASTAFile(outputFASTA, chain, residue, mutant);
	cleanup(outputFASTA);

	if (inputProtein.isMultiChain()) {
		string toScwrlFilename = outputBase + "_MID1.pdb";
		string FromScwrlFilename = outputBase + "_MID2.pdb";

		ofstream toScwrl;
		ifstream fromScwrl;
		ofstream finalPDB;

		if (!openFile(toScwrlFilename, toScwrl) || !openFile(FinalPDBFilename, finalPDB)) {
			cleanup(toScwrl);
			cleanup(finalPDB);
			return -1;
		}

		inputProtein.createPDBFile_scwrl(toScwrl, chain);
		scwrlCommand = "./external/scwrl/Scwrl4 -i " + toScwrlFilename + " -o " + FromScwrlFilename + " -s " + FastaFilename;
		system(scwrlCommand.c_str());

		if (!openFile(FromScwrlFilename, fromScwrl)) {
			return -1;
		}

		Protein scwrlChain(fromScwrl);
		inputProtein.createPDBFile(scwrlChain, finalPDB, chain);

		cleanup(toScwrl);
		cleanup(fromScwrl);
		cleanup(finalPDB);

		return 0;
	} else {
		scwrlCommand = "./external/scwrl/Scwrl4 -i " + PDBFilename + " -o " + FinalPDBFilename + " -s " + FastaFilename;
		return system(scwrlCommand.c_str());
	}
}

int Mutator::mutateLocally(Protein inputProtein, char chain, int residue, char mutant, string outputBase) {
	Chain origChain = inputProtein.getChainWithResidue(chain, residue);
	if (origChain.isChain(' ')) {
		cout << "Error getting chain." << endl;
		return -1;
	}

	vector<Atom> atoms = origChain.getAtoms();
	int chainLength = atoms.size();

	Chain mutantChain;

	for (int i = 0; i < chainLength; i++) {
		if (atoms[i].getResidueNumber() != residue) {
			mutantChain.addAtom(atoms[i]);
		} else if (!toRemove(atoms[i], mutant)) {
			mutantChain.addAtom(mutateAtom(atoms[i], mutant));
		}
	}

	inputProtein.setChainWithResidue(chain, residue, mutantChain);
	inputProtein.updateSerials(1);

	ofstream finalPDB;
	// Open the input and output files.  Quit program on failure
	if (!openFile(outputBase + ".pdb", finalPDB)) {
		return -1;
	}

	finalPDB << inputProtein.ToString();
	cleanup(finalPDB);
	return 0;
}

float Mutator::getSurfaceArea(int chainNum, int residue, string PDBFilename, bool alreadyRun) {
	string proteinOnly = PDBFilename.substr(0, 4);
		
	if (!alreadyRun) {
		chdir(SURF_RACER_FOLDER);
		
		string pipeFilepath = proteinOnly + "_command.txt";
		string command = string("./") + SURF_RACER_EXE + " < " + pipeFilepath;
		
		ofstream pipeFile;
		if (!openFile(pipeFilepath, pipeFile)) {
			chdir("../..");
			return -2;
		}
		
		pipeFile << 1 << endl;
		pipeFile << proteinOnly << ".pdb" << endl;
		pipeFile << 2 << endl;
		pipeFile << 1 << endl;
		pipeFile << endl;
		pipeFile << endl;
		
		cout << "Running Surface Racer" << endl; 
		
		cleanup(pipeFile);
		system(command.c_str());
		
		chdir("../..");
	}
	
	cout << endl << endl << "Retrieving surface area" << endl;
	string outputFilepath = SURF_RACER_FOLDER + proteinOnly + "_residue.txt";
	
	ifstream surfaceAreaFile;
	if (!openFile(outputFilepath, surfaceAreaFile)) {
		return -1;
	}
	
	int currentChain = 0;
	int oldResidueNum = -1;
	int currentResidueNum;
	string line;
	while (getline(surfaceAreaFile, line)) {
		currentResidueNum = atoi(line.substr(0,3).c_str());
		
		if (currentResidueNum < oldResidueNum) {
			currentChain += 1;
		}
		
		if (currentChain == chainNum && currentResidueNum == residue) {
			return atof(line.substr(8,8).c_str());
		}
		oldResidueNum = currentResidueNum;
	}
	cleanup(surfaceAreaFile);
	
	return -1;
}

Atom Mutator::mutateAtom(Atom originAtom, char mutant) {
	string atomType = originAtom.getAtomName();

	if (atomType == "N" || atomType == "CA" || atomType == "C" || atomType == "O") {
		originAtom.setAminoChar(mutant);
	} else if (atomType == "CB" && (mutant == 'A' || mutant == 'S')) {
		originAtom.setAminoChar(mutant);
	} else if (mutant == 'S') {
		return mutateAtomToSerine(originAtom);
	}

	return originAtom;
}

Atom Mutator::mutateAtomToSerine(Atom originAtom) {
	string atomType = originAtom.getAtomName();
	char aminoChar = originAtom.getAminoChar();

	switch (aminoChar) {
		case 'C':
			if (atomType == "SG")
				originAtom.updateAtomInformation('S', "OG", "O");
			break;
		case 'D':
		case 'E':
		case 'F':
		case 'H':
		case 'K':
		case 'L':
		case 'M':
		case 'N':
		case 'P':
		case 'Q':
		case 'R':
		case 'W':
		case 'Y':
			if (atomType == "CG")
				originAtom.updateAtomInformation('S', "OG", "O");
			break;
		case 'I':
		case 'T':
		case 'V':
			if (atomType == "CG2")
				originAtom.updateAtomInformation('S', "OG", "O");
			break;
	}

	return originAtom;
}

int Mutator::performEnergyMinimization(string outputBase) {
	string emCommand = "cd external/em/ && ./runEnergyMinimization.sh " + outputBase;
	int result = system(emCommand.c_str());

	if (result == 0) {
		cout << endl << endl << "Energy minimization complete in file: " << outputBase << "_em.pdb" << endl;
	} else {
		cout << endl << endl << "Energy minimization failed." << endl;
	}
	return result;
}

int Mutator::performMutation(Protein inputProtein, char chain, int residue, char mutant, string PDBFilename, string outputBase, char performEM, char performSR) {
	char originResidue = inputProtein.getResidueType(chain, residue);
	bool surfaceAlreadyRun = false;
	int chainNum = -1;
	
	///////////////////////////////////////////////////////////////////////////
	ofstream debug("debug.txt",ios_base::app);
	if (debug.is_open()) {
		debug << PDBFilename << "\t" << chain << "\t" << residue << "\t" << mutant << "\t";
	}
	///////////////////////////////////////////////////////////////////////////

	if (originResidue == '\0') {
		cout << "Error getting residue to mutate." << endl;
		///////////////////////////////////////////////////////////////////////////
		if (debug.is_open()) {
			debug << "Error: Getting source residue" << endl;
		}
		///////////////////////////////////////////////////////////////////////////
		return -1;
	}
	if (originResidue == mutant) {
		cout << "Original residue matches mutant.  Nothing to perform." << endl;
		///////////////////////////////////////////////////////////////////////////
		if (debug.is_open()) {
			debug << "Error: No mutation needed" << endl;
		}
		///////////////////////////////////////////////////////////////////////////
		return 0;
	}

	if (performEM == 'm') {
		string surfaceName = SURF_RACER_FOLDER + PDBFilename.substr(0, 4) + "_residue.txt";
		ifstream check(surfaceName.c_str());

		if (check.is_open()) {
			cout << "Surface area has already been calculated for the base protein." << endl;
			surfaceAlreadyRun = true;
			cleanup(check);
		} else {
			string surfacePDBName = SURF_RACER_FOLDER + PDBFilename.substr(0, 4) + ".pdb";
			cout << "Saving PDB for energy minimization check." << endl;

			ofstream surfRacerPDB;
			if (!openFile(surfacePDBName, surfRacerPDB)) {
				///////////////////////////////////////////////////////////////////////////
				if (debug.is_open()) {
					debug << "Error: SurfRacer PDB" << endl;
				}
				///////////////////////////////////////////////////////////////////////////
				return -1;
			}

			surfRacerPDB << inputProtein.ToString_ChainsOnly();
			cleanup(surfRacerPDB);
		}
		chainNum = inputProtein.getChainNumber(chain, residue);
	}

	int result;
	if (canMutateLocally(originResidue, mutant)) {
		cout << "Mutating Locally: Chain " << chain << ", Residue " << residue << ", " << Atom::getAminoString(originResidue)
			 << " -> " << Atom::getAminoString(mutant) << endl << endl;
		result = mutateLocally(inputProtein, chain, residue, mutant, outputBase);
	} else {
		cout << "Mutating with Scwrl4: Chain " << chain << ", Residue " << residue << ", " << Atom::getAminoString(originResidue)
			 << " -> " << Atom::getAminoString(mutant) << endl << endl;
		result = sendToScwrl(inputProtein, chain, residue, mutant, PDBFilename, outputBase);
	}

	if (result == -1) {
		cout << "Failed during mutation." << endl;
		///////////////////////////////////////////////////////////////////////////
		if (debug.is_open()) {
			debug << "Error: Mutation error" << endl;
		}
		///////////////////////////////////////////////////////////////////////////
		return -1;
	}

	cout << "Mutation complete in file: " << outputBase << ".pdb" << endl;
	
	
	///////////////////////////////////////////////////////////////////////////
	if (debug.is_open()) {
		// debug << originResidue << "\t" << getSurfaceArea(chainNum, residue, PDBFilename, surfaceAlreadyRun) << "\t";
		
		if (largeToSmallMutation(originResidue, mutant)) {
			debug << "L->S\t";
		} else {
			debug << "Not L->S\t";
		}
		
		if (inputProtein.isMultiChain()) {
			debug << "Multi Chain\t";
		} else {
			debug << "Single Chain\t";
		}
		
		if (canMutateLocally(originResidue, mutant)) {
			debug << "Local\t";
		} else {
			debug << "Scwrl\t";
		}
		
		if (performEM == 'n') {
			debug << outputBase << ".pdb\tForced: Skipped" << endl;
		} else if (performEM == 'y') {
			debug << outputBase << "_em.pdb\tForced: Run" << endl;
		} else if (shouldPerformEM(originResidue, mutant, chainNum, residue, PDBFilename, surfaceAlreadyRun, performSR)) {
			debug << outputBase << "_em.pdb\tRun" << endl;
		} else {
			debug << outputBase << ".pdb\tSkipped" << endl;
		}
	}
	///////////////////////////////////////////////////////////////////////////
	
	if (performEM == 'n') {
		cout << endl << "Energy minimization skipped per arguments." << endl << endl;
		return 0;
	} else if (performEM == 'y' && (performSR == 'n' || performSR == 'm')) { // if EM is to be run without SR filtering
		cout << endl << "Energy minimization being run without SR as per arguments." << endl << endl;
		return performEnergyMinimization(outputBase);
  } else if(performEM == 'y' && performSR == 'y'){ // if EM is to be filtered by SR result
		cout << endl << "Energy minimization being run with SR as per arguments." << endl << endl;
    if (shouldPerformEM(originResidue, mutant, chainNum, residue, PDBFilename, surfaceAlreadyRun, performSR)) {
	    cout << endl << "Conditions met, performing energy minimization." << endl << endl;
	    return performEnergyMinimization(outputBase);
		} else {
			cout << endl << "Conditions have not been met, skipping energy minimization." << endl << endl;
			return 0;
		}
  } else { // default option if no arguments are specified
		cout << endl << "Energy minimization being run." << endl << endl;
		return performEnergyMinimization(outputBase);
	}

	return 0;
}
