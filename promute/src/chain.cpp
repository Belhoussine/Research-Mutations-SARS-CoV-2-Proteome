#include "chain.h"

Chain::Chain() {
	chainID = ' ';
}

Chain::Chain(const Chain &other) {
	chainID = other.chainID;
	atoms = other.atoms;
}

bool Chain::isChain(char chainID) {
	return this->chainID == chainID;
}

bool Chain::containsResidue(int residueNum) {
	int length = atoms.size();
	for (int i = 0; i < length; i++) {
		if (atoms[i].getResidueNumber() == residueNum) {
			return true;
		}
	}
	return false;
}

vector<Atom> Chain::getAtoms() {
	return atoms;
}

int Chain::getNumAtoms() {
	return atoms.size();
}

int Chain::getMaxSerial() {
	int length = atoms.size();
	int maxSerial = -1;
	int atomSerial;

	for (int i = 0; i < length; i++) {
		atomSerial = atoms[i].getSerialNumber();
		if (maxSerial < atomSerial) {
			maxSerial = atomSerial;
		}
	}
	return maxSerial;
}

void Chain::addAtom(string line) {
	if (line.compare(0,4,"ATOM") == 0 || line.compare(0,6,"HETATM") == 0) {
		Atom tempAtom(line);
		char lineChain = toUpper(tempAtom.getChain());

		if (chainID == ' ') {
			chainID = lineChain;
		}
		if (chainID == lineChain) {
			atoms.push_back(tempAtom);
		}
	} else {
		printf("Not an Atom");
	}
}

void Chain::addAtom(Atom newAtom) {
	char lineChain = toUpper(newAtom.getChain());

	if (chainID == ' ') {
		chainID = lineChain;
	}
	if (chainID == lineChain) {
		atoms.push_back(newAtom);
	}
}

void Chain::clearChain() {
	chainID = ' ';
	atoms.clear();
}

int Chain::updateSerials(int minSerialNum) {
	int length = atoms.size();
	int counter = minSerialNum;

	for (int i = 0; i < length; i++) {
		atoms[i].setSerialNumber(counter);
		counter++;
	}

	return counter;
}

char Chain::getResidueType(int residueNum) {
	int length = atoms.size();
	for (int i = 0; i < length; i++) {
		if (atoms[i].getResidueNumber() == residueNum) {
			return atoms[i].getAminoChar();
		}
	}
	return '\0';
}

char Chain::getChainID() {
	return chainID;
}

int Chain::getMinResidueNum() {
	int length = atoms.size();
	int minResidue = atoms[0].getResidueNumber();
	int currentRes;
	
	for (int i = 1; i < length; i++) {
		currentRes = atoms[i].getResidueNumber();
		if (currentRes < minResidue) {
			minResidue = currentRes;
		}
	}
	return minResidue;
}

int Chain::getMaxResidueNum() {
	int length = atoms.size();
	int maxResidue = atoms[0].getResidueNumber();
	int currentRes;
	
	for (int i = 1; i < length; i++) {
		currentRes = atoms[i].getResidueNumber();
		if (currentRes > maxResidue) {
			maxResidue = currentRes;
		}
	}
	return maxResidue;
}

string Chain::fastaString(int residueNum, char destination) {
	int numAtoms = atoms.size();
	string outputString = "";
	char aminoChar;

	int currentResNum = -1;
	int atomResNum;

	for (int i = 0; i < numAtoms; i++) {
		atomResNum = atoms[i].getResidueNumber();
		if (currentResNum != atomResNum) {
			currentResNum = atomResNum;
			aminoChar = atoms[i].getAminoChar();

			if (currentResNum == residueNum) {
				outputString += toUpper(destination);
			} else if (aminoChar != '\0') {
				outputString += toLower(aminoChar);
			}
		}
	}
	return outputString;
}

string Chain::ToString() {
	int numAtoms = atoms.size();
	string outputString = "";

	for (int i = 0; i < numAtoms; i++) {
		outputString += atoms[i].ToString() + "\n";
	}
	return outputString;
}
