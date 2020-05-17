#include "protein.h"

Protein::Protein() {
	init();
}

Protein::Protein(ifstream &pdb) {
	init();
	initializeProtein(pdb);
}

Protein::Protein(const Protein &other) {
	proteinName = other.proteinName;
	changeComment = other.changeComment;
	preRemarks = other.preRemarks;
	postRemarks = other.postRemarks;
	chains = other.chains;
}

void Protein::init() {
	proteinName = "";
	preRemarks = "";
	postRemarks = "";

	changeComment  = "REMARK   File generated with ProMute 2.0\n";
	//changeComment += "REMARK   INSERT URL HERE\n";
	changeComment += "REMARK\n";
	changeComment += "REMARK   Code written and maintained by\n";
	changeComment += "REMARK\n";
	changeComment += "REMARK   Katie Hursh\n";
	changeComment += "REMARK   Filip Jagodzinski\n";
	changeComment += "REMARK   Western Washington University\n";
	changeComment += "REMARK   Last Modified: 22 May 2018\n";
	changeComment += "REMARK\n";
	//changeComment += "REMARK   Hydrogens added by Scwrl have been maintained\n";
}

bool Protein::isMultiChain() {
	return chains.size() > 1;
}

void Protein::initializeProtein(ifstream &pdb) {
	Atom currAtom;
	Chain currChain;
	string line;
	char currentChainID = '\0';
	char atomChain;
	bool firstAtom = true;
	bool pre = true;

	while (getline(pdb, line)) {
		if (line.compare(0,4,"ATOM") == 0 || line.compare(0,6,"HETATM") == 0 || line.compare(0,3,"TER") == 0) {
			pre = false;
			currAtom.initializeAtom(line);

			atomChain = currAtom.getChain();

			if (currentChainID != atomChain) {
				currentChainID = atomChain;

				if (firstAtom) {
					firstAtom = false;
				} else {
					chains.push_back(currChain);
					currChain.clearChain();
				}
			}
			currChain.addAtom(currAtom);
		} else if (line.compare(0,5,"MODEL") == 0) {
			changeComment += "REMARK   Only the first model was kept\n";
		} else if (line.compare(0,3,"END") == 0) {
			if (line.compare(0,6,"ENDMDL") != 0) {
				postRemarks += line;
			}
			chains.push_back(currChain);
			changeComment += "REMARK\n";
			return;
		} else if (pre && line.compare(0,6,"SEQRES") != 0 && line.compare(0,6,"NUMMDL") != 0) {
			preRemarks += line + "\n";
		} else if (line.compare(0,6,"SEQRES") != 0 && line.compare(0,6,"NUMMDL") != 0) {
			postRemarks += line + "\n";
		}
	}
	changeComment += "REMARK\n";
}

void Protein::createFASTAFile(ofstream &fasta, char chain, int residueNum, char mutant) {
	int numChains = chains.size();
	for (int i = 0; i < numChains; i++) {
		if (chains[i].isChain(chain)) {
			fasta << chains[i].fastaString(residueNum, mutant);
		}
	}
}

void Protein::createPDBFile_scwrl(ofstream &pdb, char chain) {
	int numChains = chains.size();

	pdb << preRemarks;
	for (int i = 0; i < numChains; i++) {
		if (chains[i].isChain(chain)) {
			pdb << chains[i].ToString();
		}
	}
	pdb << postRemarks;
}

void Protein::createPDBFile(Protein scwrlChain, ofstream &finalPDB, char chain) {
	int numChains = chains.size();
	int nextSerialNum = 1;
	bool foundMutation = false;

	finalPDB << changeComment << preRemarks;
	for (int i = 0; i < numChains; i++) {
		if (!chains[i].isChain(chain)) {
			nextSerialNum = chains[i].updateSerials(nextSerialNum);
			finalPDB << chains[i].ToString();
		} else if (!foundMutation) {
			foundMutation = true;
			nextSerialNum = scwrlChain.updateSerials(nextSerialNum);
			finalPDB << scwrlChain.ToString_ChainsOnly();
		}
	}
	finalPDB << postRemarks;
}

void Protein::addToChangeComment(string append) {
	int length = append.size();
	
	changeComment += append;
	if (append[length-1] != '\n') {
		changeComment += "\n";
	}
	
}

int Protein::getMaxSerial() {
	int length = chains.size();
	int maxSerial = -1;
	int chainMaxSerial;

	for (int i = 0; i < length; i++) {
		chainMaxSerial = chains[i].getMaxSerial();
		if (maxSerial < chainMaxSerial) {
			maxSerial = chainMaxSerial;
		}
	}
	return maxSerial;
}

int Protein::updateSerials(int minSerialNum) {
	int length = chains.size();
	int counter = minSerialNum;

	for (int i = 0; i < length; i++) {
		counter = chains[i].updateSerials(counter);
	}
	return counter;
}

int Protein::getChainNumber(char chain, int residueNum) {
	int length = chains.size();
	for (int i = 0; i < length; i++) {
		if (chains[i].isChain(chain) && chains[i].containsResidue(residueNum)) {
			return i;
		}
	}
	
	return -1;
}

char Protein::getResidueType(char chain, int residueNum) {
	int length = chains.size();
	char resType = '\0';
	for (int i = 0; i < length && resType == '\0'; i++) {
		if (chains[i].isChain(chain)) {
			resType = chains[i].getResidueType(residueNum);
		}
	}
	return resType;
}

Chain Protein::getChainWithResidue(char chain, int residueNum) {
	int length = chains.size();
	for (int i = 0; i < length; i++) {
		if (chains[i].isChain(chain) && chains[i].containsResidue(residueNum)) {
			return chains[i];
		}
	}
	
	Chain empty;
	return empty;
}

void Protein::setChainWithResidue(char chain, int residueNum, Chain newChain) {
	int length = chains.size();
	for (int i = 0; i < length; i++) {
		if (chains[i].isChain(chain) && chains[i].containsResidue(residueNum)) {
			chains[i] = newChain;
			return;
		}
	}
}

vector<char> Protein::getChainIDs() {
	vector<char> output;
	
	int numChains = chains.size();
	int chainsGrabbed = 0;
	char currentChain;
	bool found;
	
	for (int i = 0; i < numChains; i++) {
		found = false;
		currentChain = chains[i].getChainID();
		
		for (int j = 0; j < chainsGrabbed; j++) {
			if (currentChain == output[j]) {
				found = true;
				break;
			}
		}
		
		if (!found) {
			output.push_back(currentChain);
			chainsGrabbed += 1;
		}
	}
	
	return output;
}

int Protein::getMinResidueNum(char chain) {
	int minResidue = 10000;
	int localMin;
	int length = chains.size();
	
	for (int i = 0; i < length; i++) {
		if (chains[i].isChain(chain)) {
			localMin = chains[i].getMinResidueNum();
			if (localMin < minResidue) {
				minResidue = localMin;
			}
		}
	}
	
	if (minResidue == 10000) {
		return -1;
	} else {
		return minResidue;
	}
}

int Protein::getMaxResidueNum(char chain) {
	int maxResidue = -1;
	int localMax;
	int length = chains.size();
	
	for (int i = 0; i < length; i++) {
		if (chains[i].isChain(chain)) {
			localMax = chains[i].getMaxResidueNum();
			if (localMax > maxResidue) {
				maxResidue = localMax;
			}
		}
	}
	
	return maxResidue;
}

string Protein::ToString() {
	return changeComment + preRemarks + ToString_ChainsOnly() + postRemarks;
}

string Protein::ToString_ChainsOnly() {
	string outputString = "";

	int numChains = chains.size();
	for (int i = 0; i < numChains; i++) {
		outputString += chains[i].ToString();
	}

	return outputString;
}
