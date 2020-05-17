#include "atom.h"

Atom::Atom() {}
Atom::Atom(string line) {
	initializeAtom(line);
}
Atom::Atom(const Atom &other) {
	recordName = other.recordName;
	serialNum = other.serialNum;
	name = other.name;
	altLocation = other.altLocation;
	resName = other.resName;
	chain = other.chain;
	resSequence = other.resSequence;
	insertionCode = other.insertionCode;
	xcoord = other.xcoord;
	ycoord = other.ycoord;
	zcoord = other.zcoord;
	occupancy = other.occupancy;
	temperature = other.temperature;
	element = other.element;
	charge = other.charge;
}

char Atom::aminoCharacter() {
	if (resName == "ALA") return 'A';
	else if (resName == "CYS") return 'C';
	else if (resName == "ASP") return 'D';
	else if (resName == "GLU") return 'E';
	else if (resName == "PHE") return 'F';
	else if (resName == "GLY") return 'G';
	else if (resName == "HIS") return 'H';
	else if (resName == "ILE") return 'I';
	else if (resName == "LYS") return 'K';
	else if (resName == "LEU") return 'L';
	else if (resName == "MET") return 'M';
	else if (resName == "ASN") return 'N';
	else if (resName == "PRO") return 'P';
	else if (resName == "GLN") return 'Q';
	else if (resName == "ARG") return 'R';
	else if (resName == "SER") return 'S';
	else if (resName == "THR") return 'T';
	else if (resName == "VAL") return 'V';
	else if (resName == "TRP") return 'W';
	else if (resName == "TYR") return 'Y';
	else if (resName == "ASX") return 'B';
	else if (resName == "XLE") return 'J';
	else if (resName == "XAA") return 'X';
	else if (resName == "GLX") return 'Z';
	else return '\0';
}

void Atom::initializeAtom(string line) {
	int length = line.length();

	string local = line;
	if (length < 80) {
		local.append(80 - length, ' ');
	}

	recordName		= local.substr(0, 6);
	serialNum		= local.substr(6, 5);
	name			= local.substr(12, 4);
	altLocation		= local.substr(16, 1);
	resName			= local.substr(17, 3);
	chain			= local.substr(21, 1);
	resSequence		= local.substr(22, 4);
	insertionCode 	= local.substr(26, 1);
	xcoord			= local.substr(30, 8);
	ycoord			= local.substr(38, 8);
	zcoord			= local.substr(46, 8);
	occupancy		= local.substr(54, 6);
	temperature		= local.substr(60, 6);
	element			= local.substr(76, 2);
	charge			= local.substr(78, 2);
}

string Atom::getAtomName() {
	string trimName = name;
	for (int i = 3; i >= 0; i--) {
		if (trimName[i] == ' ') {
			trimName.erase(i,1);
		}
	}
	return trimName;
}

char Atom::getChain() {
	return chain[0];
}

char Atom::getAminoChar() {
	return aminoCharacter();
}

int Atom::getResidueNumber() {
	return atoi(resSequence.c_str());
}

int Atom::getSerialNumber() {
	return atoi(serialNum.c_str());
}

void Atom::setAminoChar(char newAminoChar) {
	resName = Atom::getAminoString(newAminoChar);
}

void Atom::setSerialNumber(int newSerialNum) {
	if (newSerialNum > 0 && newSerialNum < 1000000) {
		char buffer[7];
		snprintf(buffer, 7, "%5d", newSerialNum);
		serialNum = string(buffer);
	}
}

void Atom::updateAtomInformation(char newAminoChar, string newName, string newElement) {
	resName = Atom::getAminoString(newAminoChar);

	switch (newName.length()) {
		case 0:
			name = "    ";
			break;
		case 1:
			name = " " + newName + "  ";
			break;
		case 2:
			name = " " + newName + " ";
			break;
		case 3:
			name = " " + newName;
			break;
		case 4:
			name = newName;
			break;
		default:
			name = newName.substr(0,4);
			break;
	}

	switch (newElement.length()) {
		case 0:
			element = "  ";
			break;
		case 1:
			element = " " + newElement;
			break;
		case 2:
			element = newElement;
			break;
		default:
			element = newElement.substr(0,2);
			break;
	}
}

string Atom::ToString() {
	string outputString = recordName + serialNum + " " + name + altLocation + resName;
	outputString += " " + chain + resSequence + insertionCode + "   " + xcoord;
	outputString += ycoord + zcoord + occupancy + temperature + "          ";
	outputString += element + charge;
	return outputString;
}

string Atom::getAminoString(char aminoChar) {
	switch (aminoChar) {
		case 'A': return "ALA";
		case 'C': return "CYS";
		case 'D': return "ASP";
		case 'E': return "GLU";
		case 'F': return "PHE";
		case 'G': return "GLY";
		case 'H': return "HIS";
		case 'I': return "ILE";
		case 'K': return "LYS";
		case 'L': return "LEU";
		case 'M': return "MET";
		case 'N': return "ASN";
		case 'P': return "PRO";
		case 'Q': return "GLN";
		case 'R': return "ARG";
		case 'S': return "SER";
		case 'T': return "THR";
		case 'V': return "VAL";
		case 'W': return "TRP";
		case 'Y': return "TYR";
		case 'B': return "ASX";
		case 'J': return "XLE";
		case 'X': return "XAA";
		case 'Z': return "GLX";
		default:
			return "   ";
	}
}

int Atom::getAtomCount(char aminoChar) {
	switch (aminoChar) {
		case 'A': return 13;
		case 'C': return 14;
		case 'D': return 16;
		case 'E': return 19;
		case 'F': return 23;
		case 'G': return 10;
		case 'H': return 20;
		case 'I': return 22;
		case 'K': return 24;
		case 'L': return 22;
		case 'M': return 20;
		case 'N': return 17;
		case 'P': return 17;
		case 'Q': return 20;
		case 'R': return 26;
		case 'S': return 14;
		case 'T': return 17;
		case 'V': return 19;
		case 'W': return 27;
		case 'Y': return 24;
		case 'B': return 17;  //16?
		case 'J': return 22;
		case 'Z': return 20;  //19?
		default:
			return -1; //no definition...
	}
}
