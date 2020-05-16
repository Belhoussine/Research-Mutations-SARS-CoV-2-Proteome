#include "utility.h"

// Returns a lowercase version of the character passed to the function
char toLower(char character) {
	if (character >= 65 && character <= 90) {
		return character + 32;
	} else {
		return character;
	}
}

// Returns an uppercase version of the character passed to the function
char toUpper(char character) {
	if (character >= 97 && character <= 122) {
		return character - 32;
	} else {
		return character;
	}
}

// Returns an uppercase version of the string passed to the function
string toUpper(string str) {
	int length = str.length();
	for (int i = 0; i < length; i++) {
		str[i] = toUpper(str[i]);
	}
	return str;
}

bool openFile(string filename, ofstream &file) {
	file.open(filename.c_str());
	if (!file.is_open()) {
		cout << "Problem opening file for writing: " << filename << endl;
		return false;
	}
	return true;
}
bool openFile(string filename, ifstream &file) {
	file.open(filename.c_str());
	if (!file.is_open()) {
		cout << "Problem opening input file: " << filename << endl;
		return false;
	}
	return true;
}
bool openPDBFile(string filename, ifstream &pdb) {
	pdb.open(filename.c_str());
	if (!pdb.is_open() && !downloadPDBFile(filename, pdb)) {
		cout << "Problem opening Original PDB file: " << filename << endl;
		return false;
	}
	return true;
}

// Downloads and opens the given PDB file from PDB_URL. Returns 0 on success, -1 on failure
bool downloadPDBFile(string PDBFilename, ifstream &pdb) {
	string command = "curl -o \"" + PDBFilename + "\" \"" + PDB_URL + PDBFilename + "\"";
	system(command.c_str());

	pdb.open(PDBFilename.c_str());
	if (!pdb.is_open()) {
		return false;
	}
	return true;
}

void cleanup(ofstream &file) {
	if (file.is_open()) {
		file.close();
	}
}

void cleanup(ifstream &file) {
	if (file.is_open()) {
		file.close();
	}
}
