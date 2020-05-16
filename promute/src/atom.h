#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <stdlib.h>
#include <cstdio>
using namespace std;

class Atom {
	private:
		string recordName;
		string serialNum;
		string name;
		string altLocation;
		string resName;
		string chain;
		string resSequence;
		string insertionCode;
		string xcoord;
		string ycoord;
		string zcoord;
		string occupancy;
		string temperature;
		string element;
		string charge;

		char aminoCharacter();

	public:
		Atom();
		Atom(string line);
		Atom(const Atom &other);

		void initializeAtom(string line);
		string getAtomName();
		char getChain();
		char getAminoChar();
		int getResidueNumber();
		int getSerialNumber();
		void setAminoChar(char newAminoChar);
		void setSerialNumber(int newSerialNum);
		void updateAtomInformation(char newAmino, string newName, string newElement);

		string ToString();

		static string getAminoString(char aminoChar);
		static int getAtomCount(char aminoChar);
};

#endif