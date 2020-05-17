#ifndef CHAIN_H
#define CHAIN_H

#include <string>
#include <vector>
#include <stdlib.h>
#include "atom.h"
#include "utility.h"
using namespace std;

class Chain {
	private:
		char chainID;
		vector<Atom> atoms;

	public:
		Chain();
		Chain(const Chain &other);
		
		bool isChain(char chainID);
		bool containsResidue(int residueNum);
		vector<Atom> getAtoms();
		int getNumAtoms();
		int getMaxSerial();
		void addAtom(string line);
		void addAtom(Atom newAtom);
		void clearChain();
		int updateSerials(int minSerialNum);
		char getResidueType(int residueNum);
		char getChainID();
		int getMinResidueNum();
		int getMaxResidueNum();

		string fastaString(int residueNum, char destination);
		string ToString();
};

#endif