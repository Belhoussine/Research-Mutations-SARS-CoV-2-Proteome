#ifndef PROTEIN_H
#define PROTEIN_H

#include <string>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "chain.h"
using namespace std;

class Protein {
	private:
		string proteinName;
		string changeComment;
		string preRemarks;
		string postRemarks;
		vector<Chain> chains;

		void init();

	public:
		Protein();
		Protein(ifstream &pdb);
		Protein(const Protein &other);
		
		bool isMultiChain();
		void initializeProtein(ifstream &pdb);
		void createFASTAFile(ofstream &fasta, char chain, int residueNum, char mutant);
		void createPDBFile_scwrl(ofstream &pdb, char chain);
		void createPDBFile(Protein scwrlChain, ofstream &finalPDB, char chain);
		void addToChangeComment(string append);
		int getMaxSerial();
		int updateSerials(int minSerialNum);
		int getChainNumber(char chain, int residueNum);
		char getResidueType(char chain, int residueNum);
		Chain getChainWithResidue(char chain, int residueNum);
		void setChainWithResidue(char chain, int residueNum, Chain newChain);
		vector<char> getChainIDs();
		int getMinResidueNum(char chain);
		int getMaxResidueNum(char chain);

		string ToString();
		string ToString_ChainsOnly();
};

#endif