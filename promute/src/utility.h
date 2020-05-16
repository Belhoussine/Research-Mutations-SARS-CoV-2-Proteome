#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <string.h>

using namespace std;

const string PDB_URL = "https://files.rcsb.org/download/";

char toLower(char character);
char toUpper(char character);
string toUpper(string str);
bool openFile(string filename, ofstream &file);
bool openFile(string filename, ifstream &file);
bool openPDBFile(string filename, ifstream &file);
bool downloadPDBFile(string PDBFilename, ifstream &pdb);
void cleanup(ofstream &file);
void cleanup(ifstream &file);

#endif