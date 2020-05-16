#include <vector>
#include "utility.h"
#include "protein.h"

vector<char> getChains(Protein inputProtein, string chains) {
	vector<char> output;
	if (chains.length() != 3) {
		return output;
	}
	
	char startChain = chains[0];
	char endChain = chains[2];
	
	vector<char> proteinChains = inputProtein.getChainIDs();
	int numChains = proteinChains.size();
	int startPosition = -1;
	int endPosition = -1;
	
	if (startChain == 'X' && endChain == 'X') {
		return proteinChains;
	} else if (startChain == 'X') {
		startPosition = 0;
	} else if (endChain == 'X') {
		endPosition = numChains - 1;
	}
	
	for (int i = 0; i < numChains; i++) {
		if (startPosition == -1 && startChain == proteinChains[i]) {
			startPosition = i;
		}
		if (endPosition == -1 && endChain == proteinChains[i]) {
			endPosition = i;
		}
		
		if (startPosition != -1 && endPosition != -1) {
			break;
		}
	}
	
	if (startPosition == -1 || endPosition == -1) {
		return output;
	}
	
	for (int i = startPosition; i <= endPosition; i++) {
		output.push_back(proteinChains[i]);
	}
	return output;
}

vector<int> getResiduesFromArg(string residues) {
	vector<int> output;
	int stringLen = residues.length();
	
	if (stringLen == 0) {
		return output;
	}
	
	int positionOfComma = residues.find_first_of(",");
	if (positionOfComma == 0 || positionOfComma == stringLen - 1) {
		return output;
	} else if (positionOfComma != -1) {
		vector<int> left, right;
		int leftCount, rightCount;
		
		left = getResiduesFromArg(residues.substr(0,positionOfComma));
		right = getResiduesFromArg(residues.substr(positionOfComma + 1, stringLen - (positionOfComma + 1)));
		
		leftCount = left.size();
		rightCount = right.size();
		
		if (leftCount == 0 || rightCount == 0) {
			return output;
		}
		
		for (int i = 0; i < leftCount; i++) {
			output.push_back(left[i]);
		}
		for (int i = 0; i < rightCount; i++) {
			output.push_back(right[i]);
		}
	} else {
		int positionOfColon = residues.find_first_of(":");
		if (positionOfColon == 0 || positionOfColon == stringLen - 1) {
			return output;
		} else if (positionOfColon != -1) {
			int leftInt = atoi(residues.substr(0,positionOfColon).c_str());
			int rightInt = atoi(residues.substr(positionOfColon + 1, stringLen - (positionOfColon + 1)).c_str());
			
			if (residues.substr(0,positionOfColon) == "X") {
				leftInt = 1;
			}
			if (residues.substr(positionOfColon + 1, stringLen - (positionOfColon + 1)) == "X") {
				rightInt = 9999;
			}
			
			for (int i = leftInt; i <= rightInt; i++) {
				output.push_back(i);
			}
		} else {
			output.push_back(atoi(residues.c_str()));
		}
	}
	
	return output;
}

vector<int> getResidues(Protein inputProtein, char chain, string residues) {
	vector<int> output;
	
	vector<int> arguments = getResiduesFromArg(residues);
	int numToCheck = arguments.size();
	
	if (numToCheck == 0) {
		return output;
	}
	
	int minResidue = inputProtein.getMinResidueNum(chain);
	int maxResidue = inputProtein.getMaxResidueNum(chain);
	int currentArg;
	
	for (int i = 0; i < numToCheck; i++) {
		currentArg = arguments[i];
		if (currentArg >= minResidue && currentArg <= maxResidue) {
			output.push_back(currentArg);
		}
	}
	return output;
}

vector<char> getTargets(string targets) {
	vector<char> output;
	
	if (targets == "X") {
		output.push_back('A');
		output.push_back('C');
		output.push_back('D');
		output.push_back('E');
		output.push_back('F');
		output.push_back('G');
		output.push_back('H');
		output.push_back('I');
		output.push_back('K');
		output.push_back('L');
		output.push_back('M');
		output.push_back('N');
		output.push_back('P');
		output.push_back('Q');
		output.push_back('R');
		output.push_back('S');
		output.push_back('T');
		output.push_back('V');
		output.push_back('W');
		output.push_back('Y');
	} else if (targets == "POL") {
		output.push_back('C');
		output.push_back('H');
		output.push_back('M');
		output.push_back('N');
		output.push_back('Q');
		output.push_back('S');
		output.push_back('T');
		output.push_back('W');
		output.push_back('Y');
	} else if (targets == "CHAR") {
		output.push_back('D');
		output.push_back('E');
		output.push_back('K');
		output.push_back('R');
	} else if (targets == "PHOBIC") {
		output.push_back('A');
		output.push_back('F');
		output.push_back('G');
		output.push_back('I');
		output.push_back('L');
		output.push_back('P');
		output.push_back('V');
	} else if (targets.size() == 1) {
		output.push_back(targets[0]);
	}
	
	return output;
}

int main(int argc, char* argv[]) {
	if (argc != 6) {
		cout << "Invalid number of arguments" << endl;
		cout << "Usage: " << argv[0] << " protein chains residues targets outputname" << endl;
		return 1;
	}

	string proteinName = toUpper(string(argv[1]));
	string chains = toUpper(string(argv[2]));
	string residues = toUpper(string(argv[3]));
	string targets = toUpper(string(argv[4]));

	string PDBFilename = proteinName + ".pdb";
	string outputFilename = string(argv[5]);

	ifstream inputPDB;
	if (!openPDBFile(PDBFilename, inputPDB)) {
		return -1;
	}

	Protein inputProtein(inputPDB);
	cleanup(inputPDB);
	
	ofstream outputScript;
	if (!openFile(outputFilename, outputScript)) {
		return -1;
	}
	
	vector<char> chainVector = getChains(inputProtein, chains);
	vector<char> targetVector = getTargets(targets);
	vector<int> resVector;
	string chainOutput;
	
	int numChains = chainVector.size();
	int numTargets = targetVector.size();
	int numResidues;
	for (int i = 0; i < numChains; i++) {
		chainOutput = "./proMute " + proteinName + " " + chainVector[i] + " ";
		resVector = getResidues(inputProtein, chainVector[i], residues);
		numResidues = resVector.size();
	
		for (int j = 0; j < numResidues; j++) {
			for (int k = 0; k < numTargets; k++) {
				outputScript << chainOutput << resVector[j] << " " << targetVector[k] << endl;
			}
		}
	}
	
	cleanup(outputScript);
}
