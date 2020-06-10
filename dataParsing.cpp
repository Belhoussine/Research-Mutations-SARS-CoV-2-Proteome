#include "iostream"
#include <fstream>
#include <string> 
#include <stdio.h>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <iomanip> 
#include <random>
#include <errno.h>
#include <exception>
#include <cerrno>
#include <string.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <sys/types.h>
#include <dirent.h>

using namespace std;

struct fileInfo{
    string name;
    int step;
    double value;
    double diff;
};

typedef std::vector<std::string> stringvec;

const double wtEnergy = -5491.0831;

//If given file exists, returns stream of file. Else exits with error. 
ifstream getFile(string file);

void getFileInfo(ifstream &myFile, vector<fileInfo> &arr, string name);
 
void read_directory(const std::string& name, stringvec& v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}

int main() {

//Make arrayList of fileInfo

stringvec v;
read_directory(".", v);

vector<fileInfo> myArray;
vector<fileInfo> &arr = myArray;

//Go through all files
for(vector <string> :: iterator it = v.begin(); it != v.end(); ++it){
    string currTitle = *it;
    if((currTitle.compare(".")!=0)&&(currTitle.compare("..")!=0)&&(currTitle.compare(".vscode")!=0)&&
    (currTitle.compare("dataParsing")!=0)&&(currTitle.compare("dataParsing.cpp")!=0)&&(currTitle.compare("results.csv")!=0)){
        ifstream balFile = getFile(*it);
        ifstream &myFile = balFile;
        getFileInfo(myFile, arr, *it);
    }
}
int i=0;
ofstream myfile;
myfile.open ("results.csv");
for(vector <fileInfo> :: iterator it = arr.begin(); it != arr.end(); ++it){
    fileInfo currFile = *it;
    myfile<<currFile.name<< "," <<currFile.step<< "," <<currFile.value<< "," <<currFile.diff<<endl;
}



// cout<<arr.at(3).name<<endl;
// cout<<arr.at(3).step<<endl;
// cout<<arr.at(3).value<<endl;
// cout<<arr.at(3).diff<<endl;

// std::copy(v.begin(), v.end(),
// std::ostream_iterator<std::string>(std::cout, "\n"));

}

//If given file exists, returns stream of file. Else exits with error.
ifstream getFile(string file){
    ifstream myFile;
    myFile.open(file, ios::in);
    if(myFile.fail()){
        //Upon file failing to open, prints error message and exits
        cerr << "Error \"" << strerror(errno) << "\" when opening file: " << file <<endl;
        // exit(1);
    }
    return myFile;
}

//Takes an array of accounts and the file, then goes through to fill up the array with the accounts from the file.
void getFileInfo(ifstream &myFile, vector<fileInfo> &arr, string name){
    string line;
    string::size_type sz; 

    fileInfo newAddition;
    newAddition.name = name;
    // cout<<newAddition.name<<endl;
    //getting file to right section (after all the junk at the start)
    int i=0;
    while(i<153){
        getline(myFile, line);
        i++;
    }

    //Going through energys and finding smallest difference in absolute value.
    //results.at(1) == step
    //results.at(11) == value
    double minAbs = abs(wtEnergy);

    while(getline(myFile, line)){
        if(line.find("ENERGY")!=string::npos){
            istringstream s_stream(line);
            std::vector<std::string> results(std::istream_iterator<std::string>{s_stream},
                                 std::istream_iterator<std::string>());
            //If the absolute value of the wildtype energy - the value of the current step is less than the min absolute value so far found, change minAbsolute value and save
            if((abs((wtEnergy-stod(results.at(11), &sz))))<minAbs){
                minAbs = abs((wtEnergy-stod(results.at(11), &sz)));
                newAddition.step = stoi(results.at(1), &sz);
                newAddition.value =stod(results.at(11), &sz);
                newAddition.diff = abs((wtEnergy-stod(results.at(11), &sz)));
            }
        }
    }
    arr.push_back(newAddition);
}

