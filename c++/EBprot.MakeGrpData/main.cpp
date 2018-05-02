
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: Hiromi Jenna Koh
 *
 * Created on 26 Aug, 2016, 11:03 AM
 */

#include "global.hpp"

#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <regex>
#include <ctime>

using namespace std;

/*
 *
 */


int main(int argc, char** argv) {
    
    if(argc < 2){
        cerr<<"\nUSAGE: EBprot.MakeGrpData.exe <data.file> <parameter.file>\n";
        return 0;
    }

    clock_t begin = clock();
    
    string dataF = argv[1];
    string paramF = argv[2];
    Input_t *UserInput = new Input_t();
    //Peptide_t *PepStruct = new Peptide_t();
    map<string, vector<string> > ProtPepMap;
    map<string, Peptide_t> PeptideMap;
    cerr<<"-------------------------------------------------------------------------"<<endl;    
    cerr<< "\nInitiating EBprot.MakeGrpData...\n\nMaking data for group comparisons...\n" <<endl;
    readUserInput(paramF, *UserInput);
    
    vector<string> header= makeRatio(dataF,UserInput->grplabels, UserInput->groupsize, ProtPepMap, PeptideMap);

    makeComparisonDat(UserInput->minsample, UserInput->grplabels, ProtPepMap,PeptideMap);
    cerr<<"Generating output file now...\n"<<endl;
    outputWtData(UserInput->logData,header,UserInput->newlabels,UserInput->grplabels,ProtPepMap,PeptideMap);
   
    clock_t end = clock();
    double elapse_sec = double(end - begin)/(double)CLOCKS_PER_SEC;
    cerr<< "EBprot.MakeGrpData module has completed running in " << elapse_sec <<" seconds.\n" <<endl;
    cerr<<"Please run EBprotV2 software for differential expression analysis.\n";
    cerr<< "\n------------------------------ Completed ------------------------------\n\n"<<endl;
    
    return 0;
    
}



