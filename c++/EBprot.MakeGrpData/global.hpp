/* 
 * File:   global.hpp
 * Author: Hiromi Jena Koh
 *
 * Created on May 9, 2016, 2:22 PM
 */

#ifndef GLOBAL_HPP
#define	GLOBAL_HPP

#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <set>
#include <stdlib.h>
using namespace std;

struct Input_t{
    vector<string> grplabels;
    int minsample=1;
    vector<int> groupsize;
    vector<string> newlabels;
    bool logData;
};

struct Peptide_t{
    map<string, vector<double> > groupMap;
    vector<double> weightedAvg;
};


const int NULL_VAL = 0;
const double NA_VAL = 1e6;
const double pi = 3.14159265358979323846;
const double e = 2.71828182845904523536;

bool GREP(string STRING, string PATTERN);
vector<string> stringSplit(string str, char delim);
string concatenate(vector<string> vec, char delim);
vector<double> grabVec(vector<double> vec, int start, int end);
int countNonNAs(vector<double> vec);

bool readUserInput(string file, Input_t &UserInput);
vector<string> makeRatio(string file, vector<string> SampleLab, vector<int> SampleSize, map<string, vector<string> > &PROTmap, map<string, Peptide_t> &PEPmap);
void makeComparisonDat(int minsam,vector<string> labels, map<string, vector<string> > &PROTmap, map<string, Peptide_t> &PEPmap);
void outputWtData(bool transform,vector<string> header, vector<string> labels,vector<string> grplab, map<string, vector<string> >&PROTmap, map<string, Peptide_t> &PEPmap);

#endif   /* GLOBAL_HPP */
