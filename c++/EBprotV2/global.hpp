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

    string ExpDesign;
    bool log2transform;
    bool outlierRm;
    vector<string> labels;
    int minpep=1;
    int mink=5;
    double threshold = 0.05;
    double lbound=0.2;
    double rbound=0.8;
};


struct Protein_t{
    vector<double> PPscore;
    vector<double> PostOdds;
    vector<double> bfdr;
    vector<double> medratio;
    vector<int> numPep;
    vector<int> numPepRM;
};

const int NULL_VAL = 0;
const double NA_VAL = 1e6;
const double pi = 3.14159265358979323846;
const double e = 2.71828182845904523536;

bool GREP(string STRING, string PATTERN);
vector<string> stringSplit(string str, char delim);
string concatenate(vector<string> vec, char delim);
bool allMissing(vector<double> vec);
void rmMissing(vector<double> &vec);
void TruncateVal(vector<double> &vec);
vector<double> getNull(vector<double> vec, double UpperB, double LowerB);
vector<double> UniqueVec(vector<double> vec);
double TrapeziumRule(vector<double> x, vector<double> y);
vector<string> NotIn(set<string> set1, set<string> set2);
vector<string> NotIn2(set<string>set1, set<string> set2);
set<string> getKey(map<string, vector<string> > &MAP);
void cleanMap(map<string, Protein_t> &MAP);


//statistics measure
double CalSD(vector<double> vec);
double CalMeanMed(vector<double> vec,int type);
double getMode(vector<double> ratio,int Nbins);
double getModeDensity(vector<double> den, vector<double> r);
double DiesAt(vector<double> den, vector<double> r, int type);
vector<double> computeDensity(vector<double> x);
vector<double> getquantile(vector<double> vec, double lb, double rb);
double getSum(vector<double> vec);
double getMin(vector<double> vec);
double getMax(vector<double> vec);
vector<double> NORMpdf(vector<double> vec);
double pdfStdNorm(double x, double m, double ss);
vector<double> logNormTrun(vector<double> vec, double mu, double grid, double a, double b);
vector<double> logNormal(vector<double> vec, double mu, double v);
double sdTruncatedNormal(vector<double> v, double mu, double a, double b);
vector<double> NORMpdf2(vector<double> v, double m, double ss);

bool readUserInput(string file, Input_t &UserInput);
set<string> readData(string file, vector<string> SampleLab, bool transform, map<string, vector<string> >&PROTmap, map<string,vector<double> >&PEPmap);
vector<int> OutlierRm(int mink, map<string, vector<string> >&PROTmap, map<string, double> &PEPmap);
set<string> createMap(int index, int minPep, map<string,vector<string> >&PROTmap, map<string, vector<string> >&PROTmapN, map<string,vector<double> >&PEPmap, map<string, double> &PEPmapN);
set<string> createMap2(int numrep, int minPep, map<string,vector<string> >&PROTmap, map<string, vector<string> >&PROTmapN, map<string,vector<double> >&PEPmap, map<string, double> &PEPmapN);
vector<double> nonparamFit(vector<string> protNA,vector<int> pepRm,map<string, vector<string> >&PROTmap, map<string, double>&PEPmap, Input_t &UserInput, map<string, Protein_t> &PROTEIN, map<double,vector<double> >&densityMap);
void outputDensity(string lab,vector<double> ratio, map<double, vector<double> >&DENmap);
vector<double> EstProtPiEM(vector<double> pi,map<string, vector<string> >&PROTmap, map<string, double>&PEPmap, map<double,vector<double> >&DENmap);
void computePPscore(vector<double> pi, map<string, vector<string> >&PROTmap, map<string, double>&PEPmap,map<double, vector<double> >&DENmap, map<string, Protein_t>&PROTEIN);
void computeBFDR(double bfdr_thres, int index,map<string, Protein_t> &PROTEIN);
void outputResult(Input_t *UserInput, map<string, Protein_t> &PROTEIN);
void UnmergeMap(int len,map<string, vector<string> > &PROTmap, map<string, Protein_t> &PROTEINfit, map<string, Protein_t> &PROTEIN);


#endif   /* GLOBAL_HPP */
