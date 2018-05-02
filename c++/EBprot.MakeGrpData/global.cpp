/* 
 * File:   global.cpp
 * Author: Hiromi Jenna Koh
 * 
 * Created on May 9, 2015, 2:22 PM
 */

#include "global.hpp"

#include <gsl/gsl_randist.h> 
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_histogram.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <set>
#include <cmath>
#include <regex>
#include <limits>
#include <iomanip>

using namespace std;


// Function splits a string on character
vector<string> stringSplit(string str, char delim) {
    vector<string> ret;
    string buf = "";
    int i = 0;
    while((unsigned) i < str.length()) {
        if(str[i] != delim) buf += str[i];
        else if(buf.length() > 0) {
            ret.push_back(buf);
            buf = "";
        }
        i++;
    }
    if(!buf.empty()) ret.push_back(buf);
    return(ret);
}

string concatenate(vector<string> vec, char delim){
    string ret=vec.at(0);
    for(int j=1; j<vec.size(); j++){
        ret += delim + vec.at(j);
    }
    return ret;
}


bool GREP(string STRING, string PATTERN){
	bool ret=false;
	size_t found = STRING.find(PATTERN);
	if(found!= string::npos) ret = true;
	return ret;
}



//trim from start
static inline string &ltrim(string &s) {
        s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
        return s;
}

// trim from end
static inline string &rtrim(string &s) {
        s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline string &trim(string &s) {
        return ltrim(rtrim(s));
}


bool readUserInput(string file, Input_t &UserInput){
   
     bool ret = true;
     ifstream File;
     File.open(file.c_str(), ios::in);
     if(!File.is_open()){
        cerr<<"Error! Unable to open file...\n";
        exit(1);
     }
    
     string line, d;
     vector<string> v,grpLab;
     int counter=0;
     while(!File.eof()){
         getline(File,line);
         string start = line.substr(0,1);
	 regex expr("#");
	 if(regex_match(start,expr)) continue;

	 if(GREP(line,"GROUP_SIZE")){
	 	v = stringSplit(line,'=');
		d = trim(v.at(1));
		vector<string> d_new = stringSplit(d,' ');
		vector<int> dd;
		for(int j=0; j<d_new.size(); j++) dd.push_back(stoi(d_new.at(j)));
		UserInput.groupsize = dd;
	 }
	 if(GREP(line,"MIN_SAMPLE")){
		v = stringSplit(line,'=');
		d = trim(v.at(1));
	 	UserInput.minsample = stoi(d);
	 }
	 if(GREP(line,"GROUP_LABELS")){
		v = stringSplit(line,'=');
		d = trim(v.at(1));
		vector<string> lab = stringSplit(d,' ');
		UserInput.grplabels=lab;
	 }

     if(GREP(line,"IS_LOG")){
		v = stringSplit(line,'=');
		d = trim(v.at(1));
		bool ret = false;
		if(d=="TRUE"|d=="true") ret=true;
		UserInput.logData=ret;
	 }

	 int grpsize = UserInput.grplabels.size();
	 if(GREP(line,"CONTRAST_MATRIX")) {
		counter = grpsize;
                cerr<<"The contrast matrix provided is:"<<endl;
		continue;
	 }
	 if(counter>0 && (!line.empty())){
	 	cerr<<line<<endl;
		v = stringSplit(line,' ');
		if(v.size() != grpsize) {cerr<<"Problem with dimension of contrast matrix.";exit(1);}
		for(int j=0; j<v.size(); j++){
		   if(v.at(j)=="1") {
			string tmp_lab= UserInput.grplabels.at(grpsize-counter)+"/"+UserInput.grplabels.at(j);
			grpLab.push_back(tmp_lab);
		   }
		}
		counter--;
		//cerr<<counter<<endl;
		
	    	if(counter==0) UserInput.newlabels = grpLab;
	}
     }
     cerr<<"\nThere are "<<grpLab.size()<<" groups in the input data.\n"<<endl;
     cerr<<"User has specified to make the follow comparisons:\n";
     for(int k=0; k<grpLab.size(); k++) {cerr<< (k+1)<<")\t"<<grpLab.at(k)<<endl;} 
     return ret;
}


vector<string> makeRatio(string file, vector<string> SampleLab,vector<int> SampleSize,map<string, vector<string> >&PROTmap, map<string, Peptide_t> &PEPmap){

     vector<string> out;
     set<string> proteins;
     ifstream File;
     File.open(file.c_str(), ios::in);

     if(!File.is_open()){
	cerr<<"Error! Unable to open data file.\n";
	exit(1);
     }

     int ctr=1;
     string line;
     vector<string> v,v_r;
     bool skipH = true;

     while(!File.eof()){
	getline(File,line);
	if(skipH){
	    skipH=false;
	    v_r = stringSplit(line,'\t');
 	    for(int k=0; k<2;k++) out.push_back(v_r.at(k));
	    v_r.erase(v_r.begin(), v_r.begin()+2);
	    int grpsum =0;
	    for(int k=0; k<SampleSize.size(); k++) grpsum += SampleSize.at(k);
	    if(v_r.size() != grpsum){
		cerr<<"Error, the total number of sample size per group  provided in parameter file do not match the number of samples in data file.\n";
		cerr<<"Please check and try again."<<endl;
		exit(1);
	    }
	    continue;
	}

	if(line=="") break;
	v = stringSplit(line,'\t');
	string pep = v.at(0);
	string prot = v.at(1);
	vector<double> r;
        for(int i=0; i<v_r.size(); i++) {
		string tmp_d = v.at(i+2);
		double d;
		if(tmp_d=="NA"|tmp_d=="") d = NA_VAL;
		else d = stof(tmp_d);		
      		r.push_back(d);
	}

	map<string,vector<string> >::iterator m_iter;
	map<string,Peptide_t>::iterator iter;
        map<string, vector<double> > grpMap;
	
           m_iter = PROTmap.find(prot);
	   iter = PEPmap.find(pep);
	   vector<string> pep_vec;
	   Peptide_t *PepStruct = new Peptide_t();
	   if(m_iter==PROTmap.end()) {
		pep_vec.push_back(pep);
		PROTmap[prot]=pep_vec;
		proteins.insert(prot);
	   }
	   else m_iter->second.push_back(pep);

	   if(iter==PEPmap.end()) PEPmap[pep]= *PepStruct;
	   if(iter!=PEPmap.end()){
	      cerr<<"Detected duplicates of peptides in data file. Please check and try again.\n"<<endl;
	      exit(1);
	   }
	   iter = PEPmap.find(pep);
	   //vector<double> grpVec;
	   int start = 0;
	   int end;
	   for(int j=0; j<SampleLab.size(); j++){
		end = start+SampleSize.at(j);
		grpMap[SampleLab.at(j)] = grabVec(r,(start+1),end);
	        start = end;
	   }
	   iter->second.groupMap = grpMap;
	   grpMap.clear();
	   delete PepStruct;	
     }
     //map<string, Peptide_t>::iterator iter;;
     //map<string, vector<double> >::iterator g_iter;
     //int checker=0;
     //for(iter = PEPmap.begin(); iter!=PEPmap.end();iter++){
     //   if(checker<5){

	     //grpMap_tmp = iter->second.groupMap;
     //        for(int j=0; j<SampleLab.size(); j++){
     //		g_iter = iter->second.groupMap.find(SampleLab.at(j));
     //		vector<double> vec_r = g_iter->second;
     //        }
     //        cerr<<endl;
	     //cerr<<"There are "<< mapsize <<" proteins in the map of peptide "<<iter->first<<endl;
     //	     checker++;
     //   }
     //}
     cerr<<"\nThere are "<<proteins.size()<<" unique proteins and "<<PEPmap.size()<<" peptides detected from the dataset.\n"<<endl;
     return out;
}

vector<double> grabVec(vector<double> vec, int start, int end){

    vector<double> ret;
    for(int i=(start-1); i<end;i++) ret.push_back(vec.at(i));
    return ret;
}

int countNonNAs(vector<double> vec){
     int count=vec.size();
     for(int i=0; i<vec.size(); i++) {
	if(vec.at(i)==NA_VAL) count--;
     }
     return count;
}
     

void makeComparisonDat(int minsam,vector<string> labels, map<string,vector<string> > &PROTmap, map<string, Peptide_t> &PEPmap){

     map<string, vector<string> >::iterator prot_iter;
     map<string, Peptide_t>::iterator pep_iter, p_iter;
     map<string, vector<double> >::iterator g_iter, g_iter2;

     vector<string> peptides;
     vector<double> tmp_r,wtavg;
     int counter=0;
     for(prot_iter=PROTmap.begin(); prot_iter!=PROTmap.end(); prot_iter++){
 	 peptides = prot_iter->second;
	 for(int i=0; i<peptides.size(); i++){
              pep_iter = PEPmap.find(peptides.at(i));
	      //if(counter<5) cerr<<"Peptide : \t"<<peptides.at(i)<<endl;
	      for(int j=0; j<labels.size(); j++){
		   g_iter = pep_iter->second.groupMap.find(labels.at(j));
	           tmp_r = g_iter->second;
                   double den =0.0;
		   double avg=0.0;
		   double retavg;
   		   bool hasData=false;
		   if(countNonNAs(tmp_r)< minsam) retavg = NA_VAL;
		   else{
       		     for(int k=0; k<tmp_r.size(); k++){
		         double num = 0.0;
		         double curr_v = tmp_r.at(k);
		         if(curr_v!= NA_VAL){
			     hasData = true;
			     //if(counter<5) cerr<<"data: "<<curr_v<<endl;
			     for(int h=0; h<peptides.size(); h++){
				   p_iter = PEPmap.find(peptides.at(h));
			           g_iter2 = p_iter->second.groupMap.find(labels.at(j));
				   double tt = g_iter2->second.at(k);
			           if(tt!=NA_VAL) {num++;den++;}
			      }
                              avg += (num*curr_v);
			      //if(counter<5){
			      //cerr<<"Numerator: "<<num<<endl;
			      //cerr<<"Denominator: "<<den<<endl;
			      //}
		          }
                      }
                      if(!hasData) retavg = NA_VAL;
	              if(hasData) retavg = (avg/den);
	           }
	        wtavg.push_back(retavg);
		}//end labels
                pep_iter->second.weightedAvg = wtavg;
	        wtavg.clear();
       		counter++;
         }//end peptide
         //if(counter<5) cerr<<"\n----------------------------------------------------\n"<<endl;
     }//end proteins

}


void outputWtData(bool transform,vector<string> header,vector<string> labels, vector<string> grplab,map<string, vector<string> >&PROTmap, map<string, Peptide_t> &PEPmap){

   ofstream outF("weighted_data.txt");
   ofstream outF2("weighted_grpcomparisons.txt");
   string GrpA, GrpB;
   map<string, vector<string> >::iterator prot_iter;
   map<string, Peptide_t>::iterator pep_iter;
   vector<string> peptides;
   string prot,pep, grp;   
   vector<double> wtavg;
   double grpAvg, grpA, grpB;

   outF<<header.at(0)<<"\t"<<header.at(1);
   for(int j=0; j<grplab.size(); j++) outF<<"\t"<<grplab.at(j);
   outF<<endl;

   outF2<<header.at(0)<<"\t"<<header.at(1);
   for(int j=0; j<labels.size(); j++) outF2<<"\t"<<labels.at(j);
   outF2<<endl;
   int count=0;

   for(prot_iter=PROTmap.begin(); prot_iter!=PROTmap.end(); prot_iter++){
 	peptides = prot_iter->second;
  	prot = prot_iter->first;
        for(int i=0; i<peptides.size(); i++){
	      pep = peptides.at(i);
	      pep_iter = PEPmap.find(pep);
	      wtavg = pep_iter->second.weightedAvg;
	      if(countNonNAs(wtavg)>0){
	        outF<<pep<<"\t"<<prot;
	        for(int k=0; k<wtavg.size(); k++) {
		   if(wtavg.at(k)!=NA_VAL) outF<<"\t"<<wtavg.at(k);
		   else outF<<"\t"<<"NA";
	        }
	        outF<<endl;
  	        vector<double> AVG;
                for(int k=0; k<labels.size();k++){
		   vector<string> tmpL = stringSplit(labels.at(k),'/');
                   GrpA = tmpL.at(0);
		   GrpB = tmpL.at(1);
		   for(int h=0; h<grplab.size();h++){
		      if(grplab.at(h)==GrpA) grpA = wtavg.at(h);
	              if(grplab.at(h)==GrpB) grpB = wtavg.at(h);
		   }
		   if(transform) grpAvg = grpA - grpB;
		   if(!transform) grpAvg = grpA/grpB;
	           if(grpA==NA_VAL|grpB==NA_VAL) grpAvg = NA_VAL;
	           AVG.push_back(grpAvg);
		 }
		 if(countNonNAs(AVG)>0){
		      outF2<<pep<<"\t"<<prot;
		      for(int h=0;h<AVG.size();h++){
		         double val = AVG.at(h);
		         if(val==NA_VAL) outF2 <<"\t"<<"NA";
 		         if(val!=NA_VAL) outF2<< "\t"<<val;
		      }
		      outF2<<endl;	
		 }//end AVGIF
	       }//end wtavgIF
	}//end for loop
    }
    outF.close();
    outF2.close();

}
