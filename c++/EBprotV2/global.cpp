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

set<string> getKey(map<string, vector<string> >&MAP){
     set<string> ret;
     map<string, vector<string> >::iterator iter;
     for(iter=MAP.begin(); iter!=MAP.end(); iter++) ret.insert(iter->first);
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

bool allMissing(vector<double> vec){
   bool ret=true;
   for(int i=0; i<vec.size(); i++){
	if(vec.at(i)!=NA_VAL) ret=false;
   }
   return ret;
}   


double CalSD(vector<double> vec){
        
	double ret;
	double mean = CalMeanMed(vec,1);
        double temp = 0;
	double densize =0;
        for(int i = 0; i < vec.size(); i++)
        {
            if(vec.at(i)!=NA_VAL){
		temp += ((vec.at(i) - mean) * (vec.at(i) - mean)) ;
        	densize++; 
	    }
	}
        ret = sqrt(temp /(densize-1));
	return ret;
}

double CalMeanMed(vector<double> vec, int type){

    double ret;
    // type 1 calculate Mean
    if(type==1){
	double sum = 0;
	double densize =0;
        for(int i = 0; i < vec.size(); i++){
            if(vec.at(i)!=NA_VAL){
		sum += vec.at(i);
		densize++;
	    }
	}
     	ret = (sum / densize);
     }
     //if type 2 calculate median
     if(type==2){
        rmMissing(vec);
        sort(vec.begin(), vec.end());
    	if ((vec.size()%2) == 0) {ret = (vec.at((vec.size()/2)-1)+vec.at(vec.size()/2))/2;}
    	else  {
		div_t res = div(vec.size(), 2);
		ret = vec.at(res.quot);
	}
     }
     return ret;
}

void rmMissing(vector<double> &vec){

   for(int i=0; i<vec.size(); i++){
	if(vec.at(i) == NA_VAL) vec.erase(vec.begin()+i);
   }
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
     vector<string> v;

     while(!File.eof()){
         getline(File,line);
         string start = line.substr(0,1);
	 regex expr("#");
	 if(regex_match(start,expr)) continue;

         if(GREP(line,"EXPERIMENTAL_DESIGN")){
		v = stringSplit(line,'=');
	        d = trim(v.at(1));
		UserInput.ExpDesign=d;
	 }
	 if(GREP(line,"LOG2_TRANSFORM")){
	        v = stringSplit(line,'=');
		d = trim(v.at(1));
		bool d_new=false;
		if(d=="true") d_new=true;
		UserInput.log2transform=d_new;
	 }
	 if(GREP(line,"OUTLIER_RM")){
	   	v = stringSplit(line,'=');
		d = trim(v.at(1));
		bool d_new=false;
		if(d=="true") d_new=true;
		UserInput.outlierRm=d_new;
	 }
	 if(GREP(line,"MIN_PEP")){
	 	v = stringSplit(line,'=');
		d = trim(v.at(1));
		UserInput.minpep=stoi(d);
	 }
	 if(GREP(line,"MIN_K")){
		v = stringSplit(line,'=');
		d = trim(v.at(1));
	 	UserInput.mink = stoi(d);
	 }
         if(GREP(line,"BFDR_THRESHOLD")){
	        v = stringSplit(line,'=');
	        d = trim(v.at(1));
	        UserInput.threshold = stof(d);
	 }
	 if(GREP(line,"LEFT_B")){
		v = stringSplit(line,'=');
	 	d = trim(v.at(1));
	        UserInput.lbound=stof(d);
	 }
	 if(GREP(line,"RIGHT_B")){
		v = stringSplit(line,'=');
		d = trim(v.at(1));
		UserInput.rbound=stof(d);
	 }
	 if(GREP(line,"LABELS")){
		v = stringSplit(line,'=');
		d = trim(v.at(1));
		vector<string> lab = stringSplit(d,' ');
		UserInput.labels=lab;
	 }

     }
     cerr<<"The type of experimental design specified is "<<UserInput.ExpDesign<<endl;
     cerr<<"There are "<< UserInput.labels.size()<<" samples in the input data."<<endl;

     return ret;
}


set<string> readData(string file, vector<string> SampleLab,bool transform, map<string, vector<string> >&PROTmap, map<string, vector<double> >&PEPmap){

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
	    v_r.erase(v_r.begin(), v_r.begin()+2);
	    if(v_r.size() != SampleLab.size()){
		cerr<<"Error, the number of sample labels provided in parameter file do no match number of samples in data file.\n";
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
		if(transform){
		   if(tmp_d=="NA"|tmp_d=="0"|tmp_d=="") d = NA_VAL;
		   else d = log(stof(tmp_d))/log(2);
		}
		if(!transform){
		   //if data already log, zero implies a ratio of 1, different from NA.
		   if(tmp_d=="NA"|tmp_d=="") d = NA_VAL;
		   else d = stof(tmp_d);
		}
		r.push_back(d);
	}
        bool Missing = allMissing(r);

	map<string,vector<string> >::iterator m_iter;
	map<string, vector<double> >::iterator iter;
	
        if(!Missing){
	   //only insert protein if it has non missing ratios across samples 
           proteins.insert(prot);
           m_iter = PROTmap.find(prot);
	   iter = PEPmap.find(pep);
	   vector<string> pep_vec;

	   if(m_iter==PROTmap.end()) {
		pep_vec.push_back(pep);
		PROTmap[prot]=pep_vec;
	   }
	   else m_iter->second.push_back(pep);

	   if(iter==PEPmap.end()) PEPmap[pep]=r;
	   if(iter!=PEPmap.end()){
	      cerr<<"Detected duplicates of peptides in data file. Please check and try again.\n"<<endl;
	      exit(1);
	   }
	}  	
     }
     return proteins;
}


vector<int> OutlierRm(int mink, map<string, vector<string> >&PROTmap, map<string,double> &PEPmap){

    double ref_sd=0;
    map<string, vector<string> >::iterator m_iter;
    map<string, double>::iterator p_iter;
    vector<double> SDvec;
   
    for(m_iter = PROTmap.begin(); m_iter!=PROTmap.end(); m_iter++){

	string curr_prot = m_iter->first;
	vector<string> curr_pep = m_iter->second;
	// find reference SD based on protein with at least mink peptides
	if(curr_pep.size()>= mink){
           vector<double> pepR;
	   for(int i=0; i<curr_pep.size();i++){
	      p_iter = PEPmap.find(curr_pep.at(i));
	      pepR.push_back(p_iter->second);
	   }
	   double tmp_sd = CalSD(pepR);
	   SDvec.push_back(tmp_sd);
	}	
    }
    cerr<<"SDvec size: "<<SDvec.size();
    ref_sd = CalMeanMed(SDvec,2);
    cerr<<"\nIn the Outlier Filtering step, "<<SDvec.size()<<" proteins with at least "<<mink<<" peptides were used to calculate a reference SD."<<endl;
    cerr<<"The median of the standard deviations of the peptide distribution is "<< ref_sd <<endl;
        
    vector<string> protRM;
    vector<int> numPepRm;
    int counter=0;

    // filtering is only done on proteins with at least 2 peptides
    for(m_iter = PROTmap.begin();m_iter!=PROTmap.end(); m_iter++){
	
	bool pepleft=false;
        string curr_prot = m_iter->first;
	vector<string> curr_pep = m_iter->second;
        vector<double> pepR;
	double pepr;
	int peprm=0;
 
	if(curr_pep.size()>1){
	   for(int i=0; i<curr_pep.size(); i++){
	      p_iter = PEPmap.find(curr_pep.at(i));
	      pepr = p_iter->second;
	      if(pepr!=NA_VAL) pepR.push_back(pepr);
	   }
           double ubound = (CalMeanMed(pepR,2)+ 4*ref_sd);
	   double lbound = (CalMeanMed(pepR,2)- 4*ref_sd);
           vector<string> newcurr_pep;

	   for(int k=0; k<curr_pep.size(); k++){
		p_iter = PEPmap.find(curr_pep.at(k));
		pepr = p_iter->second;
  		if(pepr>ubound | pepr<lbound) {
		    PEPmap.erase(p_iter);
		    peprm++;
		}
	        else{	
		    pepleft=true;
		    newcurr_pep.push_back(curr_pep.at(k));
		}
	    }
	    m_iter->second = newcurr_pep;
            if(!pepleft) {
		protRM.push_back(curr_prot);
		//PROTmap.erase(m_iter);
	    }
        }
	if(find(protRM.begin(), protRM.end(), curr_prot)== protRM.end()) numPepRm.push_back(peprm);	
	counter++;	
    }

    if(protRM.size()>0){
	  int nprotRM = protRM.size();
	  if(nprotRM>30) {
		cerr<<"\n\nThere are "<<nprotRM<<" proteins removed in this step and the first 30 are:\n";
	  }
	  else cerr<<"\n\nThe list of "<<nprotRM<<" proteins removed in this step are:\n";
          cerr<<"--------------------------------------------------------------"<<endl;
	  for (int j=0; j<nprotRM; j++) {
		if(j>29) continue;
		cerr<<j+1<<") "<< protRM.at(j)<<endl;
	  }
	        
          cerr<<"---------------------------------------------------------------"<<endl;   

    }

    // remove from PROTmap
    for(int j=0; j<protRM.size(); j++){
	m_iter = PROTmap.find(protRM.at(j));
        PROTmap.erase(m_iter);
    }

    return numPepRm;

}



set<string> createMap(int index,int minPep, map<string,vector<string> >&PROTmap,map<string,vector<string> >&PROTmapN, map<string, vector<double> >&PEPmap, map<string, double>&PEPmapN){

     set<string> ret, protRM;
     map<string, vector<string> >::iterator m_iter;
     map<string, vector<double> >::iterator p_iter;
 
     for(m_iter =PROTmap.begin(); m_iter!=PROTmap.end(); m_iter++){

	string curr_p = m_iter->first;
	vector<string> pep_vec = m_iter->second;
        double tmp_r;
	vector<string> peptidesN;

        for(int i=0; i<pep_vec.size(); i++){
	    string curr_pep = pep_vec.at(i);
	    p_iter = PEPmap.find(curr_pep);
	    double val_r = p_iter->second.at(index);
	    if(val_r!=NA_VAL) {
		PEPmapN[curr_pep] = val_r;
	 	peptidesN.push_back(curr_pep);
	    }
        }
        if(peptidesN.size()>0 & peptidesN.size()<minPep) protRM.insert(curr_p); 
	if(peptidesN.size()>=minPep) {
	    PROTmapN[curr_p] = peptidesN;
	    ret.insert(curr_p);
	}
    }
    cerr<<"There were "<< protRM.size()<<" proteins removed as a results of having less than the minimum peptides required by user.\n"<<endl;

    return ret;
}


set<string> createMap2(int numrep,int minPep, map<string,vector<string> >&PROTmap,map<string,vector<string> >&PROTmapN, map<string, vector<double> >&PEPmap, map<string, double>&PEPmapN){

     set<string> ret,protRM;
     map<string, vector<string> >::iterator m_iter;
     map<string, vector<double> >::iterator p_iter;
     string curr_p, curr_pnew; 
     int counter=0;
     for(m_iter =PROTmap.begin(); m_iter!=PROTmap.end(); m_iter++){

	curr_p = m_iter->first;
	vector<string> pep_vec = m_iter->second;
	vector<string> peptidesN;
        for(int j=0; j<numrep; j++){
	   curr_pnew = curr_p +"______R"+ to_string(j);
           for(int i=0; i<pep_vec.size(); i++){
	       string curr_pep = pep_vec.at(i);
	       p_iter = PEPmap.find(curr_pep);
	       double val_r = p_iter->second.at(j);
	       string curr_pepnew;
	       if(val_r!=NA_VAL) {
	            string curr_pepnew = curr_pep + "______R"+ to_string(j);
		    PEPmapN[curr_pepnew] = val_r;
	 	    peptidesN.push_back(curr_pepnew);
	       }
	   }
           if(peptidesN.size()>0 & peptidesN.size()<minPep) protRM.insert(curr_pnew);     
	   if(peptidesN.size()>=minPep) {
	      PROTmapN[curr_pnew] = peptidesN;
              ret.insert(curr_pnew);
              if(counter<10){
                if(counter==0) {
			cerr<<"\nFormating data into long format...\n";
			cerr<<"---------------------------------------"<<endl;
		}
		cerr<<"Protein: "<<curr_pnew<<endl;
		cerr<<"Peptides: "<<endl;
		for(int l=0; l<peptidesN.size(); l++) cerr<<peptidesN.at(l)<<endl;
                cerr<<"----------------------------------------------"<<endl;
	      }
	      counter++;
	   }
           peptidesN.clear();
        }

    }

    cerr<<"There are "<< protRM.size()<<" proteins removed as a result of having less than the minimum peptides required by user.\n"<<endl;

    return ret;
}






vector<double> nonparamFit(vector<string> protNA,vector<int> pepRm,map<string, vector<string> >&PROTmap, map<string, double>&PEPmap, Input_t &UserInput, map<string, Protein_t>&PROTEIN, map<double, vector<double> >&densityMap){

     map<string, vector<string> >::iterator prot_iter;
     map<string, double>::iterator pep_iter;
     map<string, Protein_t>::iterator P_iter;
     vector<double> allRatio;
     int count=0, pos;

     cerr<<"PROTmap size: "<<PROTmap.size()<<endl;
     cerr<<"length of pepRM: "<<pepRm.size()<<endl;

     for(prot_iter=PROTmap.begin(); prot_iter!=PROTmap.end(); prot_iter++){
	 string curr_p = prot_iter->first;
         Protein_t *Protein = new Protein_t();
	 P_iter = PROTEIN.find(curr_p);
         if(P_iter==PROTEIN.end()) {
		PROTEIN[curr_p] = *Protein;
		P_iter = PROTEIN.find(curr_p);
	 }
	 vector<string> peptides = prot_iter->second;
	 vector<double> pepr;
	 for(int j=0; j<peptides.size();j++){
		string curr_pep = peptides.at(j);
	  	pep_iter = PEPmap.find(curr_pep);
		pepr.push_back(pep_iter->second);
	        allRatio.push_back(pep_iter->second);
	 }
	 rmMissing(pepr);
	 P_iter->second.numPep.push_back(pepr.size());
         if(pepRm.size()>0) P_iter->second.numPepRM.push_back(pepRm.at(count));
	 else P_iter->second.numPepRM.push_back(0);
         P_iter->second.medratio.push_back(CalMeanMed(pepr,2));
	 pos = P_iter->second.medratio.size();
	 count++;
         delete Protein;
         pepr.clear();
         peptides.clear();
     }

     if(protNA.size()>0){
     for(int j=0; j<protNA.size(); j++){
	string curr_p = protNA.at(j);
	Protein_t *Protein = new Protein_t();
	P_iter = PROTEIN.find(curr_p);
	if(P_iter!=PROTEIN.end()){
	  P_iter->second.numPep.push_back(NA_VAL);
	  P_iter->second.numPepRM.push_back(NA_VAL);
	  P_iter->second.medratio.push_back(NA_VAL);
	  P_iter->second.PPscore.push_back(NA_VAL);
	  P_iter->second.PostOdds.push_back(NA_VAL);
	}
	if(P_iter==PROTEIN.end()){
	   PROTEIN[curr_p] = *Protein;
	   P_iter = PROTEIN.find(curr_p);
	   P_iter->second.numPep.push_back(NA_VAL);
	   P_iter->second.numPepRM.push_back(NA_VAL);
	   P_iter->second.medratio.push_back(NA_VAL);
	   P_iter->second.PPscore.push_back(NA_VAL);
	   P_iter->second.PostOdds.push_back(NA_VAL);
        }
        delete Protein;
     }
     }
     double Ubound = (UserInput.rbound);
     double Lbound = (UserInput.lbound);
     //allRatio = UniqueVec(allRatio);
     sort(allRatio.begin(), allRatio.end());
     vector<double> density = computeDensity(allRatio);
     vector<double> quantile_r = getquantile(allRatio,Lbound,Ubound);
          
     cerr<<"\nUser specified truncation occurs at "<< Lbound*100 <<"th and "<<Ubound*100<<"th percentile."<<endl;
     cerr<<setprecision(3)<<std::fixed;
     cerr<<"This occurs at ratio: ("<<quantile_r.at(0)<<", "<<quantile_r.at(1)<<")\n"<<endl;

     vector<double> Null_r = getNull(allRatio, quantile_r.at(1), quantile_r.at(0));
     //Null_r = UniqueVec(Null_r);
     sort(Null_r.begin(), Null_r.end());

     //bool allsame = true;
     //double first = Null_r.at(0);
     //for(int h=1; h<Null_r.size();h++) if(Null_r.at(h)!=first) allsame=false;
     //if(allsame) cerr<<"Null_r vector all same!"<<endl;

     double MODE = getModeDensity(density,allRatio);
     cerr<<"The estimates for the Null components: "<<endl;
     cerr<<"Mode occurs at ratio : "<<MODE<<endl;
    
     double a = getMin(Null_r) - 1e-6;
     double b = getMax(Null_r) + 1e-6;

     cerr<<"a: " <<a<<"\tb:"<<b<<endl;

     double sdNull = sdTruncatedNormal(Null_r,MODE,a , b);
     cerr<<"sdNull is: "<<sdNull<<endl;

     vector<double> denNullnorm = NORMpdf2(Null_r,MODE, sdNull);
     //for(int k=0; k<denNullnorm.size();k++) cerr<<"fvec_G: " <<denNullnorm.at(k)<<"\t";
     //cerr<<endl;
     vector<double> denNullNP;
     bool found;
     for(int i=0; i<Null_r.size(); i++){
        found = false;		
	for(int j=0; j<allRatio.size(); j++){ 
	    if(!found){     
		if(allRatio.at(j)==Null_r.at(i)) {
			denNullNP.push_back(density.at(j));
			found=true;
	        }
	    }
	    if(found) continue;
        }
     }
     // for(int k=0; k<5; k++) cerr<<"fvec_NP: "<<denNullNP.at(k)<<endl;
     double pi_null;
     
     for(int k=0; k<denNullnorm.size(); k++){
	if(Null_r.at(k)==MODE) {
	pi_null = denNullNP.at(k)/denNullnorm.at(k);
	cerr<<"Density Null at Mode: "<<denNullNP.at(k)<<endl;
     	cerr<<"Density of guassian Null at Mode: "<< denNullnorm.at(k)<<endl;
        break;
	}
     }

     cerr<<"pi_null: "<<pi_null<<endl;
     if(pi_null>=1) pi_null=0.99;
     vector<double> denNULL_G = NORMpdf2(allRatio, MODE, sdNull);
     vector<double> denALT;
     for(int i=0; i<denNULL_G.size(); i++){
         double tmp = (density.at(i) - (pi_null*denNULL_G.at(i)))/(1-pi_null);
         if(tmp<0) tmp=0;
         denALT.push_back(tmp);
     }
     TruncateVal(denNULL_G);
     TruncateVal(denALT);

     vector<double> denUp, denDown;
     for(int h=0; h<denALT.size(); h++){
 	double tmp = denALT.at(h);
        if(allRatio.at(h) < MODE) {
	      denUp.push_back(1e-100);
 	      denDown.push_back(tmp);
        }
	if(allRatio.at(h) >=MODE){
    	      denUp.push_back(tmp);
 	      denDown.push_back(1e-100);
   	}
     }


     //Calculating area under curve to estimate pi_up and pi_down
     double Area_up = TrapeziumRule(allRatio,denUp);
     double Area_down = TrapeziumRule(allRatio,denDown);
     double pi_up = (1-pi_null)*(Area_up/(Area_up+Area_down));
     double pi_down = (1-pi_null)*(Area_down/(Area_up+Area_down));
     
     cerr<<"Area Up: "<<Area_up<<endl;
     cerr<<"Area Down: "<<Area_down<<endl;

     cerr<<"\nPeptide-level proportions are estimated as follows: "<<endl;
     cerr<<"pi_neg: "<<pi_down<<endl;
     cerr<<"pi_null: "<<pi_null<<endl;
     cerr<<"pi_pos: "<<pi_up<<"\n"<<endl;


     for(int j=0; j<allRatio.size(); j++){
     	denUp.at(j) = denUp.at(j)/Area_up;
     	denDown.at(j) = denDown.at(j)/Area_down;
     }


     //storing density information in a map
     //vector will be in the order of density, null density, up density , down density, prop
     map<double, vector<double> >::iterator d_iter;
     for(int i=0; i<allRatio.size(); i++){
	 double curr_r = allRatio.at(i);
         d_iter = densityMap.find(curr_r);
         if(d_iter==densityMap.end()) {
		double prop = (pi_null*denNULL_G.at(i))/((pi_down*denDown.at(i))+(pi_up*denUp.at(i)));
		if(curr_r==MODE) prop = 2*prop;
                vector<double> den;
		den.push_back(density.at(i));
                den.push_back(denNULL_G.at(i)*pi_null);
		den.push_back(denUp.at(i)*pi_up);
		den.push_back(denDown.at(i)*pi_down);	
		den.push_back(prop);
	        densityMap[curr_r] = den;
	 }
     }

     vector<double> PI_pep;
     PI_pep.push_back(pi_down);
     PI_pep.push_back(pi_null);	
     PI_pep.push_back(pi_up); 
     vector<double> PI_prot =EstProtPiEM(PI_pep, PROTmap, PEPmap, densityMap); 
     string label = UserInput.labels.at(pos-1);
     if(UserInput.ExpDesign=="timecourse") label = "timecourse";
     outputDensity(label,allRatio,densityMap);
     //cerr<<"finish output Density file"<<endl;
     computePPscore(PI_prot,PROTmap,PEPmap,densityMap, PROTEIN);
     //cerr<<"finish computePPscore"<<endl;
     return PI_prot;
}


void computePPscore(vector<double> pi, map<string, vector<string> >&PROTmap, map<string, double>&PEPmap, map<double, vector<double> >&DENmap, map<string, Protein_t>&PROTEIN){

     double pi_neg = pi.at(0);
     double pi_null = pi.at(1);
     double pi_pos = pi.at(2);

     map<string, vector<string> >::iterator prot_iter;
     map<string, double>::iterator pep_iter;
     map<double, vector<double> >::iterator d_iter,d_iter2;
     map<string, Protein_t>::iterator P_iter;

     string curr_prot, curr_pep;
     vector<string> peptides;
     double aP,bP,cP,mP,aN,bN,cN,mN, postodd, ppscore;
     int sumlen=0;
     int counter=0;

     double MODE;
     vector<double> fullden;
     vector<double> allr;
     vector<double> upDen;
     vector<double> downDen;
     //cerr<<"DENmap size: "<<DENmap.size()<<endl;
     for(d_iter= DENmap.begin(); d_iter!=DENmap.end(); d_iter++){
	allr.push_back(d_iter->first);
	fullden.push_back(d_iter->second.at(0));
	upDen.push_back(d_iter->second.at(2));
	downDen.push_back(d_iter->second.at(3));
        //cerr<<d_iter->second.at(3)<<endl;
     }
     MODE = getModeDensity(fullden,allr);
     //cerr<<"size of upDen: "<<upDen.size()<<endl;
     double UpDiesAt = DiesAt(upDen,allr,1);
     //cerr<<"size of downDen: "<<downDen.size()<<endl;
     double DownDiesAt = DiesAt(downDen, allr,2);
     //cerr<<"or here????"<<endl;
     cerr<<"Positive density dies at " <<UpDiesAt<<" and negative dies at "<<DownDiesAt<<"\n"<<endl;
    
     
     for(prot_iter=PROTmap.begin(); prot_iter!=PROTmap.end(); prot_iter++){
	 curr_prot = prot_iter->first;
         peptides = prot_iter->second;
         P_iter = PROTEIN.find(curr_prot);
	 sumlen+=peptides.size();
	 int NumP=0,NumN=0;
         double pep_r, pep_r2,Ppep_UP=0,Ppep_DOWN=0,Ppep_NULL=0,Npep_UP=0,Npep_DOWN=0,Npep_NULL=0;
         double minup_val, mindown_val, minnull_val, tmp_up, tmp_down, tmp_null;
	 vector<double> tmp_UP,tmp_DOWN, tmp_NULL;
         for(int i=0; i<peptides.size(); i++){
	      pep_iter = PEPmap.find(peptides.at(i));
	      pep_r2 = pep_iter->second;
              d_iter2= DENmap.find(pep_r2);
	      tmp_UP.push_back(d_iter2->second.at(2));
	      tmp_DOWN.push_back(d_iter2->second.at(3));
	      tmp_NULL.push_back(d_iter2->second.at(1));
	 }

         vector<double> vec_up;
	 vector<double> vec_down;
	 vector<double> vec_null;
	 if(peptides.size()>1){
	    for(int h=0;h<peptides.size();h++){
	      if(tmp_UP.at(h)!=getMin(tmp_UP)) vec_up.push_back(tmp_UP.at(h));
	      if(tmp_DOWN.at(h)!=getMin(tmp_DOWN)) vec_down.push_back(tmp_DOWN.at(h));
	      if(tmp_NULL.at(h)!=getMin(tmp_NULL)) vec_null.push_back(tmp_NULL.at(h));
	    }
	    if(vec_up.size()>0)  minup_val = 0.2*getMin(vec_up);
	    if(vec_down.size()>0)  mindown_val = 0.2*getMin(vec_down);
	    if(vec_null.size()>0)  minnull_val = 0.2*getMin(vec_null); 

	    if(vec_up.size()==0)  minup_val = getMin(tmp_UP);
	    if(vec_down.size()==0)  mindown_val = getMin(tmp_DOWN);
	    if(vec_null.size()==0)  minnull_val = getMin(tmp_NULL);

	 }
	 if(peptides.size()==1){
            minup_val = getMin(tmp_UP);
	    mindown_val = getMin(tmp_DOWN);
	    minnull_val = getMin(tmp_NULL); 
	 }
         if(minup_val<0.001) minup_val=0.001;
         if(mindown_val<0.001) mindown_val =0.001;
	 if(minnull_val <0.001) minnull_val = 0.001;

         tmp_UP.clear(); tmp_DOWN.clear(); tmp_NULL.clear();
	 vec_up.clear(); vec_down.clear(); vec_null.clear();

         for(int j=0; j<peptides.size(); j++){
	        pep_iter = PEPmap.find(peptides.at(j));
	        pep_r = pep_iter->second;
		d_iter = DENmap.find(pep_r);
		tmp_up = d_iter->second.at(2);
	        tmp_down = d_iter->second.at(3);
	        tmp_null = d_iter->second.at(1);
		if(pep_r>=MODE){
		   if(pep_r>UpDiesAt){
	              Ppep_UP += log(tmp_up);
		      Ppep_DOWN += log(tmp_down);
	 	      Ppep_NULL += log(tmp_null);
		   }
		   if(pep_r <= UpDiesAt){
		      Ppep_UP += log(minup_val);
		      Ppep_DOWN += log(mindown_val);
		      Ppep_NULL += log(minnull_val); 
		   }
		   NumP++;
		}
		if(pep_r<MODE){
		   if(pep_r < DownDiesAt){
	              Npep_UP += log(tmp_up);
		      Npep_DOWN += log(tmp_down);
	 	      Npep_NULL += log(tmp_null);
		   }
		   if(pep_r>= DownDiesAt){
		      Npep_UP += log(minup_val);
		      Npep_DOWN += log(mindown_val);
		      Npep_NULL += log(minnull_val);
		   }
		   NumN++;
		}
	        //if(curr_prot=="Q9Y623"|curr_prot=="P55268"){
		//   cerr<<"Ppep_UP: "<<Ppep_UP<<"\tPpep_DOWN: "<<Ppep_DOWN<<"\tPpep_NULL: "<<Ppep_NULL<<endl;
		//   cerr<<"Npep_UP: "<<Npep_UP<<"\tNpep_DOWN: "<<Npep_DOWN<<"\tNpep_NULL: "<<Npep_NULL<<endl;
	       //}
         }
	 aP = Ppep_UP;
	 bP = Ppep_NULL;
	 cP = Ppep_DOWN;
	 mP = max(aP,max(bP,cP));
	 double Pcomp_UP = pi_pos*exp(aP-mP);
	 double Pcomp_NULL = pi_null*exp(bP-mP);
	 double Pcomp_DOWN = pi_neg*exp(cP-mP);

	 aN = Npep_UP;
	 bN = Npep_NULL;
	 cN = Npep_DOWN;
	 mN = max(aN,max(bN,cN));
	 double Ncomp_UP = pi_pos*exp(aN-mN);
	 double Ncomp_NULL = pi_null*exp(bN-mN);
	 double Ncomp_DOWN = pi_neg*exp(cN-mN);


	 double Psum_comp = Pcomp_UP + Pcomp_NULL;
	 double Nsum_comp = Ncomp_NULL + Ncomp_DOWN;
	 double ZupP = Pcomp_UP/(Pcomp_UP+Pcomp_NULL);
	 double ZnullP = Pcomp_NULL/(Pcomp_UP +Pcomp_NULL+Pcomp_DOWN);
	 double ZdownP = Pcomp_DOWN/(Pcomp_DOWN+Pcomp_NULL);

	 double ZupN = Ncomp_UP/(Ncomp_UP+Ncomp_NULL);
	 double ZnullN = Ncomp_NULL/(Ncomp_UP+Ncomp_NULL+Ncomp_DOWN);
	 double ZdownN = Ncomp_DOWN/(Ncomp_DOWN+Ncomp_NULL);

         double null_LLP = log(pi_null) +(bP-mP);
	 double Up_LLP = log(pi_pos) + (aP-mP);
	 double Down_LLP = log(pi_neg) + (cP-mP);

         double null_LLN = log(pi_null) +(bN-mN);
	 double Up_LLN = log(pi_pos) + (aN-mN);
	 double Down_LLN = log(pi_neg) + (cN-mN);

	 double Zup = ((NumP*ZupP)+(NumN*ZupN))/(NumP+NumN);
	 double Zdown = ((NumP*ZdownP)+(NumN*ZdownN))/(NumP+NumN);
         if (Zup!=Zup) Zup = 0;
	 if(Zdown!=Zdown) Zdown =0;

         ppscore = Zup-Zdown;
 	 //if(ppscore!=ppscore){
	//	cerr<<"protein: "<<curr_prot<<endl;
	//	cerr<<"Zup: "<<Zup<<"\tZdown: "<<Zdown<<endl;
	//	cerr<<"Number postive, negative: "<< NumP<<"\t"<<NumN<<endl;
	//	cerr<<"Positive Zup Zdown: "<<ZupP<<"\t"<<ZdownP<<endl;
	//	cerr<<"Negative Zup Zdown: "<<ZupN<<"\t"<<ZdownN<<endl;
        // }

	 P_iter->second.PPscore.push_back(ppscore);

	 double Up_LLp = (Up_LLP-null_LLP);
	 double Up_LLn = (Up_LLN-null_LLN);
	 //double Up_LL = ((NumP*Up_LLp)+(NumN*Up_LLn))/(NumP+NumN);
	 double Down_LLp = (Down_LLP-null_LLP);
	 double Down_LLn = (Down_LLN-null_LLN);
	 //double Down_LL = ((NumP*Down_LLp)+(NumN*Down_LLn))/(NumP+NumN);

         postodd = max(Up_LLp,Down_LLn);
	 P_iter->second.PostOdds.push_back(postodd);

         // if(counter<10){
	   //  cerr<<"MODE: "<<MODE<<"\tNumPos: "<<NumP<<"\tNumNeg: "<<NumN<<endl;
	     //cerr<<"ZupP: "<<ZupP<<"\tZdownP: "<<ZdownP<<"\tZupN: "<<ZupN<<"\tZdownN: "<<ZdownN<<endl;
	     //cerr<<"Zup: "<<Zup<<"\tZdown: "<<Zdown<<endl;
	     //cerr<<"Up_LL: "<<Up_LL<<"\tDown_LL: "<<Down_LL<<"\tNull_LL: "<< null_LL<<endl;
	     //cerr<<"PPscore: "<<ppscore<<"\tPostOdds: "<<postodd<<endl;
	//}
        counter++;
     }	  

}


double TrapeziumRule(vector<double> x, vector<double> y){

     double area=0;
     for(int i=0; i<(x.size()-1); i++){
	double len = x.at(i+1)-x.at(i);
        double ht = (y.at(i+1)+y.at(i))/2;
        area+=(ht*len);
     }
     return area;
}

void TruncateVal(vector<double> &vec){
   for(int i=0; i<vec.size(); i++){
        if(vec.at(i)< 1e-100) vec.at(i) = 1e-100;
   }
}


vector<double> NORMpdf2(vector<double> v, double m, double ss){
   vector<double> ret,vec_n;
   //centering vec before finding pdf
   for(int i=0; i<v.size(); i++){
     double tmp = (v.at(i)-m);
     double tmp2 = gsl_ran_gaussian_pdf(tmp,ss);
     ret.push_back(tmp2);
   }
   return ret;
}


//function to find the variance which maximizes the log-Lik of Truncated normal
double sdTruncatedNormal(vector<double> v, double mu, double a, double b){
   double ret;
   vector<double> grid,LogLik;
   double tmp = 0;
   for(int i=0; i<50000; i++){
      tmp += 0.0001;
      grid.push_back(tmp);
   }
   for(int j=0; j<grid.size(); j++){
      double ll = getSum(logNormTrun(v,mu,grid.at(j),a,b));
      LogLik.push_back(ll);
   }
   
   double maxLL = getMax(LogLik);
   vector<int> indices;
   for(int k=0; k<LogLik.size();k++) if(LogLik.at(k)==maxLL) indices.push_back(k);
   vector<double> var_vec;
   for(int i=0; i<indices.size(); i++) var_vec.push_back(grid.at(indices.at(i)));
   ret = CalMeanMed(var_vec,1);

   return ret;
}

//log-Lik for Truncated normal
vector<double> logNormTrun(vector<double> vec, double mu, double grid, double a, double b){
   vector<double> ret;
   int count=0;

   double x1 = log(pdfStdNorm(b,mu,grid) - pdfStdNorm(a,mu,grid))/log(2);

   for(int i=0; i<vec.size(); i++){
      double val;
      double v = grid*grid;
      double x2 = - (pow((vec.at(i)-mu),2)/(2*v))/log(2) -0.5*(log(2*pi*v)/log(2));
      val = x2 - x1;
      ret.push_back(val);
   }
   return ret;
}


//log-Lik for Normal dist
vector<double> logNormal(vector<double> vec, double mu, double v){

    vector<double> ret;
    for(int i=0; i<vec.size(); i++){
        double t = - (pow((vec.at(i)-mu), 2)/(2*v))/log(2) - 0.5*(log(2*pi*v)/log(2));
        ret.push_back(t);
    }
    return ret;
}


double getModeDensity(vector<double> den, vector<double> r){

   double mode;
   vector<int> index;
   double maxden = getMax(den);
   for(int j=0; j<den.size(); j++) if(den.at(j)==maxden) index.push_back(j);

   vector<double> mode_vec;
   for(int i=0; i<index.size(); i++) mode_vec.push_back(r.at(index.at(i)));
   mode = CalMeanMed(mode_vec,1);
   return mode;
}

double DiesAt(vector<double> den, vector<double> r, int type){
   double ret;
   vector<double> rr;
   double peak = getModeDensity(den,r);
   //cerr<<"peak"<< peak<<endl;
   double min = getMin(den);
   //cerr<<"min: "<<min<<endl;
   for(int j=0; j<den.size(); j++) {
	//for positive -1 for negative -2
	if(type==1){if(den.at(j)==min & den.at(j)<=peak) rr.push_back(r.at(j));}
	if(type==2){if(den.at(j)==min & den.at(j)>=peak) rr.push_back(r.at(j));}
   }
    //cerr<<rr.size()<<endl;
    if(rr.size()==0) ret =peak;
    else{
       if(type==1) ret = getMax(rr);
       if(type==2) ret = getMin(rr);
    }	 
        
   return ret;
}   



vector<double> UniqueVec(vector<double> vec){

   sort(vec.begin(), vec.end());
   vec.erase(unique(vec.begin(), vec.end()), vec.end());
   return vec;
}


double getMin(vector<double> vec){
   double min=vec.at(0);
   if(vec.size()==1) min = vec.at(0);
   else{
     for(int i=1; i<vec.size(); i++){
	double min_n = vec.at(i);
	if(min_n<min) min = min_n;
     }
   }
   return min;
}

double getMax(vector<double> vec){
   double max=vec.at(0);
   for(int i=1; i<vec.size(); i++){
	double max_new = vec.at(i);
	if(max_new>max) max = max_new;
   }
   return max;
}

double getMode(vector<double> ratio, int Nbins){

     double mode = 0;
     double j, max_j, max_id, a, b;
     double binWidth = 0.0001;
     vector<double> v;
     vector<double>::iterator it;
     v.resize(Nbins, 0); // initialize vector to have Nbins, all at zero
     for(it = ratio.begin(); it != ratio.end(); it++) {
         for(j = 0; j < (Nbins - 1); j++) {
	       a =  (j * binWidth);
	       b = a + binWidth;
               if( (*it >= a) & (*it < b) ) {
			v[j]++;
                        break; // done, get out of this iteration
	       }
	  }
      }
      // Now find the largest value in 'v'
      max_j = getMax(v);
      //for(int i = 0; i < Nbins; i++) if( v[i] > max_j ) max_j = v[i];
      // now find the element of 'v' that holds the max value
      max_id = 0;
      for(int i = 0; i < Nbins; i++) if( v[i] == max_j ) max_id = i;
      mode = (max_id*binWidth) + (binWidth/2);
      return mode; 
}


vector<double> getquantile(vector<double> vec, double lb, double rb){
     vector<double> ret;
     sort(vec.begin(), vec.end());
     double q1 = vec.at(round(lb*vec.size())-1);
     double q2 = vec.at(round(rb*vec.size())-1);

     ret.push_back(q1);
     ret.push_back(q2);

     return ret;
}


vector<double> getNull(vector<double> vec, double UpperB, double LowerB){

     vector<double> ret;
     for(int i=0; i<vec.size(); i++) {
	if(vec.at(i)<=UpperB & vec.at(i)>=LowerB) ret.push_back(vec.at(i));
     } 
     return ret;
}
	


vector<double> computeDensity(vector<double> x){

    vector<double> den;
    double sigma = CalSD(x);
    int size = x.size();
    double h=pow((4*pow(sigma,5))/(3*size),0.2);

    for(int i=0; i<size; i++){
        vector<double> tmp_den;
	for(int j=0; j<size; j++){
	    double d = (x.at(j)-x.at(i))/h;
	    tmp_den.push_back(d);
	}
	double d_n = getSum(NORMpdf(tmp_den))/(size*h);
        den.push_back(d_n);
    }
    return den;
}


double getSum(vector<double> vec){
    double sum=0;
    for(int i=0; i<vec.size(); i++) sum+=vec.at(i);
    return sum;
}

vector<double> NORMpdf(vector<double> vec){

     vector<double> pdf;
     for(int i=0; i<vec.size(); i++){
        double p = gsl_ran_ugaussian_pdf(vec.at(i));
	pdf.push_back(p);
     }
     return pdf;
}


double pdfStdNorm(double x, double m, double ss){
    double q;
    q = (x-m)/ss;
    double ret = gsl_cdf_ugaussian_P(q);
    return ret;
}


void outputDensity(string lab,vector<double> ratio, map<double, vector<double> >&DENmap){

    string filenm;
    filenm = "densityfile_"+lab+".txt";
    ofstream outF(filenm);
    double r, density,den_null, den_up, den_down, prop;
    map<double, vector<double> >::iterator iter;
    vector<double> vec;
    ratio = UniqueVec(ratio);

    outF<<"Ratio\tDensity\tNull_density\tPos_density\tNeg_density\tPropNullvsAlt"<<endl;
    for(int i=0; i<ratio.size(); i++){
	r = ratio.at(i);
        iter = DENmap.find(r);
	vec = iter->second;
        density = vec.at(0);
        den_null = vec.at(1);
        den_up = vec.at(2);
        den_down = vec.at(3);
        prop = vec.at(4);
	outF<<r<<"\t"<<density<<"\t"<<den_null<<"\t"<<den_up<<"\t"<<den_down<<"\t"<<prop<<endl;
    }
    outF.close();

}


vector<double> EstProtPiEM(vector<double> pi,map<string, vector<string> >&PROTmap, map<string, double>&PEPmap, map<double, vector<double> >&DENmap){

        map<string, vector<string> >::iterator prot_iter;
 	map<string, double>::iterator pep_iter;
        map<double, vector<double> >::iterator d_iter;

        double tol = 100;
        double old_pneg,old_ppos, old_pnull, pneg, ppos,pnull;
	string curr_prot, curr_pep;
	vector<string> pep_vec;
        double r,a,b,c;
        vector<double> ret;
        int counter=0;

        //initialization
        pneg = pi.at(0);
	pnull = pi.at(1);
	ppos = pi.at(2);

        while(tol> 1e-5){
	     old_pneg = pneg;
	     old_pnull = pnull;
	     old_ppos = ppos;
	     int sumlen=0;
             double Zdown=0, Znull=0, Zup=0;

	     // E-Step //   
	     for(prot_iter=PROTmap.begin(); prot_iter!=PROTmap.end(); prot_iter++){
		 curr_prot = prot_iter->first;
		 pep_vec = prot_iter->second;
		 sumlen+= pep_vec.size();
		 vector<double> pep_ratio;
		 double pep_UP=0, pep_DOWN=0, pep_NULL=0;
		 for(int j=0; j<pep_vec.size(); j++){
			curr_pep = pep_vec.at(j);
		 	pep_iter = PEPmap.find(curr_pep);
			r = pep_iter->second;
			pep_ratio.push_back(r);
			d_iter = DENmap.find(r);
			pep_UP += log(d_iter->second.at(2));
			pep_DOWN += log(d_iter->second.at(3));
			pep_NULL += log(d_iter->second.at(1));
	         }
		 a = pep_UP;
		 b = pep_NULL;
		 c = pep_DOWN;
		 double m = max(a,max(b,c));
		 double comp_UP = ppos*exp(a-m);
		 double comp_DOWN = pneg*exp(c-m);
		 double comp_NULL = pnull*exp(b-m);
		 Zdown += (comp_DOWN/(comp_UP+comp_DOWN+comp_NULL))*pep_vec.size();
		 Znull += (comp_NULL/(comp_UP+comp_DOWN+comp_NULL))*pep_vec.size();
		 Zup += (comp_UP/(comp_UP+comp_DOWN+comp_NULL))*pep_vec.size();
	    }
	    // M-step //
	    pneg = Zdown/sumlen;
	    pnull = Znull/sumlen;
	    ppos = Zup/sumlen;
	    if(pneg<0.01) pneg=0.01;
	    if(pnull<0.01) pnull=0.01;
	    if(ppos <0.01) ppos =0.01;
             
            double denom = pneg+pnull+ppos;
	    pneg = pneg/denom;
	    pnull = pnull/denom;
	    ppos = ppos/denom;
	    
	    counter++;	
  	    tol = abs(max(pneg-old_pneg, max(pnull-old_pnull, ppos-old_ppos)));
	    if (counter%100==0) cerr<<"Iteration "<< counter<<"...\n";		
	}
	cerr<<"\nProtein-level proportions pi are estimated as follows:"<<endl;
	cerr<<"pi_neg: "<<pneg<<endl;
	cerr<<"pi_null: "<<pnull<<endl;
	cerr<<"pi_pos: "<<ppos<<endl;
        ret.push_back(pneg);
	ret.push_back(pnull);
	ret.push_back(ppos);
	return ret;
}

void computeBFDR(double bfdr_thres,int index, map<string, Protein_t>&PROTEIN){

     map<string, Protein_t>::iterator iter;
     Protein_t Prot;
     vector<double> scorepos,scoreneg;
     vector<double> ppscore;
     double curr_score;
     vector<string> protein,na_prot;
     int numsig=0;
     int counter1=0, counter2=0;
     
     for(iter=PROTEIN.begin(); iter!=PROTEIN.end(); iter++){
        double score = iter->second.PPscore.at(index);
	if(score==NA_VAL) na_prot.push_back(iter->first);
	if(score!=NA_VAL){
	    ppscore.push_back(score);
	    protein.push_back(iter->first);
	}
     }
     for(int i=0; i<ppscore.size(); i++){
	curr_score = ppscore.at(i);
	for(int j=0; j<ppscore.size(); j++){
	    if(ppscore.at(j)>=curr_score) scorepos.push_back(1-ppscore.at(j));
	    if(ppscore.at(j)<=curr_score) scoreneg.push_back(1+ppscore.at(j));
        }
	iter = PROTEIN.find(protein.at(i));
	if(curr_score>0) iter->second.bfdr.push_back(CalMeanMed(scorepos,1));
        if(curr_score<0) iter->second.bfdr.push_back(CalMeanMed(scoreneg,1));
	if(curr_score==0) iter->second.bfdr.push_back(min(CalMeanMed(scorepos,1),CalMeanMed(scoreneg,1)));
	if(iter->second.bfdr.at(index)<bfdr_thres) numsig++;
	scorepos.clear();
	scoreneg.clear();
     }
     //fill in the rest of the map with NA for BFDR
     for(int k=0; k<na_prot.size(); k++){
	iter = PROTEIN.find(na_prot.at(k));
	iter->second.bfdr.push_back(NA_VAL);
     }

     cerr<<"\nThere were "<<numsig<<" number of significant proteins after controlling for "<<setprecision(0)<<bfdr_thres*100<<"% BFDR."<<endl;
}


void outputResult(Input_t *UserInput , map<string, Protein_t>&PROTEIN){

     ofstream outF("EBprot_results.txt");
     vector<string> labels;
     int count =0;
     string design = UserInput->ExpDesign;
     double thres = UserInput->threshold;
     vector<string> samplelab= UserInput->labels; 
     for(int i=0; i<samplelab.size(); i++){
	labels.push_back("MedianlogRatio_"+samplelab.at(i));
	labels.push_back("NumPep_"+ samplelab.at(i));
	labels.push_back("NumPepRm_"+samplelab.at(i));
	labels.push_back("PPscore_"+samplelab.at(i));
	labels.push_back("PostOdds_"+samplelab.at(i));
	labels.push_back("BFDR_"+samplelab.at(i));
     }
     if(design =="replicate"){
	labels.push_back("NumReplicates");
        labels.push_back("NumSigReplicates");
	labels.push_back("AveragePPscore");
	labels.push_back("Min_PPscore");
	labels.push_back("Max_PPscore");
	labels.push_back("NumRep_positive");
	labels.push_back("NumRep_negative");
     }

     outF <<"Protein";
     for(int i=0; i<labels.size(); i++){
	outF<<"\t"<<labels.at(i);
     }
     outF<<endl;

   
     map<string, Protein_t>::iterator iter;
     string curr_prot;
     vector<double> ratio;
     vector<double> score;
     vector<double> bfdr;

     for (iter=PROTEIN.begin(); iter!=PROTEIN.end(); iter++){
	curr_prot = iter->first;
        ratio = iter->second.medratio;
	score = iter->second.PPscore;
	bfdr = iter->second.bfdr;
        outF <<curr_prot;

	double tmp;
	int tmp2;

        for(int k=0; k<samplelab.size();k++){
        tmp = ratio.at(k);
        if(tmp == NA_VAL) outF<<"\t"<<"";
	else outF<<"\t"<<tmp;
	tmp2 = iter->second.numPep.at(k);
	if(tmp2 == NA_VAL) outF<<"\t"<<"";
	else outF<<"\t"<<tmp2;
 	tmp2 = iter->second.numPepRM.at(k);
	if(tmp2 ==NA_VAL) outF<<"\t"<<"";
	else outF<<"\t"<<tmp2;
	tmp = score.at(k);
	if(tmp==NA_VAL) outF<<"\t"<<"";
	else outF<<"\t"<<tmp;
	tmp = iter->second.PostOdds.at(k);
	if(tmp==NA_VAL) outF<<"\t"<<"";
	else outF<<"\t"<<tmp;
	tmp = bfdr.at(k);
	if(tmp==NA_VAL) outF<<"\t"<<"";
	else outF<<"\t"<<tmp;
	}

        if(design=="replicate"){
	    rmMissing(ratio);
	    rmMissing(score);
	    rmMissing(bfdr);
	    outF<<"\t"<<ratio.size();
	    int numsig=0;
	    for(int k=0; k<score.size(); k++) if(bfdr.at(k)< thres) numsig++;	    
	    outF<<"\t"<<numsig;
	    outF<<"\t"<<CalMeanMed(score,1)<<"\t"<<getMin(score)<<"\t"<<getMax(score);
            int numpos=0;
	    int numneg=0;
	    for(int k=0; k<ratio.size(); k++){
		if(ratio.at(k)>=0)numpos++;
	        else numneg++;
	    }
	    outF<<"\t"<<numpos<<"\t"<<numneg;
        }
	outF<<endl;	
        ratio.clear();
	score.clear();
	bfdr.clear();
      }
      outF.close();
} 

vector<string> NotIn2(set<string> set1, set<string> set2){
    vector<string> v;
    set<string> SET2;
    set<string>::iterator it;
    for(it = set2.begin(); it!=set2.end(); it++){
	string tmp = stringSplit(*it, '_').at(0);
        SET2.insert(tmp);
    }
    set_difference(set1.begin(),set1.end(), SET2.begin(), SET2.end(), back_inserter(v));
    
    return v;
}

vector<string> NotIn(set<string> set1, set<string> set2){
    vector<string> v; 
    set_difference (set1.begin(),set1.end(), set2.begin(), set2.end(), back_inserter(v));
    return v;
}  

void cleanMap(map<string, Protein_t> &MAP){

   map<string, Protein_t>::iterator iter;
   vector<string> Protrm;
   vector<double> vec;
   for(iter=MAP.begin(); iter!= MAP.end(); iter++){
         vec = iter->second.PPscore;
	 if(allMissing(vec)) Protrm.push_back(iter->first);
	vec.clear();
   }

   if(Protrm.size()>0){
	for(int j=0; j<Protrm.size(); j++){
	    iter = MAP.find(Protrm.at(j));
	    MAP.erase(iter);
	 }	 
   }
}


void UnmergeMap(int len, map<string, vector<string> > &PROTmap, map<string, Protein_t> &PROTEINfit , map<string, Protein_t> &PROTEIN){

      map<string, vector<string> >::iterator p_iter;
      map<string, Protein_t>::iterator iter1, iter2;
      string tmp_p, new_prot, tmp_pep, new_pep;
      vector<string> peptides;
      double score, bfdr,ratio, odds;
      int pepnum, peprm;
      char buffer[33];

      for(p_iter= PROTmap.begin(); p_iter!=PROTmap.end(); p_iter++){
	tmp_p = p_iter->first;
	peptides = p_iter->second;
	for(int k=0; k<len; k++){
 	    new_prot = tmp_p+ "______R" + to_string(k);
	    iter1 = PROTEINfit.find(new_prot);
            
	    if(iter1!=PROTEINfit.end()){
	        ratio = iter1->second.medratio.at(0);
	        bfdr = iter1->second.bfdr.at(0);
	        score = iter1->second.PPscore.at(0);
	        odds = iter1->second.PostOdds.at(0);
		pepnum = iter1->second.numPep.at(0);
	        peprm = iter1->second.numPepRM.at(0);

		if(k==0){
		    Protein_t *Protein = new Protein_t();
		    PROTEIN[tmp_p] = *Protein;
		    iter2 = PROTEIN.find(tmp_p);
		}
		if(k!=0) {
		    iter2 = PROTEIN.find(tmp_p);
	       	    if(iter2 ==PROTEIN.end())cerr<<"Problem, should not reach the end!"<<endl;
	        }
		iter2->second.medratio.push_back(ratio);
	        iter2->second.PPscore.push_back(score);
		iter2->second.PostOdds.push_back(odds);
		iter2->second.bfdr.push_back(bfdr);
		iter2->second.numPep.push_back(pepnum);
		iter2->second.numPepRM.push_back(peprm);
	    }
	    if(iter1==PROTEINfit.end()){
		if(k==0){
		    Protein_t *Protein = new Protein_t();
		    PROTEIN[tmp_p] = *Protein;
		    iter2 = PROTEIN.find(tmp_p);
		}
		if(k!=0) iter2 = PROTEIN.find(tmp_p);
	
		iter2->second.medratio.push_back(NA_VAL);
		iter2->second.PPscore.push_back(NA_VAL);
		iter2->second.PostOdds.push_back(NA_VAL);
		iter2->second.bfdr.push_back(NA_VAL);
		iter2->second.numPep.push_back(NA_VAL);
	        iter2->second.numPepRM.push_back(NA_VAL);
	   }
      }
     
   }
}
