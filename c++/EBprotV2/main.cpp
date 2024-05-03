
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: Hiromi Koh, Huiyi Tay
 *
 * Created on 9 May, 2016, 11:03 AM
 */

#include "global.hpp"

#include <gsl/gsl_cdf.h> 
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
        cerr<<"\nUSAGE: EBprotNP.exe <data.file> <parameter.file>\n";
        return 0;
    }
    clock_t begin = clock();
    
    string dataF = argv[1];
    string paramF = argv[2];
    
    map <string, vector<string> > PROTmap;
    map <string, vector<double> > PEPmap;
    map <string, Protein_t> PROTEIN;
    Input_t *UserInput = new Input_t();
        
    cerr<<"==================================================================="<<endl;
    cerr<< "Beginning Nonparametric EBprot Analysis...\n\n" <<endl;
    readUserInput(paramF, *UserInput);
  
    cerr<< "Reading in the data file...\n"<<endl;
    bool logtrans = UserInput->log2transform;
    set<string> protList = readData(dataF, UserInput->labels, logtrans, PROTmap, PEPmap);
    cerr<< "\nThere are a total of "<<protList.size()<< " proteins and "<< PEPmap.size()<< " peptides.\n"<<endl; 

    vector<string> labels = UserInput->labels;
    double thres = UserInput->threshold;
    string design = UserInput->ExpDesign;
    
    cerr<<"Carrying out analysis for a "<<design<<" experimental design."<<endl;

    // independent and replicate design //
    if(design =="independent"|design == "replicate"){

	 for(int i=0; i<labels.size(); i++){
 
	    map<string, double> PEPmapNEW;
	    map<string, vector<string> >PROTmapNEW;
            map<double, vector<double> >DENmap;
	    cerr<<"=============================================================="<<endl;
	    cerr<<"Starting model fitting for "<<labels.at(i)<<"...\n"<<endl;
	    vector<int> PepRM;
            int minPep = UserInput->minpep;
 	    set<string> protinRep = createMap(i,minPep,PROTmap,PROTmapNEW, PEPmap,PEPmapNEW);
            
	    if(UserInput->outlierRm) {
	        cerr<<"Filtering out outlying peptides from proteins..."<<endl;
		PepRM = OutlierRm(UserInput->mink, PROTmapNEW,PEPmapNEW);
    	     }
    	     else{ cerr<<"No filtering of outlying peptides will be carried out.\n\n"<<endl;}
   
             set<string> Proteins = getKey(PROTmapNEW);
             vector<string> protNA = NotIn(protList,Proteins);
             vector<double> ALLratio =nonparamFit(protNA,PepRM, PROTmapNEW, PEPmapNEW, *UserInput, PROTEIN, DENmap);

             map<string, Protein_t>::iterator iter;
	     int count=1;
	     cerr<<"\nReporting the results for first 6 proteins:"<<endl;
	     cerr<<"----------------------------------------------------------------"<<endl;
	     for(iter = PROTEIN.begin();iter!=PROTEIN.end(); iter++){
		if(iter->second.medratio.at(i)!=NA_VAL){
		   if(count<=5){
		        cerr<<"Protein: "<< iter->first<<endl;
			cerr<<"Median log2ratio: "<< iter->second.medratio.at(i)<<endl;
			cerr<<"Number of peptide: "<< iter->second.numPep.at(i)<<endl;
			cerr<<"Number of peptide removed: "<< iter->second.numPepRM.at(i)<<endl;;
			cerr<<"PPscore: "<< iter->second.PPscore.at(i)<<endl;
			//cerr<<"Posterior Odds: "<< iter->second.PostOdds.at(i)<<endl;
			cerr<<"----------------------------------------------------------------"<<endl;
		    }
		    count++;
		}
	     }
	     computeBFDR(thres,i, PROTEIN);
	 }
	 
    }

    if(design =="timecourse"){
          map<string, double> PEPmapNEW;
	  map<string, vector<string> >PROTmapNEW;
	  map<double, vector<double> >DENmap;
	  map<string, Protein_t> PROTEINfit;
          cerr<<"\nStarting model fitting procedure...\n"<<endl;
          int minPep = UserInput->minpep;
          set<string> protinRep = createMap2(labels.size(),minPep,PROTmap, PROTmapNEW,PEPmap, PEPmapNEW);
	  vector<int> PepRM;
	  //cerr<<"PROTmap size: "<<PROTmap.size()<<endl;
	  //cerr<<"PEPmap size: "<<PEPmap.size()<<endl;
	  //cerr<<"PROTmapNEW size: "<<PROTmapNEW.size()<<endl;
          //cerr<<"PEPmapNEW size: "<<PEPmapNEW.size()<<endl;
          if(UserInput->outlierRm){
	     cerr<<"\nFiltering out outlying peptides from proteins..." <<endl;
	     PepRM = OutlierRm(UserInput->mink, PROTmapNEW, PEPmapNEW);
	  }
	  else{ cerr<<"No filtering of outlying peptides will be carried out.\n\n"<<endl;}
          set<string> Proteins = getKey(PROTmapNEW);
          vector<string> protNA = NotIn(protinRep,Proteins);
          vector<double> ALLratio =nonparamFit(protNA,PepRM, PROTmapNEW, PEPmapNEW, *UserInput, PROTEINfit, DENmap);

          map<string, Protein_t>::iterator iter;
	  int count=1;
	  cerr<<"\nReporting the results for first 6 proteins:"<<endl;
	  cerr<<"----------------------------------------------------------------"<<endl;
	  for(iter = PROTEINfit.begin();iter!=PROTEINfit.end(); iter++){
		if(iter->second.medratio.at(0)!=NA_VAL){
		   if(count<=5){
		        cerr<<"Protein: "<< iter->first<<endl;
			cerr<<"Median log2ratio: "<< iter->second.medratio.at(0)<<endl;
			cerr<<"Number of peptide: "<< iter->second.numPep.at(0)<<endl;
			cerr<<"Number of peptide removed: "<< iter->second.numPepRM.at(0)<<endl;;
			cerr<<"PPscore: "<< iter->second.PPscore.at(0)<<endl;
			cerr<<"Posterior Odds: "<< iter->second.PostOdds.at(0)<<endl;
			cerr<<"----------------------------------------------------------------"<<endl;
		    }
		    count++;
		}
	     }
	     computeBFDR(thres,0, PROTEINfit);
	     UnmergeMap(labels.size(), PROTmap,PROTEINfit, PROTEIN);
    }
    cleanMap(PROTEIN);
    outputResult(UserInput, PROTEIN);
   
    clock_t end = clock();
    double elapse_sec = double(end - begin)/(double)CLOCKS_PER_SEC;
    cerr<< "\n\nEBprotNP has completed running in " << elapse_sec <<" seconds.\n" <<endl;
    cerr<< "================================= Completed ===================================\n"<<endl;
    
    return 0;
    
}



