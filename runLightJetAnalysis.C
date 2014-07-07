#include "LightJetAnalysis.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "TString.h"

#if boostflag==1
#include "boost/program_options.hpp"
namespace po = boost::program_options;
#endif

using std::cout;
using std::endl;
using std::string;
using std::map;
using namespace std;

#define OUTPUT_NAME "out.root"
#define N_ANALYZED_EVENTS 1000 //this value is only necessary in the commented code
#define DATA_FILE_NAME "/atlas/output/pnef/skimmed.20140703.15.27_ClustersAndTruth.PythJ1to3mc12aJETMET.jetmet2012.root"

static inline void init(char* data_file_name, int argc, char* argv);

static inline void init(char* data_file_name, int argc, char* argv) 
{
      // argument parsing  ------------------------
    cout << "Called as: ";
    for(int i = 0; i < argc; i++) {
        cout << argv[i] << " ";
    }
    cout << endl;

   /*
    #if boostflag == 1 // command line parsing if boost is installed 
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("Debug",          po::value<int>(&debug)->default_value(false) ,     "Debug flag")
      ("InFile",         po::value<string>(&data_file_name)->default_value(DATA_FILE_NAME) , "input file")  //default value of file set to macro
      ("OutFile",        po::value<string>(&output_name)->default_value(OUTPUT_NAME), "output file name")
      ("NEvents",        po::value<int>(&nevents)->default_value(N_ANALYZED_EVENTS), "number of events to analyze") //default number of events for which to explicitly make TTrees set to macro
      ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help") > 0){
        cout << desc << "\n";
        return 1;
    }
    #endif
    */
    //------

   //-------------------------
}

int main(int argc, char* argv[])
{
  const char* data_file_name = DATA_FILE_NAME;

  TFile* tf = new TFile(data_file_name, "open");
  TTree* tt = (TTree*) tf->Get("EventTree");

  LightJetAnalysis* analysis = new LightJetAnalysis(tt);
  analysis->Begin();
  analysis->Loop();
  analysis->End();

  delete analysis;
  
  return 0;
}
