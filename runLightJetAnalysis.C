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

/*#define OUTPUT_NAME "out.root"
#define N_ANALYZED_EVENTS 1000 //this value is only necessary in the commented code*/
#define DATA_FILE_NAME "/u/at/pnef/Work/Data/forTodd/20140703.15.27_ClustersAndTruth.PythJ1to3mc12aJETMET.jetmet2012.root"

///atlas/output/pnef/skimmed.20140703.15.27_ClustersAndTruth.PythJ1to3mc12aJETMET.jetmet2012.root

static inline void init(char* data_file_name, int argc, char* argv);

static inline void init(char* data_file_name, int argc, char* argv) 
{
      // argument parsing  ------------------------
    cout << "Called as: ";
    for(int i = 0; i < argc; i++) {
        cout << argv[i] << " ";
    }
    cout << endl;

   

   /* #if boostflag == 1 // command line parsing if boost is installed 
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("Debug",          po::value<int>(&debug)->default_value(false) ,     "Debug flag")
      ("InFile",         po::value<string>(&data_file_name)->default_value(DATA_FILE_NAME) , "input file")  //default value of file set to macro
      ("OutFile",        po::value<string>(&output_name)->default_value(OUTPUT_NAME), "output file name")
      ("NEvents",        po::value<int>(&nevents)->default_value(N_ANALYZED_EVENTS), "number of events to analyze") //default number of events for which to explicitly make TTrees set to macro
      ("FirstEvent",     po::value<int>(&firstevent)->default_value(0), "whatever")
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

/*int eventStart = -1;
int eventEnd = -1;
#if boostflag == 1 // command line parsing if boost is installed 
    po::options_description desc("allowed options");
    desc.add_options()
      ("start",     po::value<int>(&eventStart)->default_value(-1), "first event to start analyzing on this run") //first event
      ("end",     po::value<int>(&eventEnd)->default_value(-1), "last event to analyze this run of program") //last event
      ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
#endif*/


int main(int argc, char* argv[])
{
  int eventStart = -1;
  int eventEnd = -1;
  #if boostflag == 1 // command line parsing if boost is installed 
    po::options_description desc("allowed options");
    desc.add_options()
      ("begin",     po::value<int>(&eventStart)->default_value(-1), "first event to start analyzing on this run") //first event
      ("end",     po::value<int>(&eventEnd)->default_value(-1), "last event to analyze this run of program") //last event
      ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  #endif

  //if(argc != 3) throw "must provide exactly two arguments";

  //throw "testing that this is recent code";


  const char* data_file_name = DATA_FILE_NAME;
  cout << data_file_name << endl;
  
   cout << "event start = " << eventStart << " event end = " << eventEnd << endl;
  if(eventStart < 0 || eventEnd < 0) cout << "Default Values in use!!!!!!!!" << endl;
  /*if(eventStart < 0 || eventEnd < 0) {
    throw "Must provide two integers as arguments that are greater than zero";
  }
  if(eventStart > eventEnd) {
    throw "The starting event must be greater than the ending event";
  } */

  //cout << "first arg = " << atoi(argv[1]) << " second arg = " << atoi(argv[2]) << endl;

  TFile* tf = new TFile(data_file_name, "open");
  TTree* tt = (TTree*) tf->Get("EventTree");
  cout << tf->GetName() << endl;

  LightJetAnalysis* analysis = new LightJetAnalysis(tt);
  analysis->Begin(eventEnd / 8000);
  //analysis->Loop(atoi(argv[1]), atoi(argv[2]));
  analysis->Loop(eventStart, eventEnd);
  analysis->End();

  delete analysis;
  
  return 0;
}
