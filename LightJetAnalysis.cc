#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <algorithm>
#include <fstream>



//#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2F.h"


#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "LightJetAnalysis.h" //new edit

#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif

#ifndef ROOT_TProfile
#include "TProfile.h"
#endif

using namespace std;
using namespace fastjet;
using namespace ROOT;


#define N_BINS 100
#define N_ANALYZED_EVENTS 1000    //must be synched
#define ARBITRARY_RHO_MAX 35  
#define ARBITRARY_SIGMA_MAX 15  
#define RAPIDITY_RANGE 2
#define RAPIDITY_SPECTRUM 5
#define JET_SEARCH_RADIUS 0.4   //why not 0.7 ???
#define ARBITRARY_PT_MAX 250
#define ARBITRARY_DELT_PT_MAX 15

#define PT_MIN_A 20
#define PT_MIN_B 30
#define PT_MIN_C 50

#define TRUTH_DIF_A 20
#define TRUTH_DIF_B 30
#define TRUTH_DIF_C 40

#define N_JETS_MAX 10
#define NPV_MAX 100

#define OUTPUT_NAME "0_sigma.root"  //these two constants must be change simultaneously
#define N_SIGMA 0
#define PT_THRESHOLD 20  //minimum momentum (in GeV) for a group of clusters to be considered a jet
#define TRUTH_PT_THRESHOLD 10

#define ETA_ACCURACY_RANGE 2.4 //should be 2.4
//event specific rho, area of each jet, sigma per event

//pt - (rho + sigma)sqrt(area)

// Constructor 
LightJetAnalysis::LightJetAnalysis(TTree* tree)
{
    //if(debug) cout << "LightJetAnalysis::LightJetAnalysis Start " << endl;
    debug = false;
    Init(tree);
    if(debug) throw "Init tree method unsuccessful";
    //if(debug) cout << "LightJetAnalysis::LightJetAnalysis End " << endl;
}

// Destructor 
LightJetAnalysis::~LightJetAnalysis()
{
}

inline void LightJetAnalysis::init_rho_hist() 
{
  rho_histo = new TH1F("Rho", "Rho of Particle Clusters Across Collision Events", N_BINS, 0, ARBITRARY_RHO_MAX);
  rho_histo->SetXTitle("Rho Values of Events (GeV)");
  rho_histo->SetYTitle("Frequency of Rho Values");
}

inline void LightJetAnalysis::init_sigma_hist()
{
  sigma_histo = new TH1F("Sigma", "Sigma of Particle Clusters Across Collision Events", N_BINS, 0, ARBITRARY_SIGMA_MAX);
  sigma_histo->SetXTitle("Sigma Values of Events (GeV)");
  sigma_histo->SetYTitle("Frequency of Sigma Values");
}

TH1F* LightJetAnalysis::init_jet_hist(int min_pt) //the distribution of pt across jets
{
  const char* s = "Distribution of Jets Above Pt Threshold of ";
  char cstr[strlen(s) + 2];
  sprintf(cstr, "%s%d\n", s, min_pt);
  TH1F* histo = new TH1F("Jets", cstr, N_BINS, min_pt, ARBITRARY_PT_MAX);
  histo->SetXTitle("Pt of Jets (GeV)");
  histo->SetYTitle("Frequency of Pt Values");
  return histo;
}

TH1F* LightJetAnalysis::init_l_hist(string data_type, int min, int max, string units)
{
  TH1F* histo = new TH1F("Lead_Jets", (data_type + " of Leading Jets Across Collision Events").c_str(), N_BINS, min, max);
  histo->SetXTitle((data_type + " of Leading Jets " + units).c_str());
  histo->SetYTitle(("Frequency of " + data_type + " Values").c_str());
  return histo;
}

 inline void LightJetAnalysis::init_2D_npv_hist()
 {
  npv_histo = new TH2F("NPV_2D", "Number of Jets and Primary Vertices Across Collision Events", N_BINS, 0, 60, N_BINS, 0, 8); //recent change
  npv_histo->SetOption("COL");
  npv_histo->SetXTitle("Number of Primary Vertices");
  npv_histo->SetYTitle("Number of Jets");
 }

inline void LightJetAnalysis::init_truth_hist()
{
  truth_2d_histo = new TH2F("Truth_2D", "Density of Pt and NPV Values of Reconstructed-Truth Jet Matches", N_BINS, 0, 60, N_BINS, -30, 30); //recent change
  truth_2d_histo->SetOption("COL");
  truth_2d_histo->SetXTitle("Number of Primary Vertices");
  truth_2d_histo->SetYTitle("Difference Between Pt Truth and Pt Reconstructed");
}

TH1F* LightJetAnalysis::get_truth_hist(int min_pt, int max_pt)
{
  char buffer[80];
  sprintf(buffer, "Reconstructed-Truth Jet Matches For Truth Jets Between %d GeV and %d GeV\n", min_pt, max_pt);
  TH1F* histo = new TH1F("Jet Matches", buffer, N_BINS, -25, 25); //recent change
  histo->SetXTitle("Pt Recon. - Pt Truth (GeV)");
  histo->SetYTitle("Frequency of Pt Values");
  return histo;
}


// Begin method
void LightJetAnalysis::Begin() 
{   
    // create output file 
  fp_out = new TFile(OUTPUT_NAME, "RECREATE");   
  init_rho_hist();
  init_sigma_hist();
  jet_dist_a = init_jet_hist(PT_MIN_A);
  jet_dist_b = init_jet_hist(PT_MIN_B);
  jet_dist_c = init_jet_hist(PT_MIN_C);
  leading_jet_pt = init_l_hist("Pt", 15, 65, "(GeV)");  //recent change
  leading_jet_eta = init_l_hist("Eta", (-1) * RAPIDITY_SPECTRUM, RAPIDITY_SPECTRUM);  //recent change
  init_2D_npv_hist();
  init_truth_hist();
  truth_histo_a = get_truth_hist(TRUTH_DIF_A, TRUTH_DIF_B);
  truth_histo_b = get_truth_hist(TRUTH_DIF_B, TRUTH_DIF_C);
  arr = (bin*)calloc(N_BINS, sizeof(bin));


  jet_areas = new TH1F("Jet_Areas", "Jet areas Across All Events", N_BINS, 0, 3);
  jet_areas->SetXTitle("Area of individual jet");
  jet_areas->SetYTitle("Frequency of Area Value");
}

/*inline bool in_rapidity_range(int eta)
{
  return (eta >= (-1) * RAPIDITY_RANGE) && (eta <= RAPIDITY_RANGE);
}*/

// End

inline void LightJetAnalysis::write_npv_avg()
{
  TH1F graph("NPV", "Average Primary Vertex Count For a Selection of NPV", N_BINS, 0, NPV_MAX);
  double avg_jets;
  for(int i = 0; i < N_BINS; i++) {
    if(arr[i].nevents != 0) {
      avg_jets = arr[i].total_jets / arr[i].nevents;
    } else {
      avg_jets = 0;
    }
    graph.Fill((double)NPV_MAX * i / N_BINS, avg_jets); 
  }
  graph.SetXTitle("Number of Primary Vertices");
  graph.SetYTitle("Average Number of Jets");
  graph.Write("NPV vs <NJets>");
  delete[] arr;
}


void LightJetAnalysis::End()
{
    //rho_histo->Fit("gauss");
    rho_histo->Write("Rho");
    sigma_histo->Write("Sigma");
    jet_dist_a->Write("Jet Distribution A");
    jet_dist_b->Write("Jet Distribution B");
    jet_dist_c->Write("Jet Distribution C");
    leading_jet_pt->Write("Leading_Jet_Pt");
    leading_jet_eta->Write("Leading_Jet_Eta");
    npv_histo->Write("NPV_vs_Number_Jets");
    truth_2d_histo->Write("Pt_Recon_Matches");
    truth_histo_a->Write("Truth_A");
    truth_histo_b->Write("Truth_B");
    
    write_npv_avg();

    jet_areas->Write("Jet areas");


    delete jet_areas;

     //TCanvas* c1 = new TCanvas("c1", "X Profile", 200, 200, 200 , 200);

    npv_histo->ProfileX()->Write("X_Profile_NPV_Norm");
    truth_2d_histo->ProfileX()->Write("X_Profile_Truth");


    delete rho_histo;
    delete sigma_histo;
    delete jet_dist_a;
    delete jet_dist_b;
    delete jet_dist_c;
    delete leading_jet_eta;
    delete leading_jet_pt;
    delete npv_histo;
    delete truth_2d_histo;
    delete truth_histo_a;
    delete truth_histo_b;
    delete fp_out;
}

// loop method
void LightJetAnalysis::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0;
   
   for (Long64_t j = 0; j < min(nentries, (Long64_t)N_ANALYZED_EVENTS); j++) {
      if (LoadTree(j) < 0) break;
      nbytes += fChain->GetEntry(j);
      Analyze();
   }
}

 inline void LightJetAnalysis::get_cluster_vec(vector<PseudoJet>& clusters, bool isTruth)
 {
    float* cl_px;
    float* cl_py;
    float* cl_pz;
    float* cl_E;
    int num_clust;
    if(isTruth) {
      cl_px = truth_px;
      cl_py = truth_py;
      cl_pz = truth_pz;
      num_clust = truth_n;
      cl_E = truth_E;
    } else {
      cl_px = cl_lc_px;
      cl_py = cl_lc_py;
      cl_pz = cl_lc_pz;
      num_clust = cl_lc_n;
      cl_E = cl_lc_E;
    }
    for(int i = 0; i < num_clust; i++) {  //iterate through the total number of jets
        if(cl_E[i] >= 0){
          PseudoJet cl(cl_px[i], cl_py[i], cl_pz[i], cl_E[i]);
           clusters.push_back(cl);
        }
    }
 }

inline void LightJetAnalysis::calc_rho_sigma(double& rho, double& sigma, vector<PseudoJet>& clusters) 
{
  JetMedianBackgroundEstimator cluster_estimator(SelectorAbsRapMax(RAPIDITY_RANGE), JetDefinition(kt_algorithm, JET_SEARCH_RADIUS), active_area_explicit_ghosts);  
  cluster_estimator.set_particles(clusters); //vector
  rho = cluster_estimator.rho();
  sigma = cluster_estimator.sigma();
  rho_histo->Fill(rho);
  sigma_histo->Fill(sigma);
}

void LightJetAnalysis::calc_sub_jets(double& rho, double& sigma, vector<PseudoJet>& clusters, vector<PseudoJet>& sub_jets, bool isTruth)
{ 
  double correction = rho + (N_SIGMA * sigma);
  int min_pt = PT_THRESHOLD;
  if(isTruth) {
    correction = 0;
    min_pt = TRUTH_PT_THRESHOLD;
  }
  vector<PseudoJet> jets;
  AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(5.0));  //is this value supposed to be 5.0 or 0.5?
  ClusterSequenceArea cs(clusters, JetDefinition(antikt_algorithm, JET_SEARCH_RADIUS), area_def);
  jets = sorted_by_pt(cs.inclusive_jets());   //can add argument here for min pt 10 for truth 
  int jets_sz = jets.size();
  for(int i = 0; i < jets_sz; i++) {
    //double area = cs.area(jets[i]);    //change
    if(!jets[i].has_area()) throw "BUG: A jet does not contain an area";
    PseudoJet sub_jet;
    double subtracted_pt = jets[i].pt() - (correction) * jets[i].area();  
    if((signed)subtracted_pt > min_pt) {
      sub_jet.reset_momentum_PtYPhiM(subtracted_pt, jets[i].rapidity(), jets[i].phi_std(), jets[i].m()); 
      sub_jets.push_back(sub_jet);  //subtract from original jet 
      jet_areas->Fill(jets[i].area());   //RECENT
    }
  }

}

void inline LightJetAnalysis::graph_jet_dist(vector<PseudoJet>& sub_jets)
{
  int sz = sub_jets.size();
  for(int i = 0; i < sz; i++) {
    PseudoJet j = sub_jets[i];
    //ifarea(j.pt() > PT_THRESHOLD) {
      jet_dist_a->Fill(j.pt());
      jet_dist_b->Fill(j.pt());
      jet_dist_c->Fill(j.pt());
    //}
  }
}

void inline LightJetAnalysis::graph_leading_jets(vector<PseudoJet>& sub_jets)
{
  if(sub_jets.size() != 0 /*&& sub_jets[0].pt() > PT_THRESHOLD*/) {    //recently added second condition
    leading_jet_pt->Fill(sub_jets[0].pt());
    leading_jet_eta->Fill(sub_jets[0].rapidity());
    //cout << "Pt = " << sub_jets[0].pt() << " Eta = " << sub_jets[0].rapidity() << endl;
  }
}

int inline LightJetAnalysis::n_nontrivial_jets(vector<PseudoJet>& sub_jets)
{
  int count = 0;
  size_t sz = sub_jets.size();
  for(unsigned int i = 0; i < sz; i++) {
    if(sub_jets[i].rapidity() > -ETA_ACCURACY_RANGE && sub_jets[i].rapidity() < ETA_ACCURACY_RANGE) {
      count++;
    }
  }
  return count;
}

void inline LightJetAnalysis::graph_truth(const vector<PseudoJet>& sub_jets, const vector<PseudoJet>& truth_jets)
{
  for(unsigned int i = 0; i < sub_jets.size(); i++) {
    double min = -1;
    signed int minIndex = -1;
    if(sub_jets[i].rapidity() > -ETA_ACCURACY_RANGE && sub_jets[i].rapidity() < ETA_ACCURACY_RANGE) {
        for(unsigned int j = 0; j < truth_jets.size(); j++) {
          double deltR = sub_jets[i].delta_R(truth_jets[j]);
          //cout << "DeltR = " << deltR << endl;
          if((min == -1) || (deltR < min)) {
            min = deltR;
            minIndex = j;
          }
        }
        if(min < JET_SEARCH_RADIUS && minIndex >= 0) {
          //cout << "--------------Match found---------------" << endl;
          double truth_pt = truth_jets[minIndex].pt();
          truth_2d_histo->Fill(NPV, sub_jets[i].pt() - truth_pt);
          if(truth_pt > TRUTH_DIF_A && truth_pt < TRUTH_DIF_B) {
            truth_histo_a->Fill(sub_jets[i].pt() - truth_pt);
          }
          if(truth_pt > TRUTH_DIF_B && truth_pt < TRUTH_DIF_C) {
            truth_histo_b->Fill(sub_jets[i].pt() - truth_pt);
          }
          //cout << "difference is " << sub_jets[i].pt() - truth_pt << endl;
        }
      }
  }
}


// Analyze method, called for each event
void LightJetAnalysis::Analyze()
{
  /*cout << "{ ";
  for(int i = 0; i < cl_lc_n; i++) {
      cout << cl_lc_E[i] << ", ";
  }
  cout << "}" << endl;

  exit(1);*/

  if(NPV > 20) cout << "NPV = " << NPV << endl;

  double rho;
  double sigma;
  vector<PseudoJet> clusters; 
  vector<PseudoJet> sub_jets;
  vector<PseudoJet> truth_clusters;
  vector<PseudoJet> truth_jets;

  get_cluster_vec(clusters, false);
  get_cluster_vec(truth_clusters, true);
  calc_rho_sigma(rho, sigma, clusters);

  calc_sub_jets(rho, sigma, clusters, sub_jets, false);
  calc_sub_jets(rho, sigma, truth_clusters, truth_jets, true);

  graph_jet_dist(sub_jets);
  graph_leading_jets(sub_jets);


  int njets_x = n_nontrivial_jets(sub_jets);

  npv_histo->Fill(NPV, njets_x);  //here

  graph_truth(sub_jets, truth_jets);

  int index = (double)NPV / NPV_MAX * N_BINS;
  if(index < N_BINS) {
    arr[index].nevents++;
    arr[index].total_jets += njets_x;  //here
  }


}