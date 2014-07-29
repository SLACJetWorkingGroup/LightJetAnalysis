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
#define N_ANALYZED_EVENTS 1000000  //must be synched
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
#define TRUTH_DIF_D 50

#define N_JETS_MAX 10
#define NPV_MAX 100

#define OUTPUT_NAME "proj.root"  //these two constants must be change simultaneously
#define N_SIGMA 0 //number of gaussian sigma deviations away from the median cluster pt
#define COLOR 4 //color corresponds to line color of graphs
#define PT_THRESHOLD 4  //minimum momentum (in GeV) for a group of clusters to be considered a jet
#define TRUTH_PT_THRESHOLD 10  //minimum pt for truth jet

#define ETA_ACCURACY_RANGE 2.4 //should be 2.4
#define ETA_EXTREME_ACCURATE 0.8 //rapidity range for most reliable area of calorimeter
//event specific rho, area of each jet, sigma per event

#define R_RANGE_A 20
#define R_RANGE_B 30
#define R_RANGE_C 50
#define R_RANGE_D 100

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


// Begin method
void LightJetAnalysis::Begin() 
{   
  ofile = new ofstream;
  ofile->open("matches.txt", ios::trunc);
    // create output file 
  fp_out = new TFile(OUTPUT_NAME, "UPDATE");   //"UPDATE" or "RECREATE"  

  matched = new TH1F("Matched", "Matched Jets", N_BINS, 0, 100);
  matched->SetXTitle("Momentum (GeV)");
  matched->SetYTitle("Matched Jets");

  truth = new TH1F("Truth_Jets", "Truth Jets", N_BINS, 0, 100);
  truth->SetXTitle("Momentum (GeV)");
  truth->SetYTitle("Truth Jets");
  truth->SetLineColor(8);

  truth_proximity = new TH1F("Truth_Jets_Proximity", "Truth Jets in Proximity To Recon Jets", N_BINS, 0, 60);
  truth_proximity->SetXTitle("Momentum (GeV)");
  truth_proximity->SetYTitle("Truth Jets in Proximity");
  truth_proximity->SetLineColor(2);

  sigma_mean = new TH2F("Sigma_Mean", "R vs. Sigma(R) / Mean(R) For 20 < Pt Truth < 30", N_BINS, 0, 2, N_BINS, 0, 2);
  sigma_mean->SetYTitle("Sigma(R) / Mean(R)");
  sigma_mean->SetXTitle("R");

  R_dist = new TH1F("R_dist", "Frequency of R Values For Matched Pt Truth [20, 30] GeV", N_BINS, 0, 1.8);
  R_dist->SetXTitle("R Value (Pt Reco / Pt Truth)");
  R_dist->SetYTitle("Number of Jets");

  R_dist_b = new TH1F("R_dist", "Frequency of R Values For Matched Pt Truth [30, 50] GeV", N_BINS, 0, 1.8);
  R_dist_b->SetXTitle("R Value (Pt Reco / Pt Truth)");
  R_dist_b->SetYTitle("Number of Jets");

  R_dist_c = new TH1F("R_dist", "Frequency of R Values For Matched Pt Truth [50, 100] GeV", N_BINS, 0, 1.8);
  R_dist_c->SetXTitle("R Value (Pt Reco / Pt Truth)");
  R_dist_c->SetYTitle("Number of Jets");

  true_vs_R = new TH2F("true_vs_R", "R For Varying Values of Pt Truth", N_BINS, 0, 60, N_BINS, 0, 2);
  true_vs_R->SetXTitle("Pt Truth");
  true_vs_R->SetYTitle("R Value (Recon Pt / Truth Pt)");

  npv_vs_resolution = new TH2F("rms_vs_npv", "Statistical Variance of Pt Reco / Pt Truth Across Different NPV", 30, 0, 60, N_BINS / 5, 0, 1.5);
  npv_vs_resolution->SetYTitle("RMS of R / Mean of R");
  npv_vs_resolution->SetXTitle("NPV");

  rms_offset = new TH2F("rms_offset", "RMS of Offset Across Different NPV", 11, 5, 60, N_BINS/2, 0, 10);
  rms_offset->SetYTitle("RMS of Offset (Pt Reco. - Pt Truth)");
  rms_offset->SetXTitle("NPV");

  mean_offset = new TH2F("mean_offset", "Mean of Offset Across Different NPV", 11, 5, 60, N_BINS/2, -30, 20);
  mean_offset->SetYTitle("Mean of Offset (Pt Reco. - Pt Truth)");
  mean_offset->SetXTitle("NPV");

  pt_resolution = new TH2F("pt_resolution", "Resolution For Varying Pt Truth", 6, 0, 60, N_BINS/2, 0, 1.5);
  pt_resolution->SetYTitle("RMS(Response) / <Response>");
  pt_resolution->SetXTitle("Pt Truth (GeV)");

  pt_response = new TH2F("pt_response", "Response For Matched Pt Truth Greater 10 GeV", 5, 10, 60, N_BINS/2, 0, 1.5);
  pt_response->SetYTitle("Pt Recon / Pt Truth");
  pt_response->SetXTitle("Truth Pt (GeV)");

  npv_response = new TH2F("npv_response", "Response of Matched Jets For Varying NPV", 6, 0, 60, N_BINS/2, 0, 1.5);
  npv_response->SetYTitle("Pt Recon / Pt Truth");
  npv_response->SetXTitle("NPV");
}

/*inline bool in_rapidity_range(int eta)
{
  return (eta >= (-1) * RAPIDITY_RANGE) && (eta <= RAPIDITY_RANGE);
}*/

// End


void LightJetAnalysis::End()
{
    //rho_histo->Fit("gauss");
    //rho_histo->Write("Rho");       //
    ofile->close();
    delete ofile;



    truth_proximity->Divide(truth);
    truth_proximity->SetYTitle("(Frequency Matched Truth Pt) / (Frequency All Truth Pt)");
    truth_proximity->SetTitle("Efficiency of Identifying Pt Truth Matches");
    truth_proximity->SetXTitle("Truth Pt (GeV)");
    //truth_proximity->Write("Efficiency_0"); 

    /*matched->Divide(truth_proximity);
    matched->SetYTitle("(Frequency Recon Pt) / (Frequency Matched Truth Pt)");
    matched->SetTitle("Ratio Distribution of Recon Pt / Matched Truth Pt");
    matched->SetXTitle("GeV");
    matched->Write("Ratio_Recon_Truth"); */

    sigma_mean->SetOption("COL");
    //sigma_mean->Write("Sigma_Mean_R");
    //sigma_mean->ProfileX()->Write("Sigma_Mean_Profile"); 
    delete sigma_mean;

   // R_dist->Write("r_dist_a_0");
    //R_dist_b->Write("r_dist_b_0");
    //R_dist_c->Write("r_dist_c_0");
    delete R_dist;
    delete R_dist_b;
    delete R_dist_c;

    true_vs_R->SetOption("COL");
    //true_vs_R->Write("Truth_Vs_R_1");
    //true_vs_R->ProfileX()->Write("True_Vs_R_Profile_");
    delete true_vs_R;

    //truth->Write("Truth_Jets");
    //matched->Write("Matched_Jets");
    //truth_proximity->Write("Truth_Jets_Proximity");
    delete truth;
    delete matched;
    delete truth_proximity;

    npv_vs_resolution->SetLineColor(COLOR);
    //npv_vs_resolution->Write("npv_res_0");
    npv_vs_resolution->ProfileX()->SetYTitle("RMS of R / Mean of R");
    //npv_vs_resolution->ProfileX()->Write("npv_resolution_0");
    delete npv_vs_resolution;

    //-------------------------------------------------
    //active code as of 7/24/14

    mean_offset->SetLineColor(COLOR);
    mean_offset->ProfileX()->SetYTitle("Mean of Offset (Pt Reco. - Pt Truth)");
    //mean_offset->/*ProfileX()->*/Write("mean_offset_nr");                          //commented out 7/28/14 2:09pm

    //TH1D* offset_hist = mean_offset->ProfileX()->DrawCopy();
    /*for(int i = 0; i < 11; i++) {
      TH1D* proj = mean_offset->ProjectionY("temp", i, i);
      proj->SetLineColor(COLOR);
      mean_offset->FitSlicesY(0, i, i);
      char fileName[30];
      sprintf(fileName, "%d_to_%d_nr", i * 5, (i + 1) * 5);
      mean_offset_2->Write(fileName);
      proj->Write(fileName);
      rms_offset->Fill(2.5 + (5 * i), proj->GetRMS());
    } */

    rms_offset->SetLineColor(COLOR);
    //rms_offset->ProfileX()->SetYTitle("RMS of Offset (Pt Reco. - Pt Truth)");
    //rms_offset->ProfileX()->Write("rms_offset_nr");
    //rms_offset->SetOption("COL");
    rms_offset->SetMarkerStyle(21);
    rms_offset->SetMarkerColor(COLOR);
    //rms_offset->Write("rms_offset_nr");



    //end of active code
    //------------------------------------------------------

    pt_resolution->SetLineColor(COLOR);
    //pt_resolution->Write("pt_res_0");
    pt_resolution->ProfileX()->SetYTitle("RMS(Response) / <Response>");
    //pt_resolution->ProfileX()->Write("pt_resolution_0");

    //-----------------------------------
    //start of new graphs 7/22/14
    //-----------------------------------

    /*TH1F resolution_avg("resolution_avg", "Resolution For Varying Pt Truth", 5, 10, 60);
    resolution_avg.SetYTitle("RMS(Response) / <Response>");
    resolution_avg.SetXTitle("Pt Truth (GeV)");

    TH1F npv_avg("resolution_avg", "Resolution For Varying NPV", 6, 0, 60);
    npv_avg.SetYTitle("RMS(Response) / <Response>");
    npv_avg.SetXTitle("NPV");

    pt_response->SetOption("COL");
    pt_response->SetLineColor(COLOR);
    //pt_response->Write("pt_resp_4");
    
    TH1D* p1 = pt_response->ProjectionY("resp_10_20_4", 0, 1);  //do I need to delete projections and profiles? 
    //p1->Write("resp_10_20_4");
    
    TH1D* p2 = pt_response->ProjectionY("resp_20_30_4", 1, 2);
    //p2->Write("resp_20_30_4");

    TH1D* p3 = pt_response->ProjectionY("resp_30_40_4", 2, 3);
    //p3->Write("resp_30_40_4");

    TH1D* p4 = pt_response->ProjectionY("resp_40_50_4", 3, 4);
    //p4->Write("resp_40_50_4");
    
    TH1D* p5 = pt_response->ProjectionY("resp_50_60_4", 4, 5); */
    //p5->Write("resp_50_60_4");

    /*resolution_avg.Fill(15, p1->GetRMS() / p1->GetMean());
    resolution_avg.Fill(25, p2->GetRMS() / p2->GetMean());
    resolution_avg.Fill(35, p3->GetRMS() / p3->GetMean());
    resolution_avg.Fill(45, p4->GetRMS() / p4->GetMean());
    resolution_avg.Fill(55, p5->GetRMS() / p5->GetMean());
    //resolution_avg.SetOption("COL");
    resolution_avg.SetFillColor(COLOR); */
    //resolution_avg.Write("resolution_avg_4");

    /*TH1D* proj1 = npv_response->ProjectionY("0_10", 0, 1);
    TH1D* proj2 = npv_response->ProjectionY("10_20", 1, 2);
    TH1D* proj3 = npv_response->ProjectionY("20_30", 2, 3);
    TH1D* proj4 = npv_response->ProjectionY("30_40", 3, 4);
    TH1D* proj5 = npv_response->ProjectionY("40_50", 4, 5);
    TH1D* proj6 = npv_response->ProjectionY("50_60", 5, 6);
    npv_avg.Fill(5,  proj1->GetMean() == 0 ? 0 : (proj1->GetRMS() / proj1->GetMean()));
    npv_avg.Fill(15, proj2->GetMean() == 0 ? 0 : (proj2->GetRMS() / proj2->GetMean()));
    npv_avg.Fill(25, proj3->GetMean() == 0 ? 0 : (proj3->GetRMS() / proj3->GetMean()));
    npv_avg.Fill(35, proj4->GetMean() == 0 ? 0 : (proj4->GetRMS() / proj4->GetMean()));
    npv_avg.Fill(45, proj5->GetMean() == 0 ? 0 : (proj5->GetRMS() / proj5->GetMean()));
    npv_avg.Fill(55, proj6->GetMean() == 0 ? 0 : (proj6->GetRMS() / proj6->GetMean()));
    npv_avg.SetFillColor(COLOR); */
    //npv_avg.Write("npv_avg_4");

    /*if(proj1->GetMean() == 0 || proj2->GetMean() == 0 || proj3->GetMean() == 0 || proj4->GetMean() == 0 || proj5->GetMean() == 0 || proj6->GetMean() == 0) {
      cout << "Divide by zero error!!!!!!!!" << endl;
    } */

    delete pt_response;
    delete npv_response;

    delete rms_offset;
    delete mean_offset;
    delete pt_resolution;


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
}

void LightJetAnalysis::calc_sub_jets(double& rho, double& sigma, vector<PseudoJet>& clusters, vector<PseudoJet>& sub_jets, bool isTruth, vector<double>& areas)
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
    if((signed)subtracted_pt > min_pt && abs(jets[i].rapidity()) < 0.8) {
      sub_jet.reset_momentum_PtYPhiM(subtracted_pt, jets[i].rapidity(), jets[i].phi(), jets[i].m()); 
      sub_jets.push_back(sub_jet);  //subtract from original jet 
      areas.push_back(jets[i].area());
    }
  }

}



int inline LightJetAnalysis::n_nontrivial_jets(vector<PseudoJet>& sub_jets)
{

  int count = 0;
  size_t sz = sub_jets.size();
  for(unsigned int i = 0; i < sz; i++) {
    if(sub_jets[i].rapidity() > -ETA_ACCURACY_RANGE && sub_jets[i].rapidity() < ETA_ACCURACY_RANGE) {
      //if(sub_jets[i].pt() > 30) {
        count++;
      //}
    }
  }
  return count;
}

void inline LightJetAnalysis::graph_truth(const vector<PseudoJet>& sub_jets, const vector<PseudoJet>& truth_jets, vector<double>& R_vec, vector<double>& pt_true_matches, vector<double>& offset_vec, vector<double>& areas, double& sigma)
{
  //vector<double> recon;
  for(unsigned int i = 0; i < truth_jets.size(); i++) {
    double min = -1;
    signed int minIndex = -1;
    double truth_pt = truth_jets[i].pt();

    for(unsigned int j = 0; j < sub_jets.size(); j++) {
      if(sub_jets[j].rapidity() > -ETA_ACCURACY_RANGE && sub_jets[j].rapidity() < ETA_ACCURACY_RANGE) {
        double deltR = truth_jets[i].delta_R(sub_jets[j]);
        if(deltR < JET_SEARCH_RADIUS) {
          min = deltR;
          minIndex = j;
          break;
        }
      }
    }
    if(minIndex >= 0) {

      //cout << "truth pt = " << truth_pt << endl;
      //!!!!!!!!!!!!!!!!!
      matched->Fill(sub_jets[minIndex].pt());
      truth_proximity->Fill(truth_pt);
      //!!!!!!!!!!!!!!
      //ofile << "adding here" << endl;
      float sigma_times_area = sigma * areas[minIndex];
      *ofile << sub_jets[minIndex].pt() << " " << truth_pt << " " << NPV << " " << sigma_times_area << "&" << endl;  //no pt truth bin!!!!!!

      if(truth_pt > TRUTH_DIF_B && truth_pt < TRUTH_DIF_C) {
        offset_vec.push_back(sub_jets[minIndex].pt() - truth_pt);
      }

      R_vec.push_back(sub_jets[minIndex].pt() / truth_pt);
      pt_true_matches.push_back(truth_pt);
      //recon.push_back(sub_jets[minIndex].pt());

    }
      
  }
  //cout << "number of matches = " << pt_true_matches.size() << endl;
  //ofile >> recon;
  //ofile >> pt_true_matches;
}

 void LightJetAnalysis::rms(vector<double>& vec, double& mean, double& rms)
 {
  mean = 0;
  rms = 0;
  for(size_t i = 0; i < vec.size(); i++) {
    mean += vec[i];
  }
  mean /= vec.size();
  for(size_t i = 0; i < vec.size(); i++) {
    rms += pow(vec[i] - mean, 2);
  }
  rms /= vec.size();
  rms = sqrt(rms < 0 ? 0 : rms); 
 }

 inline void LightJetAnalysis::graph_offset(const vector<PseudoJet>& sub_jets, const vector<PseudoJet>& truth_jets, vector<double>& offset_vec, vector<double>& pt_true_matches)
 {
    /*double mean_off;
    double rms_off;
    if(offset_vec.size() != 0) {
      rms(offset_vec, mean_off, rms_off);
      rms_offset->Fill(NPV, rms_off);
      mean_offset->Fill(NPV, mean_off);
    }*/

   // TH1D offset_hist("offset", "temporary", N_BINS/2, -30, 20);
    //cout << "offset vec size = " << offset_vec.size() << endl;
    for(size_t i = 0; i < offset_vec.size(); i++) {
      mean_offset->Fill(NPV, offset_vec[i]);
      //cout << "NPV = " << NPV << "offset = " << offset_vec[i] << endl;
    }

 }

void LightJetAnalysis::graph_resolution(const vector<PseudoJet>& sub_jets, const vector<PseudoJet>& truth_jets, vector<double>& R_vec, vector<double>& pt_true_matches) 
{
    //for graph on 7/11
  double mean;
  double rms_R;
  double mean_constr;  //for resolution
  double rms_constr;  //for resolution
  vector<double> response_constr;
  //vector<double> pt_constr;


  if(R_vec.size() != 0) {
    rms(R_vec, mean, rms_R);
  }
  //use for loop to make vector<double> response_constr
  
  if(pt_true_matches.size() != R_vec.size()) throw "ERROR: # R calculations is not equal to number of Pt Truth Matches";
  
  for(size_t i = 0; i < R_vec.size(); i++) {

    //cout << "rms = " << rms_R << "mean = " << mean << endl;
    
    //cout << "pt = " << pt_true_matches[i] << "resolution = " << rms_R / mean << endl;
    pt_response->Fill(pt_true_matches[i], R_vec[i]);
    //npv_response->Fill(NPV, R_vec[i]);
    pt_resolution->Fill(pt_true_matches[i], rms_R / mean);

    if(pt_true_matches[i] > TRUTH_DIF_B && pt_true_matches[i] < TRUTH_DIF_D) {
      npv_response->Fill(NPV, R_vec[i]);
      response_constr.push_back(R_vec[i]);
      //pt_constr.push_back(pt_true_matches[i]);
    }

    /*if(pt_true_matches[i] > R_RANGE_A && pt_true_matches[i] < R_RANGE_B) {
      rms_mean->Fill(R_vec[i], rms_R / mean);
      R_dist->Fill(R_vec[i]);
    }
    if(pt_true_matches[i] > R_RANGE_B && pt_true_matches[i] < R_RANGE_C) {
      R_dist_b->Fill(R_vec[i]);
    }
    if(pt_true_matches[i] > R_RANGE_C && pt_true_matches[i] < R_RANGE_D) {
      R_dist_c->Fill(R_vec[i]);
    } */

    //true_vs_R->Fill(pt_true_matches[i], rms_R / mean);  
    //true_vs_R->Fill(pt_true_matches[i], R_vec[i]);
  }


  //note: right now you do not have a restriction on the 2d plot below, though it 
  //will need a restriction once the bug in the resolution calculation is identified
  //npv_vs_resolution->Fill(NPV, rms_R / mean);  //DELETE THIS LINE IMMEDIATELY
  if(response_constr.size() != 0) {
    rms(response_constr, mean_constr, rms_constr);
    npv_vs_resolution->Fill(NPV, rms_constr / mean_constr);  //CHANGE THIS BACK IMMEDIATELY
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

  double rho;
  double sigma;
  vector<PseudoJet> clusters; 
  vector<PseudoJet> sub_jets;
  vector<PseudoJet> truth_clusters;
  vector<PseudoJet> truth_jets;
  vector<double> areas;
  vector<double> does_absolutely_nothing;
  vector<double> R_vec;
  vector<double> pt_true_matches;
  vector<double> offset_vec;
  vector<double> pt_constr; //all pt truth matches constrained by bin 

  get_cluster_vec(clusters, false);
  get_cluster_vec(truth_clusters, true);
  calc_rho_sigma(rho, sigma, clusters);

  calc_sub_jets(rho, sigma, clusters, sub_jets, false, areas);
  calc_sub_jets(rho, sigma, truth_clusters, truth_jets, true, does_absolutely_nothing);


  //int njets_x = n_nontrivial_jets(sub_jets);

  //npv_histo->Fill(NPV, njets_x);  //here

  graph_truth(sub_jets, truth_jets, R_vec, pt_true_matches, offset_vec, areas, sigma);

  for(size_t i = 0; i < truth_jets.size(); i++) {
    truth->Fill(truth_jets[i].pt());
  }

  //graph_resolution(sub_jets, truth_jets, R_vec, pt_true_matches);
  graph_offset(sub_jets, truth_jets, offset_vec, pt_true_matches);
 //end of code for 7/11

 // cout << "NPV = " << NPV << "; NJets = " << sub_jets.size() << "; rho = " << rho << endl;


}