#ifndef  LightJetAnalysis_H
#define  LightJetAnalysis_H

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"

#include "TFile.h"
#include "TTree.h"

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"


#include "LightJetAnalysisBase.h"

using namespace std;
using namespace fastjet;

typedef struct {
    int nevents;
    double total_jets;
} bin;

class LightJetAnalysis : public LightJetAnalysisBase {
    public:
        LightJetAnalysis (TTree *tree);
        ~LightJetAnalysis ();
        void Begin();
        void End();
        void Loop();
        void Analyze();
    private:
        bool  debug;
        bin* arr;
        TFile* fp_out;
        TH1F* rho_histo;
        TH1F* sigma_histo;
        TH1F* jet_dist_a;
        TH1F* jet_dist_b;
        TH1F* jet_dist_c;
        TH1F* leading_jet_pt;
        TH1F* leading_jet_eta;
        TH2F* npv_histo;
        TH2F* truth_2d_histo;
        TH1F* truth_histo_a;
        TH1F* truth_histo_b;
        inline void init_rho_hist();
        inline void init_sigma_hist();
        TH1F* init_jet_hist(int min_pt);
        TH1F* init_l_hist(string data_type, int min, int max, string units = "");
        inline void init_2D_npv_hist();
        inline void init_1D_npv_hist();
        inline void get_cluster_vec(vector<PseudoJet>& clusters, bool isTruth);
        inline void calc_rho_sigma(double& rho, double& sigma, vector<PseudoJet>& clusters);
        inline void graph_jet_dist(vector<PseudoJet>& sub_jets);
        inline void graph_leading_jets(vector<PseudoJet>& sub_jets);
        void calc_sub_jets(double& rho, double& sigma, vector<PseudoJet>& clusters, vector<PseudoJet>& sub_jets, bool isTruth);
        inline void write_npv_avg();
        int n_nontrivial_jets(vector<PseudoJet>& sub_jets);
        inline void init_truth_hist();
        TH1F* get_truth_hist(int min_pt, int max_pt);
        inline void graph_truth(const vector<PseudoJet>& sub_jets, const vector<PseudoJet>& truth_jets);
        //static void inline create_jet_hist(vector<PseudoJet>& sub_jets);
};



#endif

     /*inline void Debug(int debug){
            debug_ = debug;
        }
        inline void SetOutName(string outname){
            outname_ = outname;
        }
        inline void SetNEvents(int nevents){
            n_events = nevents;
        }*/

           //Long64_t n_events;
