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


class LightJetAnalysis : public LightJetAnalysisBase {
    public:
        LightJetAnalysis (TTree *tree);
        ~LightJetAnalysis ();
        void Begin();
        void End();
        void Loop();
        void Analyze();
    private:
        ofstream* ofile; 
        bool  debug;
        TFile* fp_out;

        TH1F* matched;
        TH1F* truth;
        TH1F* truth_proximity;

        TH2F* sigma_mean;
        TH1F* R_dist;
        TH1F* R_dist_b;
        TH1F* R_dist_c;
        TH2F* true_vs_R;
        TH2F* npv_vs_resolution;
        TH2F* rms_offset;
        TH2F* mean_offset;
        TH2F* pt_resolution;
        TH2F* pt_response;
        TH2F* npv_response;

        inline void get_cluster_vec(vector<PseudoJet>& clusters, bool isTruth);
        inline void calc_rho_sigma(double& rho, double& sigma, vector<PseudoJet>& clusters);

        void rms(vector<double>& vec, double& mean, double& rms);
        void calc_sub_jets(double& rho, double& sigma, vector<PseudoJet>& clusters, vector<PseudoJet>& sub_jets, bool isTruth, vector<double>& areas);
        int n_nontrivial_jets(vector<PseudoJet>& sub_jets);
        TH2F* init_truth_hist(int min, int max);
        TH1F* get_truth_hist(int min_pt, int max_pt);
        inline void graph_truth(const vector<PseudoJet>& sub_jets, const vector<PseudoJet>& truth_jets, vector<double>& R_vec, vector<double>& pt_true_matches, vector<double>& offset_vec, vector<double>& areas, double& sigma);
        inline void graph_resolution(const vector<PseudoJet>& sub_jets, const vector<PseudoJet>& truth_jets, vector<double>& R_vec, vector<double>& pt_true_matches);
        inline void graph_offset(const vector<PseudoJet>& sub_jets, const vector<PseudoJet>& truth_jets, vector<double>& offset_vec, vector<double>& pt_true_matches);
};



#endif

