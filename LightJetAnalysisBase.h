//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  2 11:07:09 2014 by ROOT version 5.34/17
// from TTree EventTree/Tree with event-by-event variables
// found on file: /atlas/output/pnef/skimmed.20140701.11.56_ClustersAndTruth.PythJ1to3mc12aJETMET.jetmet2012.root
//////////////////////////////////////////////////////////

#ifndef LightJetAnalysisBase_h
#define LightJetAnalysisBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class LightJetAnalysisBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EventNumber;
   Int_t           RunNumber;
   Float_t         Weight;
   Float_t         Mu;
   Int_t           NPVtruth;
   Int_t           NPV;
   Int_t           cl_lc_n;
   Float_t         cl_lc_px[2000];   //[cl_lc_n]       //changed all arrays capacity (1530 or 1000) -> (2000)
   Float_t         cl_lc_py[2000];   //[cl_lc_n]
   Float_t         cl_lc_pz[2000];   //[cl_lc_n]
   Float_t         cl_lc_E[2000];   //[cl_lc_n]
   Int_t           cl_em_n;
   Float_t         cl_em_px[2000];   //[cl_em_n]
   Float_t         cl_em_py[2000];   //[cl_em_n]
   Float_t         cl_em_pz[2000];   //[cl_em_n]
   Float_t         cl_em_E[2000];   //[cl_em_n]
   Int_t           truth_n;
   Float_t         truth_px[2000];   //[truth_n]
   Float_t         truth_py[2000];   //[truth_n]
   Float_t         truth_pz[2000];   //[truth_n]
   Float_t         truth_E[2000];   //[truth_n]
   Float_t         truth_id[2000];   //[truth_n]

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_Mu;   //!
   TBranch        *b_NPVtruth;   //!
   TBranch        *b_NPV;   //!
   TBranch        *b_cl_lc_n;   //!
   TBranch        *b_cl_lc_px;   //!
   TBranch        *b_cl_lc_py;   //!
   TBranch        *b_cl_lc_pz;   //!
   TBranch        *b_cl_lc_E;   //!
   TBranch        *b_cl_em_n;   //!
   TBranch        *b_cl_em_px;   //!
   TBranch        *b_cl_em_py;   //!
   TBranch        *b_cl_em_pz;   //!
   TBranch        *b_cl_em_E;   //!
   TBranch        *b_truth_n;   //!
   TBranch        *b_truth_px;   //!
   TBranch        *b_truth_py;   //!
   TBranch        *b_truth_pz;   //!
   TBranch        *b_truth_E;   //!
   TBranch        *b_truth_id;   //!

   LightJetAnalysisBase(TTree *tree=0);
   virtual ~LightJetAnalysisBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef LightJetAnalysisBase_cxx
LightJetAnalysisBase::LightJetAnalysisBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/atlas/output/pnef/skimmed.20140701.11.56_ClustersAndTruth.PythJ1to3mc12aJETMET.jetmet2012.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/atlas/output/pnef/skimmed.20140701.11.56_ClustersAndTruth.PythJ1to3mc12aJETMET.jetmet2012.root");
      }
      f->GetObject("EventTree",tree);

   }
   Init(tree);
}

LightJetAnalysisBase::~LightJetAnalysisBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t LightJetAnalysisBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t LightJetAnalysisBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void LightJetAnalysisBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("Mu", &Mu, &b_Mu);
   fChain->SetBranchAddress("NPVtruth", &NPVtruth, &b_NPVtruth);
   fChain->SetBranchAddress("NPV", &NPV, &b_NPV);
   fChain->SetBranchAddress("cl_lc_n", &cl_lc_n, &b_cl_lc_n);
   fChain->SetBranchAddress("cl_lc_px", cl_lc_px, &b_cl_lc_px);
   fChain->SetBranchAddress("cl_lc_py", cl_lc_py, &b_cl_lc_py);
   fChain->SetBranchAddress("cl_lc_pz", cl_lc_pz, &b_cl_lc_pz);
   fChain->SetBranchAddress("cl_lc_E", cl_lc_E, &b_cl_lc_E);
   fChain->SetBranchAddress("cl_em_n", &cl_em_n, &b_cl_em_n);
   fChain->SetBranchAddress("cl_em_px", cl_em_px, &b_cl_em_px);
   fChain->SetBranchAddress("cl_em_py", cl_em_py, &b_cl_em_py);
   fChain->SetBranchAddress("cl_em_pz", cl_em_pz, &b_cl_em_pz);
   fChain->SetBranchAddress("cl_em_E", cl_em_E, &b_cl_em_E);
   fChain->SetBranchAddress("truth_n", &truth_n, &b_truth_n);
   fChain->SetBranchAddress("truth_px", truth_px, &b_truth_px);
   fChain->SetBranchAddress("truth_py", truth_py, &b_truth_py);
   fChain->SetBranchAddress("truth_pz", truth_pz, &b_truth_pz);
   fChain->SetBranchAddress("truth_E", truth_E, &b_truth_E);
   fChain->SetBranchAddress("truth_id", truth_id, &b_truth_id);
   Notify();
}

Bool_t LightJetAnalysisBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LightJetAnalysisBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t LightJetAnalysisBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef LightJetAnalysisBase_cxx
