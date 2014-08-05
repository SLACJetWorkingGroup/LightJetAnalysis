/* 
Written By: Todd Macdonald
Version A: Completed 7/29/14

This script reads text files that contain information related to reconstructed-truth jet matches and 
then presents this information via histograms. The histograms reflect how stable the correction method 
is across events with different levels of pileup, as well as the statistically relevance of the correction
method.

 */ 


#include <iostream.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "TLegend.h"
#include <math.h>

#define PT_RECON_MIN 4 //minimum threshold for a reconstructed jet to be considered as a real jet
#define PT_TRUTH_LOW 30 //low end pt bin for most consistent truth jet data   (recommended= 30)
#define PT_TRUTH_HIGH 60 //high end of bin for most consistent truth jet data (recommended= 50)
#define NPV_MIN 5
#define N_FILES 32

namespace std;

/*Prototypes*/
void graphInput(TH2F* offset_npv, TH2F* sigma_slope, float numSigma); //change back to float
inline int getLineColor(int numSigma);

/*Main script function*/
int slope() {

	//-------added 7/30/14 -----------
	TH2F* offset_npv = new TH2F("offset_npv", "d/dNPV Mean(Offset) At Different Corrections", 100, -2, 4, 100, -1, 1);     
	offset_npv->SetXTitle("Number of Sigma(Rho) Used for Correction");
	offset_npv->SetYTitle("Slope For Linearizion Fit Line of NPV vs. Mean(Offset)");
	offset_npv->SetOption("COL");
	//--------------------------------

	TH2F* sigma_slope = new TH2F("sigma_slope", "d/dNPV Sigma(Offset) At Different Corrections", 100, -2, 4, 40, 0, 0.1);
	for(float i = -2; i < 4; i += 0.5) { //iterate through all correction permutations 
		graphInput(offset_npv, sigma_slope, i); //graph histograms per sigma correction
	}
	/*for(int i = 0; i <= 3; i ++) { //iterate through all correction permutations 
		graphInput(offset_npv, sigma_slope, i); //graph histograms per sigma correction
	}*/
	sigma_slope->SetXTitle("Number of Sigma(Rho) Used for Correction");
	sigma_slope->SetYTitle("Slope For Linearizion Fit Line of NPV vs. Sigma(Offset)");
	TFile fp("sigma_slope.root", "RECREATE");
	sigma_slope->SetOption("COL");
	sigma_slope->Write("sigma_slope.root");

	//--------added 7/30/14---------
	TF1* foo = new TF1("linear_fit", "[0] + [1] * x", -1, 4);
	offset_npv->Fit("linear_fit");
	Double_t params[2];
	foo->GetParameters(params);
	//TLegend legend(20, 20, 100, 100, "Critical Points", to_string(-1 * params[0] / params[1]));
	//legend.AddEntry("X Intercept", )
	cout << "ideal sigma = " << -1 * params[0] / params[1] << endl;
	//delete foo;

	//--------added 8/04/14----------
	TF1* foo2 = new TF1("linear_fit_2", "[0] + [1] * x", -1, 4);
	sigma_slope->Fit("linear_fit_2");
	Double_t arr[2];
	foo2->GetParameters(arr);
	cout << "d/d(sigma of rho) = " << arr[1];
	//--------------------------------

	offset_npv->Write("offset_npv");
	delete foo;
	delete offset_npv;
	delete sigma_slope;
	//------------------------------

	return 0;
}

/*creates a set of histograms to a file for a given correction*/
void graphInput(TH2F* offset_npv, TH2F* sigma_slope, float numSigma)
{
	TH2F hist("hist", "Mean Offset Across Varying NPV", 11, 5, 60, 50, -30, 20); //temporary histogram for finding mean offset
	char buffer[100];
	float recon_pt;
	float truth_pt;
	int npv;
	float corr_quantized; //rho times area
	ifstream ifile;
	for(int i = 1; i <= N_FILES; i++) {
		char fileName[30];
		sprintf(fileName, "testing_%d.txt", i);
		ifile.open(fileName, fstream::in);
		while(ifile.getline(buffer, 100, '&')) { //the '&' symbol signifies the end of a line in a properly formatted file
			if(sscanf(buffer, "%f %f %d %f", &recon_pt, &truth_pt, &npv, &corr_quantized) != 4) { //the recon pt, truth pt, and event npv are derived from each line of text
				break;
			}
			float corr_pt = recon_pt - (numSigma * corr_quantized);
			if(corr_pt > PT_RECON_MIN && truth_pt > PT_TRUTH_LOW && truth_pt < PT_TRUTH_HIGH) { //applying pt truth bin and ensuring pt recon have sufficient momenta
				if(npv > NPV_MIN) hist.Fill(npv, corr_pt - truth_pt);
			}
		}
		//cout << "number of jets in vec = " << recon_vec.size() << endl;
		ifile.close(); //closing input file
	}
	TFile* fout = new TFile("temp.root", "RECREATE"); //creating output file; previous contents erased on first iteration
	hist.FitSlicesY();
	TH1D* sigma_hist;
    gDirectory->GetObject("hist_2", sigma_hist);  //obtaining histogram of npv vs. sigma(offset)
    TF1* linear_fit = new TF1("linear_fit", "[0] + [1] * sqrt(x)", 5, 60);  //line used to calculated linear fit
    sigma_hist->Fit("linear_fit");

    Double_t params[2];
    linear_fit->GetParameters(params); //obtaining slope and y intercept for fit lines    //UNCOMMENT CODE 8/05/14
    //cout << "Param 1 = " << params[1] << endl;
    sigma_slope->Fill(numSigma, params[1]); 

    //-------------added 7/30/14------------------
    TF1* linear_fit_b = new TF1("linear_fit_b", "[0] + [1] * x", 5, 60);  //line used to calculated linear fit
    hist.ProfileX()->Fit("linear_fit_b");
    Double_t params_b[2];
    linear_fit_b->GetParameters(params_b);
    offset_npv->Fill(numSigma, params_b[1]);
    //--------------------------------------------

    //--temporary addition 8/04/14----------------
    //remember to add x and y titles
    /*char s1[50];   
    char s2[50];
    sprintf(s1, "mean_%d", numSigma);
    sprintf(s2, "sigma_%d", numSigma); */
    /*float zero = 0;
    float one = 1;
    float two = 2;
    float three = 3;*/
    /*sigma_hist->SetTitle("NPV vs. Sigma(Offset)");
    sigma_hist->SetYTitle("Sigma(Offset)");
    sigma_hist->SetXTitle("NPV");
    hist.SetXTitle("NPV");
    hist.ProfileX()->SetYTitle("Mean Offset (GeV)");
    hist.ProfileX()->Write(s1);
    sigma_hist->Write(s2);*/
    /*if(numSigma == zero) {
    	hist.ProfileX()->Write("mean_zero");
    	sigma_hist->Write("sigma_zero");
    } else if(numSigma == one) {
    	hist.ProfileX()->Write("mean_one");
    	sigma_hist->Write("sigma_one");
    } else if(numSigma == two) {
    	hist.ProfileX()->Write("mean_two");
    	sigma_hist->Write("sigma_two");
    } else if((float)numSigma == three) {
    	hist.ProfileX()->Write("mean_three");
    	sigma_hist->Write("sigma_three");
    }*/ 
    	/*hist.ProfileX()->Write(s1);
    	sigma_hist->Write(s2);*/
    //}
    //--------------------------------------------

    delete fout;
}
