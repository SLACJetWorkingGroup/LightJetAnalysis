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

#define PT_RECON_MIN 4 //minimum threshold for a reconstructed jet to be considered as a real jet
#define PT_TRUTH_LOW 30 //low end pt bin most consistent truth jet data
#define PT_TRUTH_HIGH 50 //high end of bin for most consistent truth jet data


/*Prototypes*/
void graphInput(TH2F* sigma_slope, float numSigma);
inline int getLineColor(int numSigma);

/*Main script function*/
int slope() {
	TH2F* sigma_slope = new TH2F("sigma_slope", "d/dNPV Sigma(Offset) At Different Corrections", 100, 0, 4, 40, 0, 0.1);
	for(float i = 0; i < 4; i += 0.2) { //iterate through all correction permutations 
		graphInput(sigma_slope, i); //graph histograms for single sigma correction
	}
	sigma_slope->SetXTitle("Number of Sigma(Rho) Used for Correction");
	sigma_slope->SetYTitle("Slope For Linearizion Fit Line of NPV vs. Sigma(Offset)");
	TFile fp("sigma_slope.root", "RECREATE");
	sigma_slope->SetOption("COL");
	sigma_slope->Write("sigma_slope.root");
	delete sigma_slope;
	return 0;
}

/*creates a set of histograms to a file for a given correction*/
void graphInput(TH2F* sigma_slope, float numSigma)
{
	TH2F hist("hist", "Mean Offset Across Varying NPV", 11, 5, 60, 50, -30, 20); //temporary histogram for finding mean offset
	char buffer[100];
	float recon_pt;
	float truth_pt;
	int npv;
	float corr_quantized; //rho times area
	ifstream ifile;
	ifile.open("matches.txt", fstream::in);
	while(ifile.getline(buffer, 100, '&')) { //the '&' symbol signifies the end of a line in a properly formatted file
		if(sscanf(buffer, "%f %f %d %f", &recon_pt, &truth_pt, &npv, &corr_quantized) != 4) { //the recon pt, truth pt, and event npv are derived from each line of text
			break;
		}
		float corr_pt = recon_pt - (numSigma * corr_quantized);
		if(corr_pt > PT_RECON_MIN && truth_pt > PT_TRUTH_LOW && truth_pt < PT_TRUTH_HIGH) { //applying pt truth bin and ensuring pt recon have sufficient momenta
			hist.Fill(npv, corr_pt - truth_pt);
		}
	}
	//cout << "number of jets in vec = " << recon_vec.size() << endl;
	ifile.close(); //closing input file

	TFile* fout = new TFile("temp.root", "RECREATE"); //creating output file; previous contents erased on first iteration
	hist.FitSlicesY();
	TH1D* sigma_hist;
    gDirectory->GetObject("hist_2", sigma_hist);  //obtaining histogram of npv vs. sigma(offset)
    TF1* linear_fit = new TF1("linear_fit", "[0] + [1] * x", 5, 60);  //line used to calculated linear fit
    sigma_hist->Fit("linear_fit");

    Double_t params[2];
    linear_fit->GetParameters(params); //obtaining slope and y intercept for fit lines
    //cout << "Param 1 = " << params[1] << endl;
    sigma_slope->Fill(numSigma, params[1]);
    delete fout;
}
