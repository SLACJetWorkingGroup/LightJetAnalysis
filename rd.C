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
//#include <vector.h>
#include <string.h>

#define PT_RECON_MIN 4 //minimum threshold for a reconstructed jet to be considered as a real jet
#define PT_TRUTH_LOW	30 //low end pt bin most consistent truth jet data
#define PT_TRUTH_HIGH	50 //high end of bin for most consistent truth jet data


/*Prototypes*/
void graphInput(char* fileName, int numSigma);
inline int getLineColor(int numSigma);

/*Main script function*/
int rd() {
	for(int i = 0; i < 3; i++) { //iterate through all correction permutations 
		char buff[15];
		sprintf(buff, "matches_%d.txt", i); //parsing file name for differing sigma corrections
		graphInput(buff, i); //graph histograms for single sigma correction
	}
}

/*creates a set of histograms to a file for a given correction*/
void graphInput(char* fileName, int numSigma)
{
	TH2F hist("hist", "Mean Offset Across Varying NPV", 11, 5, 60, 50, -30, 20); //temporary histogram for finding mean offset
	hist.SetXTitle("Number of Primary Vertices");
	char buffer[100];
	float recon_pt;
	float truth_pt;
	int npv;
	ifstream ifile;
	ifile.open(fileName, fstream::in);
	while(ifile.getline(buffer, 100, '&')) { //the '&' symbol signifies the end of a line in a properly formatted file
		if(sscanf(buffer, "%f %f %d", &recon_pt, &truth_pt, &npv) != 3) { //the recon pt, truth pt, and event npv are derived from each line of text
			break;
		}
		if(truth_pt > PT_TRUTH_LOW && truth_pt < PT_TRUTH_HIGH) { //applying pt truth bin
			hist.Fill(npv, recon_pt - truth_pt);
		}
	}
	//cout << "number of jets in vec = " << recon_vec.size() << endl;
	ifile.close(); //closing input file

	TFile* fout = new TFile("output.root", getFileState(numSigma).c_str()); //creating output file; previous contents erased on first iteration
	hist.FitSlicesY();
	TH1D* sigma_hist;
    gDirectory->GetObject("hist_2", sigma_hist);  //obtaining histogram of npv vs. sigma(offset)
    sigma_hist->SetLineColor(getLineColor(numSigma));
    sigma_hist->SetYTitle("Sigma of Mean Offset");
    sigma_hist->SetXTitle("Number of Primary Vertices");
    TF1* linear_fit = new TF1("linear_fit", "[0] + [1] * x", 5, 60);  //line used to calculated linear fit
    sigma_hist->SetTitle("Gaussian Sigma of Offset Across NPV");
    sigma_hist->Fit("linear_fit");

    char outputName[10];
	sprintf(outputName, "sigma_%d", numSigma);
    sigma_hist->Write(outputName); //npv vs. sigma(offset) written to file (with fit line)

    char profileName[12];
    sprintf(profileName, "profile_%d", numSigma);
    hist.SetLineColor(getLineColor(numSigma));
    hist.ProfileX()->SetYTitle("Mean Offset (Pt Recon - Pt Truth)");
    hist.ProfileX()->Write(profileName); //npv vs. mean(offset) written to file
    
    //-----------
    //Double_t params[2];
    //linear_fit->GetParameters(params); //obtaining slope and y intercept for fit lines
    //-----------

	delete fout;
}

/*determines the line color of each set of histograms depending on the sigma correction*/
inline int getLineColor(int numSigma) 
{
	switch(numSigma) {
		case 0:
			return 4; //blue
		case 1:
			return 2; //red
		case 2:
			return 1; //black
		case 4:
			return 13; //gray
		default:
			return 9; //off blue
	}
}

/*ensures that the script will erase any previously existing contents in the output file prior to running the script*/
inline string getFileState(int numSigma) 
{
	if(numSigma == 0) {
		return "RECREATE"; //start from scratch
	}
	return "UPDATE"; //add to file
}

