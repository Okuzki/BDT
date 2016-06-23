//
// Code by Oskari Saarimaki
//
//
// triggr count on exit.
// Maybe different questions for different layers?
// Drift time?
//
//

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <math.h>

#include "TString.h"
#include "TMath.h"
#include "TRandom.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBasket.h"
#include "TH1D.h"
#include "TH2D.h"
//#include <boost/any.hpp>

const double edep_cut_forTest     = 5e-6; // GeV
const int from                    = 0;
const int to                      = from + 500;
const int nCells                  = 300;
const int trees                   = 20;
const int questions               = 4;
const double purity_check         = 0.5;
const double beta                 = 1.0;
const double eps                  = 1e-7;
const double maxTest              = 1e-4;
//const int updates                 = 2;

/*
class Hit {             // This turned out to be little unnecessary because it didn't work like I wanted.
	private:
    	double edep;
    	int layerID;
    	int cellID;
};

Hit::Hit() {};

Hit::setHit(double ed, int layer, int cell) {
	edep = ed;
	layerID = layer;
	cellID = cell;
};

double Hit::edep() { return edep; };
int Hit::layerID() { return layerID; };
int Hit::cellID()  { return cellID; };
*/

class Tester {
	public:
		void setCellEdep(const std::vector<double>& c);
		void setEdep(const std::vector<double>& c);
		void setHittype(const std::vector<int>& h);
		void setLayerIDs(const std::vector<int>& i);
		void setCellIDs(const std::vector<int>& j);
		void setTStart(const std::vector<double>& t);
		void notDone();
		bool test (const int hitNum, const int q, const std::vector<double>& WBranch, const double giniF, std::vector<double>& chosen_cut);
		double purity(const std::vector<double>& WBranch, const std::vector<int>& hittype);
		double gini(const std::vector<double>& WBranch, const std::vector<int>& hittype);
		double criterion(const double giniF, const double giniL, const double giniR);
		double dVecSum(const std::vector<double>& S);
		double err(const std::vector<double>& WBranch, const std::vector<int>& sig, const std::vector<int>& hittype);
		Tester();

	private:
		// For test
		std::vector<double>  cells_edep;
		std::vector<double>  hit_edep;
		std::vector<double>  tStart;
		std::vector<int>     hittype;
		std::vector<int>     done;
		std::vector<int>     thebest;
		std::vector<int>     layerIDs;
		std::vector<int>     cellIDs;
		std::vector<double>  holder;
		std::vector<double>  WL;
		std::vector<double>  WR;
};

Tester::Tester()                                  { done = std::vector<int>(questions,0);  thebest = std::vector<int>(questions,0); };
void Tester::setCellEdep(const std::vector<double>& c) { cells_edep = c; };
void Tester::setEdep(const std::vector<double>& c)     { hit_edep = c; };
void Tester::setHittype(const std::vector<int>& h)     { hittype = h; };
void Tester::setLayerIDs(const std::vector<int>& i)    { layerIDs = i; };
void Tester::setCellIDs(const std::vector<int>& j)     { cellIDs = j; };
void Tester::setTStart(const std::vector<double>& t)   { tStart = t; };
void Tester::notDone()                                 { for ( unsigned c1 = 0; c1 < done.size(); ++c1) done.at(c1) = 0; };

bool Tester::test(const int hitNum, const int q, const std::vector<double>& WBranch, const double giniF, std::vector<double>& chosen_cut) { // remove layerCut and add vector W and edep? OR edep could also be added to the Tester (EDIT: DONE AND DONE)
	if ( done.at(q) == 0 ){  // Remove the latter when adding other questions
		int c1               = 0;
		int c2               = 0;		
		double bestNum       = 0.0;
		double testCriterion = 0.0;
		int NBINS            = 0;
		double logBW         = 0.0;
		thebest.at(q) = -1;
		WL.clear();
		WR.clear();

		
		if ( q == 0 ) {
			NBINS         = (int)(maxTest/eps);
			logBW         = (log(maxTest)-log(eps))/(double)NBINS;

			for (c1 = 1; c1 < NBINS; ++c1) {
				for (c2 = 0; c2 < (int)WBranch.size(); ++c2) {
					if ( WBranch.at(c2) > 0.0 ) {
	
						if ( hit_edep.at(c2) < eps*exp((double)c1*logBW) ) {   // This does the actual comparing
							WL.push_back(WBranch.at(c2));
							WR.push_back(-1);
						}
						else {
							WL.push_back(-1);
							WR.push_back(WBranch.at(c2));
						}
					}
					else {
						WL.push_back(-1);
						WR.push_back(-1);
					}
				}
				testCriterion = criterion(giniF, gini(WL, hittype), gini(WR, hittype));
				/*
				if ( c1%1 == 0 ) std::cout << "eps*c1: " << setw(11) << eps*exp((double)c1*logBW) << " ";          // With these prints you can see exactly what is going through the tests
				if ( c1%1 == 0 ) std::cout << "testCriterion: " << setw(11) << testCriterion << " ";
				if ( c1%1 == 0 ) std::cout << "giniL: " << setw(11) << gini(WL, hittype) << " ";
				if ( c1%1 == 0 ) std::cout << "giniR: " << setw(11) << gini(WR, hittype) << std::endl;
				*/

				if ( testCriterion > bestNum ) {
					bestNum = testCriterion;
					thebest.at(q) = c1;
					//std::cout << "the best:         " << thebest << std::endl;
				}

				WL.clear();
				WR.clear();
			}




			done.at(q) = 1; 
			chosen_cut.at(q) = eps * exp((double)thebest.at(q) * logBW);
			//std::cout << "chosen edep cut:         " << chosen_cut.at(q) << std::endl;
			//std::cout << "chosen edep cut c1: " << thebest << std::endl;
		}
		else if ( q == 1) {

			for (c1 = 650; c1 < 1650; ++c1) {                  // From t = 650 till t = 1650 will be checked.
				for (c2 = 0; c2 < (int)WBranch.size(); ++c2) {
					if ( WBranch.at(c2) > 0.0 ) {
						if ( tStart.at(c2) < c1) {   // This does the actual comparing
							WL.push_back(WBranch.at(c2));
							WR.push_back(-1);
						}
						else {
							WL.push_back(-1);
							WR.push_back(WBranch.at(c2));
						}
					}
					else {
						WL.push_back(-1);
						WR.push_back(-1);
					}
				}
				testCriterion = criterion(giniF, gini(WL, hittype), gini(WR, hittype));
				/*
				if ( c1%1 == 0 ) std::cout << "eps*c1: " << setw(11) << eps*exp((double)c1*logBW) << " ";          // With these prints you can see exactly what is going through the tests
				if ( c1%1 == 0 ) std::cout << "testCriterion: " << setw(11) << testCriterion << " ";
				if ( c1%1 == 0 ) std::cout << "giniL: " << setw(11) << gini(WL, hittype) << " ";
				if ( c1%1 == 0 ) std::cout << "giniR: " << setw(11) << gini(WR, hittype) << endl;
				*/

				if ( testCriterion > bestNum ) {
					bestNum = testCriterion;
					thebest.at(q) = c1;
					//std::cout << "the best:         " << thebest << endl;
				}

				WL.clear();
				WR.clear();
			}

			done.at(q) = 1; 
			chosen_cut.at(q) = thebest.at(q);
		}
		else if ( q == 2) {
			NBINS         = (int)(maxTest/eps);
			logBW         = (log(maxTest)-log(eps))/(double)NBINS;
			int help          = 0;

			for (c1 = 1; c1 < NBINS; ++c1) {
				for (c2 = 0; c2 < (int)WBranch.size(); ++c2) {
					if ( WBranch.at(c2) > 0.0 ) {
						help = cellIDs.at(c2);
	
						if ( cells_edep.at( help ) < eps*exp((double)c1*logBW) ) {   // This does the actual comparing
							WL.push_back(WBranch.at(c2));
							WR.push_back(-1);
						}
						else {
							WL.push_back(-1);
							WR.push_back(WBranch.at(c2));
						}
					}
					else {
						WL.push_back(-1);
						WR.push_back(-1);
					}
				}
				testCriterion = criterion(giniF, gini(WL, hittype), gini(WR, hittype));
				/*
				if ( c1%1 == 0 ) std::cout << "eps*c1: " << setw(11) << eps*exp((double)c1*logBW) << " ";          // With these prints you can see exactly what is going through the tests
				if ( c1%1 == 0 ) std::cout << "testCriterion: " << setw(11) << testCriterion << " ";
				if ( c1%1 == 0 ) std::cout << "giniL: " << setw(11) << gini(WL, hittype) << " ";
				if ( c1%1 == 0 ) std::cout << "giniR: " << setw(11) << gini(WR, hittype) << endl;
				*/

				if ( testCriterion > bestNum ) {
					bestNum = testCriterion;
					thebest.at(q) = c1;
					//std::cout << "the best:         " << thebest << endl;
				}

				WL.clear();
				WR.clear();


			}

			done.at(q) = 1; 
			chosen_cut.at(q) = eps * exp((double)thebest.at(q) * logBW);
		}
		else {
			for (c1 = 1; c1 < 18; ++c1) {  // 18 because there are 18 layers
				for (c2 = 0; c2 < (int)layerIDs.size(); ++c2) {
					if ( WBranch.at(c2) > 0.0 ) {
		
						if ( layerIDs.at(c2) < c1 ) {   // This does the actual comparing
							WL.push_back(WBranch.at(c2));
							WR.push_back(-1);
						}
						else {
							WL.push_back(-1);
							WR.push_back(WBranch.at(c2));
						}
					}
					else {
						WL.push_back(-1);
						WR.push_back(-1);
					}
				}
				testCriterion = criterion(giniF, gini(WL, hittype), gini(WR, hittype));
				if ( testCriterion > bestNum ) {
					bestNum = testCriterion;
					thebest.at(q) = c1;
				}
				WL.clear();
				WR.clear();
			}
			done.at(q) = 1;
			chosen_cut.at(q) = thebest.at(q);
			//std::cout << "chosen layer cut:  " << chosen_cut.at(q) << std::endl;
			//std::cout << "chosen thebest:    " << thebest.at(q) << std::endl;
			
		}
	}

	if      (q == 0) return (hit_edep.at(hitNum) < chosen_cut.at(0) && hit_edep.at(hitNum) > 0);                                 // Edep cut check
    else if (q == 1) return (tStart.at(hitNum) < chosen_cut.at(1) && tStart.at(hitNum) > 0);                                     // tStart comparison
    else if (q == 2) return (cells_edep.at(cellIDs.at(hitNum)) < chosen_cut.at(2) && cells_edep.at(cellIDs.at(hitNum)) > 0);     // Cell edep cut check
    else if (q == 3) return (layerIDs.at(hitNum) < chosen_cut.at(3) && layerIDs.at(hitNum) > -1);                                // Layer cut check
	return false;

	// Maybe distance from one string hit to another?
	// I should definitely try counting all the hits on one layer on one event and use those!
}

/* 
double Tester::testResult(const std::vector<double>& T, TH1D hSignalRetentionE, TH1D (hBackgroundRejectionE) { // This had to be moved to main code, because the function couldn't recieve TH1Ds for reasons unknown.
	std::vector<double> WEqual((int)T.size(), 1.0/(int)T.size());
	double giniF = gini(WEqual, hittype);
	double resultCriterion = 0.0;
	double bestResultCrit = 0.0;
	double bestCutPoint = -1.0;
	int nSignalTestCorrect = 0;
	int nBackgroundTestCorrect = 0;

	int nSignal = 0;
	int nBackground = 0;
	for (int cp1 = 0; cp1 < hittype.size(); ++cp1) {
		if (hittype.at(cp1) > 0) nSignal++;
		else nBackground++;
	}



	WL.clear();
	WR.clear();

	for (double i = -1.0; i < 1.0; i += 0.001) {
		for (int c = 0; c < (int)T.size(); ++c) {
			if ( T.at(c) < i ) {   // This does the actual comparing
				WL.push_back(WEqual.at(c));
				WR.push_back(-1);
			}
			else {
				WL.push_back(-1);
				WR.push_back(WEqual.at(c));
			}
		}
		resultCriterion = criterion(giniF, gini(WL, hittype), gini(WR, hittype));
		if ( resultCriterion > bestResultCrit ) {
			bestResultCrit = resultCriterion;
			bestCutPoint = i;
		}
		for (int j = 0; j < WL.size(); ++j) {
			if (WL.at(j) > 0.0 && hittype.at(j) < 1) nBackgroundTestCorrect++;
			else if (WR.at(j) > 0.0 && hittype.at(j) > 0) nSignalTestCorrect++;
			else std::cout << "W = 0 in testResult." << endl;
		}
		hSignalRetentionE->Fill((double)nSignalTestCorrect/(double)nSignal);
		hBackgroundRejectionE->Fill((double)nBackgroundTestCorrect/(double)nBackground);

		nSignalTestCorrect = 0;
		nBackgroundTestCorrect = 0;
		WL.clear();
		WR.clear();
	}
	return bestCutPoint;
}
*/


/*
void Tester::cutting(const Tester tester, const std::vector<int> chosenQ, const int usedQuestions, const std::vector<double>& WTemp, const double giniF, std::vector<double> giniL, std::vector<double> giniR, std::vector<double> criterion, std::vector<double> PL, std::vector<double> PR, std::vector<double> chosen_cut, std::vector<double> WR0, std::vector<double> WL0, std::vector<double> WR1, std::vector<double> WL1, std::vector<double> WR2, std::vector<double> WL2, std::vector<double> WR3, std::vector<double> WL3) {
	testWL.clear();
	testWR.clear();
	passMainLoop = false;
	//std::cout << "DesdingSize " << chosenQ.size() << endl;
	for (c1C = 0; c1C < questions; ++c1C ){
		for (c2C = 0; c2C < usedQuestions; ++c2C) if ( c1C == chosenQ[c2C] ) passMainLoop = true;        // Because we don't want to check the same question twice.
		
		if ( !passMainLoop ) {

			for (c2C = 0; c2C < hit_(int)edep.size(); ++c2C) {
				if ( WTemp.at(c2C) > 0.0 ) {
					//hitHolder.setHit(hit_edep.at(c2C), layerIDs.at(c2C), cellIDs.at(c2C));
					if ( tester.test(c2C, c1C, WTemp, giniF, chosen_cut )) {
						testWL.push_back(WTemp.at(c2C));
						testWR.push_back(-WTemp.at(c2C));
						
					}
					else {
						testWR.push_back(WTemp.at(c2C));
						testWL.push_back(-WTemp.at(c2C));
					}
				}
				else {
					testWL.push_back(WTemp.at(c2C));  // Note that these are negative numbers which are pushed back.
					testWR.push_back(WTemp.at(c2C));
				}
			}
			
			giniL.at(c1C) = tester.gini(testWL, hittype);
			giniR.at(c1C) = tester.gini(testWR, hittype);
			criterion.at(c1C) = tester.criterion(giniF, giniL.at(c1C), giniR.at(c1C));
			PL.at(c1C) = tester.purity(testWL, hittype);
			PR.at(c1C) = tester.purity(testWR, hittype);
			if ( c1C == 0 ) {
				WR0.swap(testWR);
				WL0.swap(testWL);
			} else if ( c1C == 1 ) {
				WR1.swap(testWR);
				WL1.swap(testWL);
			} else if ( c1C == 2 ) {
				WR2.swap(testWR);
				WL2.swap(testWL);
			} else if ( c1C == 3 ) {
				WR3.swap(testWR);
				WL3.swap(testWL);
			}
			testWL.clear();
			testWR.clear();	
		}
		passMainLoop = false;
	}
}
*/

double Tester::purity(const std::vector<double>& WBranch, const std::vector<int>& hittype) {
	double Ws = 0.0;
	double Wb = 0.0;
	double P = 0.0;
	double helpP = 0.0;

	for (int c1P = 0; c1P < (int)WBranch.size(); ++c1P ) {
		helpP = WBranch.at(c1P);
		if ( helpP > 0.0 ) {
			if (hittype.at(c1P) > 0) Ws += helpP;
			else Wb += helpP;
		}
	}
	if (Ws == 0 && Wb == 0 ) P = 0.5;   // This is for safety.
	else P = Ws / (Ws + Wb);

	return P;
}

double Tester::gini(const std::vector<double>& WBranch, const std::vector<int>& hittype) {
	double Pgini = purity(WBranch, hittype);

	return (dVecSum(WBranch)*Pgini*(1-Pgini));

}

double Tester::criterion(const double giniF, const double giniL, const double giniR) {
	return giniF - giniL - giniR;
}

double Tester::dVecSum(const std::vector<double>& S) {
	double result = 0.0;
	double helpS = 0.0;
	for (int c1S = 0; c1S < (int)S.size(); ++c1S ) {
		helpS = S.at(c1S);
		if ( helpS > 0.0 ) result += helpS;
	}
	return result;
}

double Tester::err(const std::vector<double>& WBranch, const std::vector<int>& sig, const std::vector<int>& hittype) {
	double resultE = 0.0;
	double helpE = 0.0;
	bool question = true;
	for (int c1E = 0; c1E < (int)WBranch.size(); ++c1E ) {	
		if ( hittype.at(c1E) > 0 && sig.at(c1E) > 0 ) question = false;
		else if ( hittype.at(c1E) < 1 && !(sig.at(c1E) > 0) ) question = false; 
		else question = true;

		helpE = WBranch.at(c1E);
		if ( helpE > 0.0 && question ) resultE += helpE;
	}

	return resultE / dVecSum(WBranch);
}

int main() {

	TFile * sfile = new TFile("file:/home/oskari/MyThing/Rootfiles/signal.root", "READ");

	TTree * stree = (TTree*) sfile->Get("tree");
	//TFile * of = new TFile("output.root","RECREATE");
	const int binTest = 1;
	const int resolution = 200;
	TH1D *hResult = new TH1D("hResult","The results for background particles", binTest*resolution, -binTest, binTest);
	TH1D *hResultSignal  = new TH1D("hResultSignal","The results for signal particles", binTest*resolution, -binTest, binTest);
	TH1D *hPurity  = new TH1D("hPurity","The purity of the leaves", binTest*resolution, 0, binTest);
	//TH1D *hSignalRetentionE  = new TH1D("hSignalRetentionE","Signal retention efficiency", binTest*resolution, 0, binTest);
	//TH1D *hBackgroundRejectionE  = new TH1D("hBackgroundRejectionE","Background rejection efficiency", binTest*resolution, 0, binTest);

	TH2D *hEfficiency  = new TH2D("hEfficiency","Efficiency", binTest*resolution, 0, binTest, binTest*resolution, 0, binTest);

	int CdcCell_nHits = 0;
	std::vector<int> * CdcCell_layerID = 0;
	std::vector<int> * CdcCell_cellID = 0;
	std::vector<double> * CdcCell_edep = 0;
	std::vector<double> * CdcCell_stepL = 0;
	std::vector<double> * CdcCell_driftD = 0;
	std::vector<double> * CdcCell_driftT = 0;
	std::vector<double> * CdcCell_tstart = 0;
	std::vector<int> * CdcCell_posflag = 0;
	std::vector<int> * CdcCell_nPair = 0;
	std::vector<double> * CdcCell_t = 0;
	std::vector<double> * CdcCell_px = 0;
	std::vector<double> * CdcCell_py = 0;
	std::vector<double> * CdcCell_pz = 0;
	std::vector<double> * CdcCell_x = 0;
	std::vector<double> * CdcCell_y = 0;
	std::vector<double> * CdcCell_z = 0;
	std::vector<double> * CdcCell_wx = 0;
	std::vector<double> * CdcCell_wy = 0;
	std::vector<double> * CdcCell_wz = 0;
	std::vector<int> * CdcCell_hittype = 0;
	int M_nHits = 0;
	std::vector<std::string> * M_volName = 0;
	std::vector<int> * M_volID = 0;
	std::vector<double> * M_edep = 0;
	std::vector<double> * M_stepL = 0;
	std::vector<double> * M_t = 0;
	std::vector<double> * M_px = 0;
	std::vector<double> * M_py = 0;
	std::vector<double> * M_pz = 0;
	std::vector<double> * M_x = 0;
	std::vector<double> * M_y = 0;
	std::vector<double> * M_z = 0;
	std::vector<int> * M_hittype = 0;

	stree->SetBranchAddress("CdcCell_nHits",&CdcCell_nHits);
	stree->SetBranchAddress("CdcCell_layerID",&CdcCell_layerID);
	stree->SetBranchAddress("CdcCell_cellID",&CdcCell_cellID);
	stree->SetBranchAddress("CdcCell_edep",&CdcCell_edep);
	stree->SetBranchAddress("CdcCell_stepL",&CdcCell_stepL);
	stree->SetBranchAddress("CdcCell_driftD",&CdcCell_driftD);
	stree->SetBranchAddress("CdcCell_driftT",&CdcCell_driftT);
	stree->SetBranchAddress("CdcCell_tstart",&CdcCell_tstart);
	stree->SetBranchAddress("CdcCell_posflag",&CdcCell_posflag);
	stree->SetBranchAddress("CdcCell_nPair",&CdcCell_nPair);
	stree->SetBranchAddress("CdcCell_t",&CdcCell_t);
	stree->SetBranchAddress("CdcCell_px",&CdcCell_px);
	stree->SetBranchAddress("CdcCell_py",&CdcCell_py);
	stree->SetBranchAddress("CdcCell_pz",&CdcCell_pz);
	stree->SetBranchAddress("CdcCell_x",&CdcCell_x);
	stree->SetBranchAddress("CdcCell_y",&CdcCell_y);
	stree->SetBranchAddress("CdcCell_z",&CdcCell_z);
	stree->SetBranchAddress("CdcCell_wx",&CdcCell_wx);
	stree->SetBranchAddress("CdcCell_wy",&CdcCell_wy);
	stree->SetBranchAddress("CdcCell_wz",&CdcCell_wz);
	stree->SetBranchAddress("CdcCell_hittype",&CdcCell_hittype);
	stree->SetBranchAddress("M_nHits",&M_nHits);
	stree->SetBranchAddress("M_volName",&M_volName);
	stree->SetBranchAddress("M_volID",&M_volID);
	stree->SetBranchAddress("M_edep",&M_edep);
	stree->SetBranchAddress("M_stepL",&M_stepL);
	stree->SetBranchAddress("M_t",&M_t);
	stree->SetBranchAddress("M_px",&M_px);
	stree->SetBranchAddress("M_py",&M_py);
	stree->SetBranchAddress("M_pz",&M_pz);
	stree->SetBranchAddress("M_x",&M_x);
	stree->SetBranchAddress("M_y",&M_y);
	stree->SetBranchAddress("M_z",&M_z);
	stree->SetBranchAddress("M_hittype",&M_hittype);

	long c1 = 0;
	long c2 = 0;
	long c3 = 0;

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Testing sequence begin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	int nhitnoise_total = 0;
	int nhitnoise_edepcut = 0;
	int nhit_total = 0;
	int nhit_edepcut = 0;
	int ntrack_total = 0;
	//int nall_tracks = stree->GetEntries();

	for ( c1 = from; c1 < to; ++c1){
		stree->GetEntry(c1);
		for ( c2 = 0; c2 < CdcCell_nHits; ++c2){
			if ((*CdcCell_hittype).at(c2) !=1){ 
				nhitnoise_total++;
				if ((*CdcCell_edep).at(c2)<edep_cut_forTest) nhitnoise_edepcut++;
			}
			else {
				nhit_total++;
				if ((*CdcCell_edep).at(c2)<edep_cut_forTest) nhit_edepcut++;
			}
		}
		ntrack_total++;
	}

	std::cout << "Total tree tracks:          " << ntrack_total << std::endl;
	std::cout << "Total noise hits:           " << nhitnoise_total << std::endl;
	std::cout << "Above edep cut noice hits:  " << nhitnoise_edepcut << std::endl;
	std::cout << "Total signal hits:          " << nhit_total << std::endl;
	std::cout << "Above edep cut signal hits: " << nhit_edepcut << std::endl;




	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Testing sequence end ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//std::stringstream name;
	std::vector<double> cells_edep(nCells,0);
	std::vector<double> edep;
	std::vector<int> hittype;
	std::vector<int> layerID;
	std::vector<int> cellID;
	std::vector<double> tStart;
	//int nall_tracks = stree->GetEntries(); 
	std::vector<double> W;                            // Weights
	std::vector<double> alpha;                        // Alpha for each tree
	std::vector<double> T;                            // Score for each hit
	std::vector<double> PR;                           // Purity for the right son
	std::vector<double> PL;                           // Purity for the left son
	std::vector<double> giniL;                        // ginis for the left son
	std::vector<double> giniR;                        // ginis for the right son
	std::vector<double> WL0;                          // weights for the left son Q0
	std::vector<double> WR0;                          // weights for the right son Q0
	std::vector<double> WL1;                          // weights for the left son Q1
	std::vector<double> WR1;                          // weights for the right son Q1
	std::vector<double> WL2;                          // weights for the left son Q2
	std::vector<double> WR2;                          // weights for the right son Q2
	std::vector<double> WL3;                          // weights for the left son Q3
	std::vector<double> WR3;                          // weights for the right son Q3
	std::vector<double> WLC0;                         // weights for the chosen left son round zero
	std::vector<double> WRC0;                         // weights for the chosen right son round zero
	std::vector<double> WLC1;                         // weights for the chosen left son round one
	std::vector<double> WRC1;                         // weights for the chosen right son round one
	std::vector<double> WLC2;                         // weights for the chosen left son round two
	std::vector<double> WRC2;                         // weights for the chosen right son round two
	std::vector<double> WLC3;                         // weights for the chosen left son round three
	std::vector<double> WRC3;                         // weights for the chosen right son round three
	std::vector<double> WCTemporary;                  // temporary weights
	std::vector<double> testWL;                       // temporary weights for the left son
	std::vector<double> testWR;                       // temporary weights for the right son
	std::vector<double> criterion;                    // criterion for a question
	std::vector<int>    chosenQ(questions,-1);        // For remembering which questions were used and on what order
	std::vector<int>    qWhere;                       // For remembering where the questions were asked
	std::vector<int>    treeStructure(31, -100);      // For saving up the structure of the tree
	std::vector<int>    sig;                          // For marking the hit as signal or not
	std::vector<double> chosen_cut(questions,-1.0);   // To remember what limits were chosen for each question
	//std::vector<std::vector<double>> hits(4, std::vector<double>); // Information for hits
	double giniF                = 0.0;                       // gini for the father
	double weightBegin          = 0.0;                       // Weight at the very beginning
	double err                  = 0.0;
	double score                = 0.0;
	double errResult            = 0.0;
	double dVecSumResult        = 0.0;
	double trueSigRatio         = 0.0;
	double alphaSummed          = 0.0;
	double theBestCutForResults = 0.0;
	double alphaSave            = 0.0;
	int resultQ                 = 0;
	int holder                  = 0;
	int hitWasCorrect           = 0;
	int signalCount             = 0;
	int signalTrueCount         = 0;
	int usedQ                   = 0;
	bool holdThisTruth          = false;
	bool showMessage            = false;
	bool passMainLoop           = false;
	Tester tester;
	//Hit hitHolder;
	//int entries = stree->GetEntries();

	// For progress bar
	int progressBarCount = 0;
	double progressUpdate = 1.0/16.0;
	int iftoosmallcheck = 0;
	std::stringstream buff;
	// ----------------

	TTree *tree  = new TTree("tree","tree");
	//TBranch *bTreeStructure = TTree::Branch("treeStructure",&treeStructure);
	//tree->AddBranchToCache(bTreeStructure);
	tree->Branch("treeStructure",&treeStructure);
	tree->Branch("alphaSave",&alphaSave);
	tree->Branch("chosen_cut",&chosen_cut);
	tree->Branch("PL",&PL);
	tree->Branch("PR",&PR);
	tree->Branch("chosenQ",&chosenQ);
	//tree->Branch("theBestCutForResults",&theBestCutForResults);


	for ( c1 = from; c1 < to; ++c1) {
		stree->GetEntry(c1);

		for ( c2 = 0; c2 < (int)(*CdcCell_edep).size(); ++c2)
		{
			edep.push_back((*CdcCell_edep).at(c2));
			hittype.push_back((*CdcCell_hittype).at(c2));
			layerID.push_back((*CdcCell_layerID).at(c2));
			cellID.push_back((*CdcCell_cellID).at(c2));
			tStart.push_back((*CdcCell_tstart).at(c2));

		}



		for ( int cellIDc = 0; cellIDc < nCells; ++cellIDc) {
			for (c2 = 0; c2 < (int)(*CdcCell_cellID).size(); ++c2) {
				if ((*CdcCell_cellID).at(c2) == cellIDc) {
					cells_edep.at(cellIDc) += (*CdcCell_edep).at(c2);
				}
			}
		}


		// ----------------- Progress bar -----------------
		iftoosmallcheck = (int)(progressUpdate*(to - from));
		if ( iftoosmallcheck < 1 ) iftoosmallcheck = 1;
		if ( c1%iftoosmallcheck == 0 ) {
			buff.str("");
			buff.clear();
			buff << "\33[2K |";  // \33[2K removes the current line.
			for ( c2 = 0; c2 < progressBarCount; ++c2) buff << "=";
			//if ( c1 != to - 1 ) {
			//	buff << ">";
			//	c2++;
			//}
			for ( ; c2 < 1.0/progressUpdate ; ++c2) buff << "-";
			buff << "|\r";
			std::cout << buff.str();
			std::cout.flush();

			progressBarCount++;
		}
		// ------------------------------------------------



	}
	progressBarCount = 0;	



	// ~~~~~~~~~~~~~~~~~~~ Here begins the boosted decision tree ~~~~~~~~~~~~~~~~~~~

	weightBegin = 1.0/(int)edep.size();
	for ( c1 = 0; c1 < (int)edep.size(); ++c1) { // This loop is for initial weighting.
		W.push_back(weightBegin);
	}

	std::cout << "\33[2KPurity:                     "<< tester.purity(W, hittype) << std::endl;
	std::cout << "Weight sum:                 "<< tester.dVecSum(W) << std::endl; 
	std::cout << "Gini father:                "<< tester.gini(W, hittype) << std::endl; 

	tester.setCellEdep(cells_edep);
	tester.setEdep(edep);
	tester.setHittype(hittype);
	tester.setLayerIDs(layerID);
	tester.setCellIDs(cellID);
	tester.setTStart(tStart);
	T.resize(hittype.size(), -1);


	for (c1 = 0; c1 < trees; ++c1) {
		PR.clear();
		PL.clear();
		giniL.clear();
		giniR.clear();
		WL0.clear();
		WR0.clear();
		WL1.clear();
		WR1.clear();
		WL2.clear();
		WR2.clear();
		WL3.clear();
		WR3.clear();
		WLC0.clear();
		WRC0.clear();
		WLC1.clear();
		WRC1.clear();
		WLC2.clear();
		WRC2.clear();
		WLC3.clear();
		WRC3.clear();
		WCTemporary.clear();
		testWL.clear();
		testWR.clear();
		criterion.clear();
		qWhere.clear();
		sig.clear();
		for ( c2 = 0; c2 < (int)treeStructure.size(); ++c2) treeStructure.at(c2) = -100;
		giniF            = 0.0;
		weightBegin      = 0.0;
		err              = 0.0;
		score            = 0.0;
		errResult        = 0.0;
		dVecSumResult    = 0.0;
		trueSigRatio     = 0.0;
		resultQ          = 0;
		holder           = 0;
		hitWasCorrect    = 0;
		signalCount      = 0;
		signalTrueCount  = 0;
		usedQ            = 0;
		holdThisTruth    = false;
		showMessage      = ((c1 == 0) || (c1 == trees - 1));


	 
		for (c2 = 0; c2 < questions; ++c2) {
			chosenQ.at(c2) = -1;
			chosen_cut.at(c2) = -1;
		}

		
		tester.notDone();
			
		giniF = tester.gini(W, hittype);

		WCTemporary = W;

		// ----------------- Progress bar -----------------
		buff.str("");
		buff.clear();
		buff << "\33[2K |";  // \33[2K removes the current line.
		for ( c3 = 0; c3 < progressBarCount; ++c3) buff << "=";
		//if ( counter != entries - 1 ) {
		//	buff << ">";
		//	c3++;
		//}
		for ( ; c3 < 29 ; ++c3) buff << "-";  // 29
		buff << "|\r";
		std::cout << buff.str();
		std::cout.flush();
		progressBarCount++;
		// ------------------------------------------------

		for (c2 = 0; c2 < questions; ++c2) {
			for (c3 = 0; c3 < (int)edep.size(); ++c3) {
				//hitHolder.setHit(edep.at(c3), layerID.at(c3), cellID.at(c3));
				if ( tester.test(c3, c2, WCTemporary, giniF, chosen_cut )) {
					testWL.push_back(W.at(c3));
					testWR.push_back(-W.at(c3));
				}
				else {
					testWR.push_back(W.at(c3));
					testWL.push_back(-W.at(c3));
				}
			}
			giniL.push_back(tester.gini(testWL, hittype));
			giniR.push_back(tester.gini(testWR, hittype));
			criterion.push_back(tester.criterion(giniF, giniL.at(c2), giniR.at(c2)));
			PL.push_back(tester.purity(testWL, hittype));
			PR.push_back(tester.purity(testWR, hittype));
			if ( c2 == 0 ) {
				WR0.swap(testWR);
				WL0.swap(testWL);
			} else if ( c2 == 1 ) {
				WR1.swap(testWR);
				WL1.swap(testWL);
			} else if ( c2 == 2 ) {
				WR2.swap(testWR);
				WL2.swap(testWL);
			} else if ( c2 == 3 ) {
				WR3.swap(testWR);
				WL3.swap(testWL);
			}
			testWL.clear();
			testWR.clear();

			// ----------------- Progress bar -----------------
			progressBarCount++;
			buff.str("");
			buff.clear();
			buff << "\33[2K |";  // \33[2K removes the current line.
			for ( c3 = 0; c3 < progressBarCount; ++c3) buff << "=";
			//if ( counter != entries - 1 ) {
			//	buff << ">";
			//	c3++;
			//}
			for ( ; c3 < 29 ; ++c3) buff << "-";
			buff << "|\r";
			std::cout << buff.str();
			std::cout.flush();
			progressBarCount++;
			// ------------------------------------------------
		}


		if ( showMessage ) {                         // For testing purposes
				std::cout << "\33[2KThe first splitting:" << std::endl;
				std::cout << "Criterion" << 0 << ": " << std::setw(11) << criterion.at(0) << " | ";
				std::cout << "PurityL"   << 0 << ": " << std::setw(11) << PL.at(0)        << " | ";
				std::cout << "PurityR"   << 0 << ": " << std::setw(11) << PR.at(0)        << " | ";
				std::cout << "GiniL"     << 0 << ": " << std::setw(11) << giniL.at(0)     << " | ";
				std::cout << "GiniR"     << 0 << ": " << std::setw(11) << giniR.at(0)     << std::endl;
		}
		chosenQ.at(0) = 0;			
		for ( c2 = 1; c2 < questions; ++c2 ) {
			if ( showMessage ) {                         // For testing purposes
				std::cout << "Criterion" << c2 << ": " << std::setw(11) << criterion.at(c2) << " | ";
				std::cout << "PurityL"   << c2 << ": " << std::setw(11) << PL.at(c2)        << " | ";
				std::cout << "PurityR"   << c2 << ": " << std::setw(11) << PR.at(c2)        << " | ";
				std::cout << "GiniL"     << c2 << ": " << std::setw(11) << giniL.at(c2)     << " | ";
				std::cout << "GiniR"     << c2 << ": " << std::setw(11) << giniR.at(c2)     << std::endl;
			}
			if (criterion.at(c2) > criterion.at(chosenQ.at(0))) chosenQ.at(0) = c2;
		}
			//std::cout << std::endl;
		treeStructure.at(0) = chosenQ.at(0);

		criterion.at(chosenQ.at(0)) = -1;
		if ( giniL.at(chosenQ.at(0)) > giniR.at(chosenQ.at(0)) ) {
			giniF = giniL.at(chosenQ.at(0));
			giniL.at(chosenQ.at(0)) = -1;    // Because these numbers are unused later on.
			PL.at(chosenQ.at(0)) = -1;
			qWhere.push_back(1);
		}
		else {
			giniF = giniR.at(chosenQ.at(0));
			giniR.at(chosenQ.at(0)) = -1;    
			PR.at(chosenQ.at(0)) = -1;
			qWhere.push_back(0);
		}

		if ( chosenQ.at(0) == 0 ) {        // This might not be needed if replaced with something like below.
			WRC0 = WR0;
			WLC0 = WL0;
		} else if ( chosenQ.at(0) == 1 ) {
			WRC0 = WR1;
			WLC0 = WL1;				
		} else if ( chosenQ.at(0) == 2 ) {
			WRC0 = WR2;
			WLC0 = WL2;
		} else if ( chosenQ.at(0) == 3 ) {
			WRC0 = WR3;
			WLC0 = WL3;
		}
		/*
		if ( chosenQ.at(0) == 0 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR0;
			else WCTemporary = WL0;
		} else if ( chosenQ.at(0) == 1 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR1;
			else WCTemporary = WL1;				
		} else if ( chosenQ.at(0) == 2 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR2;
			else WCTemporary = WL2;
		} else if ( chosenQ.at(0) == 3 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR3;
			else WCTemporary = WL3;
		}
		*/
		if ( c1 == 0 || c1 == trees - 1 ) {
			std::cout << "ChosenQ.at(0)" << ":    " << chosenQ.at(0) << std::endl;
			std::cout << "The new giniF" << ": " << giniF << std::endl;
		}

		// ----------------- Progress bar -----------------
		progressBarCount++;
		buff.str("");
		buff.clear();
		buff << "\33[2K |";  // \33[2K removes the current line.
		for ( c3 = 0; c3 < progressBarCount; ++c3) buff << "=";
		//if ( counter != entries - 1 ) {
		//	buff << ">";
		//	c3++;
		//}
		for ( ; c3 < 29 ; ++c3) buff << "-";
		buff << "|\r";
		std::cout << buff.str();
		std::cout.flush();
		// ------------------------------------------------

		//giniL.clear();
		//giniR.clear();
		//PL.clear();
		//PR.clear();
		//criterion.clear();
		tester.notDone();

		if ( qWhere.at(0) == 0 ) WCTemporary.swap(WRC0);
		else WCTemporary.swap(WLC0);
		for (c2 = 0; c2 < questions ; ++c2 ){
			if (c2 == chosenQ.at(0)) continue;
			for (c3 = 0; c3 < (int)edep.size(); ++c3) {
				if ( WCTemporary.at(c3) > 0.0 ) {
					//hitHolder.setHit(edep.at(c3), layerID.at(c3), cellID.at(c3));
					if ( tester.test(c3, c2, WCTemporary, giniF, chosen_cut )) {
						testWL.push_back(W.at(c3));
						testWR.push_back(-W.at(c3));
						
					}
					else {
						testWR.push_back(W.at(c3));
						testWL.push_back(-W.at(c3));
					}
				}
				else {
					testWL.push_back(-W.at(c3));
					testWR.push_back(-W.at(c3));
				}
			}
				
			giniL.at(c2) = tester.gini(testWL, hittype);
			giniR.at(c2) = tester.gini(testWR, hittype);
			criterion.at(c2) = tester.criterion(giniF, giniL.at(c2), giniR.at(c2));
			PL.at(c2) = tester.purity(testWL, hittype);
			PR.at(c2) = tester.purity(testWR, hittype);
			if ( c2 == 0 ) {
				WR0.swap(testWR);
				WL0.swap(testWL);
			} else if ( c2 == 1 ) {
				WR1.swap(testWR);
				WL1.swap(testWL);
			} else if ( c2 == 2 ) {
				WR2.swap(testWR);
				WL2.swap(testWL);
			} else if ( c2 == 3 ) {
				WR3.swap(testWR);
				WL3.swap(testWL);
			}
			testWL.clear();
			testWR.clear();	

			// ----------------- Progress bar -----------------
			progressBarCount++;
			buff.str("");
			buff.clear();
			buff << "\33[2K |";  // \33[2K removes the current line.
			for ( c3 = 0; c3 < progressBarCount; ++c3) buff << "=";
			//if ( counter != entries - 1 ) {
			//	buff << ">";
			//	c3++;
			//}
			for ( ; c3 < 29 ; ++c3) buff << "-";
			buff << "|\r";
			std::cout << buff.str();
			std::cout.flush();
			progressBarCount++;
			// ------------------------------------------------
		}

		if ( showMessage ) {        // For testing purposes
				std::cout << "\33[2KThe second splitting:" << std::endl;
				std::cout << "Criterion" << 0 << ": " << std::setw(11) << criterion.at(0) << " | ";
				std::cout << "PurityL"   << 0 << ": " << std::setw(11) << PL.at(0)        << " | ";
				std::cout << "PurityR"   << 0 << ": " << std::setw(11) << PR.at(0)        << " | ";
				std::cout << "GiniL"     << 0 << ": " << std::setw(11) << giniL.at(0)     << " | ";
				std::cout << "GiniR"     << 0 << ": " << std::setw(11) << giniR.at(0)     << std::endl;
		}
		chosenQ.at(1) = 0;			
		for ( c2 = 1; c2 < questions; ++c2 ) {
			if ( showMessage ) {               // For testing purposes
				std::cout << "Criterion" << c2 << ": " << std::setw(11) << criterion.at(c2) << " | ";
				std::cout << "PurityL"   << c2 << ": " << std::setw(11) << PL.at(c2)        << " | ";
				std::cout << "PurityR"   << c2 << ": " << std::setw(11) << PR.at(c2)        << " | ";
				std::cout << "GiniL"     << c2 << ": " << std::setw(11) << giniL.at(c2)     << " | ";
				std::cout << "GiniR"     << c2 << ": " << std::setw(11) << giniR.at(c2)     << std::endl;
			}
			if (criterion.at(c2) > criterion.at(chosenQ.at(1))) chosenQ.at(1) = c2;
		}
		criterion.at(chosenQ.at(1)) = -1;

		/*
		holdThisTruth = (giniL.at(chosenQ.at(1)) > giniR.at(chosenQ.at(1)))
		
		for (c2 = 0; c2 < qWhere.size(); ++c2) {   // Something like this could be planned if one wants to optimise
			if ( holdThisTruth ) {
				if ( qWhere.at(c2)  )
			}
		}
		*/
		if ( giniL.at(chosenQ.at(0)) == -1.0 ) treeStructure.at(1) = chosenQ.at(1);
		else if ( giniR.at(chosenQ.at(0)) == -1.0 ) treeStructure.at(2) = chosenQ.at(1);
		else std::cout << "Error------------------------------------" << std::endl;

		if ( giniL.at(chosenQ.at(1)) > giniR.at(chosenQ.at(1)) ) {
			//treeStructure.at(1) = chosenQ.at(1);
			if ( qWhere.at(0) == 0 ) {
				if ( giniL.at(chosenQ.at(1)) > giniL.at(chosenQ.at(0)) ) {
					giniF = giniL.at(chosenQ.at(1)); 
					giniL.at(chosenQ.at(1)) = -1;    
					PL.at(chosenQ.at(1)) = -1;
					qWhere.push_back(1);
				}
				else {
					giniF = giniL.at(chosenQ.at(0));
					giniL.at(chosenQ.at(0)) = -1;    
					PL.at(chosenQ.at(0)) = -1;
					qWhere.push_back(2);
				}
			}
			else {
				if ( giniL.at(chosenQ.at(1)) > giniR.at(chosenQ.at(0)) ) {
					giniF = giniL.at(chosenQ.at(1));
					giniL.at(chosenQ.at(1)) = -1;    
					PL.at(chosenQ.at(1)) = -1;
					qWhere.push_back(1);
				}
				else {
					giniF = giniR.at(chosenQ.at(0));
					giniR.at(chosenQ.at(0)) = -1;    
					PR.at(chosenQ.at(0)) = -1;
					qWhere.push_back(2);
				}
			}
		}
		else {
			//treeStructure.at(2) = chosenQ.at(1);
			if ( qWhere.at(0) == 0 ) {
				if ( giniR.at(chosenQ.at(1)) > giniL.at(chosenQ.at(0)) ) {
					giniF = giniR.at(chosenQ.at(1)); 
					giniR.at(chosenQ.at(1)) = -1;    
					PR.at(chosenQ.at(1)) = -1;
					qWhere.push_back(0);
				}
				else {
					giniF = giniL.at(chosenQ.at(0));
					giniL.at(chosenQ.at(0)) = -1;    
					PL.at(chosenQ.at(0)) = -1;
					qWhere.push_back(2);
				}
			}
			else {
				if ( giniR.at(chosenQ.at(1)) > giniR.at(chosenQ.at(0)) ) {
					giniF = giniR.at(chosenQ.at(1)); 
					giniR.at(chosenQ.at(1)) = -1;    
					PR.at(chosenQ.at(1)) = -1;
					qWhere.push_back(0);
				}
				else {
					giniF = giniR.at(chosenQ.at(0));
					giniR.at(chosenQ.at(0)) = -1;    
					PR.at(chosenQ.at(0)) = -1;
					qWhere.push_back(2);
				}
			}
		}

		if ( c1 == 0 || c1 == trees - 1 ) {
			std::cout << "ChosenQ.at(1)" << ":    " << chosenQ.at(1) << std::endl;
			std::cout << "The new giniF" << ": " << giniF << std::endl;
		}

		// ----------------- Progress bar -----------------
		progressBarCount++;
		buff.str("");
		buff.clear();
		buff << "\33[2K |";  // \33[2K removes the current line.
		for ( c3 = 0; c3 < progressBarCount; ++c3) buff << "=";
		//if ( counter != entries - 1 ) {
		//	buff << ">";
		//	c3++;
		//}
		for ( ; c3 < 29 ; ++c3) buff << "-";
		buff << "|\r";
		std::cout << buff.str();
		std::cout.flush();
		// ------------------------------------------------




		if ( chosenQ.at(1) == 0 ) {            // This might not be needed if relplaced with something like below
			WRC1 = WR0;
			WLC1 = WL0;
		} else if ( chosenQ.at(1) == 1 ) {
			WRC1 = WR1;
			WLC1 = WL1;				
		} else if ( chosenQ.at(1) == 2 ) {
			WRC1 = WR2;
			WLC1 = WL2;
		} else if ( chosenQ.at(1) == 3 ) {
			WRC1 = WR3;
			WLC1 = WL3;
		}
		/*
		if ( chosenQ.at(0) == 0 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR0;
			else WCTemporary = WL0;
		} else if ( chosenQ.at(0) == 1 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR1;
			else WCTemporary = WL1;
		} else if ( chosenQ.at(0) == 2 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR2;
			else WCTemporary = WL2;
		} else if ( chosenQ.at(0) == 3 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR3;
			else WCTemporary = WL3;
		}
		*/

		tester.notDone();

		if ( qWhere.at(1) == 1 ) WCTemporary.swap(WLC1);
		else if ( qWhere.at(1) == 0 ) WCTemporary.swap(WRC1);
		else if ( qWhere.at(1) == 2 ) {
			if ( qWhere.at(0) == 0 ) WCTemporary.swap(WLC0);    // These have to be done in other way than on main level because if one branch has been chosen once it cannot be chosen again.
			else WCTemporary.swap(WRC0);
		}
		for (c2 = 0; c2 < questions; ++c2 ){
			if (c2 == chosenQ.at(0) || c2 == chosenQ.at(1)) continue;        // Because we don't want to check the same question twice.

			for (c3 = 0; c3 < (int)edep.size(); ++c3) {
				if ( WCTemporary.at(c3) > 0.0 ) {
					//hitHolder.setHit(edep.at(c3), layerID.at(c3), cellID.at(c3));
					if ( tester.test(c3, c2, WCTemporary, giniF, chosen_cut )) {
						testWL.push_back(W.at(c3));
						testWR.push_back(-W.at(c3));
						
					}
					else {
						testWR.push_back(W.at(c3));
						testWL.push_back(-W.at(c3));
					}
				}
				else {
					testWL.push_back(-W.at(c3));
					testWR.push_back(-W.at(c3));
				}
			}
			
			giniL.at(c2) = tester.gini(testWL, hittype);
			giniR.at(c2) = tester.gini(testWR, hittype);
			criterion.at(c2) = tester.criterion(giniF, giniL.at(c2), giniR.at(c2));
			PL.at(c2) = tester.purity(testWL, hittype);
			PR.at(c2) = tester.purity(testWR, hittype);
			if ( c2 == 0 ) {
				WR0.swap(testWR);
				WL0.swap(testWL);
			} else if ( c2 == 1 ) {
				WR1.swap(testWR);
				WL1.swap(testWL);
			} else if ( c2 == 2 ) {
				WR2.swap(testWR);
				WL2.swap(testWL);
			} else if ( c2 == 3 ) {
				WR3.swap(testWR);
				WL3.swap(testWL);
			}
			testWL.clear();
			testWR.clear();	

			// ----------------- Progress bar -----------------
			progressBarCount++;
			buff.str("");
			buff.clear();
			buff << "\33[2K |";  // \33[2K removes the current line.
			for ( c3 = 0; c3 < progressBarCount; ++c3) buff << "=";
			//if ( counter != entries - 1 ) {
			//	buff << ">";
			//	c3++;
			//}
			for ( ; c3 < 29 ; ++c3) buff << "-";
			buff << "|\r";
			std::cout << buff.str();
			std::cout.flush();
			progressBarCount++;
			// ------------------------------------------------
		}

		if ( showMessage ) {     // For testing purposes
				std::cout << "\33[2KThe third splitting:" << std::endl;
				std::cout << "Criterion" << 0 << ": " << std::setw(11) << criterion.at(0) << " | ";
				std::cout << "PurityL"   << 0 << ": " << std::setw(11) << PL.at(0)        << " | ";
				std::cout << "PurityR"   << 0 << ": " << std::setw(11) << PR.at(0)        << " | ";
				std::cout << "GiniL"     << 0 << ": " << std::setw(11) << giniL.at(0)     << " | ";
				std::cout << "GiniR"     << 0 << ": " << std::setw(11) << giniR.at(0)     << std::endl;
		}
		chosenQ.at(2) = 0;			
		for ( c2 = 1; c2 < questions; ++c2 ) {
			if ( showMessage ) {     // For testing purposes
				std::cout << "Criterion" << c2 << ": " << std::setw(11) << criterion.at(c2) << " | ";
				std::cout << "PurityL"   << c2 << ": " << std::setw(11) << PL.at(c2)        << " | ";
				std::cout << "PurityR"   << c2 << ": " << std::setw(11) << PR.at(c2)        << " | ";
				std::cout << "GiniL"     << c2 << ": " << std::setw(11) << giniL.at(c2)     << " | ";
				std::cout << "GiniR"     << c2 << ": " << std::setw(11) << giniR.at(c2)     << std::endl;
			}
			if (criterion.at(c2) > criterion.at(chosenQ.at(2))) chosenQ.at(2) = c2;
		}
		criterion.at(chosenQ.at(2)) = -1;

		if ( qWhere.at(1) == 1 ) {
			if ( treeStructure.at(1) >= 0 ) treeStructure.at(3) = chosenQ.at(2);
			else treeStructure.at(5) = chosenQ.at(2);
		}
		else if ( qWhere.at(1) == 0 ) {
			if ( treeStructure.at(1) >= 0 ) treeStructure.at(4) = chosenQ.at(2);
			else treeStructure.at(6) = chosenQ.at(2);
		}
		else if ( qWhere.at(1) == 2 ) {
			if ( treeStructure.at(1) >= 0 ) treeStructure.at(2) = chosenQ.at(2);
			else treeStructure.at(1) = chosenQ.at(2);
		}


			
		giniF = giniR.at(chosenQ.at(2));
		qWhere.push_back(0);
		holder = 2;
		holdThisTruth = false;

		for (c2 = qWhere.size() - 1; c2 > 0; --c2) {
			if( giniL.at(chosenQ.at(c2)) > giniF ) {
				giniF = giniL.at(chosenQ.at(c2));
				holder = c2;
				holdThisTruth = true;
			}
			if( giniR.at(chosenQ.at(c2 - 1)) > giniF ) {
				giniF = giniR.at(chosenQ.at(c2 - 1));
				holder = c2 - 1;
				holdThisTruth = false;
			}
		}
		if( giniL.at(chosenQ.at(0)) > giniF ) {
			giniF = giniL.at(chosenQ.at(0));
			holder = 0;
			holdThisTruth = true;
		}

		if ( holder == 2 ) {
			if ( holdThisTruth ) {
				giniL.at(chosenQ.at(holder)) = -1;
				PL.at(chosenQ.at(holder))    = -1;
				qWhere.at(qWhere.size() - 1) = 1;
			} else {
				giniR.at(chosenQ.at(holder)) = -1;
				PR.at(chosenQ.at(holder))    = -1;
				qWhere.at(qWhere.size() - 1) = 0;
			}
		} else if ( holder == 1 ) {
			if ( holdThisTruth ) {
				giniL.at(chosenQ.at(holder)) = -1;
				PL.at(chosenQ.at(holder))    = -1;
				qWhere.at(qWhere.size() - 1) = 2;
			} else {
				giniR.at(chosenQ.at(holder)) = -1;
				PR.at(chosenQ.at(holder))    = -1;
				qWhere.at(qWhere.size() - 1) = 2;
			}
		} else if ( holder == 0 ) {
			if ( holdThisTruth ) {
				giniL.at(chosenQ.at(holder)) = -1;
				PL.at(chosenQ.at(holder))    = -1;
				qWhere.at(qWhere.size() - 1) = 3;
			} else {
				giniR.at(chosenQ.at(holder)) = -1;
				PR.at(chosenQ.at(holder))    = -1;
				qWhere.at(qWhere.size() - 1) = 3;
			}
		}
			

		/* // Abandoned idea for now.
		giniF = giniL[chosenQ.at(2)];
		if ( giniL[chosenQ.at(2)] > giniR[chosenQ.at(2)] ) {  
			if ( qWhere.at(1) == 0 ) {
				if ( giniL[chosenQ.at(2)] > giniL.at(chosenQ.at(1)) ) {
					if ( qWhere.at(0) == 0 ) {
						if ( giniL[chosenQ.at(2)] > giniL.at(chosenQ.at(0)) ) {
							giniF = giniL[chosenQ.at(2)]; 
							qWhere.push_back(1);
						}
						else {
							qiniF = giniL.at(chosenQ.at(0));
							qWhere.push_back(3);
						}
					}
					else {
						if ( giniL[chosenQ.at(2)] > giniR.at(chosenQ.at(0)) ) {
							giniF = giniL[chosenQ.at(2)]; 
							qWhere.push_back(1);
						}
						else {
							qiniF = giniR.at(chosenQ.at(0));
							qWhere.push_back(3);
						}
					}
				}
				else {
					giniF = giniL.at(chosenQ.at(0));
					itIsLeft1 = 2;
				}
			}
			else {
				if ( giniL.at(chosenQ.at(1)) > giniR.at(chosenQ.at(0)) ) {
					giniF = giniL.at(chosenQ.at(1));
					itIsLeft1 = 1;
				}
				else {
					giniF = giniR.at(chosenQ.at(0));
					itIsLeft1 = 2;
				}
			}
		}
		else {
			giniF = giniR[chosenQ.at(2)]
			if ( qWhere.at(1) == 0 ) {
				if ( giniR[chosenQ.at(2)] > giniL.at(chosenQ.at(1)) ) {
					if ( qWhere.at(0) == 0 ) {
						if ( giniR[chosenQ.at(2)] > giniL.at(chosenQ.at(0)) ) {
							giniF = giniR[chosenQ.at(2)]; 
							qWhere.push_back(0);
						}
						else {
							qiniF = giniL.at(chosenQ.at(0));
							qWhere.push_back(3);
						}
					}
					else {
						if ( giniR[chosenQ.at(2)] > giniR.at(chosenQ.at(0)) ) {
							giniF = giniR[chosenQ.at(2)]; 
							qWhere.push_back(1);
						}
						else {
							qiniF = giniR.at(chosenQ.at(0));
							qWhere.push_back(3);
						}
					}
				}
				else {
					giniF = giniL.at(chosenQ.at(0));
					itIsLeft1 = 2;
				}
			}
			else {
				if ( giniR.at(chosenQ.at(1)) > giniR.at(chosenQ.at(0)) ) {
					giniF = giniR.at(chosenQ.at(1)); 
					itIsLeft1 = 0;
				}
				else {
					giniF = giniR.at(chosenQ.at(0));
					itIsLeft1 = 2;
				}
			}
		}
		*/

		if ( c1 == 0 || c1 == trees - 1 ) {
			std::cout << "ChosenQ.at(2)" << ":    " << chosenQ.at(2) << std::endl;
			std::cout << "The new giniF" << ": " << giniF << std::endl;
		}


		// ----------------- Progress bar -----------------
		progressBarCount++;
		buff.str("");
		buff.clear();
		buff << "\33[2K |";  // \33[2K removes the current line.
		for ( c3 = 0; c3 < progressBarCount; ++c3) buff << "=";
		//if ( counter != entries - 1 ) {
		//	buff << ">";
		//	c3++;
		//}
		for ( ; c3 < 29 ; ++c3) buff << "-";
		buff << "|\r";
		std::cout << buff.str();
		std::cout.flush();
		// ------------------------------------------------




		if ( chosenQ.at(2) == 0 ) {            // This might not be needed if relplaced with something like below
			WRC2 = WR0;
			WLC2 = WL0;
		} else if ( chosenQ.at(2) == 1 ) {
			WRC2 = WR1;
			WLC2 = WL1;				
		} else if ( chosenQ.at(2) == 2 ) {
			WRC2 = WR2;
			WLC2 = WL2;
		} else if ( chosenQ.at(2) == 3 ) {
			WRC2 = WR3;
			WLC2 = WL3;
		}

		/*
		if ( chosenQ.at(0) == 0 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR0;
			else WCTemporary = WL0;
		} else if ( chosenQ.at(0) == 1 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR1;
			else WCTemporary = WL1;				
		} else if ( chosenQ.at(0) == 2 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR2;
			else WCTemporary = WL2;
		} else if ( chosenQ.at(0) == 3 ) {
			if ( qWhere.at(0) > 0 ) WCTemporary = WR3;
			else WCTemporary = WL3;
		}
		*/

		tester.notDone();
			
		if ( qWhere.at(2) == 1 )  WCTemporary.swap(WLC2);
		else if ( qWhere.at(2) == 0 )  WCTemporary.swap(WRC2);
		else if ( qWhere.at(2) == 2 ) {
			if ( qWhere.at(1) == 0 )  WCTemporary.swap(WLC1);
			else if ( qWhere.at(1) == 1 )  WCTemporary.swap(WRC1);
			else if ( qWhere.at(1) == 2 ) {
				if ( giniL.at(chosenQ.at(1)) > giniR.at(chosenQ.at(1)) )  WCTemporary.swap(WLC1);
				else  WCTemporary.swap(WRC1); 
			}
		} 
		else if ( qWhere.at(2) == 3 ) {
			if ( qWhere.at(0) == 0 )  WCTemporary.swap(WLC0);
			else  WCTemporary.swap(WRC0);
		}
			
		//tester.cutting(tester, chosenQ, 3, WCTemporary, giniF, giniL, giniR, criterion, PL, PR, chosen_cut, WR0, WL0, WR1, WL1, WR2, WL2, WR3, WL3);

		
		testWL.clear();
		testWR.clear();
		for (c2 = 0; c2 < questions; ++c2 ) {
			for (c3 = 0; c3 < 3; ++c3) if ( chosenQ.at(c3) == c2) passMainLoop = true;        // Because we don't want to check the same question twice.
			
			if ( !passMainLoop ) {
				
				for (c3 = 0; c3 < (int)edep.size(); ++c3) {
					if ( WCTemporary.at(c3) > 0.0 ) {
						//hitHolder.setHit(edep.at(c3), hittype.at(c3), layerID.at(c3));
						if ( tester.test(c3, c2, WCTemporary, giniF, chosen_cut )) {
							testWL.push_back(WCTemporary.at(c3));
							testWR.push_back(-WCTemporary.at(c3));
							
						}
						else {
							testWR.push_back(WCTemporary.at(c3));
							testWL.push_back(-WCTemporary.at(c3));
						}
					}
					else {
						testWL.push_back(WCTemporary.at(c3));  // Note that these are negative numbers which are pushed back.
						testWR.push_back(WCTemporary.at(c3));
					}
				}
				
				giniL.at(c2) = tester.gini(testWL, hittype);
				giniR.at(c2) = tester.gini(testWR, hittype);
				criterion.at(c2) = tester.criterion(giniF, giniL.at(c2), giniR.at(c2));
				PL.at(c2) = tester.purity(testWL, hittype);
				PR.at(c2) = tester.purity(testWR, hittype);
				if ( c2 == 0 ) {
					WR0.swap(testWR);
					WL0.swap(testWL);
				} else if ( c2 == 1 ) {
					WR1.swap(testWR);
					WL1.swap(testWL);
				} else if ( c2 == 2 ) {
					WR2.swap(testWR);
					WL2.swap(testWL);
				} else if ( c2 == 3 ) {
					WR3.swap(testWR);
					WL3.swap(testWL);
				}
				testWL.clear();
				testWR.clear();	
			}
			passMainLoop = false;

			// ----------------- Progress bar -----------------
			progressBarCount++;
			buff.str("");
			buff.clear();
			buff << "\33[2K |";  // \33[2K removes the current line.
			for ( c3 = 0; c3 < progressBarCount; ++c3) buff << "=";
			//if ( counter != entries - 1 ) {
			//	buff << ">";
			//	c3++;
			//}
			for ( ; c3 < 29 ; ++c3) buff << "-";
			buff << "|\r";
			std::cout << buff.str();
			std::cout.flush();
			progressBarCount++;
			// ------------------------------------------------
		}
		progressBarCount = 0;
		usedQ++;
			

		if ( showMessage ) {                         // For testing purposes
				std::cout << "\33[2KThe fourth splitting:" << std::endl;
				std::cout << "Criterion"  << 0 << ": " << std::setw(11) << criterion.at(0) << " | ";
				std::cout << "PurityL"    << 0 << ": " << std::setw(11) << PL.at(0)        << " | ";
				std::cout << "PurityR"    << 0 << ": " << std::setw(11) << PR.at(0)        << " | ";
				std::cout << "GiniL"      << 0 << ": " << std::setw(11) << giniL.at(0)     << " | ";
				std::cout << "GiniR"      << 0 << ": " << std::setw(11) << giniR.at(0)     << std::endl;
		}
		chosenQ.at(3) = 0;
		for ( c2 = 1; c2 < questions; ++c2 ) {
			if ( showMessage ) {                         // For testing purposes
				std::cout << "Criterion"  << c2 << ": " << std::setw(11) << criterion.at(c2) << " | ";
				std::cout << "PurityL"    << c2 << ": " << std::setw(11) << PL.at(c2)        << " | ";
				std::cout << "PurityR"    << c2 << ": " << std::setw(11) << PR.at(c2)        << " | ";
				std::cout << "GiniL"      << c2 << ": " << std::setw(11) << giniL.at(c2)     << " | ";
				std::cout << "GiniR"      << c2 << ": " << std::setw(11) << giniR.at(c2)     << std::endl;
			}
			if (criterion.at(c2) > criterion.at(chosenQ.at(3))) chosenQ.at(3) = c2;
		}



		if ( qWhere.at(2) == 1 ) {
			if (qWhere.at(1) == 0 || qWhere.at(1) == 1) {
				if ( treeStructure.at(3) >= 0 ) treeStructure.at(7) = chosenQ.at(3);
				else if ( treeStructure.at(4) >= 0 ) treeStructure.at(9) = chosenQ.at(3);
				else if ( treeStructure.at(5) >= 0 ) treeStructure.at(11) = chosenQ.at(3);
				else if ( treeStructure.at(6) >= 0 ) treeStructure.at(13) = chosenQ.at(3);
				else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
			}
			else if (qWhere.at(1) == 2) {
				if ( qWhere.at(0) == 1 ) {
					if ( giniL.at(chosenQ.at(2)) > giniR.at(chosenQ.at(2)) ) treeStructure.at(11) = chosenQ.at(3);
					else if ( giniL.at(chosenQ.at(2)) < giniR.at(chosenQ.at(2)) ) treeStructure.at(13) = chosenQ.at(3);
					else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
				}
				else {
					if ( giniL.at(chosenQ.at(2)) > giniR.at(chosenQ.at(2)) ) treeStructure.at(7) = chosenQ.at(3);
					else if ( giniL.at(chosenQ.at(2)) < giniR.at(chosenQ.at(2)) ) treeStructure.at(9) = chosenQ.at(3);
					else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
				}
			}
			else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
		}
		else if ( qWhere.at(2) == 0 ) {
			if (qWhere.at(1) == 0 || qWhere.at(1) == 1) {
				if ( treeStructure.at(3) >= 0 ) treeStructure.at(8) = chosenQ.at(3);
				else if ( treeStructure.at(4) >= 0 ) treeStructure.at(10) = chosenQ.at(3);
				else if ( treeStructure.at(5) >= 0 ) treeStructure.at(12) = chosenQ.at(3);
				else if ( treeStructure.at(6) >= 0 ) treeStructure.at(14) = chosenQ.at(3);
				else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
			}
			else if (qWhere.at(1) == 2) {
				if ( qWhere.at(0) == 1 ) {
					if ( giniL.at(chosenQ.at(2)) > giniR.at(chosenQ.at(2)) ) treeStructure.at(12) = chosenQ.at(3);
					else if ( giniL.at(chosenQ.at(2)) < giniR.at(chosenQ.at(2)) ) treeStructure.at(14) = chosenQ.at(3);
					else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
				}
				else {
					if ( giniL.at(chosenQ.at(2)) > giniR.at(chosenQ.at(2)) ) treeStructure.at(8) = chosenQ.at(3);
					else if ( giniL.at(chosenQ.at(2)) < giniR.at(chosenQ.at(2)) ) treeStructure.at(10) = chosenQ.at(3);
					else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
				}
			}
			else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
		}
		else if ( qWhere.at(2) == 2 ) {
			if ( treeStructure.at(3) >= 0 ) treeStructure.at(4) = chosenQ.at(3);
			else if ( treeStructure.at(4) >= 0 ) treeStructure.at(3) = chosenQ.at(3);
			else if ( treeStructure.at(5) >= 0 ) treeStructure.at(6) = chosenQ.at(3);
			else if ( treeStructure.at(6) >= 0 ) treeStructure.at(5) = chosenQ.at(3);
			else if ( qWhere.at(1) == 2 ) {
				if ( qWhere.at(0) == 1 ) {
					if ( giniL.at(chosenQ.at(1)) > giniR.at(chosenQ.at(1)) ) {
						if ( giniL.at(chosenQ.at(3)) > giniR.at(chosenQ.at(3)) ) treeStructure.at(7) = chosenQ.at(3);
						else if ( giniL.at(chosenQ.at(3)) < giniR.at(chosenQ.at(3)) ) treeStructure.at(8) = chosenQ.at(3);
						else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
					}
					else {
						if ( giniL.at(chosenQ.at(3)) > giniR.at(chosenQ.at(3)) ) treeStructure.at(9) = chosenQ.at(3);
						else if ( giniL.at(chosenQ.at(3)) < giniR.at(chosenQ.at(3)) ) treeStructure.at(10) = chosenQ.at(3);
						else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
					}
				}
				else {
					if ( giniL.at(chosenQ.at(1)) > giniR.at(chosenQ.at(1)) ) {
						if ( giniL.at(chosenQ.at(3)) > giniR.at(chosenQ.at(3)) ) treeStructure.at(11) = chosenQ.at(3);
						else if ( giniL.at(chosenQ.at(3)) < giniR.at(chosenQ.at(3)) ) treeStructure.at(12) = chosenQ.at(3);
						else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
					}
					else {
						if ( giniL.at(chosenQ.at(3)) > giniR.at(chosenQ.at(3)) ) treeStructure.at(13) = chosenQ.at(3);
						else if ( giniL.at(chosenQ.at(3)) < giniR.at(chosenQ.at(3)) ) treeStructure.at(14) = chosenQ.at(3);
						else std::cout << "There is actually a mistake here. Please repair me." << std::endl;
					}
				}
			}
		}
		else if ( qWhere.at(2) == 3 ) {
			if ( treeStructure.at(1) >= 0 ) treeStructure.at(2) = chosenQ.at(3);
			else treeStructure.at(1) = chosenQ.at(3);
		}


		// ================================= Decide if hits are signal or not =================================


		for (c2 = 0; c2 < (int)W.size(); ++c2) sig.push_back(0);


		for (c2 = 0, c3 = 1; c2 < 15; ++c2, c3 += 2) {
			if ( ( treeStructure.at(c3) < 0 || treeStructure.at(c3) > 3 ) && treeStructure.at(c2) > -1 && treeStructure.at(c2) < 4) {
				if ( PL.at(treeStructure.at(c2)) > purity_check ) treeStructure.at(c3) = 66;
				else treeStructure.at(c3) = -66;
			}
			if ( ( treeStructure.at(c3 + 1) < 0 || treeStructure.at(c3 + 1) > 3 ) && treeStructure.at(c2) > -1 && treeStructure.at(c2) < 4) {
				if ( PR.at(treeStructure.at(c2)) > purity_check ) treeStructure.at(c3 + 1) = 66;
				else treeStructure.at(c3 + 1) = -66;
			}
		}

			
		for (c2 = 0; c2 < (int)W.size(); ++c2) {
			if      ( PL.at(0) > purity_check && WL0.at(c2) > 0 ) sig.at(c2) = 1;
			else if ( PR.at(0) > purity_check && WR0.at(c2) > 0 ) sig.at(c2) = 1;
			else if ( PL.at(1) > purity_check && WL1.at(c2) > 0 ) sig.at(c2) = 1;
			else if ( PR.at(1) > purity_check && WR1.at(c2) > 0 ) sig.at(c2) = 1;
			else if ( PL.at(2) > purity_check && WL2.at(c2) > 0 ) sig.at(c2) = 1;
			else if ( PR.at(2) > purity_check && WR2.at(c2) > 0 ) sig.at(c2) = 1;
			else if ( PL.at(3) > purity_check && WL3.at(c2) > 0 ) sig.at(c2) = 1;
			else if ( PR.at(3) > purity_check && WR3.at(c2) > 0 ) sig.at(c2) = 1;
		}
			
		errResult = tester.err(W, sig, hittype);

		alphaSave = beta * log((1 - errResult) / errResult );
			
		alpha.push_back(alphaSave);

		if (alpha.at(c1) > 0.0) alphaSummed += alpha.at(c1);
		else alphaSummed -= alpha.at(c1);

		for ( c2 = 0; c2 < (int)PR.size(); ++c2) {
			if (PR.at(c2) >= 0.0) hPurity->Fill(PR.at(c2));
			if (PL.at(c2) >= 0.0) hPurity->Fill(PL.at(c2));
		}


			
		for (c2 = 0; c2 < (int)W.size(); ++c2){
			if      ( PL.at(0) > purity_check && WL0.at(c2) > 0 ) T.at(c2) += alpha.at(c1);
			else if ( PR.at(0) > purity_check && WR0.at(c2) > 0 ) T.at(c2) += alpha.at(c1);
			else if ( PL.at(1) > purity_check && WL1.at(c2) > 0 ) T.at(c2) += alpha.at(c1);
			else if ( PR.at(1) > purity_check && WR1.at(c2) > 0 ) T.at(c2) += alpha.at(c1);
			else if ( PL.at(2) > purity_check && WL2.at(c2) > 0 ) T.at(c2) += alpha.at(c1);
			else if ( PR.at(2) > purity_check && WR2.at(c2) > 0 ) T.at(c2) += alpha.at(c1);
			else if ( PL.at(3) > purity_check && WL3.at(c2) > 0 ) T.at(c2) += alpha.at(c1);
			else if ( PR.at(3) > purity_check && WR3.at(c2) > 0 ) T.at(c2) += alpha.at(c1);
			else T.at(c2) -= alpha.at(c1);


			if (hittype.at(c2) > 0) signalTrueCount++;

			if ( (sig.at(c2) < 1 && hittype.at(c2) > 0 ) || (sig.at(c2) > 0 && hittype.at(c2) < 1 ) ) {
				W.at(c2) *= exp(alpha.at(c1));
				//T.at(c2) += -alpha.at(c1); // This was actually wrong!
			}
			else { 
				hitWasCorrect++;
				//T.at(c2) +=  alpha.at(c1);
				if ( hittype.at(c2) > 0) signalCount++;
			}
		}
			
		/*
		if ( c1 == 0) {                                             // With this you can see all the weights 53.00 -> 03.27
			for (c2 = 0; c2 < (int)W.size() -8; c2 += 9) {
					std::cout << std::setw(11) << W.at(c2) << " ";
					std::cout << std::setw(11) << W[c2 + 1] << " ";
					std::cout << std::setw(11) << W[c2 + 2] << " ";
					std::cout << std::setw(11) << W[c2 + 3] << " ";
					std::cout << std::setw(11) << W[c2 + 4] << " ";
					std::cout << std::setw(11) << W[c2 + 5] << " ";
					std::cout << std::setw(11) << W[c2 + 6] << " ";
					std::cout << std::setw(11) << W[c2 + 7] << " ";
					std::cout << std::setw(11) << W[c2 + 8] << " ";
					std::cout << std::endl;
				}
			for (; c2 < (int)W.size(); ++c2) {
				std::cout << std::setw(11) << W.at(c2) << " ";
			}
			std::cout << std::endl;
		}
		*/


		if ( false /* c1 == trees - 1 || c1 == 0 */) {                                              // With this you can see all the results and actual signal information as pairs
			for (c2 = 0; c2 < (int)sig.size() -18; c2 += 19) {
					std::cout << std::setw(4) << sig.at(c2)      << " " << std::setw(2) << hittype.at(c2)      << " ";
					std::cout << std::setw(4) << sig.at(c2 + 1 ) << " " << std::setw(2) << hittype.at(c2 + 1 ) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 2 ) << " " << std::setw(2) << hittype.at(c2 + 2 ) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 3 ) << " " << std::setw(2) << hittype.at(c2 + 3 ) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 4 ) << " " << std::setw(2) << hittype.at(c2 + 4 ) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 5 ) << " " << std::setw(2) << hittype.at(c2 + 5 ) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 6 ) << " " << std::setw(2) << hittype.at(c2 + 6 ) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 7 ) << " " << std::setw(2) << hittype.at(c2 + 7 ) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 8 ) << " " << std::setw(2) << hittype.at(c2 + 8 ) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 9 ) << " " << std::setw(2) << hittype.at(c2 + 9 ) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 10) << " " << std::setw(2) << hittype.at(c2 + 10) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 11) << " " << std::setw(2) << hittype.at(c2 + 11) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 12) << " " << std::setw(2) << hittype.at(c2 + 12) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 13) << " " << std::setw(2) << hittype.at(c2 + 13) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 14) << " " << std::setw(2) << hittype.at(c2 + 14) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 15) << " " << std::setw(2) << hittype.at(c2 + 15) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 16) << " " << std::setw(2) << hittype.at(c2 + 16) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 17) << " " << std::setw(2) << hittype.at(c2 + 17) << " ";
					std::cout << std::setw(4) << sig.at(c2 + 18) << " " << std::setw(2) << hittype.at(c2 + 18) << " ";
					std::cout << std::endl;
				}
			for (; c2 < (int)sig.size(); ++c2) {
				std::cout << std::setw(4) << sig.at(c2) << " " << std::setw(2) << hittype.at(c2) << " ";
			}
			std::cout << std::endl;
		}
			


		dVecSumResult = tester.dVecSum(W);
		if (  c1 == 0 || c1 == trees - 1 ) std::cout << "W sum before normalisation: " << dVecSumResult << std::endl;

		for (c2 = 0; c2 < (int)W.size(); ++c2){
			W.at(c2) /= dVecSumResult;
		}

		if (  c1 == 0 || c1 == trees - 1 ) std::cout << "W sum after normalisation:  " << tester.dVecSum(W) << std::endl;



		if ( signalTrueCount != 0 ) trueSigRatio = (double)signalCount / (double)signalTrueCount;
		else if ( signalCount > 0 ) std::cout << "============ AN ERROR HAS OCCURED ============";
		else trueSigRatio = -1.0; 
		if ( true /* c1%updates == 0 || c1 == trees - 1 */ ) {
			std::cout << "~~~         err["  << std::setw(4) << c1 << "]: " << std::setw(11) << errResult << " ";
			std::cout << " ~~        alpha[" << std::setw(4) << c1 << "]: " << std::setw(12) << alpha.at(c1) << " ";
			std::cout << " ~~        ratio[" << std::setw(4) << c1 << "]: " << std::setw(11) << (double)hitWasCorrect / (double)(int)W.size()  << " ";
			std::cout << " ~~       sratio[" << std::setw(4) << c1 << "]: " << std::setw(10) << trueSigRatio << std::endl;
			std::cout <<  "    chosen_cut0[" << std::setw(4) << c1 << "]: " << std::setw(11) << chosen_cut.at(0) << " ";
			std::cout << " ~~  chosen_cut1[" << std::setw(4) << c1 << "]: " << std::setw(12) << chosen_cut.at(1) << " ";
			std::cout << " ~~  chosen_cut2[" << std::setw(4) << c1 << "]: " << std::setw(11) << chosen_cut.at(2) << " ";
			std::cout << " ~~  chosen_cut3[" << std::setw(4) << c1 << "]: " << std::setw(10) << chosen_cut.at(3) << " ~~~" << std::endl;
		}

		tree->Fill();
	}

	for ( c1 = 0; c1 < (int)T.size(); ++c1) {
		T.at(c1) /= alphaSummed;
		if ( hittype.at(c1) > 0 )
			hResultSignal->Fill( T.at(c1) );
		else
			hResult->Fill( T.at(c1) );
	}



	// =============================================================================================================

	std::vector<double> WEqual((int)T.size(), 1.0/(double)T.size());
	giniF = tester.gini(WEqual, hittype);
	double resultCriterion = 0.0;
	double bestResultCrit = 0.0;
	theBestCutForResults = -1.0;
	int nSignalTestCorrect = 0;
	int nBackgroundTestCorrect = 0;
	int nSignal = 0;
	int nBackground = 0;
	for (int cp1 = 0; cp1 < (int)hittype.size(); ++cp1) {
		if (hittype.at(cp1) > 0) nSignal++;
		else nBackground++;
	}
	testWL.clear();
	testWR.clear();

	for (double i = -1.0; i < 1.0; i += 0.001) {
		for (int c = 0; c < (int)T.size(); ++c) {
			if ( T.at(c) < i ) {   // This does the actual comparing
				testWL.push_back(WEqual.at(c));
				testWR.push_back(-1);
			}
			else {
				testWL.push_back(-1);
				testWR.push_back(WEqual.at(c));
			}
		}
		resultCriterion = tester.criterion(giniF, tester.gini(testWL, hittype), tester.gini(testWR, hittype));
		if ( resultCriterion > bestResultCrit ) {
			bestResultCrit = resultCriterion;
			theBestCutForResults = i;
		}
		for (int j = 0; j < (int)testWL.size(); ++j) {
			if (testWL.at(j) > 0.0 && hittype.at(j) < 1) nBackgroundTestCorrect++;
			else if (testWR.at(j) > 0.0 && hittype.at(j) > 0) nSignalTestCorrect++;
			//else std::cout << "W = 0 in testResult." << std::endl;
		}
		//hSignalRetentionE->Fill((double)nSignalTestCorrect/(double)nSignal);
		//hBackgroundRejectionE->Fill((double)nBackgroundTestCorrect/(double)nBackground);
		hEfficiency->Fill((double)nSignalTestCorrect/(double)nSignal, (double)nBackgroundTestCorrect/(double)nBackground);

		nSignalTestCorrect = 0;
		nBackgroundTestCorrect = 0;
		testWL.clear();
		testWR.clear();
	}

	std::cout << "The best cut for results: " << theBestCutForResults << std::endl;

	for (c2 = 0; c2 < (int)T.size(); ++c2) {
		if      ( T.at(c2) > theBestCutForResults ) sig.at(c2) = 1;
		else sig.at(c2) = 0;
	}

	// =============================================================================================================



	hitWasCorrect = 0;
	signalCount = 0;
	for (c2 = 0; c2 < (int)W.size(); ++c2){

		//if (hittype.at(c2) > 0) signalTrueCount++;

		if ( (sig.at(c2) < 1 && hittype.at(c2) > 0 ) || (sig.at(c2) > 0 && hittype.at(c2) < 1 ) ) {
			//W.at(c2) *= exp(alpha.at(c1));
			//T.at(c2) += -alpha.at(c1); // This was actually wrong!
		}
		else { 
			hitWasCorrect++;
			//T.at(c2) +=  alpha.at(c1);
			if ( hittype.at(c2) > 0) signalCount++;
		}
	}

	std::cout << "The final results: " << std::endl;
	if ( c1 == trees - 1 || c1 == 0) {                                              // With this you can see all the results and actual signal information as pairs
		for (c2 = 0; c2 < (int)sig.size() -18; c2 += 19) {
				std::cout << std::setw(4) << sig.at(c2)      << " " << std::setw(2) << hittype.at(c2)      << " ";
				std::cout << std::setw(4) << sig.at(c2 + 1 ) << " " << std::setw(2) << hittype.at(c2 + 1 ) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 2 ) << " " << std::setw(2) << hittype.at(c2 + 2 ) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 3 ) << " " << std::setw(2) << hittype.at(c2 + 3 ) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 4 ) << " " << std::setw(2) << hittype.at(c2 + 4 ) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 5 ) << " " << std::setw(2) << hittype.at(c2 + 5 ) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 6 ) << " " << std::setw(2) << hittype.at(c2 + 6 ) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 7 ) << " " << std::setw(2) << hittype.at(c2 + 7 ) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 8 ) << " " << std::setw(2) << hittype.at(c2 + 8 ) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 9 ) << " " << std::setw(2) << hittype.at(c2 + 9 ) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 10) << " " << std::setw(2) << hittype.at(c2 + 10) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 11) << " " << std::setw(2) << hittype.at(c2 + 11) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 12) << " " << std::setw(2) << hittype.at(c2 + 12) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 13) << " " << std::setw(2) << hittype.at(c2 + 13) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 14) << " " << std::setw(2) << hittype.at(c2 + 14) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 15) << " " << std::setw(2) << hittype.at(c2 + 15) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 16) << " " << std::setw(2) << hittype.at(c2 + 16) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 17) << " " << std::setw(2) << hittype.at(c2 + 17) << " ";
				std::cout << std::setw(4) << sig.at(c2 + 18) << " " << std::setw(2) << hittype.at(c2 + 18) << " ";
				std::cout << std::endl;
			}
		for (; c2 < (int)sig.size(); ++c2) {
			std::cout << std::setw(4) << sig.at(c2) << " " << std::setw(2) << hittype.at(c2) << " ";
		}
		std::cout << std::endl;
	}

	if ( signalTrueCount != 0 ) trueSigRatio = (double)signalCount / (double)signalTrueCount;
	else if ( signalCount > 0 ) std::cout << "============ AN ERROR HAS OCCURED ============";
	else trueSigRatio = -1.0; 
	std::cout << "~~~ ratio: " << std::setw(11) << (double)hitWasCorrect / (double)W.size()  << " ";
	std::cout << " ~~ sratio: " << std::setw(10) << trueSigRatio << std::endl;




	// Create a root file for the histos.

	TString fileName = "BDTOutput.root";
	TFile *fout = new TFile(fileName,"RECREATE"); 
	fout->cd();

	// Write all histograms
	hResult->Write();
	hResultSignal->Write();
	hPurity->Write();
	//hSignalRetentionE->Write();
	//hBackgroundRejectionE->Write();
	hEfficiency->Write();
	tree->Write();

	fout->Close();

	return 0;
    



    

}









