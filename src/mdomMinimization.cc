/** @file mdomMinimization.cc
 *  @brief Minimizer. Only used for concrete studies.
 * 
 *  @author Cristian Jesus Lozano Mariscal (c.lozano@wwu.de)
 * 
 *  @version Geant4 10.7
 */

#include "mdomMinimization.hh"
#include "G4SystemOfUnits.hh"
#include "G4SimplexDownhill.hh"
#include "mdomAnalysisManager.hh"

#include "G4ios.hh"


extern std::vector<std::tuple<G4int,G4int,G4int> > AcumulateHits;
extern MdomMinimization Minimization;
extern G4double LogL;
extern G4String	gEvent2Reconstruct;

extern std::vector<double> readColumnDouble (G4String fn, int col);


MdomMinimization::MdomMinimization(){
	control = true;
}


MdomMinimization::~MdomMinimization(){}


void MdomMinimization::ReadEventToReconstruct()
{
	for (int i=0; i < 24; i++) {
		readingEvent =  readColumnDouble(gEvent2Reconstruct,i+1);
		EventToReconstruct[i] = readingEvent.at(0);
	}
	control = false;
}


int MdomMinimization::Factorial(int a) 
{
	int factorial = 1;
	if (a==0 || a==1) {
		return factorial;
	}
	else {
		for (int b=1 ; b<=a ; b++) {
			factorial=b*factorial;
		}
	return factorial;
	}
}


G4double MdomMinimization::BuildPoisson(G4double hits,G4double expected)
{
	// Build the log of the poisson functions
	G4double calculation = exp(-expected)*pow(expected,hits)/Factorial(hits);
	G4double function;
	if (calculation == 0) {
		function = log(0.0000001); //because log(0)=-inf!!!
	}
	else {
		function = log(calculation);
	}
	//G4cout << "#hits -> " << hits << " # Expected -> " << expected << "# Calculation -> " << calculation <<" # Poisson -> " << function << G4endl;
	if (hits > 0 ) { 
		return function;
	}
	else {
		return 0;
	}
}


void MdomMinimization::BuildAcumulateHits(struct EvtStat evtStat)
{
	// Save the hits from mdomAnalysisManager
	//AcumulateHits.clear();
	for ( int j=0; j<(G4int)evtStat.hitsPMTs.size(); j++ ) {
		AcumulateHits.push_back(evtStat.hitsPMTs[j]);
	}
}


G4double MdomMinimization::BuildLikelihood(G4double TotalHits[24])
{
	if (control){
		ReadEventToReconstruct(); // To read the file just once
	}
	G4double logL= 0;

	for ( int j=0; j<24; j++ ) {
		Poiss = BuildPoisson(EventToReconstruct[j],TotalHits[j]);
		logL = logL + Poiss;
	}
	LogL = -logL; //Global Variable
}


G4double MdomMinimization::Minimize()
{
	return LogL; // Send LogL to mdom.cc (or to whoever calls it)
}

 

   
   
   





