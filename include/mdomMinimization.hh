#ifndef mdomMinimizationr_h
#define mdomMinimization_h 1

#include "G4Types.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"


#include <vector>
#include <fstream>
#include <tuple>

class MdomMinimization
{
	public:
	 MdomMinimization();
	 ~MdomMinimization();
	 
	 
	void ReadEventToReconstruct();
	G4double BuildPoisson(G4double hits, G4double expected);
	void BuildAcumulateHits(struct EvtStat evtStat);
	G4double BuildLikelihood(G4double TotalHits[24]);
	G4double Minimize();
	bool control;
	
	private:
	int Factorial(int a);
	std::vector<double> readingEvent;
	G4double EventToReconstruct[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	G4double Poiss;
};

#endif
