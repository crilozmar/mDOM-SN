#ifndef mdomSNTools_h
#define mdomSNTools_h 1

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

class mdomSNTools
{
public:
  	mdomSNTools();    
	~mdomSNTools();
	 void RandomPosition(G4ThreeVector& Position);
	G4double InverseCumul(std::vector<G4double>  xvals, std::vector<G4double>  yvals, G4int nPoints);
	G4double InverseCumulAlgorithm(std::vector<G4double>  x, std::vector<G4double>  f, std::vector<G4double>  a, std::vector<G4double>  Fc, G4int  nPoints);
	G4int findtime(G4double time, std::vector<G4double> timearray);
	G4double EnergyDistribution(G4double Emean, G4double Emean2, G4double& alpha);
	G4double linealinterpolation(G4double realX,G4double lowerX, G4double upperX, G4double lowerY,G4double upperY);
	G4double GetAlpha(G4double Emean,G4double Emean2);
	G4double WeighMe(G4double energy, G4double sigma, G4double NTargets);

private:     
	std::vector <G4double> mdompos;
	G4double Rmin;
	bool CheckVolumeFormDOMs(G4ThreeVector position);
};

void GetSlopes(std::vector<G4double>  xvals, std::vector<G4double>  yvals, G4int nPoints, std::vector<G4double>&  x, std::vector<G4double>&  f, std::vector<G4double>&  a, std::vector<G4double>&  Fc);
G4double NumberOfTargets(G4int targetPerMolecule);
void MakeEnergyDistribution(G4double Emean, G4double alpha, G4int nPoints, std::vector<G4double>& x, std::vector<G4double>& f);

#endif
