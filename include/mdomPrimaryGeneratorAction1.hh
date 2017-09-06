#ifndef mdomPrimaryGeneratorAction1_h
#define mdomPrimaryGeneratorAction1_h 1
 
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <vector>


class G4ParticleGun;
class G4Event;

class mdomPrimaryGeneratorAction1 : public G4VUserPrimaryGeneratorAction
{
public:
	mdomPrimaryGeneratorAction1(G4ParticleGun*);    
	~mdomPrimaryGeneratorAction1();

public:
	void GeneratePrimaries(G4Event* anEvent);

  public:        
    G4double InverseCumul(int ControlParameter);   
    void RandomPosition();
    G4int                  nPoints0;     //tabulated function
    std::vector<G4double>  x0;
    std::vector<G4double>  f0;           //f(x)
    std::vector<G4double>  a0;           //slopes
    std::vector<G4double>  Fc0;          //cumulative of f
    G4int                  nPoints1;     //tabulated function
    std::vector<G4double>  x1;
    std::vector<G4double>  f1;           //f(x)
    std::vector<G4double>  a1;           //slopes
    std::vector<G4double>  Fc1;          //cumulative of f
    G4int                  nPoints2;     //tabulated function
    std::vector<G4double>  x2;
    std::vector<G4double>  f2;           //f(x)
    std::vector<G4double>  a2;           //slopes
    std::vector<G4double>  Fc2;          //cumulative of f
    G4double NTargets;
    G4double me;
    G4double Gf;
    
	std::vector<double> nu_time;
	std::vector<double> nu_luminosity;
	std::vector<double> nu_meanenergy;
	std::vector<double> nu_meanenergysquare;
	  
private:     
    G4ParticleGun*         ParticleGun;

    G4int                  ControlParameter;
    G4int                  nPoints;     //tabulated function
    std::vector<G4double>  x;
    std::vector<G4double>  f;           //f(x)
    std::vector<G4double>  a;           //slopes
    std::vector<G4double>  Fc;          //cumulative of f

  private:
    void LuminosityDist(); 
    G4double GetAlpha(G4double Emean,G4double Emean2);
    G4int findtime(G4double time);
    G4double linealinterpolation(G4double realX,G4double lowerX, G4double upperX, G4double lowerY,G4double upperY);
    void Fe_nu(G4double Emean, G4double Emean2);     
    void DistFunction(G4double Enu);      
    G4double ElectronEnergy(G4double nu_energy, G4double costheta);
    G4double NumberOfTargets(G4int targetPerMolecule);
    G4double TotalCrossSection(G4double energy);
    G4double WeighMe(G4double energy);
    G4double fixenergy;
    G4double fixenergy2;
    G4double alpha;

};

#endif
