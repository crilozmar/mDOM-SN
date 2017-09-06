#ifndef mdomPrimaryGeneratorAction3_h
#define mdomPrimaryGeneratorAction3_h 1
 
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <vector>


class G4ParticleGun;
class G4Event;

class mdomPrimaryGeneratorAction3 : public G4VUserPrimaryGeneratorAction
{
public:
	mdomPrimaryGeneratorAction3(G4ParticleGun*);    
	~mdomPrimaryGeneratorAction3();

public:
	void GeneratePrimaries(G4Event* anEvent);

  public:        
    G4double InverseCumul(int ControlParameter);   
    void RandomPosition();
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
    
	std::vector<double> data_energy;
	std::vector<double> data_fe;
	  
private:     
    G4ParticleGun*         ParticleGun;

    G4int                  ControlParameter;
    G4int                  nPoints;     //tabulated function
    std::vector<G4double>  x;
    std::vector<G4double>  f;           //f(x)
    std::vector<G4double>  a;           //slopes
    std::vector<G4double>  Fc;          //cumulative of f

  private:

    void Fe_build();     
    void DistFunction(G4double Enu);      
    G4double ElectronEnergy(G4double nu_energy, G4double costheta);
    G4double NumberOfTargets(G4int targetPerMolecule);
    G4double TotalCrossSection(G4double energy);
    G4double WeighMe(G4double energy);

};

#endif
