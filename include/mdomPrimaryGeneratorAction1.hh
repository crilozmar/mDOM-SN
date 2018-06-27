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

	void GeneratePrimaries(G4Event* anEvent);
     
    G4double NTargets;
    G4double me;
    G4double Gf;
    G4int                  nPoints_lum;     //tabulated function
    std::vector<G4double>  x_lum;
    std::vector<G4double>  f_lum;           //f(x)
    std::vector<G4double>  a_lum;           //slopes
    std::vector<G4double>  Fc_lum;          //cumulative of f
    G4int                  nPoints1;     //tabulated function
    std::vector<G4double>  fixFe_X;
    std::vector<G4double>  fixFe_Y;           //f(x)
    std::vector<G4double>  fixFe_a;           //slopes
    std::vector<G4double>  fixFe_Fc;          //cumulative of f
    G4int fixE_nPoints;
    G4int                  angdist_nPoints;    
    std::vector<G4double>  angdist_x;
    std::vector<G4double>  angdist_y;           

	std::vector<double> nu_time;
	std::vector<double> nu_luminosity;
	std::vector<double> nu_meanenergy;
	std::vector<double> nu_meanenergysquare;
	  
private:     
    G4ParticleGun*         ParticleGun;
	void DistFunction(G4double Enu);
    G4double ElectronEnergy(G4double nu_energy, G4double costheta);
    G4double TotalCrossSection(G4double energy);
    G4double fixenergy;
    G4double fixenergy2;
    G4double alpha;

};

#endif
