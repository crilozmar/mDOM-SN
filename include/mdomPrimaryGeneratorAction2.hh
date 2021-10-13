#ifndef mdomPrimaryGeneratorAction2_h
#define mdomPrimaryGeneratorAction2_h 1
 
#include "G4VUserPrimaryGeneratorAction.hh"

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;

class mdomPrimaryGeneratorAction2 : public G4VUserPrimaryGeneratorAction
{
public:
	mdomPrimaryGeneratorAction2(G4ParticleGun*);    
	~mdomPrimaryGeneratorAction2();

	void GeneratePrimaries(G4Event* anEvent);
    G4double InverseCumul(int ControlParameter);  
    G4int                  nPoints0;     //tabulated function
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
    
    G4double ThresholdEnergy;
    G4double Delta;
    G4double y2;
    G4double NTargets;
    G4double me;
    G4double mp;
    G4double mn;
    G4double Gf;
    G4double consg;
    
	std::vector<double> nubar_time;
	std::vector<double> nubar_luminosity;
	std::vector<double> nubar_meanenergy;
	std::vector<double> nubar_meanenergysquare;

private:     
    G4ParticleGun*         ParticleGun;
	void DistFunction(G4double Enu);
    G4double PositronEnergy(G4double nubar_energy, G4double costheta);
    G4double TotalCrossSection(G4double energy);
    G4double fixenergy;
    G4double fixenergy2;
    G4double alpha;
};

#endif
