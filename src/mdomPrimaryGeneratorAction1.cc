/** @file mdomPrimaryGeneratorAction1.cc
 *  @brief Primary generator. Produces electrons coming from elastic scattering of SN neutrinos
 * 
 *  @author Cristian Jesus Lozano Mariscal (c.lozano@wwu.de)
 * 
 *  @version Geant4 10.7
 */

#include "mdomPrimaryGeneratorAction.hh"
#include "mdomPrimaryGeneratorAction1.hh"
#include "mdomDetectorConstruction.hh"
#include "mdomAnalysisManager.hh"


#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTypes.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "mdomSNTools.hh"

extern G4String	gnufluxname;
extern std::vector<double> readColumnDouble (G4String fn, int col);
extern  G4double gRadius;
extern  G4double gHeight;
extern MdomAnalysisManager gAnalysisManager;
extern G4double	gSNmeanEnergy;
extern G4bool	gfixmeanenergy;
extern G4double 	gfixalpha;

mdomPrimaryGeneratorAction1::mdomPrimaryGeneratorAction1(G4ParticleGun* gun)
: ParticleGun(gun)
{ 
  // building energy distribution of electronic neutrinos...
  //
  if (gfixmeanenergy == false) {
	nu_time = readColumnDouble(gnufluxname, 1);
	nu_luminosity = readColumnDouble(gnufluxname, 2); //this is not luminosity but flux, I do not want to change it everywhere...
	nu_meanenergy = readColumnDouble(gnufluxname, 3);
	nu_meanenergysquare = readColumnDouble(gnufluxname, 4);
	for (unsigned int u = 0; u <nu_time.size(); u++) {
		nu_time[u] = nu_time.at(u)*s;
		nu_meanenergy[u] = nu_meanenergy.at(u)*MeV;
		nu_meanenergysquare[u] = nu_meanenergysquare.at(u)*MeV*MeV;
	}
	// Since the luminosity spectrum is not gonna change, it is worthy to compute already the slopes and store them
	nPoints_lum =  nu_time.size();
	GetSlopes(nu_time,  nu_luminosity, nPoints_lum, x_lum, f_lum, a_lum, Fc_lum);

  } else { //save all slopes and stuff since it always gonna be the same
	fixenergy = gSNmeanEnergy*MeV;
	alpha = gfixalpha; 
	fixenergy2 = fixenergy*fixenergy*(2+alpha)/(1+alpha); //Only for crosscheck
	std::vector<G4double> x1;
	std::vector<G4double> f1;
	fixE_nPoints = 500;
	MakeEnergyDistribution(fixenergy, alpha, fixE_nPoints, x1, f1);
	GetSlopes(x1, f1, fixE_nPoints, fixFe_X, fixFe_Y, fixFe_a, fixFe_Fc);
  }
  
  Gf = 1.166e-5*1e-6/(MeV*MeV);
  me = electron_mass_c2;
  
  NTargets = NumberOfTargets(10); //10 electrons per molecule
  
}



mdomPrimaryGeneratorAction1::~mdomPrimaryGeneratorAction1()
{ }



void mdomPrimaryGeneratorAction1::GeneratePrimaries(G4Event* anEvent)
{
  // Particle and position
  G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  ParticleGun->SetParticleDefinition(particle);
  
  mdomSNTools SNToolBox;
  
  G4ThreeVector Position;
  SNToolBox.RandomPosition(Position);
  ParticleGun->SetParticlePosition(Position);

  
  beggining:
  
  G4double timeofspectrum;
  G4double Emean;
  G4double Emean2;
  G4double nu_energy;
  if (gfixmeanenergy == false) {
  //set energy from a tabulated distribution
  //
	timeofspectrum = SNToolBox.InverseCumulAlgorithm(x_lum, f_lum, a_lum, Fc_lum,nPoints_lum);
	
	G4int timepos = SNToolBox.findtime(timeofspectrum, nu_time);
	Emean = SNToolBox.linealinterpolation(timeofspectrum,nu_time.at(timepos-1), nu_time.at(timepos), nu_meanenergy.at(timepos-1),nu_meanenergy.at(timepos));
	Emean2 = SNToolBox.linealinterpolation(timeofspectrum,nu_time.at(timepos-1), nu_time.at(timepos), nu_meanenergysquare.at(timepos-1),nu_meanenergysquare.at(timepos));
	
	nu_energy = 0;
	G4int count = 0;
	while (nu_energy <= 0.1*MeV) {
	nu_energy = SNToolBox.EnergyDistribution(Emean, Emean2, alpha);
	count += 1;
	if (count > 10) {
		goto beggining;
      // It would not have sense to have a energy below that for the neutrino, and angular distribution can fail. The counter 	is to avoid entering in a infinite bucle
		}
	}
  } else {
	timeofspectrum = 0.0;
	Emean = fixenergy;
	Emean2 = fixenergy2;
	nu_energy =  SNToolBox.InverseCumulAlgorithm(fixFe_X, fixFe_Y, fixFe_a, fixFe_Fc, fixE_nPoints);
  }
  
  // THIS IS IMPORTANT: IF CHANGE THE DIRECTION OF THE SIMULATED SUPERNOVA, YOU HAVE TO CHANGE ALSO THE WEIGH FUNCTION!
  DistFunction(nu_energy);
  
  G4double costheta = SNToolBox.InverseCumul(angdist_x, angdist_y, angdist_nPoints);
  G4double sintheta = std::sqrt(1. - costheta*costheta);
  G4double phi = twopi*G4UniformRand();
  
  G4double zdir = -costheta; 
  G4double xdir = -sintheta*std::cos(phi);
  G4double ydir = -sintheta*std::sin(phi);
  
  // G4cout << xdir << "  " << ydir << "  " << zdir << G4endl;
  
  // from nu_energy and costheta, we get e- energy
  G4double e_energy = ElectronEnergy(nu_energy, costheta);
  
  ParticleGun->SetParticleEnergy(e_energy); 
  ParticleGun->SetParticleMomentumDirection(G4ThreeVector(xdir,ydir,zdir));
  
  // now calculate the weight
  G4double sigma = TotalCrossSection(nu_energy);
  G4double Weigh = SNToolBox.WeighMe(sigma, NTargets);
  
  //sending stuff to analysismanager
  gAnalysisManager.nuTime = timeofspectrum;
  gAnalysisManager.nuMeanEnergy = Emean;
  gAnalysisManager.nuEnergy = nu_energy;
  gAnalysisManager.cosTheta = costheta;
  gAnalysisManager.primaryEnergy = e_energy; 
  gAnalysisManager.weigh = Weigh;

  
  //G4cout << "?????????????????????????????????????" << Weigh << G4endl;
  // G4cout << timeofspectrum<< "      " << nu_energy/MeV << "        " << e_energy/MeV << "        " << costheta << G4endl;
  //create vertex
  // 
  ParticleGun->GeneratePrimaryVertex(anEvent);
}




void mdomPrimaryGeneratorAction1::DistFunction(G4double Enu)
{
  // Angular distribution
  //
  // Carlo Giunti and Chung W.Kim (2007), Fundamentals of Neutrino Physics and Astrophysics, Oxford University Press
  // Chapter 5, eq. 5.29

  angdist_nPoints = 500;
  
  G4double g1=0.73;
  G4double g2=0.23;
  G4double sigma0=2.*pow(Gf,2)*pow(me,2)/pi*pow(197.326e-15,2); 
  
  angdist_x.resize(angdist_nPoints); 
  angdist_y.resize(angdist_nPoints);
  
  G4double min = 0.; G4double max = 1.; 
  G4double delta = (max-min)/G4double(angdist_nPoints-1);
  for (G4int i=0; i<angdist_nPoints; i++) {
    angdist_x[i] = min + i*delta; //costheta
  }
  
  for (G4int j=0; j<angdist_nPoints; j++) {
    G4double dem = pow((pow((me+Enu),2)-pow(Enu,2)*pow(angdist_x[j],2)),2);
    G4double factor1 = 4.*pow(Enu,2)*pow((me+Enu),2)*angdist_x[j]/dem;
    G4double factor2 = 2.*me*Enu*pow(angdist_x[j],2)/dem;
    G4double factor3 = 2.*pow(me,2)*pow(angdist_x[j],2)/dem;
    angdist_y[j] = sigma0*factor1*(pow(g1,2)+pow(g2,2)*pow((1.-factor2),2)-g1*g2*factor3); // dsigma/dcos(theta)
  }
}



G4double mdomPrimaryGeneratorAction1::ElectronEnergy(G4double nu_energy, G4double costheta)
{
  // Get electron energy from elastic scattering as a function of incident neutrino energy and scatter angle
  //
  // Carlo Giunti and Chung W.Kim (2007), Fundamentals of Neutrino Physics and Astrophysics, Oxford University Press
  // Chapter 5, eq. 5.27
  G4double energy=2*me*pow(nu_energy,2)*pow(costheta,2)/(pow((me+nu_energy),2)-pow(nu_energy,2)*pow(costheta,2));
  return energy;
}



G4double mdomPrimaryGeneratorAction1::TotalCrossSection(G4double energy) {
  // Returns value of the TotalCrossSection for certain energy to use it in WeighMe
  //
  //M. Buchkremer, Electroweak Interactions: Neutral currents in neutrino-lepton elastic
  // scattering experiments, Universit ??e Catholique de Louvain /CP3, 2011.
  G4double sin2thetaw = 0.231;
  G4double sigma = pow(Gf,2)*me*energy/(2.*pi)*(1.+4.*sin2thetaw+16./3.*pow(sin2thetaw,2))*pow(197.326e-15,2)*m*m;
  //G4cout << "CRRROOOOOOOOSSSSSEEEEECCCCTTTIIIIIOOOOONNNN   " << sigma/(m*m) << " of energy "<< energy/MeV << G4endl;
  return sigma;
}

