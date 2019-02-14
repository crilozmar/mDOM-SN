#include "mdomPrimaryGeneratorAction.hh"
#include "mdomPrimaryGeneratorAction3.hh"
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


extern double gworldsize;
extern G4String	gSunspectrum;
extern std::vector<double> readColumnDouble (G4String fn, int col);
extern  G4double gRadius;
extern  G4double gHeight;
extern G4bool gSun_nue;
extern G4bool gSun_numutau;
extern MdomAnalysisManager gAnalysisManager;

//This class is for the solar neutrinos, which are only elastic scattering (my assumption) but the model does not depend on time

mdomPrimaryGeneratorAction3::mdomPrimaryGeneratorAction3(G4ParticleGun* gun)
: ParticleGun(gun)
{ 
  // building energy distribution of electronic neutrinos...
  //
  
  Gf = 1.166e-5*1e-6/(MeV*MeV);
  me = electron_mass_c2;
  
  NTargets = NumberOfTargets(10); //10 electrons per molecule
  
  Fe_build();

}



mdomPrimaryGeneratorAction3::~mdomPrimaryGeneratorAction3()
{ }



void mdomPrimaryGeneratorAction3::GeneratePrimaries(G4Event* anEvent)
{
  // Particle and position
  G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  ParticleGun->SetParticleDefinition(particle);

  mdomSNTools SNToolBox;
  G4ThreeVector Position;
  SNToolBox.RandomPosition(Position);
  ParticleGun->SetParticlePosition(Position);
  //set energy from a tabulated distribution
  //
  beggining:
  
  ControlParameter = 1;
  G4double nu_energy = InverseCumul(ControlParameter);  
  
  G4int count = 0;
  while (nu_energy <= 0.1*MeV) {
         //G4cout << "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ" << G4endl;
  nu_energy = InverseCumul(ControlParameter); 
  count += 1;
  if (count > 10) {
      goto beggining;
      // It would not have sense to have a energy below that for the neutrino, and angular distribution can fail. The counter is to avoid entering in a infinite bucle
    }
  } 
  // angle distribution. We suppose the incident neutrino would come with momentum direction (0,0,-1) 
  // THIS IS IMPORTANT: IF CHANGE THE DIRECTION OF THE SIMULATED SUPERNOVA, YOU HAVE TO CHANGE ALSO THE WEIGH FUNCTION!
  ControlParameter = 2;
  DistFunction(nu_energy);
  
  G4double costheta = InverseCumul(ControlParameter);
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
  
  G4double Weigh = WeighMe(nu_energy);
  
  //sending stuff to analysismanager
  gAnalysisManager.nuTime = -999.99*s; //This is just to use the same AnalysisManager and my same python script than for the SNe....
  gAnalysisManager.nuMeanEnergy = -999.99*MeV; //This is just to use the same AnalysisManager and my same python script than for the SNe....
  gAnalysisManager.nuEnergy = nu_energy;
  gAnalysisManager.cosTheta = costheta;
  gAnalysisManager.primaryEnergy = e_energy; 
  gAnalysisManager.weigh = Weigh;

  
  //G4cout << "?????????????????????????????????????" << Weigh << G4endl;
  //G4cout << nu_energy/MeV << "        " << e_energy/MeV << "        " << costheta << G4endl;
  //create vertex
  // 
  ParticleGun->GeneratePrimaryVertex(anEvent);
}


void mdomPrimaryGeneratorAction3::Fe_build()
{// Energy distribution
  // tabulated function 
  // f is assumed positive, linear per segment, continuous

    // Fe(E) corresponding to solar neutrino energies
    
    // Data and formulas from:
    //  Phys. Rev. C, 54, 411 (1996)
   //  http://www.sns.ias.edu/~jnb/
    
  data_energy = readColumnDouble(gSunspectrum, 1);
  data_fe = readColumnDouble(gSunspectrum, 2);
  nPoints1 = data_energy.size();
  x1.resize(nPoints1); f1.resize(nPoints1);

  for (unsigned int u = 0; u <data_energy.size(); u++) {
    x1[u] = data_energy.at(u)*MeV;
    f1[u] = data_fe.at(u);
  }
  //compute slopes
  //
  a1.resize(nPoints1);
  for (G4int j=0; j<nPoints1-1; j++) { 
    a1[j] = (f1[j+1] - f1[j])/(x1[j+1] - x1[j]);
  };
  //compute cumulative function
  //
  Fc1.resize(nPoints1);  
  Fc1[0] = 0.;
  for (G4int j=1; j<nPoints1; j++) {
    Fc1[j] = Fc1[j-1] + 0.5*(f1[j] + f1[j-1])*(x1[j] - x1[j-1]);
  };     
}




void mdomPrimaryGeneratorAction3::DistFunction(G4double Enu)
{
  // Angular distribution
  //
  // Carlo Giunti and Chung W.Kim (2007), Fundamentals of Neutrino Physics and Astrophysics, Oxford University Press
  // Chapter 5, eq. 5.29

  nPoints2 = 500;
  G4double g1;
  if (gSun_nue == true) {
	g1=0.73;
  }else if (gSun_numutau == true) {
	g1=-0.27;
  }
  G4double g2=0.23;
  G4double sigma0=2.*pow(Gf,2)*pow(me,2)/pi*pow(197.326e-15,2); 
  
  x2.resize(nPoints2); 
  f2.resize(nPoints2);
  
  G4double min = 0.; G4double max = 1.; 
  G4double delta = (max-min)/G4double(nPoints2-1);
  for (G4int i=0; i<nPoints2; i++) {
    x2[i] = min + i*delta; //costheta
  }
  
  for (G4int j=0; j<nPoints2; j++) {
    G4double dem = pow((pow((me+Enu),2)-pow(Enu,2)*pow(x2[j],2)),2);
    G4double factor1 = 4.*pow(Enu,2)*pow((me+Enu),2)*x2[j]/dem;
    G4double factor2 = 2.*me*Enu*pow(x2[j],2)/dem;
    G4double factor3 = 2.*pow(me,2)*pow(x2[j],2)/dem;
    f2[j] = sigma0*factor1*(pow(g1,2)+pow(g2,2)*pow((1.-factor2),2)-g1*g2*factor3); // dsigma/dcos(theta)
  }
  //compute slopes
  //
  a2.resize(nPoints2);
  for (G4int j=0; j<nPoints2-1; j++) { 
    a2[j] = (f2[j+1] - f2[j])/(x2[j+1] - x2[j]);
  };
  //compute cumulative function
  //
  Fc2.resize(nPoints2);  
  Fc2[0] = 0.;
  for (G4int j=1; j<nPoints2; j++) {
    Fc2[j] = Fc2[j-1] + 0.5*(f2[j] + f2[j-1])*(x2[j] - x2[j-1]);
  };     
}




G4double mdomPrimaryGeneratorAction3::ElectronEnergy(G4double nu_energy, G4double costheta)
{
  // Get electron energy from elastic scattering as a function of incident neutrino energy and scatter angle
  //
  // Carlo Giunti and Chung W.Kim (2007), Fundamentals of Neutrino Physics and Astrophysics, Oxford University Press
  // Chapter 5, eq. 5.27
  G4double energy=2*me*pow(nu_energy,2)*pow(costheta,2)/(pow((me+nu_energy),2)-pow(nu_energy,2)*pow(costheta,2));
  return energy;
}



G4double mdomPrimaryGeneratorAction3::NumberOfTargets(G4int targetPerMolecule) {
  //Just to calculate number of targets in ice per cubic meter, assuming ice as H2O pure
  G4double Density = 921.6*kg/m3; //Density of ice at -50 celsius degrees
  G4double MolarMass = 18.01528e-3*kg; //kg per mol
  G4double Na = 6.022140857e23;
  G4double Nm = Density/MolarMass*Na;//molecules/m^3 of ice
  G4double Nt = Nm*targetPerMolecule;
  return Nt;
}


G4double mdomPrimaryGeneratorAction3::TotalCrossSection(G4double energy) {
  // Returns value of the TotalCrossSection for certain energy to use it in WeighMe
  //
  //M. Buchkremer, Electroweak Interactions: Neutral currents in neutrino-lepton elastic
  // scattering experiments, Universit ́e Catholique de Louvain /CP3, 2011.
  G4double sigma;
  G4double sin2thetaw = 0.231;
  if (gSun_nue == true) {
	sigma = pow(Gf,2)*me*energy/(2.*pi)*(1.+4.*sin2thetaw+16./3.*pow(sin2thetaw,2))*pow(197.326e-15,2)*m*m;
  } else if (gSun_numutau == true)  {
	sigma = pow(Gf,2)*me*energy/(2.*pi)*(1.-4.*sin2thetaw+16./3.*pow(sin2thetaw,2))*pow(197.326e-15,2)*m*m;
  }
  //G4cout << "CRRROOOOOOOOSSSSSEEEEECCCCTTTIIIIIOOOOONNNN   " << sigma/(m*m) << " of energy "<< energy/MeV << G4endl;
  return sigma;
}



G4double mdomPrimaryGeneratorAction3::WeighMe(G4double energy) {
  // GIve the weigh of the interaction because of the cross section as:
  // Weigh = sigma * Nt * r, where r is the distance of the neutrino in the ice 
  // and Nt is the number of target in ice (electrons in this case)
  //
  // THIS IS THOUGH FOR NEUTRINOS COMING FROM Z AXIS WITH A WORLD AS A CYLINDER, if you changed some of that you will also have to change this function
  
  G4double Sigma = TotalCrossSection(energy);
  double weigh = Sigma*NTargets*(2*gHeight);
  //G4cout << "Peso -> " << weigh << "  con Ntargets  " << NTargets/(1/m3) << "  y gHeight " << 2*gHeight/m << G4endl;
  return weigh;
}



G4double mdomPrimaryGeneratorAction3::InverseCumul(int controlparameter)
{
  // InverseCumul gives a random value based in a given distribution
  // tabulated function
  // f is assumed positive, linear per segment, continuous 
  // --> cumulative function is second order polynomial
  // --> ControlParameter == 1 -->> Energy distribution
  // --> ControlParameter == 2 -->> Angle distribution
 if (controlparameter == 1) {
    x = x1;
    f = f1;
    a = a1;
    nPoints = nPoints1;
    Fc = Fc1;
  } else if (controlparameter == 2) {
    x = x2;
    f = f2;
    a = a2;
    nPoints = nPoints2;
    Fc = Fc2;
  } else {
    G4cout << "ERROR --> INVALID CONTROL PARAMETER!" << G4endl;
    return 0;
  }
  
  //choose y randomly
  G4double y_rndm = G4UniformRand()*Fc[nPoints-1];
  //find bin
  G4int j = nPoints-2;
  while ((Fc[j] > y_rndm) && (j > 0)) j--;
  //y_rndm --> x_rndm :  Fc(x) is second order polynomial
  G4double x_rndm = x[j];
  G4double aa = a[j];
  if (aa != 0.) {
    G4double b = f[j]/aa, c = 2*(y_rndm - Fc[j])/aa;
    G4double delta = b*b + c;
    G4int sign = 1; if (aa < 0.) sign = -1;
    x_rndm += sign*std::sqrt(delta) - b;    
  } else if (f[j] > 0.) {
    x_rndm += (y_rndm - Fc[j])/f[j];
  };
  return x_rndm;
}


