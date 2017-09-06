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

extern double gworldsize;
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
	nu_luminosity = readColumnDouble(gnufluxname, 2);
	nu_meanenergy = readColumnDouble(gnufluxname, 3);
	nu_meanenergysquare = readColumnDouble(gnufluxname, 4);

	for (unsigned int u = 0; u <nu_time.size(); u++) {
		nu_meanenergy[u] = nu_meanenergy.at(u)*MeV;
		nu_meanenergysquare[u] = nu_meanenergysquare.at(u)*MeV*MeV;
	}
  
	LuminosityDist();
  } else {
	fixenergy = gSNmeanEnergy*MeV;
	alpha = gfixalpha; 
	fixenergy2 = fixenergy*fixenergy*(2+alpha)/(1+alpha); //Only for crosscheck
	ControlParameter = 1;
	Fe_nu(fixenergy, fixenergy2);
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

  RandomPosition();
  
  
  beggining:
  
  G4double timeofspectrum;
  G4double Emean;
  G4double Emean2;
  if (gfixmeanenergy == false) {
  //set energy from a tabulated distribution
  //
	ControlParameter = 0;
	timeofspectrum = InverseCumul(ControlParameter);
	
	G4int timepos = findtime(timeofspectrum);
	Emean = linealinterpolation(timeofspectrum,nu_time.at(timepos-1), nu_time.at(timepos), nu_meanenergy.at(timepos-1),nu_meanenergy.at(timepos));
	Emean2 = linealinterpolation(timeofspectrum,nu_time.at(timepos-1), nu_time.at(timepos), nu_meanenergysquare.at(timepos-1),nu_meanenergysquare.at(timepos));
	
	Fe_nu(Emean, Emean2);
  } else {
	timeofspectrum = 0.0;
	Emean = fixenergy;
	Emean2 = fixenergy2;
  }
  
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




void mdomPrimaryGeneratorAction1::RandomPosition() {
  // Just to give a random coordinates for the primaries inside the Ice
    /*
  G4double Pos_cosTheta = 2*G4UniformRand() - 1;  //cosTheta uniform in [0, pi]
  G4double Pos_sinTheta = std::sqrt(1. - Pos_cosTheta*Pos_cosTheta);
  G4double Pos_phi      = twopi*G4UniformRand();  //phi uniform in [0, 2*pi]
  G4ThreeVector ur(Pos_sinTheta*std::cos(Pos_phi),Pos_sinTheta*std::sin(Pos_phi),Pos_cosTheta);
  */
  G4double Rmin = (356./2.+27.5+1)*mm; 
  G4double Rmax = pow(3,1./2.)*gworldsize*m;
  G4double Rmax2 = gworldsize*m;
  /*G4double fRmin3 = Rmin*Rmin*Rmin;
  G4double fRmax3 = Rmax*Rmax*Rmax;
  
  G4double R3 = fRmin3 + G4UniformRand()*(fRmax3 - fRmin3);
  G4double R  = std::pow(R3, 1./3);  
  
  G4ThreeVector Position = R*ur;

  if (Position[0] > gRadius || Position[1] > gRadius || Position[2] > gHeight) {
    RandomPosition();
  } else {

  //G4cout << gRadius/m << ", " << gHeight/m << G4endl;
  //G4cout << Position[0]/m << ", " << Position[1]/m << ", " << Position[2]/m << G4endl;
    */
  G4double posornegX = 1;
  if (G4UniformRand()<0.5) { posornegX = -1;}
    G4double posornegY = 1;
  if (G4UniformRand()<0.5) { posornegY = -1;}
    G4double posornegZ = 1;
  if (G4UniformRand()<0.5) { posornegZ = -1;}
  G4double posz = posornegZ*(G4UniformRand()*gHeight);
  G4double posx = posornegX*(G4UniformRand()*gRadius);
  G4double posy = posornegY*(G4UniformRand()*gRadius);
  G4ThreeVector Position(posx,posy,posz);
  G4double R3 = pow(pow(Position[0],2)+pow(Position[1],2)+pow(Position[2],2),1./2.);
  G4double R2 = pow(pow(Position[0],2)+pow(Position[1],2),1./2.);
  if (( R3 <= Rmin ) || (R3 >= Rmax) || (R2 >= Rmax2)) {
     RandomPosition();
  } else {
  ParticleGun->SetParticlePosition(Position);
  }
}



void mdomPrimaryGeneratorAction1::LuminosityDist()
{
  // Prob Dist to generate neutrinos during the burst time
  nPoints0 = nu_time.size();
  
  x0.resize(nPoints0); f0.resize(nPoints0);
  for (G4int j=0; j<nPoints0; j++) {
    x0[j] = nu_time.at(j); f0[j] = nu_luminosity.at(j);
  };
  //compute slopes
  //
  a0.resize(nPoints0);
  for (G4int j=0; j<nPoints0-1; j++) { 
    a0[j] = (f0[j+1] - f0[j])/(x0[j+1] - x0[j]);
  };
  //compute cumulative function
  //
  Fc0.resize(nPoints0);  
  Fc0[0] = 0.;
  for (G4int j=1; j<nPoints0; j++) {
    Fc0[j] = Fc0[j-1] + 0.5*(f0[j] + f0[j-1])*(x0[j] - x0[j-1]);
  };     
}


void mdomPrimaryGeneratorAction1::Fe_nu(G4double Emean, G4double Emean2)
{// Energy distribution
  // tabulated function 
  // f is assumed positive, linear per segment, continuous
    nPoints1 = 500;

    // Fe(E) corresponding to neutrino energies from a heavy neutron star, ls220
    
    // Data and formulas from:
    //  Irene Tamborra et al., High-resolution supernova neutrino spectra represented by a simple fit, PHYSICAL REVIEW D 86, 125031 (2012)
   //
   if (gfixmeanenergy == false) {
	alpha = GetAlpha(Emean, Emean2);
   } 
   
  G4double min = 0.; G4double max = 50.; 
  G4double delta = (max-min)/G4double(nPoints1-1);
  x1.resize(nPoints1); f1.resize(nPoints1);

  for (G4int i=0; i<nPoints1; i++) {
    x1[i] = (min + i*delta)*MeV; //Energy
  }

  for (G4int j=0; j<nPoints1; j++) {
  f1[j] = pow(x1[j],alpha)*exp(-(alpha+1.)*x1[j]/Emean); // F(e), energy dist. function
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

/// This is a helper function
G4int mdomPrimaryGeneratorAction1::findtime(G4double time)
{
  for (unsigned int j=0; j<nu_time.size(); j++) {
    if (time <= nu_time.at(j)) {
      return j;
    };
  };
 G4cout << "FATAL ERROR -> Not posible to find time of spectrum!!!" << G4endl;
 return 0;
}


G4double mdomPrimaryGeneratorAction1::linealinterpolation(G4double realX,G4double lowerX, G4double upperX, G4double lowerY,G4double upperY) {
	G4double slope = (upperY-lowerY)/(upperX-lowerX);
	G4double result = (slope*(realX-lowerX)+lowerY);
	return result;
}



G4double mdomPrimaryGeneratorAction1::GetAlpha(G4double Emean,G4double Emean2)
{
  // Get Alpha Parameter for the energy distribution
  G4double factor1 = Emean2/pow(Emean,2);
  G4double factor2;
  std::vector<G4double> div;
  G4int sizeloop = 100000;
  div.resize(sizeloop);
  
  G4double min = 0.; G4double max = 40.; 
  G4double delta = (max-min)/G4double(sizeloop);
  for (G4int i=0; i<sizeloop-1; i++) {
    div[i] = min + i*delta;
  };
  for (G4int k=0; k<sizeloop; k++) {
	  factor2 = (2.+div[k])/(1.+div[k]);
	  if (factor2 <= factor1) {
		return div[k];
	  }
  }
  G4cout << "FATAL ERROR -> Alpha Parameter. Not posible to get the energy distribution!!!" << G4endl;
  return 0;
}




void mdomPrimaryGeneratorAction1::DistFunction(G4double Enu)
{
  // Angular distribution
  //
  // Carlo Giunti and Chung W.Kim (2007), Fundamentals of Neutrino Physics and Astrophysics, Oxford University Press
  // Chapter 5, eq. 5.29

  nPoints2 = 500;
  
  G4double g1=0.73;
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




G4double mdomPrimaryGeneratorAction1::ElectronEnergy(G4double nu_energy, G4double costheta)
{
  // Get electron energy from elastic scattering as a function of incident neutrino energy and scatter angle
  //
  // Carlo Giunti and Chung W.Kim (2007), Fundamentals of Neutrino Physics and Astrophysics, Oxford University Press
  // Chapter 5, eq. 5.27
  G4double energy=2*me*pow(nu_energy,2)*pow(costheta,2)/(pow((me+nu_energy),2)-pow(nu_energy,2)*pow(costheta,2));
  return energy;
}



G4double mdomPrimaryGeneratorAction1::NumberOfTargets(G4int targetPerMolecule) {
  //Just to calculate number of targets in ice per cubic meter, assuming ice as H2O pure
  G4double Density = 921.6*kg/m3; //Density of ice at -50 celsius degrees
  G4double MolarMass = 18.01528e-3*kg; //kg per mol
  G4double Na = 6.022140857e23;
  G4double Nm = Density/MolarMass*Na;//molecules/m^3 of ice
  G4double Nt = Nm*targetPerMolecule;
  return Nt;
}


G4double mdomPrimaryGeneratorAction1::TotalCrossSection(G4double energy) {
  // Returns value of the TotalCrossSection for certain energy to use it in WeighMe
  //
  //M. Buchkremer, Electroweak Interactions: Neutral currents in neutrino-lepton elastic
  // scattering experiments, Universit ́e Catholique de Louvain /CP3, 2011.
  G4double sin2thetaw = 0.231;
  G4double sigma = pow(Gf,2)*me*energy/(2.*pi)*(1.+4.*sin2thetaw+16./3.*pow(sin2thetaw,2))*pow(197.326e-15,2)*m*m;
  //G4cout << "CRRROOOOOOOOSSSSSEEEEECCCCTTTIIIIIOOOOONNNN   " << sigma/(m*m) << " of energy "<< energy/MeV << G4endl;
  return sigma;
}



G4double mdomPrimaryGeneratorAction1::WeighMe(G4double energy) {
  // GIve the weigh of the interaction because of the cross section as:
  // Weigh = sigma * Nt * r, where r is the distance of the neutrino in the ice 
  // and Nt is the number of target in ice (electrons in this case)
  //
  // THIS IS THOUGH FOR SN NEUTRINOS COMING FROM Z AXIS WITH A WORLD AS A CYLINDER, if you changed some of that you will also have to change this function
  
  G4double Sigma = TotalCrossSection(energy);
  double weigh = Sigma*NTargets*(2*gHeight);
  //G4cout << "Peso -> " << weigh << "  con Ntargets  " << NTargets/(1/m3) << "  y gHeight " << 2*gHeight/m << G4endl;
  return weigh;
}



G4double mdomPrimaryGeneratorAction1::InverseCumul(int controlparameter)
{
  // InverseCumul gives a random value based in a given distribution
  // tabulated function
  // f is assumed positive, linear per segment, continuous 
  // --> cumulative function is second order polynomial
  // --> ControlParameter == 0 -->> LuminosityDIst
  // --> ControlParameter == 1 -->> Energy distribution
  // --> ControlParameter == 2 -->> Angle distribution
  if (controlparameter == 0) {
    x = x0;
    f = f0;
    a = a0;
    nPoints = nPoints0;
    Fc = Fc0;
  } else if (controlparameter == 1) {
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


