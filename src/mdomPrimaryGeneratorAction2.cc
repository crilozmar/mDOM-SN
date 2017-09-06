#include "mdomPrimaryGeneratorAction.hh"
#include "mdomPrimaryGeneratorAction2.hh"
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
extern G4String	gnubarfluxname;
extern std::vector<double> readColumnDouble (G4String fn, int col);
extern  G4double gRadius;
extern  G4double gHeight;
extern MdomAnalysisManager gAnalysisManager;

extern G4double	gSNmeanEnergy;
extern G4bool		gfixmeanenergy;
extern G4double 	gfixalpha;

mdomPrimaryGeneratorAction2::mdomPrimaryGeneratorAction2(G4ParticleGun* gun)
: ParticleGun(gun)
{        
  // building energy distribution of electronic antineutrinos...
  //
  if (gfixmeanenergy == false) {
	nubar_time = readColumnDouble(gnubarfluxname, 1);
	nubar_luminosity = readColumnDouble(gnubarfluxname, 2);
	nubar_meanenergy = readColumnDouble(gnubarfluxname, 3);
	nubar_meanenergysquare = readColumnDouble(gnubarfluxname, 4);

	for (unsigned int u = 0; u <nubar_time.size(); u++) {
		nubar_meanenergy[u] = nubar_meanenergy.at(u)*MeV;
		nubar_meanenergysquare[u] = nubar_meanenergysquare.at(u)*MeV*MeV;
		}
		
	  LuminosityDist();
  } else {
	fixenergy = gSNmeanEnergy*MeV;
	alpha = gfixalpha; 
	fixenergy2 = fixenergy*fixenergy*(2+alpha)/(1+alpha); //Only for crosscheck
	ControlParameter = 1;
	Fe_nubar(fixenergy, fixenergy2);
  }
      
  Gf = 1.166e-5*1e-6/(MeV*MeV);
  me = electron_mass_c2;
  mp = proton_mass_c2;
  mn = neutron_mass_c2;
  consg = 1.26;

  
  NTargets = NumberOfTargets(2); //2 protons (hydrogen) per molecule

}



mdomPrimaryGeneratorAction2::~mdomPrimaryGeneratorAction2()
{ }



void mdomPrimaryGeneratorAction2::GeneratePrimaries(G4Event* anEvent)
{
    // Particle and position
  G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("e+");
  ParticleGun->SetParticleDefinition(particle);
  RandomPosition();
  
  beggining:


  //set energy from a tabulated distribution
  //
  G4double timeofspectrum;
  G4double Emean;
  G4double Emean2;
  if (gfixmeanenergy == false) {
	ControlParameter = 0;
	timeofspectrum = InverseCumul(ControlParameter);
	
	G4int timepos = findtime(timeofspectrum);
	Emean = linealinterpolation(timeofspectrum,nubar_time.at(timepos-1), nubar_time.at(timepos), nubar_meanenergy.at(timepos-1),nubar_meanenergy.at(timepos));
	Emean2 = linealinterpolation(timeofspectrum,nubar_time.at(timepos-1), nubar_time.at(timepos), nubar_meanenergysquare.at(timepos-1),nubar_meanenergysquare.at(timepos));
	
	Fe_nubar(Emean, Emean2);
 } else {
	timeofspectrum = 0.0;
	Emean = fixenergy;
	Emean2 = fixenergy2;
  }

  G4double nubar_energy = InverseCumul(ControlParameter);  
  ThresholdEnergy = neutron_mass_c2 + electron_mass_c2 - proton_mass_c2+0.1*MeV; //+0.1MeV because angularcrosssection fails if the energy is too close to the threshold.

  ControlParameter = 1;
  G4int count = 0;
  while (nubar_energy <= ThresholdEnergy) {
    nubar_energy = InverseCumul(ControlParameter); 
    count += 1;
    if (count > 10) {
      goto beggining;
      // To avoid entering in an almost infinity bucle if energy of the chosen time is very unlikely to be higher than 
      // the threshold, go to the beggining and choose a different time of the burst
    }
  }
   
  // angle distribution. We suppose the incident antineutrino would come with momentum direction (0,0,-1)
  ControlParameter = 2;
  DistFunction(nubar_energy);
  
  G4double costheta = InverseCumul(ControlParameter);
  G4double sintheta = std::sqrt(1. - costheta*costheta);
  G4double phi = twopi*G4UniformRand();
  
  G4double zdir = -costheta; 
  G4double xdir = -sintheta*std::cos(phi);
  G4double ydir = -sintheta*std::sin(phi);
  
  // from nu_energy and costheta, we get e- energy
  G4double p_energy = PositronEnergy(nubar_energy, costheta);
  
  ParticleGun->SetParticleEnergy(p_energy); 
  ParticleGun->SetParticleMomentumDirection(G4ThreeVector(xdir,ydir,zdir));
  
  G4double Weigh = WeighMe(nubar_energy);
  
  //sending stuff to analysismanager
  gAnalysisManager.nuTime = timeofspectrum;
  gAnalysisManager.nuMeanEnergy = Emean;
  gAnalysisManager.nuEnergy = nubar_energy;
  gAnalysisManager.cosTheta = costheta;
  gAnalysisManager.primaryEnergy = p_energy; 
  gAnalysisManager.weigh = Weigh;

  
  // G4cout << timeofspectrum<< "      " << nubar_energy/MeV << "        " << p_energy/MeV << "        " << costheta << G4endl;
  //create vertex
  //   
  ParticleGun->GeneratePrimaryVertex(anEvent);
}





void mdomPrimaryGeneratorAction2::RandomPosition() {
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



void mdomPrimaryGeneratorAction2::LuminosityDist()
{
  // Prob Dist to generate neutrinos during the burst time
  nPoints0 = nubar_time.size();
  x0.resize(nPoints0); f0.resize(nPoints0);
  for (G4int j=0; j<nPoints0; j++) {
    x0[j] = nubar_time.at(j); f0[j] = nubar_luminosity.at(j);
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



void mdomPrimaryGeneratorAction2::Fe_nubar(G4double Emean, G4double Emean2)
{
  // Energy distribution
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

G4int mdomPrimaryGeneratorAction2::findtime(G4double time)
{
  // This is a helper function
  for (unsigned int j=0; j<nubar_time.size(); j++) {
    if (time <= nubar_time.at(j)) {
      return j;
    };
  };
 G4cout << "FATAL ERROR -> Not posible to find time of spectrum!!!" << G4endl;
 return 0;
}



G4double mdomPrimaryGeneratorAction2::linealinterpolation(G4double realX,G4double lowerX, G4double upperX, G4double lowerY,G4double upperY) {
	G4double slope = (upperY-lowerY)/(upperX-lowerX);
	G4double result = (slope*(realX-lowerX)+lowerY);
	return result;
}



G4double mdomPrimaryGeneratorAction2::GetAlpha(G4double Emean,G4double Emean2)
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



void mdomPrimaryGeneratorAction2::DistFunction(G4double Enubar)
{
  // Angular distribution
  //
  // P. Vogel, J. F. Beacom. (1999). The angular distribution of the reaction nubar_e + p -> e+ + n. Phys.Rev., D60, 053003 
  // Eq. 14
  //
  nPoints2 = 100;
  

  G4double costhetac = 0.974;
  G4double consf = 1.;
  G4double consf2 = 3.706;
  Delta = mn-mp;
  y2 = (pow(Delta,2)-pow(me,2))/2.;
  G4double M = mp;
  G4double sigma0 = pow(Gf,2)*pow(costhetac,2)/pi*1.024*pow(197.326e-15,2);

  x2.resize(nPoints2); 
  f2.resize(nPoints2);
  
  G4double min = -1.; G4double max = 1.; 
  G4double delta = (max-min)/G4double(nPoints2-1);
  for (G4int i=0; i<nPoints2; i++) {
    x2[i] = min + i*delta; //costheta
  }
  
    G4double Ee0 = Enubar - Delta;
    G4double pe0 = pow((pow(Ee0,2)-pow(me,2)),1./2.);
    G4double ve0 = pe0/Ee0;
  
  for (G4int j=0; j<nPoints2; j++) {
    G4double Ee1 = Ee0*(1.-Enubar/M*(1.-ve0*x2[j]))-y2/M;
    G4double pe1 = pow((pow(Ee1,2)-pow(me,2)),1./2.);
    G4double ve1 = pe1/Ee1;
    G4double part1 = 2.*(consf+consf2)*consg*((2.*Ee0+Delta)*(1.-ve0*x2[j])-pow(me,2)/Ee0);
    G4double part2 = (pow(consf,2)+pow(consg,2))*(Delta*(1.+ve0*x2[j])+pow(me,2)/Ee0);
    G4double part3 = (pow(consf,2)+3.*pow(consg,2))*((Ee0+Delta)*(1.- x2[j]/ve0)-Delta);
    G4double part4 = (pow(consf,2)+pow(consg,2))*((Ee0+Delta)*(1.-x2[j]/ve0)-Delta)*ve0*x2[j];
    G4double Tau = part1 + part2 + part3 + part4;

    f2[j] = sigma0/2.*((pow(consf,2)+3.*pow(consg,2))+(pow(consf,2)-pow(consg,2))*ve1*x2[j])*Ee1*pe1-sigma0/2.*(Tau/M)*Ee0*pe0;// dsigma/dcos(theta)
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




G4double mdomPrimaryGeneratorAction2::PositronEnergy(G4double Enubar, G4double costheta)
{
  // Get positron energy from inverse beta decay as a function of incident antineutrino energy and scatter angle
  //
  // P. Vogel, J. F. Beacom. (1999). The angular distribution of the reaction nu_e + p -> e+ + n. Phys.Rev., D60, 053003 
  // Eq. 13
  //
  G4double Ee0 = Enubar - Delta;
  G4double pe0 = pow((pow(Ee0,2)-pow(electron_mass_c2,2)),1./2.);
  G4double ve0 = pe0/Ee0;
  G4double Energy = Ee0*(1.-Enubar/proton_mass_c2*(1.-ve0*costheta))-y2/proton_mass_c2;
  return Energy;
}





G4double mdomPrimaryGeneratorAction2::NumberOfTargets(G4int targetPerMolecule) {
  //Just to calculate number of targets in ice per cubic meter, assuming ice as H2O pure
  G4double Density = 921.6*kg/(m*m*m); //Density of ice at -50 celsius degrees
  G4double MolarMass = 18.01528e-3*kg; //kg per mol
  G4double Na = 6.022140857e23;
  G4double Nm = Density/MolarMass*Na;//molecules/m^3 of ice
  G4double Nt = Nm*targetPerMolecule;
  return Nt;
}


G4double mdomPrimaryGeneratorAction2::TotalCrossSection(G4double energy) {
  // Returns value of the TotalCrossSection for certain energy to use it in WeighMe
  //
  // T. Totani, K. Sato, H. E. Dalhed, J. R. Wilson,Future detection of supernova neutrino burst and explosion mechanism, Astrophys. J. 496 ,1998 216–225, Preprint astro-ph/9710203, 1998, Equation 9
  G4double hbar = 6.58211899e-16*1e-6*MeV*s;
  G4double c = 3e8*m/s;
  G4double sigma0 = pow(2.*Gf*me*hbar*c,2)/pi;
  G4double constante = -0.00325/(MeV*MeV);
  G4double deltaWM = constante*(energy-(mn-mp)/2.);
  G4double Ee = energy+mp-mn;
  G4double pec = pow((pow(Ee,2)-pow(me,2)),1./2.);
  G4double sigma = 1./4.*sigma0*(1.+3.*pow(consg,2))*(1.+deltaWM)*Ee*pec/pow(me,2);
  //G4cout << "CRRROOOOOOOOSSSSSEEEEECCCCTTTIIIIIOOOOONNNN   " << sigma/(m*m) << " of energy "<< energy/MeV << G4endl;
  return sigma;
}



G4double mdomPrimaryGeneratorAction2::WeighMe(G4double energy) {
  // GIve the weigh of the interaction because of the cross section as:
  // Weigh = sigma * Nt * r, where r is the distance of the neutrino in the ice 
  // and Nt is the number of target in ice (H nucleus in this case)
  //
  // THIS IS THOUGH FOR SN ANTINEUTRINOS COMING FROM Z AXIS WITH A WORLD AS A CYLINDER, if you changed some of that you will also have to change this function
  
  G4double Sigma = TotalCrossSection(energy);
  G4double weigh = Sigma*NTargets*(2*gHeight);
  //G4cout << "Peso -> " << weigh << "  con Ntargets  " << NTargets/(1/m3) << "  y gHeight " << 2*gHeight/m << G4endl;
  return weigh;
}


// InverseCumul gives a random value based in a given distribution
G4double mdomPrimaryGeneratorAction2::InverseCumul(int controlparameter)
{
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
