#include "mdomSNTools.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4VTouchable.hh"
#include <stdlib.h>
#include "G4Navigator.hh"
#include "G4TouchableHistoryHandle.hh"

extern double gworldsize;
extern  G4double gRadius;
extern  G4double gHeight;
extern G4bool gfixmeanenergy;
extern G4Navigator* aNavigator;
extern G4double	gmdomseparation;
extern G4int	gn_mDOMs;

mdomSNTools::mdomSNTools(){
	/*
	Rmin = (356./2.+27.5+1)*mm;  //mdom size along its bigger part
	mdompos.resize(gn_mDOMs);
	for (int k = 0; k < gn_mDOMs; k++){
		G4double zpos;
		if (gn_mDOMs % 2 == 0) {
			zpos = gmdomseparation*(gn_mDOMs/2-k-1./2.);
		} else {
			zpos = gmdomseparation*(gn_mDOMs/2-k);
		}
		mdompos.at(k) = zpos;
	}*/
}

mdomSNTools::~mdomSNTools(){
}

bool mdomSNTools::CheckVolumeFormDOMs(G4ThreeVector position){
	// check whether the given position is inside one of the modules
	// Translation of frame of reference and then check R in spherical coordinates. Module will be assume to be spherical with Rmin
	//It returns true when the particle is inside the module!
	//G4cout << position.getX()/m << " "<<position.getY()/m << " "<< position.getZ()/m <<G4endl;
	//G4ThreeVector ProblematicPoint = G4ThreeVector(-161.517*mm,85.6623*mm,22.7964*mm);
  	aNavigator->LocateGlobalPointAndSetup(position);
  	G4TouchableHistoryHandle aTouchable = aNavigator->CreateTouchableHistoryHandle();
	G4int HistoryDepth = aTouchable->GetHistoryDepth();
	//G4cout << HistoryDepth << G4endl;
	if (HistoryDepth > 0) {return true;}/* 
	for ( int j=0; j<(G4int)mdompos.size(); j++ ) {
		// First part: check sphere volume
		G4ThreeVector thismodulepos;
		thismodulepos = G4ThreeVector(0*m,0*m,mdompos[j]);
		G4ThreeVector translation = thismodulepos - position ;
		G4double distance = translation.getR();
		if (distance < Rmin) {
			return true;
		}
		//Second: check holder volume, simulated as cylinder
		// Harness has approximately this size... there is no global variables for it
		G4double radius = 200 * mm;
		G4double height = 23*mm;
		G4double Rho = pow(pow(translation.getX(),2)+pow(translation.getY(),2),1./2.);
		G4double z = abs(translation.getZ());
		//G4cout << Rho/m << " " << radius/m << " " << z/m << " " << height/m << G4endl;
		if ( (Rho < radius) && (z < height) ) {
			return true;
		}
	}*/
	return false;
}

void mdomSNTools::RandomPosition(G4ThreeVector& Position) {
  // Just to give a random coordinates for the primaries inside the Ice
  G4double Rmax = pow(3,1./2.)*gworldsize*m;
  G4double Rmax2 = gworldsize*m;
  
  G4double posz;
  G4double posx;
  G4double posy;
  G4ThreeVector tempPosition;
  G4double R3;
  G4double R2;
  R3 = 0*m;
  R2 = 0*m;
  G4bool boolparameter = true;
  while ( ( boolparameter==true ) || (R3 >= Rmax) || (R2 >= Rmax2)) {
	G4double posornegX = 1;
	if (G4UniformRand()<0.5) { posornegX = -1;}
		G4double posornegY = 1;
	if (G4UniformRand()<0.5) { posornegY = -1;}
		G4double posornegZ = 1;
	if (G4UniformRand()<0.5) { posornegZ = -1;}
	posz = posornegZ*(G4UniformRand()*gHeight);
	posx = posornegX*(G4UniformRand()*gRadius);
	posy = posornegY*(G4UniformRand()*gRadius);
	tempPosition = G4ThreeVector(posx,posy,posz);
	R3 = pow(pow(tempPosition[0],2)+pow(tempPosition[1],2)+pow(tempPosition[2],2),1./2.);
	R2 = pow(pow(tempPosition[0],2)+pow(tempPosition[1],2),1./2.);
	boolparameter = CheckVolumeFormDOMs(tempPosition);
  } 
	Position = G4ThreeVector(tempPosition);
}


G4double mdomSNTools::linealinterpolation(G4double realX,G4double lowerX, G4double upperX, G4double lowerY,G4double upperY) {
	G4double slope = (upperY-lowerY)/(upperX-lowerX);
	G4double result = (slope*(realX-lowerX)+lowerY);
	return result;
}


G4double mdomSNTools::EnergyDistribution(G4double Emean, G4double Emean2, G4double& alpha)
{// Energy distribution, returns the value of the energy from it using the inverse cumulative algorithm
  // tabulated function 
  // f is assumed positive, linear per segment, continuous
    G4int nPoints1 = 500;

    // Fe(E) corresponding to neutrino energies from a heavy neutron star, ls220
    
    // Data and formulas from:
    //  Irene Tamborra et al., High-resolution supernova neutrino spectra represented by a simple fit, PHYSICAL REVIEW D 86, 125031 (2012)
   //
   if (gfixmeanenergy == false) {
	alpha = GetAlpha(Emean, Emean2);
   } 
  std::vector<G4double> x1;
  std::vector<G4double> f1;
  MakeEnergyDistribution(Emean, alpha, nPoints1, x1 , f1);
	G4double choosenenergy = InverseCumul(x1, f1, nPoints1);
	return choosenenergy;
}

G4double mdomSNTools::GetAlpha(G4double Emean,G4double Emean2)
{
  // Get Alpha Parameter for the energy distribution
	G4double alpha = (2*pow(Emean,2)-Emean2)/(Emean2-pow(Emean,2));
	return alpha;
}

G4double mdomSNTools::InverseCumul(std::vector<G4double>  xvals, std::vector<G4double>  yvals, G4int nPoints)
{
  // InverseCumul gives a random value based in a given distribution
  // tabulated function
  // f is assumed positive, linear per segment, continuous 
    std::vector<G4double>  x_g;
    std::vector<G4double>  f_g;           //f(x)
    std::vector<G4double>  a_g;           //slopes
    std::vector<G4double>  Fc_g;          //cumulative of f

  GetSlopes(xvals, yvals, nPoints, x_g, f_g, a_g, Fc_g);
  G4double x_rndm =  InverseCumulAlgorithm(x_g, f_g, a_g, Fc_g, nPoints);
  return x_rndm;
}


G4int mdomSNTools::findtime(G4double time, std::vector<G4double> timearray)
{
  for (unsigned int j=0; j<timearray.size(); j++) {
    if (time <= timearray.at(j)) {
      return j;
    };
  };
 G4cout << "FATAL ERROR -> Not posible to find time of spectrum!!!" << G4endl;
 return 0;
}
   
G4double mdomSNTools::InverseCumulAlgorithm(std::vector<G4double>  x, std::vector<G4double>  f, std::vector<G4double>  a, std::vector<G4double>  Fc, G4int  nPoints)
{
  // InverseCumul calculation taking:
  //  std::vector<G4double>  x;
  //  std::vector<G4double>  f;           //f(x)
  //  std::vector<G4double>  a;           //slopes
  //  std::vector<G4double>  Fc;          //cumulative of f
	
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



//Helper functions 
void GetSlopes(std::vector<G4double>  xvals, std::vector<G4double>  yvals, G4int nPoints, std::vector<G4double>&  x, std::vector<G4double>&  f, std::vector<G4double>&  a, std::vector<G4double>&  Fc)
{
	//helper function which gets slopes and prepare the function to be use for the inverse cumulative algorithm
	// x -> new x to store the stuff, not necessary but convenient to do not overwrite into the same variable
	// f -> same but for y
	// a -> where to store the slopes
	// Fc -> cumulative function 
	
  // create a copy of the array. Really not necessary...
  x.resize(nPoints); f.resize(nPoints);
  for (G4int j=0; j<nPoints; j++) {
    x[j] = xvals.at(j); f[j] = yvals.at(j);
  };
  //compute slopes
  //
  a.resize(nPoints);
  for (G4int j=0; j<nPoints-1; j++) { 
    a[j] = (f[j+1] - f[j])/(x[j+1] - x[j]);
  };
  //compute cumulative function
  //
  Fc.resize(nPoints);  
  Fc[0] = 0.;
  for (G4int j=1; j<nPoints; j++) {
    Fc[j] = Fc[j-1] + 0.5*(f[j] + f[j-1])*(x[j] - x[j-1]);
  };     
}


G4double NumberOfTargets(G4int targetPerMolecule) {
  //Just to calculate number of targets in ice per cubic meter, assuming ice as H2O pure
  G4double Density = 921.6*kg/m3; //Density of ice at -50 celsius degrees
  G4double MolarMass = 18.01528e-3*kg; //kg per mol
  G4double Na = 6.022140857e23;
  G4double Nm = Density/MolarMass*Na;//molecules/m^3 of ice
  G4double Nt = Nm*targetPerMolecule;
  return Nt;
}


void MakeEnergyDistribution(G4double Emean, G4double alpha, G4int nPoints, std::vector<G4double>& x, std::vector<G4double>& f)
{
  G4double min = 0.; G4double max = 50.; 
  G4double delta = (max-min)/G4double(nPoints-1);
  x.resize(nPoints); f.resize(nPoints);

  for (G4int i=0; i<nPoints; i++) {
    x[i] = (min + i*delta)*MeV; //Energy
  }

  for (G4int j=0; j<nPoints; j++) {
  f[j] = pow(x[j],alpha)*exp(-(alpha+1.)*x[j]/Emean); // F(e), energy dist. function
  }
}


G4double mdomSNTools::WeighMe(G4double sigma, G4double NTargets) {
  // GIve the weigh of the interaction because of the cross section as:
  // Weigh = sigma * Nt * r, where r is the distance of the neutrino in the ice 
  // and Nt is the number of target in ice (electrons in this case)
  //
  // THIS IS THOUGH FOR SN NEUTRINOS COMING FROM Z AXIS WITH A WORLD AS A CYLINDER, if you changed some of that you will also have to change this function
  
  double weigh = sigma*NTargets*(2*gHeight);
  //G4cout << "Peso -> " << weigh << "  con Ntargets  " << NTargets/(1/m3) << "  y gHeight " << 2*gHeight/m << G4endl;
  return weigh;
}


