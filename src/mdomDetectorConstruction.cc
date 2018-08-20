#include "mdomDetectorConstruction.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Polycone.hh"
#include "G4Ellipsoid.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4NistManager.hh"
#include <sstream>
#include <iostream>
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4Transform3D.hh"
//since Geant4.10: include units manually
#include "G4SystemOfUnits.hh"
// #include "G4Cache.hh"

// extern G4String ghitsfilename;
// extern G4String	gabsfilename;
// extern G4String	greffilename;
// extern G4String	gscatfilename;
// extern G4String	gdetectortype;

extern G4int gPMT;		// PMT type
extern G4int gEnvironment;	// lab filling
extern G4bool gmdomharness;
extern G4bool gropes;
extern G4double	gmdomseparation;
extern G4int	gn_mDOMs;
extern G4Navigator* aNavigator;

extern G4double gRefCone_angle; // opening semi-angle of RefCone
extern G4int gGlass;
extern G4int gGel;
extern G4int gConeMat;
extern G4int gHolderColor;
extern G4int gDOM;
extern G4int	gDepthpos;

extern G4double gworldsize;
extern  G4double gRadius;
extern  G4double gHeight;

extern G4double gscintYield; 
extern G4double gscintTimeConst;
extern G4double gscintSpectrum;
extern G4double gTemperature;

extern std::vector<double> readColumnDouble (G4String fn, int col);
extern std::vector<G4String> explode (G4String s, char d);

double hc_eVnm = 1239.84193; // h*c in eV * nm
G4double hc_eVnm_w_units = 1239.84193*eV*nm; // h*c in eV * nm

// defining mDOM dimensions:
G4double GlasOutRad = 0.5*356*mm;	// outer radius of galss cylinder (pressure vessel); roughly 0.5*13"
G4double GlasThick = 13*mm;			// maximum glass thickness
G4double GlasInRad = GlasOutRad - GlasThick;
G4double CylHigh = 27.5*mm;			// height of cylindrical part of glass half-vessel
G4double GelThick = 2*mm;			// distance between inner glass surface and holding structure, filled with gel
G4double GelPMT = 2*mm;				// distance between inner glass surface and tip of PMTs
G4double RefConeDZ = 15*mm;			// half-height of reflector
G4double MPMTzoffset = 10*mm;		// z-offset of middle PMT circles w.r.t. center of glass sphere
G4double MPMTroffset = 2.6*mm;		// middle PMT circles are slightly further out due to MPMTzoffset

//	hard-coded limits
const double PHOTON_NRG_MAX = (hc_eVnm / 150)*eV;
const double PHOTON_NRG_MIN = (hc_eVnm / 730)*eV;

mdomDetectorConstruction::mdomDetectorConstruction()
:World_solid(0), World_logical(0), World_physical(0)
{}

mdomDetectorConstruction::~mdomDetectorConstruction()
{}


// --------------Mie scattering in ice, Ice absorption and Ice refraction-----------------
// Formula and data from :
//	M.G.Aartsen et al., Nuclear Instruments and Methods in Physics Research A 711 (2013) 73–89.
//	M. Ackermann, et al., Journal of Geophysical Research 111 (2006) D13203.

const G4int NUMENTRIES_ICE = 61;
const G4int NUMENTRIES_DEPTH = 110;

G4double ENERGY_spice[NUMENTRIES_ICE] = {
  1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
  1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
  1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
  1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
  1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
  2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
  2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
  2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
  2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
  2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
  3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
  3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
  3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
  4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
  5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV,
  hc_eVnm/150*eV
}; 

G4double be400inv_spice[NUMENTRIES_DEPTH]= {
  13.2*m, 14.0*m, 14.7*m, 17.0*m, 16.0*m, 14.4*m, 16.0*m,
  20.8*m, 26.7*m, 34.7*m, 39.7*m ,38.7*m, 27.8*m, 16.6*m, 
  13.7*m, 13.5*m ,15.7*m ,15.7*m ,14.7*m ,17.6*m ,21.6*m,
  24.0*m ,20.0*m ,17.8*m ,28.9*m ,36.9*m ,42.1*m ,46.5*m,
  45.4*m, 39.1*m ,30.6*m ,26.5*m ,19.3*m ,20.8*m ,20.1*m,
  20.3*m ,24.5*m ,33.5*m ,36.2*m ,35.4*m ,32.3*m ,40.2*m,
  44.7*m ,34.5*m ,30.6*m ,27.5*m ,19.7*m ,21.4*m ,28.8*m,
  38.3*m ,38.4*m ,44.2*m ,50.5*m ,46.6*m ,36.8*m ,26.7*m,
  20.3*m ,17.4*m ,16.1*m ,9.4*m ,10.6*m ,13.2*m ,10.9*m,
  6.8*m ,5.5*m ,5.0*m ,7.2*m ,9.8*m ,12.2*m ,21.1*m,
  54.3*m ,50.5*m ,33.5*m ,34.6*m ,48.4*m ,53.2*m ,46.3*m,
  32.9*m ,27.4*m ,30.5*m ,28.9*m ,35.1*m ,39.9*m ,48.0*m,
  53.3*m ,54.8*m ,57.9*m ,61.1*m ,76.8*m ,79.0*m ,75.6*m,
  75.3*m ,78.0*m ,59.4*m ,51.8*m ,32.9*m ,23.9*m ,28.6*m,
  32.5*m ,44.5*m ,56.9*m ,57.5*m ,54.3*m ,61.3*m ,68.8*m,
  77.6*m ,79.8*m ,89.4*m ,80.7*m ,56.7*m
};

G4double a400inv_spice[NUMENTRIES_DEPTH]= {
  45.1*m, 48.6*m, 53.2*m, 57.6*m, 57.6*m, 52.2*m, 60.1*m, 
  74.6*m, 96.6*m, 110.5*m, 135.6*m, 134.7*m, 98.2 *m, 64.7 *m, 
  48.5 *m, 44.3 *m, 54.4 *m, 56.7 *m, 52.1 *m, 60.7 *m, 72.7 *m, 
  78.9 *m, 68.7 *m, 66.6 *m, 100.0*m, 128.6*m, 148.2*m, 165.7*m, 
  156.0*m, 138.5*m, 113.9*m, 90.2 *m, 73.5 *m, 75.9 *m, 67.8 *m, 
  68.6 *m, 83.8 *m, 119.5*m, 121.6*m, 108.3*m, 113.4*m, 139.1*m, 
  148.1*m, 122.8*m, 113.8*m, 89.9 *m, 71.7 *m, 70.6 *m, 95.9 *m, 
  116.5*m, 143.6*m, 169.4*m, 178.0*m, 156.5*m, 135.3*m, 103.9*m, 
  75.2 *m, 66.2 *m, 53.7 *m, 33.6 *m, 36.2 *m, 44.0 *m, 40.4 *m, 
  24.9 *m, 20.1 *m, 17.9 *m, 28.4 *m, 34.4 *m, 41.6 *m, 84.4 *m, 
  173.1*m, 180.8*m, 116.7*m, 120.4*m, 164.4*m, 172.8*m, 149.2*m, 
  108.4*m, 91.1 *m, 98.9 *m, 94.0 *m, 113.1*m, 134.8*m, 154.1*m, 
  157.6*m, 180.5*m, 179.7*m, 185.2*m, 227.2*m, 220.8*m, 223.9*m, 
  256.6*m, 264.4*m, 193.7*m, 159.1*m, 118.7*m, 86.2 *m, 104.0*m, 
  119.7*m, 140.6*m, 203.5*m, 201.8*m, 178.2*m, 206.0*m, 205.2*m, 
  232.1*m, 259.4*m, 276.1*m, 244.3*m, 185.2*m
};


G4double Depth_spice[NUMENTRIES_DEPTH] = {
  1398.4*m, 1408.4*m, 1418.4*m, 1428.4*m, 1438.4*m, 1448.4*m, 1458.4*m, 
  1468.4*m, 1478.4*m, 1488.4*m, 1498.4*m, 1508.5*m, 1518.6*m, 1528.7*m, 
  1538.8*m, 1548.7*m, 1558.7*m, 1568.5*m, 1578.5*m, 1588.5*m, 1598.5*m, 
  1608.5*m, 1618.5*m, 1628.5*m, 1638.5*m, 1648.4*m, 1658.4*m, 1668.4*m, 
  1678.5*m, 1688.5*m, 1698.5*m, 1708.5*m, 1718.5*m, 1728.5*m, 1738.5*m, 
  1748.5*m, 1758.5*m, 1768.5*m, 1778.5*m, 1788.5*m, 1798.5*m, 1808.5*m, 
  1818.4*m, 1828.4*m, 1838.4*m, 1848.4*m, 1858.4*m, 1868.5*m, 1878.5*m, 
  1888.5*m, 1898.5*m, 1908.5*m, 1918.5*m, 1928.5*m, 1938.5*m, 1948.5*m, 
  1958.5*m, 1968.5*m, 1978.5*m, 1988.4*m, 1998.4*m, 2008.4*m, 2018.5*m, 
  2028.5*m, 2038.5*m, 2048.5*m, 2058.5*m, 2068.5*m, 2078.5*m, 2088.5*m, 
  2098.5*m, 2108.5*m, 2118.4*m, 2128.4*m, 2138.4*m, 2148.4*m, 2158.3*m, 
  2168.3*m, 2178.3*m, 2188.2*m, 2198.2*m, 2208.2*m, 2218.2*m, 2228.2*m, 
  2238.3*m, 2248.3*m, 2258.3*m, 2268.2*m, 2278.2*m, 2288.1*m, 2298.0*m, 
  2308.0*m, 2318.0*m, 2328.0*m, 2338.0*m, 2348.0*m, 2357.9*m, 2367.8*m, 
  2377.8*m, 2387.8*m, 2397.9*m, 2408.0*m, 2418.0*m, 2428.1*m, 2438.1*m, 
  2448.2*m, 2458.2*m, 2468.3*m, 2478.4*m, 2488.4*m
};


// gforward, gbackward, forward backward ratio
G4double MIE_spice_const[3]={0.972,0, 1};

G4double mdomDetectorConstruction::GetWaveLenght(float energy){
  // 	G4double h = 4.13566733e-15*eV*s;
  // 	G4double c = 3e17*nm/s;
  // 	G4double WaveLenght = h*c/energy;
  G4double WaveLenght = hc_eVnm_w_units / energy;
  return WaveLenght;
}

G4double mdomDetectorConstruction::Mie_Scattering(int u, int depth_pos){
  // depth_pos is the coordinate for the chosen depth in Depth_spice. For example to choose
  // depth=2278.2 m, we use depth_pos = 88
  G4double alpha = 0.90;
  G4double av_costheta = 0.9;
  G4double lambd = GetWaveLenght(ENERGY_spice[u]);
  G4double be_inv = 1./(1./(be400inv_spice[depth_pos])*pow((lambd/(400.*nm)), -alpha));
  G4double b_inv = be_inv*(1.-av_costheta);
  return b_inv;
}

G4double mdomDetectorConstruction::Spice_Temperature(G4double depth){
  G4double spice_temp = 221.5-0.00045319/m*depth+5.822e-6/m2*pow(depth, 2.);
  return spice_temp;
}

G4double mdomDetectorConstruction::Spice_Absorption(int u, int depth_pos){
  G4double kappa = 1.08;
  G4double paramA = 6954./m;
  G4double paramB = 6618*nm;
  G4double lambd = GetWaveLenght(ENERGY_spice[u]);
  G4double adust = 1./(a400inv_spice[depth_pos])*pow(lambd/(400.*nm), -kappa);
  G4double deltatau = Spice_Temperature(Depth_spice[depth_pos])-Spice_Temperature(1730.);
  G4double a_inv=1./(adust+paramA*exp(-paramB/lambd)*(1.+0.01*deltatau));
  return a_inv;
}

G4double mdomDetectorConstruction::Spice_Refraction(int u){
  // unknown depth. Parametrization by Thomas Kittler.
  G4double lambd = GetWaveLenght(ENERGY_spice[u])*1e-3;
  G4double n_phase = 1.55749 - 1.57988/nm * lambd + 3.99993/(nm*nm) * pow(lambd,2) - 4.68271/(nm*nm*nm) * pow(lambd,3) + 2.09354/(nm*nm*nm*nm) * pow(lambd,4);
  // 	G4double n_group = n_phase * ( 1. + 0.227106 - 0.954648/nm * lambd + 1.42568/(nm*nm) * pow(lambd,2) - 0.711832/(nm*nm*nm) * pow(lambd,3));
  return n_phase;	// using this now after discussion with Timo
}

//
// building the actual detector --------------------------------------------------------------------------------------------------------
//

G4VPhysicalVolume* mdomDetectorConstruction::Construct() {
  
  // some decisions concering the setup to be simulated (need to be inside function to work...):
  G4String PMT_types[] = {"12199s", "etel", "12199e"};
  G4String PMT_type = PMT_types[gPMT];
  
  // 	double gRefCone_angle = 45;		// XXX remove !
  G4double RefCone_angle;
  RefCone_angle = gRefCone_angle*deg;
  
  G4String RefCone_materials[] = {"V95", "V98", "Aluminium", "Total98"};
  G4String RefCone_material = RefCone_materials[gConeMat];	
  G4String Glass_types[] = {"VitroVex", "Chiba", "Kopp", "myVitroVex", "myChiba", "WOMQuartz", "fusedSilica"};
  G4String Glass_type = Glass_types[gGlass];	
  G4String Gel_types[] = {"Wacker", "Chiba", "IceCube", "Wacker_company"};
  G4String Gel_type = Gel_types[gGel];
  G4String OM_types[] = {"mDOM", "PDOM","halfmDOM","halfmDOM2","PMTwithSample"};
  G4String OM_type = OM_types[gDOM];	
  G4String Holder_colors[] = {"Black", "White"};
  G4String Holder_color = Holder_colors[gHolderColor];
  G4String World_types[] = {"air", "ice", "spice"};
  G4String World_type = World_types[gEnvironment];	
  
  //	declare variables, including visualisation attributes
  int i,k;
  double pi = M_PI;
  
  std::vector<G4String> n;
  G4Transform3D transformers;
  
  //	visualisation properties only, no effect on actual optics	
  G4VisAttributes* Glass_vis= new G4VisAttributes(G4Colour(0.7,0.7,0.8,0.2));
  G4VisAttributes* Gel_vis= new G4VisAttributes(G4Colour(0.45,0.5,0.35,0.2));
  G4VisAttributes* Alu_vis= new G4VisAttributes(G4Colour(0.8,0.8,0.9,1.0));
  G4VisAttributes* stahl_vis= new G4VisAttributes(G4Colour(0.75,0.75,0.85,1.0));
  G4VisAttributes* Absorber_vis= new G4VisAttributes(G4Colour(0.2,0.2,0.2,1.0));
  G4VisAttributes* PhotoCathode_vis= new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.3));
  G4VisAttributes* World_vis= new G4VisAttributes(G4Colour(0.0,0.0,0.0,0.05));
  // 	G4VisAttributes* World_vis= new G4VisAttributes(G4Colour(1.0,1.0,1.0,1.0));
  G4VisAttributes* Board_vis= new G4VisAttributes(G4Colour(0,1,0,1));
  
  //////////////////////////////
  // Optical surface definitions
  //////////////////////////////
  
  // --------------------------------- reflectivity of different cone materials --------------------------------------------------
  // "old" Aluminium with constant R
  // 	G4double RefCone_nrg[2] = {PHOTON_NRG_MIN, PHOTON_NRG_MAX} ; 
  // 	G4double RefCone_refl[2]      = { 0.92, 0.92 };
  
  // ideal polished un-oxidated Aluminium
  G4double AluPhotonEnergy[5]	= {hc_eVnm / 248*eV, hc_eVnm / 400*eV, hc_eVnm / 532*eV, hc_eVnm / 633*eV, hc_eVnm / 800*eV}; 
  G4double AluReflectivity[5]	= {0.9260, 0.9200, 0.9160, 0.9070, 0.8660};
  
  // Almeco V95
  // G4double V95PhotonEnergy[42] = {
  // hc_eVnm / 250*eV,
  // hc_eVnm / 255*eV,
  // hc_eVnm / 260*eV,
  // hc_eVnm / 265*eV,
  // hc_eVnm / 270*eV,
  // hc_eVnm / 275*eV,
  // hc_eVnm / 280*eV,
  // hc_eVnm / 285*eV,
  // hc_eVnm / 290*eV,
  // hc_eVnm / 295*eV,
  // hc_eVnm / 300*eV,
  // hc_eVnm / 305*eV,
  // hc_eVnm / 310*eV,
  // hc_eVnm / 315*eV,
  // hc_eVnm / 320*eV,
  // hc_eVnm / 325*eV,
  // hc_eVnm / 330*eV,
  // hc_eVnm / 335*eV,
  // hc_eVnm / 340*eV,
  // hc_eVnm / 345*eV,
  // hc_eVnm / 350*eV,
  // hc_eVnm / 355*eV,
  // hc_eVnm / 360*eV,
  // hc_eVnm / 365*eV,
  // hc_eVnm / 370*eV,
  // hc_eVnm / 375*eV,
  // hc_eVnm / 380*eV,
  // hc_eVnm / 385*eV,
  // hc_eVnm / 390*eV,
  // hc_eVnm / 395*eV,
  // hc_eVnm / 400*eV,
  // hc_eVnm / 430*eV,
  // hc_eVnm / 460*eV,
  // hc_eVnm / 490*eV,
  // hc_eVnm / 520*eV,
  // hc_eVnm / 550*eV,
  // hc_eVnm / 580*eV,
  // hc_eVnm / 610*eV,
  // hc_eVnm / 640*eV,
  // hc_eVnm / 670*eV,
  // hc_eVnm / 700*eV,
  // hc_eVnm / 730*eV
  // };
  
  // G4double V95Reflectivity[42] = {
  // 0.1710,
  // 0.1840,
  // 0.1969,
  // 0.2099,
  // 0.2228,
  // 0.2358,
  // 0.2253,
  // 0.2147,
  // 0.2042,
  // 0.1936,
  // 0.1836,
  // 0.1676,
  // 0.1219,
  // 0.0873,
  // 0.0679,
  // 0.0891,
  // 0.1892,
  // 0.3590,
  // 0.5413,
  // 0.6834,
  // 0.7765,
  // 0.8302,
  // 0.8313,
  // 0.8650,
  // 0.8988,
  // 0.9104,
  // 0.9221,
  // 0.9299,
  // 0.9377,
  // 0.9425,
  // 0.9472,
  // 0.9571,
  // 0.9582,
  // 0.9579,
  // 0.9564,
  // 0.9536,
  // 0.9464,
  // 0.9389,
  // 0.9293,
  // 0.9168,
  // 0.8993,
  // 0.8763
  // };
  
  // // Almeco V98
  // G4double V98Reflectivity[42] = {
  // 0.1408,
  // 0.1458,
  // 0.1508,
  // 0.1557,
  // 0.1607,
  // 0.1657,
  // 0.1670,
  // 0.1682,
  // 0.1695,
  // 0.1708,
  // 0.1722,
  // 0.2009,
  // 0.2700,
  // 0.2951,
  // 0.2807,
  // 0.2287,
  // 0.1256,
  // 0.0360,
  // 0.1619,
  // 0.4265,
  // 0.6006,
  // 0.6897,
  // 0.7687,
  // 0.8206,
  // 0.8724,
  // 0.8910,
  // 0.9097,
  // 0.9233,
  // 0.9369,
  // 0.9463,
  // 0.9558,
  // 0.9775,
  // 0.9807,
  // 0.9839,
  // 0.9850,
  // 0.9837,
  // 0.9843,
  // 0.9843,
  // 0.9840,
  // 0.9819,
  // 0.9810,
  // 0.9799
  // };
  
  // latest simulations provided by Almeco, 
  // coating of Vega material enhances refelctivity at desired wavelength, effect depends on n of surroundig material
  // V95 in air
  G4double V95AirPhotonEnergy[28] = {
    
    hc_eVnm / 150.0*eV,
    hc_eVnm / 200.0*eV,
    hc_eVnm / 225.0*eV,
    hc_eVnm / 250.0*eV,
    hc_eVnm / 259.0*eV,
    hc_eVnm / 270.5*eV,
    hc_eVnm / 281.9*eV,
    hc_eVnm / 295.0*eV,
    hc_eVnm / 303.2*eV,
    hc_eVnm / 317.1*eV,
    hc_eVnm / 323.6*eV,
    hc_eVnm / 329.4*eV,
    hc_eVnm / 332.7*eV,
    hc_eVnm / 335.1*eV,
    hc_eVnm / 336.7*eV,
    hc_eVnm / 340.8*eV,
    hc_eVnm / 344.9*eV,
    hc_eVnm / 348.2*eV,
    hc_eVnm / 350.7*eV,
    hc_eVnm / 358.0*eV,
    hc_eVnm / 372.8*eV,
    hc_eVnm / 399.8*eV,
    hc_eVnm / 449.7*eV,
    hc_eVnm / 499.6*eV,
    hc_eVnm / 598.6*eV,
    hc_eVnm / 648.5*eV,
    hc_eVnm / 698.4*eV,
    hc_eVnm / 730*eV
  };
  
  G4double V95AirReflectivity[28] = {
    
    0.000,
    0.071,
    0.117,
    0.265,
    0.277,
    0.277,
    0.269,
    0.263,
    0.245,
    0.180,
    0.166,
    0.182,
    0.215,
    0.261,
    0.315,
    0.488,
    0.695,
    0.813,
    0.865,
    0.908,
    0.946,
    0.964,
    0.968,
    0.966,
    0.944,
    0.928,
    0.898,
    0.849
  };
  
  // V95 in gel
  G4double V95GelPhotonEnergy[27] = {
    
    hc_eVnm / 150.0*eV,
    hc_eVnm / 200.0*eV,
    hc_eVnm / 225.0*eV,
    hc_eVnm / 250.0*eV,
    hc_eVnm / 259.0*eV,
    hc_eVnm / 269.6*eV,
    hc_eVnm / 281.9*eV,
    hc_eVnm / 290.9*eV,
    hc_eVnm / 309.7*eV,
    hc_eVnm / 314.6*eV,
    hc_eVnm / 322.0*eV,
    hc_eVnm / 328.6*eV,
    hc_eVnm / 333.5*eV,
    hc_eVnm / 336.7*eV,
    hc_eVnm / 341.7*eV,
    hc_eVnm / 349.0*eV,
    hc_eVnm / 353.1*eV,
    hc_eVnm / 363.3*eV,
    hc_eVnm / 380.9*eV,
    hc_eVnm / 398.9*eV,
    hc_eVnm / 449.7*eV,
    hc_eVnm / 499.6*eV,
    hc_eVnm / 549.5*eV,
    hc_eVnm / 599.0*eV,
    hc_eVnm / 648.9*eV,
    hc_eVnm / 698.9*eV,
    hc_eVnm / 730*eV
  };
  
  G4double V95GelReflectivity[27] = {
    
    0.000,
    0.071,
    0.117,
    0.206,
    0.215,
    0.205,
    0.184,
    0.180,
    0.132,
    0.104,
    0.095,
    0.116,
    0.184,
    0.301,
    0.500,
    0.801,
    0.855,
    0.906,
    0.940,
    0.952,
    0.958,
    0.954,
    0.946,
    0.932,
    0.916,
    0.893,
    0.857
  };
  
  // V98 in air
  G4double V98AirPhotonEnergy[26] = {
    
    hc_eVnm / 150.0*eV,
    hc_eVnm / 200.0*eV,
    hc_eVnm / 225.0*eV,
    hc_eVnm / 249.6*eV,
    hc_eVnm / 261.9*eV,
    hc_eVnm / 278.3*eV,
    hc_eVnm / 294.7*eV,
    hc_eVnm / 308.7*eV,
    hc_eVnm / 319.8*eV,
    hc_eVnm / 329.2*eV,
    hc_eVnm / 333.7*eV,
    hc_eVnm / 339.1*eV,
    hc_eVnm / 341.5*eV,
    hc_eVnm / 345.6*eV,
    hc_eVnm / 348.9*eV,
    hc_eVnm / 356.3*eV,
    hc_eVnm / 362.9*eV,
    hc_eVnm / 376.8*eV,
    hc_eVnm / 399.4*eV,
    hc_eVnm / 436.0*eV,
    hc_eVnm / 468.8*eV,
    hc_eVnm / 530.4*eV,
    hc_eVnm / 600.2*eV,
    hc_eVnm / 649.4*eV,
    hc_eVnm / 716.7*eV,
    hc_eVnm / 730*eV
  };
  
  G4double V98AirReflectivity[26] = {
    
    0.000,
    0.102,
    0.118,
    0.306,
    0.320,
    0.316,
    0.316,
    0.294,
    0.246,
    0.151,
    0.058,
    0.024,
    0.099,
    0.302,
    0.500,
    0.703,
    0.815,
    0.877,
    0.941,
    0.965,
    0.975,
    0.982,
    0.982,
    0.978,
    0.972,
    0.971
  };
  
  // V98 in gel
  G4double V98GelPhotonEnergy[28] = {
    
    hc_eVnm / 150.0*eV,
    hc_eVnm / 200.0*eV,
    hc_eVnm / 225.0*eV,
    hc_eVnm / 248.8*eV,
    hc_eVnm / 259.4*eV,
    hc_eVnm / 270.1*eV,
    hc_eVnm / 286.5*eV,
    hc_eVnm / 304.6*eV,
    hc_eVnm / 319.4*eV,
    hc_eVnm / 325.1*eV,
    hc_eVnm / 331.3*eV,
    hc_eVnm / 335.8*eV,
    hc_eVnm / 340.3*eV,
    hc_eVnm / 343.6*eV,
    hc_eVnm / 347.7*eV,
    hc_eVnm / 352.2*eV,
    hc_eVnm / 356.3*eV,
    hc_eVnm / 365.4*eV,
    hc_eVnm / 391.6*eV,
    hc_eVnm / 401.5*eV,
    hc_eVnm / 413.0*eV,
    hc_eVnm / 440.1*eV,
    hc_eVnm / 477.8*eV,
    hc_eVnm / 522.2*eV,
    hc_eVnm / 599.3*eV,
    hc_eVnm / 650.2*eV,
    hc_eVnm / 699.5*eV,
    hc_eVnm / 730*eV
  };
  
  G4double V98GelReflectivity[28] = {
    
    0.000,
    0.102,
    0.118,
    0.243,
    0.250,
    0.239,
    0.222,
    0.193,
    0.139,
    0.085,
    0.035,
    0.014,
    0.073,
    0.160,
    0.422,
    0.598,
    0.705,
    0.799,
    0.901,
    0.927,
    0.933,
    0.953,
    0.967,
    0.972,
    0.974,
    0.974,
    0.972,
    0.970
  };
  
  // Total98
  // hypopthetical material with R = 98% for all wavelengths, inpored by CTA ligth concentrators with additinal coating
  G4double T98PhotonEnergy[2] = {PHOTON_NRG_MIN, PHOTON_NRG_MAX};
  G4double T98Reflectivity[2] = {0.98, 0.98};
  
  // -------------------------------------- chosing actual cone material -----------------------------------------------------------
  
  G4OpticalSurface* RefCone_optical= new G4OpticalSurface("RefCone optical");
  RefCone_optical->SetModel(unified);
  RefCone_optical->SetType(dielectric_metal);
  RefCone_optical->SetFinish(polished);
  G4MaterialPropertiesTable *RefConeOpticalProperties = new G4MaterialPropertiesTable();
  if (RefCone_material == "Aluminium") {RefConeOpticalProperties -> AddProperty("REFLECTIVITY", AluPhotonEnergy, AluReflectivity, 5);}
  if (RefCone_material == "V95") {RefConeOpticalProperties -> AddProperty("REFLECTIVITY", V95GelPhotonEnergy, V95GelReflectivity, 27);}
  if (RefCone_material == "V98") {RefConeOpticalProperties -> AddProperty("REFLECTIVITY", V98GelPhotonEnergy, V98GelReflectivity, 28);}
  if (RefCone_material == "Total98") {RefConeOpticalProperties -> AddProperty("REFLECTIVITY", T98PhotonEnergy, T98Reflectivity, 2);}
  RefCone_optical->SetMaterialPropertiesTable(RefConeOpticalProperties);
  
  
  // ------------------------------------- fancy white holder -------------------------------------------------------------------------
  
  G4OpticalSurface* Holder_optical = new G4OpticalSurface("Holder optical");
  Holder_optical->SetModel(unified);
  Holder_optical->SetType(dielectric_dielectric);
  Holder_optical->SetFinish(groundfrontpainted);	// assures totally Lambertian (=diffuse) reflection
  G4double sigma_alpha = 0.1;	// roughness of surface (via microfacets), 0->rough, 1->nearly polished
  Holder_optical->SetSigmaAlpha(sigma_alpha);	
  G4MaterialPropertiesTable *HolderOpticalProperties = new G4MaterialPropertiesTable();
  G4double HolderPhotonEnergy[2] = {PHOTON_NRG_MIN, PHOTON_NRG_MAX};
  G4double HolderReflectivity[2] = {0.98, 0.98};	// inspired by OptoIndex catalogue
  HolderOpticalProperties -> AddProperty("REFLECTIVITY", HolderPhotonEnergy, HolderReflectivity, 2);
  Holder_optical->SetMaterialPropertiesTable(HolderOpticalProperties);
  
  // ------------------------------------- PDOM_Harness -------------------------------------------------------------------------
  G4OpticalSurface* PDOM_Harness_optical = new G4OpticalSurface("steel ring surface");
  PDOM_Harness_optical->SetModel(unified);
  PDOM_Harness_optical->SetType(dielectric_metal);
  PDOM_Harness_optical->SetFinish(ground);
  G4double PDOM_HarnessPhotonEnergy[2] = {PHOTON_NRG_MIN, PHOTON_NRG_MAX} ; 
  G4double PDOM_HarnessReflectivity[2]      = { 0.60, 0.60 };
  G4MaterialPropertiesTable *PDOM_Harness_MPT = new G4MaterialPropertiesTable();
  PDOM_Harness_MPT->AddProperty("REFLECTIVITY", PDOM_HarnessPhotonEnergy, PDOM_HarnessReflectivity,2);
  PDOM_Harness_optical->SetMaterialPropertiesTable(PDOM_Harness_MPT);
  
  
  //--------- Material definition ---------
  
  G4NistManager* MatDatBase = G4NistManager::Instance();
  double ambient_temperature = (-35+273.15)*kelvin;
  double ambient_temperature_air = (20+273.15)*kelvin;
  double ambient_pressure = 200*bar;
  double ambient_pressure_air = 1*bar;
  
  G4Material* Mat_LabAir = new G4Material("Air", 1.290*mg/cm3, MatDatBase->FindOrBuildMaterial("G4_AIR"), kStateGas,  ambient_temperature_air, ambient_pressure_air);
  
  G4Material* Mat_Ice = new G4Material("Ice", 0.917*g/cm3, MatDatBase->FindOrBuildMaterial("G4_WATER"), kStateSolid, ambient_temperature, ambient_pressure);
  
  G4Material* Mat_Spice = new G4Material("Spice", 0.917*g/cm3, MatDatBase->FindOrBuildMaterial("G4_WATER"), kStateSolid, ambient_temperature, ambient_pressure);
  
  G4Material* Mat_Vacuum = new G4Material("Vacuum", 0.3*1.290*mg/cm3, 1, kStateGas, ambient_temperature, 0.3*bar);
  Mat_Vacuum->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_AIR"), 100.0*perCent);
  
  G4Material* Mat_HighVacuum = new G4Material("High Vacuum", 0.3*1.290*mg/cm3, 1, kStateGas, ambient_temperature, pow(10,-9)*bar);
  Mat_HighVacuum->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_AIR"), 100.0*perCent);
  
  G4Material* Mat_Absorber = new G4Material("Absorber Black Paint", 1.0*g/cm3, 1, kStateSolid, ambient_temperature);
  Mat_Absorber->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_C"), 100.0*perCent);
  
  G4Material* Mat_Reflector = new G4Material("Reflective Material", 2.7*g/cm3, 1, kStateSolid, ambient_temperature);
  Mat_Reflector->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_Al"), 100.0*perCent);
  
  G4Material* Mat_Stahl = new G4Material("Stahl", 8*g/cm3, 1, kStateSolid, ambient_temperature, ambient_pressure);
  Mat_Stahl->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_Fe"), 100.0*perCent);
  
  G4Material* Mat_Vessel_Glass = new G4Material("Vessel Glass", 2.302*g/cm3, 2, kStateSolid, ambient_temperature, ambient_pressure);
  Mat_Vessel_Glass->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_Si"), 0.4674349);
  Mat_Vessel_Glass->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_O"), 0.5325651);
  
  G4Material* Mat_Tube_Glass = new G4Material("Tube Glass", 2.302*g/cm3, 2, kStateSolid, ambient_temperature);
  Mat_Tube_Glass->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_Si"), 0.4674349);
  Mat_Tube_Glass->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_O"), 0.5325651);
  
  G4Material* Mat_Gel = new G4Material("Optical Gel", 0.97*g/cm3, 4, kStateSolid, ambient_temperature);
  Mat_Gel->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_Si"), 0.3787);
  Mat_Gel->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_O"), 0.2158);
  Mat_Gel->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_C"), 0.3239);
  Mat_Gel->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_H"), 0.0816);
  
  G4Material* Mat_BiAlkali = new G4Material("Photocathode", 2.0*g/cm3, 3, kStateSolid, ambient_temperature);
  Mat_BiAlkali->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_Sb"), 0.5);
  Mat_BiAlkali->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_K"), 0.35);
  Mat_BiAlkali->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_Cs"), 0.15);
  
  // ------------ Generate & Add Material Properties Table ------------
  // ... insert possibility to read from textfiles here if needed ...
  
  //	Generate static (hard-wired) material properties
  
  //	Mie Scattering
  //G4int Depth_pos = 88; //depth = 2278.2 m
  G4double MIE_spice[NUMENTRIES_ICE]={};
  for (unsigned int u = 0; u < NUMENTRIES_ICE; u++) {
    MIE_spice[u] = Mie_Scattering(u , gDepthpos); 
  }
  //	Ice Absorption
  G4double ABS_spice[NUMENTRIES_ICE]={};
  for (unsigned int u = 0; u < NUMENTRIES_ICE; u++) {
    ABS_spice[u] = Spice_Absorption(u , gDepthpos); 
  }
  
  //	Ice refraction
  // 	G4cout << " Energy   " << "  n_group " << G4endl;  // printing for checking ...
  G4double refval_spice[NUMENTRIES_ICE]={};
  for (unsigned int u = 0; u < NUMENTRIES_ICE; u++) {
    // 	  G4cout << ENERGY_spice[u]/eV << "  " << Spice_Refraction(u) << G4endl;
    refval_spice[u] = Spice_Refraction(u);
  }
  
  G4MaterialPropertiesTable* proptable_ice = new G4MaterialPropertiesTable();
  proptable_ice->AddProperty("RINDEX", ENERGY_spice, refval_spice, NUMENTRIES_ICE)->SetSpline(true);
  Mat_Ice->SetMaterialPropertiesTable(proptable_ice);
  
  G4MaterialPropertiesTable* proptable_spice = new G4MaterialPropertiesTable();
  proptable_spice->AddProperty("RINDEX", ENERGY_spice, refval_spice, NUMENTRIES_ICE)->SetSpline(true);
  proptable_spice->AddProperty("ABSLENGTH", ENERGY_spice, ABS_spice, NUMENTRIES_ICE)->SetSpline(true);
  proptable_spice->AddProperty("MIEHG",ENERGY_spice,MIE_spice,NUMENTRIES_ICE)->SetSpline(true);
  proptable_spice->AddConstProperty("MIEHG_FORWARD",MIE_spice_const[0]);
  proptable_spice->AddConstProperty("MIEHG_BACKWARD",MIE_spice_const[1]);
  proptable_spice->AddConstProperty("MIEHG_FORWARD_RATIO",MIE_spice_const[2]);
  Mat_Spice->SetMaterialPropertiesTable(proptable_spice);
  
  G4double BiAlkaliPhotonEnergy[2] = {PHOTON_NRG_MIN, PHOTON_NRG_MAX};
  G4double BiAlkaliAbsLen[2] = {1*mm, 1*mm};
  G4double BiAlkaliRefIndex[2] = {4.,4.};
  G4MaterialPropertiesTable* proptable_BiAlkali = new G4MaterialPropertiesTable();
  proptable_BiAlkali->AddProperty("ABSLENGTH", BiAlkaliPhotonEnergy, BiAlkaliAbsLen, 2);
  proptable_BiAlkali->AddProperty("RINDEX", BiAlkaliPhotonEnergy, BiAlkaliRefIndex, 2);
  Mat_BiAlkali->SetMaterialPropertiesTable(proptable_BiAlkali);
  
  // ------------------------------ tube glass -------------------------------------------------------------
  
  G4double TubeGlassRIndPhotonEnergy[59] = {
    
    hc_eVnm / 150*eV,
    hc_eVnm / 160*eV,
    hc_eVnm / 170*eV,
    hc_eVnm / 180*eV,
    hc_eVnm / 190*eV,
    hc_eVnm / 200*eV,
    hc_eVnm / 210*eV,
    hc_eVnm / 220*eV,
    hc_eVnm / 230*eV,
    hc_eVnm / 240*eV,
    hc_eVnm / 250*eV,
    hc_eVnm / 260*eV,
    hc_eVnm / 270*eV,
    hc_eVnm / 280*eV,
    hc_eVnm / 290*eV,
    hc_eVnm / 300*eV,
    hc_eVnm / 310*eV,
    hc_eVnm / 320*eV,
    hc_eVnm / 330*eV,
    hc_eVnm / 340*eV,
    hc_eVnm / 350*eV,
    hc_eVnm / 360*eV,
    hc_eVnm / 370*eV,
    hc_eVnm / 380*eV,
    hc_eVnm / 390*eV,
    hc_eVnm / 400*eV,
    hc_eVnm / 410*eV,
    hc_eVnm / 420*eV,
    hc_eVnm / 430*eV,
    hc_eVnm / 440*eV,
    hc_eVnm / 450*eV,
    hc_eVnm / 460*eV,
    hc_eVnm / 470*eV,
    hc_eVnm / 480*eV,
    hc_eVnm / 490*eV,
    hc_eVnm / 500*eV,
    hc_eVnm / 510*eV,
    hc_eVnm / 520*eV,
    hc_eVnm / 530*eV,
    hc_eVnm / 540*eV,
    hc_eVnm / 550*eV,
    hc_eVnm / 560*eV,
    hc_eVnm / 570*eV,
    hc_eVnm / 580*eV,
    hc_eVnm / 590*eV,
    hc_eVnm / 600*eV,
    hc_eVnm / 610*eV,
    hc_eVnm / 620*eV,
    hc_eVnm / 630*eV,
    hc_eVnm / 640*eV,
    hc_eVnm / 650*eV,
    hc_eVnm / 660*eV,
    hc_eVnm / 670*eV,
    hc_eVnm / 680*eV,
    hc_eVnm / 690*eV,
    hc_eVnm / 700*eV,
    hc_eVnm / 710*eV,
    hc_eVnm / 720*eV,
    hc_eVnm / 730*eV
  };
  
  G4double TubeGlassRInd[59] = {
    
    1.7371,
    1.6932,
    1.6603,
    1.6349,
    1.6149,
    1.5987,
    1.5855,
    1.5744,
    1.5651,
    1.5572,
    1.5504,
    1.5445,
    1.5394,
    1.5348,
    1.5308,
    1.5273,
    1.5241,
    1.5213,
    1.5187,
    1.5164,
    1.5143,
    1.5123,
    1.5106,
    1.5090,
    1.5075,
    1.5061,
    1.5049,
    1.5037,
    1.5026,
    1.5016,
    1.5007,
    1.4998,
    1.4990,
    1.4983,
    1.4976,
    1.4969,
    1.4963,
    1.4957,
    1.4951,
    1.4946,
    1.4941,
    1.4937,
    1.4932,
    1.4928,
    1.4924,
    1.4920,
    1.4917,
    1.4913,
    1.4910,
    1.4907,
    1.4904,
    1.4901,
    1.4899,
    1.4896,
    1.4894,
    1.4891,
    1.4889,
    1.4887,
    1.4885
  };
  
  // --------------------- general energies for refractive indices -----------------------------------------		
  G4double GeneralRIndPhotonEnergy[49] = {
    
    hc_eVnm / 250*eV,
    hc_eVnm / 260*eV,
    hc_eVnm / 270*eV,
    hc_eVnm / 280*eV,
    hc_eVnm / 290*eV,
    hc_eVnm / 300*eV,
    hc_eVnm / 310*eV,
    hc_eVnm / 320*eV,
    hc_eVnm / 330*eV,
    hc_eVnm / 340*eV,
    hc_eVnm / 350*eV,
    hc_eVnm / 360*eV,
    hc_eVnm / 370*eV,
    hc_eVnm / 380*eV,
    hc_eVnm / 390*eV,
    hc_eVnm / 400*eV,
    hc_eVnm / 410*eV,
    hc_eVnm / 420*eV,
    hc_eVnm / 430*eV,
    hc_eVnm / 440*eV,
    hc_eVnm / 450*eV,
    hc_eVnm / 460*eV,
    hc_eVnm / 470*eV,
    hc_eVnm / 480*eV,
    hc_eVnm / 490*eV,
    hc_eVnm / 500*eV,
    hc_eVnm / 510*eV,
    hc_eVnm / 520*eV,
    hc_eVnm / 530*eV,
    hc_eVnm / 540*eV,
    hc_eVnm / 550*eV,
    hc_eVnm / 560*eV,
    hc_eVnm / 570*eV,
    hc_eVnm / 580*eV,
    hc_eVnm / 590*eV,
    hc_eVnm / 600*eV,
    hc_eVnm / 610*eV,
    hc_eVnm / 620*eV,
    hc_eVnm / 630*eV,
    hc_eVnm / 640*eV,
    hc_eVnm / 650*eV,
    hc_eVnm / 660*eV,
    hc_eVnm / 670*eV,
    hc_eVnm / 680*eV,
    hc_eVnm / 690*eV,
    hc_eVnm / 700*eV,
    hc_eVnm / 710*eV,
    hc_eVnm / 720*eV,
    hc_eVnm / 730*eV
  };
  
  // --------------------- VitroVex-----------------------------------------		
  
  G4double VitroVexGlassRInd[49] = {
    
    1.5229,
    1.5179,
    1.5135,
    1.5095,
    1.5060,
    1.5029,
    1.5001,
    1.4975,
    1.4952,
    1.4931,
    1.4912,
    1.4895,
    1.4879,
    1.4864,
    1.4851,
    1.4838,
    1.4827,
    1.4816,
    1.4806,
    1.4797,
    1.4789,
    1.4781,
    1.4773,
    1.4766,
    1.4759,
    1.4753,
    1.4747,
    1.4742,
    1.4737,
    1.4732,
    1.4727,
    1.4723,
    1.4719,
    1.4715,
    1.4711,
    1.4708,
    1.4704,
    1.4701,
    1.4698,
    1.4695,
    1.4693,
    1.4690,
    1.4687,
    1.4685,
    1.4683,
    1.4681,
    1.4678,
    1.4676,
    1.4674
  };
  
  G4double VitroVexGlassPhotonEnergy[35] = {
    
    hc_eVnm / 730.0*eV,
    hc_eVnm / 609.6*eV,
    hc_eVnm / 599.5*eV,
    hc_eVnm / 589.6*eV,
    hc_eVnm / 579.6*eV,
    hc_eVnm / 569.5*eV,
    hc_eVnm / 559.5*eV,
    hc_eVnm / 549.6*eV,
    hc_eVnm / 539.5*eV,
    hc_eVnm / 529.6*eV,
    hc_eVnm / 519.6*eV,
    hc_eVnm / 509.6*eV,
    hc_eVnm / 499.7*eV,
    hc_eVnm / 489.7*eV,
    hc_eVnm / 479.6*eV,
    hc_eVnm / 469.6*eV,
    hc_eVnm / 459.7*eV,
    hc_eVnm / 449.7*eV,
    hc_eVnm / 439.7*eV,
    hc_eVnm / 429.8*eV,
    hc_eVnm / 419.7*eV,
    hc_eVnm / 409.7*eV,
    hc_eVnm / 399.7*eV,
    hc_eVnm / 389.8*eV,
    hc_eVnm / 379.7*eV,
    hc_eVnm / 369.8*eV,
    hc_eVnm / 359.8*eV,
    hc_eVnm / 349.7*eV,
    hc_eVnm / 339.8*eV,
    hc_eVnm / 329.7*eV,
    hc_eVnm / 319.8*eV,
    hc_eVnm / 309.8*eV,
    hc_eVnm / 299.8*eV,
    hc_eVnm / 296.6*eV,
    hc_eVnm / 230.0*eV
  };
  
  G4double VitroVexGlassAbsLen[35] = {
    
    530.0*mm,
    530.0*mm,
    540.0*mm,
    580.0*mm,
    650.0*mm,
    750.0*mm,
    730.0*mm,
    650.0*mm,
    630.0*mm,
    600.0*mm,
    600.0*mm,
    580.0*mm,
    580.0*mm,
    500.0*mm,
    420.0*mm,
    400.0*mm,
    380.0*mm,
    350.0*mm,
    360.0*mm,
    360.0*mm,
    350.0*mm,
    450.0*mm,
    590.0*mm,
    610.0*mm,
    600.0*mm,
    300.0*mm,
    550.0*mm,
    400.0*mm,
    200.0*mm,
    100.0*mm,
    50.0*mm,
    30.0*mm,
    10.0*mm,
    0.0*mm,
    0.0*mm
  };
  
  
  
  
  // -------------------------- Chiba glass -----------------------------------------------------------
  
  G4double ChibaGlassRInd[49] = {
    
    1.5504,
    1.5445,
    1.5394,
    1.5348,
    1.5308,
    1.5273,
    1.5241,
    1.5213,
    1.5187,
    1.5164,
    1.5143,
    1.5123,
    1.5106,
    1.5090,
    1.5075,
    1.5061,
    1.5049,
    1.5037,
    1.5026,
    1.5016,
    1.5007,
    1.4998,
    1.4990,
    1.4983,
    1.4976,
    1.4969,
    1.4963,
    1.4957,
    1.4951,
    1.4946,
    1.4941,
    1.4937,
    1.4932,
    1.4928,
    1.4924,
    1.4920,
    1.4917,
    1.4913,
    1.4910,
    1.4907,
    1.4904,
    1.4901,
    1.4899,
    1.4896,
    1.4894,
    1.4891,
    1.4889,
    1.4887,
    1.4885
  };
  
  G4double ChibaGlassAbsPhotonEnergy[37] = {
    
    hc_eVnm / 250.0*eV,
    hc_eVnm / 254.0*eV,
    hc_eVnm / 258.0*eV,
    hc_eVnm / 262.0*eV,
    hc_eVnm / 266.0*eV,
    hc_eVnm / 270.0*eV,
    hc_eVnm / 274.0*eV,
    hc_eVnm / 278.0*eV,
    hc_eVnm / 282.0*eV,
    hc_eVnm / 286.0*eV,
    hc_eVnm / 290.0*eV,
    hc_eVnm / 294.0*eV,
    hc_eVnm / 298.0*eV,
    hc_eVnm / 302.0*eV,
    hc_eVnm / 306.0*eV,
    hc_eVnm / 310.0*eV,
    hc_eVnm / 314.0*eV,
    hc_eVnm / 318.0*eV,
    hc_eVnm / 322.0*eV,
    hc_eVnm / 326.0*eV,
    hc_eVnm / 330.0*eV,
    hc_eVnm / 335.0*eV,
    hc_eVnm / 340.0*eV,
    hc_eVnm / 347.0*eV,
    hc_eVnm / 351.0*eV,
    hc_eVnm / 355.0*eV,
    hc_eVnm / 368.0*eV,
    hc_eVnm / 372.0*eV,
    hc_eVnm / 376.0*eV,
    hc_eVnm / 380.0*eV,
    hc_eVnm / 384.0*eV,
    hc_eVnm / 388.0*eV,
    hc_eVnm / 392.0*eV,
    hc_eVnm / 397.0*eV,
    hc_eVnm / 401.0*eV,
    hc_eVnm / 405.0*eV,
    hc_eVnm / 730.0*eV
  };
  
  G4double ChibaGlassAbsLen[37] = {
    
    0*mm,
    0*mm,
    0*mm,
    0*mm,		
    
    // 		1.7*mm,		// probably unphysical values
    // 		1.7*mm,
    // 		1.7*mm,
    // 		1.7*mm,
    
    1.7*mm,
    1.8*mm,
    2.0*mm,
    2.2*mm,
    2.6*mm,
    3.2*mm,
    4.1*mm,
    5.4*mm,
    7.2*mm,
    9.4*mm,
    12.6*mm,
    16.8*mm,
    22.9*mm,
    29.2*mm,
    40.9*mm,
    54.5*mm,
    78.4*mm,
    99.5*mm,
    172.3*mm,
    250.7*mm,
    320.1*mm,
    480.7*mm,
    1434.1*mm,
    1434.1*mm,
    1434.1*mm,
    1434.1*mm,
    1434.1*mm,
    1434.1*mm,
    1434.1*mm,
    1434.1*mm,
    1434.1*mm,
    1434.1*mm,
    1434.1*mm
  };
  
  // ------------------------- IceCube glass --------------------------------------------------
  // values taken from DOMINANT simulation code from Chiba
  G4double IceCubeGlassPhotonEnergy[37] = { 
    
    hc_eVnm / 250.0*eV,
    hc_eVnm / 270.0*eV,
    hc_eVnm / 280.0*eV,
    hc_eVnm / 290.0*eV,
    hc_eVnm / 300.0*eV,
    hc_eVnm / 310.0*eV,
    hc_eVnm / 320.0*eV,
    hc_eVnm / 330.0*eV,
    hc_eVnm / 340.0*eV,
    hc_eVnm / 350.0*eV,
    hc_eVnm / 360.0*eV,
    hc_eVnm / 370.0*eV,
    hc_eVnm / 380.0*eV,
    hc_eVnm / 390.0*eV,
    hc_eVnm / 400.0*eV,
    hc_eVnm / 410.0*eV,
    hc_eVnm / 420.0*eV,
    hc_eVnm / 430.0*eV,
    hc_eVnm / 440.0*eV,
    hc_eVnm / 450.0*eV,
    hc_eVnm / 460.0*eV,
    hc_eVnm / 470.0*eV,
    hc_eVnm / 480.0*eV,
    hc_eVnm / 490.0*eV,
    hc_eVnm / 500.0*eV,
    hc_eVnm / 510.0*eV,
    hc_eVnm / 520.0*eV,
    hc_eVnm / 530.0*eV,
    hc_eVnm / 540.0*eV,
    hc_eVnm / 550.0*eV,
    hc_eVnm / 560.0*eV,
    hc_eVnm / 570.0*eV,
    hc_eVnm / 580.0*eV,
    hc_eVnm / 590.0*eV,
    hc_eVnm / 600.0*eV,
    hc_eVnm / 730.0*eV
  };
  
  G4double IceCubeGlassAbsLen[37] = {
    
    0.0*mm,
    0.0*mm,
    1.2*mm,
    1.4*mm,
    1.5*mm,
    3.8*mm,
    5.2*mm,
    7.8*mm,
    13.2*mm,
    24.2*mm,
    33.9*mm,
    52.8*mm,
    80.3*mm,
    107.5*mm,
    160.2*mm,
    181.2*mm,
    208.2*mm,
    236.5*mm,
    263.1*mm,
    297.1*mm,
    341.0*mm,
    401.5*mm,
    460.1*mm,
    513.1*mm,
    576.4*mm,
    648.6*mm,
    735.7*mm,
    821.4*mm,
    904.1*mm,
    995.0*mm,
    1047.6*mm,
    1106.1*mm,
    1157.8*mm,
    1229.6*mm,
    1293.7*mm,
    1293.7*mm
  };
  
  
  // --------------------- my VitroVex-----------------------------------------		
  
  G4double myVitroVexGlassPhotonEnergy[82] = {
    
    hc_eVnm / 250.0*eV,
    hc_eVnm / 305.0*eV,
    hc_eVnm / 310.0*eV,
    hc_eVnm / 315.0*eV,
    hc_eVnm / 320.0*eV,
    hc_eVnm / 325.0*eV,
    hc_eVnm / 330.0*eV,
    hc_eVnm / 335.0*eV,
    hc_eVnm / 340.0*eV,
    hc_eVnm / 345.0*eV,
    hc_eVnm / 350.0*eV,
    hc_eVnm / 355.0*eV,
    hc_eVnm / 360.0*eV,
    hc_eVnm / 365.0*eV,
    hc_eVnm / 370.0*eV,
    hc_eVnm / 375.0*eV,
    hc_eVnm / 380.0*eV,
    hc_eVnm / 385.0*eV,
    hc_eVnm / 390.0*eV,
    hc_eVnm / 395.0*eV,
    hc_eVnm / 400.0*eV,
    hc_eVnm / 405.0*eV,
    hc_eVnm / 410.0*eV,
    hc_eVnm / 415.0*eV,
    hc_eVnm / 420.0*eV,
    hc_eVnm / 425.0*eV,
    hc_eVnm / 430.0*eV,
    hc_eVnm / 435.0*eV,
    hc_eVnm / 440.0*eV,
    hc_eVnm / 445.0*eV,
    hc_eVnm / 450.0*eV,
    hc_eVnm / 455.0*eV,
    hc_eVnm / 460.0*eV,
    hc_eVnm / 465.0*eV,
    hc_eVnm / 470.0*eV,
    hc_eVnm / 475.0*eV,
    hc_eVnm / 480.0*eV,
    hc_eVnm / 485.0*eV,
    hc_eVnm / 490.0*eV,
    hc_eVnm / 495.0*eV,
    hc_eVnm / 500.0*eV,
    hc_eVnm / 505.0*eV,
    hc_eVnm / 510.0*eV,
    hc_eVnm / 515.0*eV,
    hc_eVnm / 520.0*eV,
    hc_eVnm / 525.0*eV,
    hc_eVnm / 530.0*eV,
    hc_eVnm / 535.0*eV,
    hc_eVnm / 540.0*eV,
    hc_eVnm / 545.0*eV,
    hc_eVnm / 550.0*eV,
    hc_eVnm / 555.0*eV,
    hc_eVnm / 560.0*eV,
    hc_eVnm / 565.0*eV,
    hc_eVnm / 570.0*eV,
    hc_eVnm / 575.0*eV,
    hc_eVnm / 580.0*eV,
    hc_eVnm / 585.0*eV,
    hc_eVnm / 590.0*eV,
    hc_eVnm / 595.0*eV,
    hc_eVnm / 600.0*eV,
    hc_eVnm / 605.0*eV,
    hc_eVnm / 610.0*eV,
    hc_eVnm / 615.0*eV,
    hc_eVnm / 620.0*eV,
    hc_eVnm / 625.0*eV,
    hc_eVnm / 630.0*eV,
    hc_eVnm / 635.0*eV,
    hc_eVnm / 640.0*eV,
    hc_eVnm / 645.0*eV,
    hc_eVnm / 650.0*eV,
    hc_eVnm / 655.0*eV,
    hc_eVnm / 660.0*eV,
    hc_eVnm / 665.0*eV,
    hc_eVnm / 670.0*eV,
    hc_eVnm / 675.0*eV,
    hc_eVnm / 680.0*eV,
    hc_eVnm / 685.0*eV,
    hc_eVnm / 690.0*eV,
    hc_eVnm / 695.0*eV,
    hc_eVnm / 700.0*eV,
    hc_eVnm / 730.0*eV
  };
  
  G4double myVitroVexGlassAbsLen[82] = {
    
    0.0*mm,
    0.7*mm,
    1.7*mm,
    3.8*mm,
    7.1*mm,
    12.0*mm,
    19.2*mm,
    29.6*mm,
    44.6*mm,
    66.5*mm,
    97.8*mm,
    145.7*mm,
    208.1*mm,
    297.9*mm,
    415.7*mm,
    494.2*mm,
    374.1*mm,
    327.9*mm,
    447.0*mm,
    639.4*mm,
    747.5*mm,
    742.1*mm,
    678.9*mm,
    566.4*mm,
    514.6*mm,
    478.4*mm,
    457.2*mm,
    429.9*mm,
    431.1*mm,
    416.9*mm,
    427.7*mm,
    459.2*mm,
    491.0*mm,
    525.2*mm,
    578.6*mm,
    612.5*mm,
    654.2*mm,
    685.2*mm,
    711.9*mm,
    776.5*mm,
    828.8*mm,
    889.2*mm,
    1023.0*mm,
    1113.9*mm,
    1278.5*mm,
    1410.6*mm,
    1468.5*mm,
    1615.6*mm,
    1525.9*mm,
    1565.4*mm,
    1524.8*mm,
    1471.3*mm,
    1326.7*mm,
    1455.7*mm,
    1260.1*mm,
    1136.3*mm,
    996.1*mm,
    934.6*mm,
    861.4*mm,
    803.8*mm,
    748.6*mm,
    720.0*mm,
    697.1*mm,
    684.9*mm,
    688.1*mm,
    660.8*mm,
    646.3*mm,
    674.0*mm,
    684.6*mm,
    703.3*mm,
    740.3*mm,
    774.0*mm,
    790.7*mm,
    864.4*mm,
    945.9*mm,
    1016.3*mm,
    1111.2*mm,
    1094.0*mm,
    1167.4*mm,
    1160.3*mm,
    1162.2*mm,
    1162.2*mm
  };
  
  
  // -------------------------- my Chiba glass -----------------------------------------------------------
  
  G4double myChibaGlassAbsPhotonEnergy[23] = {
    
    hc_eVnm / 250.0*eV,
    hc_eVnm / 295.0*eV,
    hc_eVnm / 300.0*eV,
    hc_eVnm / 305.0*eV,
    hc_eVnm / 310.0*eV,
    hc_eVnm / 315.0*eV,
    hc_eVnm / 320.0*eV,
    hc_eVnm / 325.0*eV,
    hc_eVnm / 330.0*eV,
    hc_eVnm / 335.0*eV,
    hc_eVnm / 340.0*eV,
    hc_eVnm / 345.0*eV,
    hc_eVnm / 350.0*eV,
    hc_eVnm / 355.0*eV,
    hc_eVnm / 360.0*eV,
    hc_eVnm / 365.0*eV,
    hc_eVnm / 370.0*eV,
    hc_eVnm / 375.0*eV,
    hc_eVnm / 380.0*eV,
    hc_eVnm / 385.0*eV,
    hc_eVnm / 390.0*eV,
    hc_eVnm / 395.0*eV,
    hc_eVnm / 730.0*eV
  };
  
  G4double myChibaGlassAbsLen[23] = {
    
    0*mm,
    0.3*mm,
    0.7*mm,
    1.6*mm,
    3.4*mm,
    6.1*mm,
    10.4*mm,
    16.2*mm,
    25.0*mm,
    37.2*mm,
    53.1*mm,
    77.9*mm,
    115.2*mm,
    173.1*mm,
    268.2*mm,
    409.1*mm,
    651.7*mm,
    946.0*mm,
    698.6*mm,
    625.1*mm,
    922.3*mm,
    1500*mm,
    1500*mm
  };
  
  
  
  // ------------------------- WOM quartz glass --------------------------------------------------
  // refractive index from Refractive index info
  // from I. H. Malitson. Interspecimen Comparison of the Refractive Index of Fused Silica, J. Opt. Soc. Am. 55, 1205-1208 (1965)
  
  G4double QuartzRIndPhotonEnergy[31] = {
    
    hc_eVnm / 150*eV,
    hc_eVnm / 160*eV,
    hc_eVnm / 170*eV,
    hc_eVnm / 190*eV,
    hc_eVnm / 210*eV,
    hc_eVnm / 230*eV,
    hc_eVnm / 250*eV,
    hc_eVnm / 270*eV,
    hc_eVnm / 290*eV,
    hc_eVnm / 310*eV,
    hc_eVnm / 330*eV,
    hc_eVnm / 350*eV,
    hc_eVnm / 370*eV,
    hc_eVnm / 390*eV,
    hc_eVnm / 410*eV,
    hc_eVnm / 430*eV,
    hc_eVnm / 450*eV,
    hc_eVnm / 470*eV,
    hc_eVnm / 490*eV,
    hc_eVnm / 510*eV,
    hc_eVnm / 530*eV,
    hc_eVnm / 550*eV,
    hc_eVnm / 570*eV,
    hc_eVnm / 590*eV,
    hc_eVnm / 610*eV,
    hc_eVnm / 630*eV,
    hc_eVnm / 650*eV,
    hc_eVnm / 670*eV,
    hc_eVnm / 690*eV,
    hc_eVnm / 710*eV,
    hc_eVnm / 730*eV
  };
  
  G4double QuartzGlassRInd[31] = {
    
    1.7029,
    1.64790241147,
    1.6114,
    1.5657,
    1.5384,
    1.5202,
    1.5074,
    1.4980,
    1.4908,
    1.4851,
    1.4806,
    1.4769,
    1.4738,
    1.4713,
    1.4691,
    1.4672,
    1.4656,
    1.4641,
    1.4629,
    1.4618,
    1.4608,
    1.4599,
    1.4591,
    1.4584,
    1.4577,
    1.4571,
    1.4565,
    1.4560,
    1.4555,
    1.4551,
    1.4546
  };
  
  G4double WOMGlassAbsPhotonEnergy[32] = {
    
    hc_eVnm / 150.0*eV,
    hc_eVnm / 182.0*eV,
    hc_eVnm / 190.0*eV,
    hc_eVnm / 200.0*eV,
    hc_eVnm / 206.0*eV,
    hc_eVnm / 215.7*eV,
    hc_eVnm / 219.6*eV,
    hc_eVnm / 223.6*eV,
    hc_eVnm / 228.4*eV,
    hc_eVnm / 233.8*eV,
    hc_eVnm / 242.4*eV,
    hc_eVnm / 253.0*eV,
    hc_eVnm / 261.4*eV,
    hc_eVnm / 267.6*eV,
    hc_eVnm / 277.9*eV,
    hc_eVnm / 285.6*eV,
    hc_eVnm / 295.1*eV,
    hc_eVnm / 312.5*eV,
    hc_eVnm / 335.2*eV,
    hc_eVnm / 368.5*eV,
    hc_eVnm / 396.5*eV,
    hc_eVnm / 438.3*eV,
    hc_eVnm / 458.9*eV,
    hc_eVnm / 487.5*eV,
    hc_eVnm / 522.4*eV,
    hc_eVnm / 559.4*eV,
    hc_eVnm / 584.8*eV,
    hc_eVnm / 603.3*eV,
    hc_eVnm / 622.3*eV,
    hc_eVnm / 639.2*eV,
    hc_eVnm / 648.8*eV,
    hc_eVnm / 730.0*eV
  };
  
  G4double WOMGlassAbsLen[32] = {
    0.0*mm,
    0.0*mm,
    3.1*mm,
    5.5*mm,
    9.3*mm,
    13.5*mm,
    15.0*mm,
    13.9*mm,
    10.3*mm,
    7.4*mm,
    5.7*mm,
    5.9*mm,
    7.3*mm,
    10.0*mm,
    22.7*mm,
    57.1*mm,
    127.0*mm,
    346.7*mm,
    1049.7*mm,
    1400.0*mm,
    522.5*mm,
    258.8*mm,
    522.5*mm,
    258.8*mm,
    258.8*mm,
    522.5*mm,
    1049.7*mm,
    346.7*mm,
    258.8*mm,
    522.5*mm,
    171.0*mm,
    346.7*mm
  };
  
  // ------------------------- Fused Silica --------------------------------------------------
  // data from Thorlabs https://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=3983
  // refractive index identical to WOM values above
  G4double FusedSilicaAbsPhotonEnergy[26] = {
    
    hc_eVnm / 150.0*eV,
    hc_eVnm / 180.0*eV,
    hc_eVnm / 190.0*eV,
    hc_eVnm / 195.0*eV,
    hc_eVnm / 200.0*eV,
    hc_eVnm / 201.0*eV,
    hc_eVnm / 202.0*eV,
    hc_eVnm / 203.0*eV,
    hc_eVnm / 204.0*eV,
    hc_eVnm / 209.0*eV,
    hc_eVnm / 210.0*eV,
    hc_eVnm / 211.0*eV,
    hc_eVnm / 212.0*eV,
    hc_eVnm / 213.0*eV,
    hc_eVnm / 214.0*eV,
    hc_eVnm / 215.0*eV,
    hc_eVnm / 219.0*eV,
    hc_eVnm / 229.0*eV,
    hc_eVnm / 230.0*eV,
    hc_eVnm / 231.0*eV,
    hc_eVnm / 238.0*eV,
    hc_eVnm / 240.0*eV,
    hc_eVnm / 242.0*eV,
    hc_eVnm / 245.0*eV,
    hc_eVnm / 246.0*eV,
    hc_eVnm / 730.0*eV
  };
  
  G4double FusedSilicaGlassAbsLen[26] = {
    
    0.*mm,
    0.*mm,
    20.*mm,
    65.0*mm,
    243.5*mm,
    267.1*mm,
    291.3*mm,
    315.3*mm,
    347.4*mm,
    518.5*mm,
    563.2*mm,
    596.2*mm,
    633.5*mm,
    667.2*mm,
    700.8*mm,
    703.6*mm,
    750.3*mm,
    761.6*mm,
    744.1*mm,
    800.5*mm,
    1125.4*mm,
    1286.0*mm,
    1483.8*mm,
    1899.6*mm,
    2000.0*mm,
    2000.0*mm
  };
  
  
  //----------------_Scintillation-----------
  // distribution of produced optical photons
  G4double TestScintVitro[9] = {
    0.000134,
    0.004432,
    0.053991,
    0.241971,
    0.398942,
    0.000134,
    0.004432,
    0.053991,
    0.241971
  };
  //----------------_Scintillation-----------
  G4double Scnt_SLOW[32] = { 
    0.000069,
    0.013156,
    0.022110,
    0.029136,
    0.035059,
    0.044152,
    0.051315,
    0.057515,
    0.061234,
    0.062612,
    0.063714,
    0.063714,
    0.062612,
    0.060959,
    0.058341,
    0.055035,
    0.049938,
    0.042223,
    0.034095,
    0.029274,
    0.024866,
    0.021559,
    0.017289,
    0.013156,
    0.010401,
    0.006268,
    0.004064,
    0.002687,
    0.001723,
    0.000896,
    0.000483,
    0.000345
  };
  /* G4double Scnt_FAST[9] = { 
   *	    0.000010, 
   *	    0.000020, 
   *	    0.000030, 
   *	    0.004000,                 
   *	    0.008000, 
   *	    0.005000, 
   *	    0.020000, 
   *	    0.001000,
   *	    0.000010 };*/
  G4double Scnt_PP[32] = { 
    hc_eVnm / 337.3*eV,
    hc_eVnm / (350.1+gscintSpectrum)*eV,
    hc_eVnm / (356.6+gscintSpectrum)*eV,
    hc_eVnm / (360.7+gscintSpectrum)*eV,
    hc_eVnm / (363.9+gscintSpectrum)*eV,
    hc_eVnm / (369.9+gscintSpectrum)*eV,
    hc_eVnm / (375.0+gscintSpectrum)*eV,
    hc_eVnm / (380.9+gscintSpectrum)*eV,
    hc_eVnm / (385.9+gscintSpectrum)*eV,
    hc_eVnm / (388.7+gscintSpectrum)*eV,
    hc_eVnm / (392.3+gscintSpectrum)*eV,
    hc_eVnm / (396.9+gscintSpectrum)*eV,
    hc_eVnm / (401.4+gscintSpectrum)*eV,
    hc_eVnm / (405.0+gscintSpectrum)*eV,
    hc_eVnm / (409.5+gscintSpectrum)*eV,
    hc_eVnm / (414.0+gscintSpectrum)*eV,
    hc_eVnm / (419.9+gscintSpectrum)*eV,
    hc_eVnm / (429.3+gscintSpectrum)*eV,
    hc_eVnm / (439.7+gscintSpectrum)*eV,
    hc_eVnm / (446.9+gscintSpectrum)*eV,
    hc_eVnm / (454.1+gscintSpectrum)*eV,
    hc_eVnm / (460.0+gscintSpectrum)*eV,
    hc_eVnm / (468.5+gscintSpectrum)*eV,
    hc_eVnm / (478.0+gscintSpectrum)*eV,
    hc_eVnm / (486.6+gscintSpectrum)*eV,
    hc_eVnm / (505.6+gscintSpectrum)*eV,
    hc_eVnm / (519.7+gscintSpectrum)*eV,
    hc_eVnm / (534.2+gscintSpectrum)*eV,
    hc_eVnm / (549.6+gscintSpectrum)*eV,
    hc_eVnm / (568.1+gscintSpectrum)*eV,
    hc_eVnm / (585.8+gscintSpectrum)*eV,
    hc_eVnm / (599.4+gscintSpectrum)*eV
  };
  
  G4double Temperature[5] = {
    -15,
    -25,
    -35,
    -45,
    -50
  };
  
  
  
  G4double FirstCompomentAmplitude[5] = {
    7.44264474e-03,
    7.09497579e-03,
    6.88619615e-03,
    6.71650386e-03,
    6.26193414e-03 
  };
  
    G4double SecondCompomentAmplitude[5] = {
    3.72712597e-03,
    3.44366607e-03,
    3.22187477e-03,
    2.79396418e-03,
    3.11992762e-03
  };
  
    G4double ThirdCompomentAmplitude[5] = {
    9.37985502e-04,
    8.57409419e-04,
    7.53362142e-04,
    6.70626881e-04,
    7.18236379e-04 
  };
  
    G4double FirstTime[5] = {
    200.235250835*ns,
    207.368858907*ns,
    221.330366343*ns,
    246.15527870*ns,
    204.9615*ns 
  };
  
    G4double SecondTime[5] = {
    1543.61044827*ns,
    1681.26113657*ns,
    1806.1209032*ns,
    1979.83710164*ns,
    1647.6259*ns
  };
  
    G4double ThirdTime[5] = {
    11851.7691222*ns,
    12832.1823583*ns,
    14405.5681616*ns,
    15307.2132678*ns,
    14335.53*ns
  };
  
  
  
  
  // ------------------------- choosing glass for simulation ------------------------------------
  // PMT glass
  G4MaterialPropertiesTable* proptable_TubeGlass = new G4MaterialPropertiesTable();
  proptable_TubeGlass->AddProperty("RINDEX", TubeGlassRIndPhotonEnergy, TubeGlassRInd, 56);
  // absorption by PMT glass taken care of by QE!
  
  // VitroVex glass
  G4MaterialPropertiesTable* proptable_VitrovexGlass = new G4MaterialPropertiesTable();
  
  //----------------_Scintillation-----------
  G4double scintYield=gscintYield/MeV;
  G4double sctintTimeConst = gscintTimeConst*ns;
  G4double sTemperature = gTemperature;
  G4int tempIndex;
  if (sTemperature>-15 || sTemperature < -50){
    G4cout << "Selected Temperature out of range. Data goes from -50C to -15C. Simulation will use default temperature -35C." << G4endl;
    sTemperature = -35;
  }
  
  for (int i = 0; i < sizeof(Temperature); ++i)
  {
    
    if (Temperature[i] == sTemperature){
      tempIndex = i;
      
    }
    
  };

  proptable_VitrovexGlass->AddConstProperty("SCINTILLATIONYIELD",scintYield);
  proptable_VitrovexGlass->AddConstProperty("FIRSTAMPLITUDE",FirstCompomentAmplitude[tempIndex]);
  proptable_VitrovexGlass->AddConstProperty("SECONDAMPLITUDE",SecondCompomentAmplitude[tempIndex]);
  proptable_VitrovexGlass->AddConstProperty("THIRDAMPLITUDE",ThirdCompomentAmplitude[tempIndex]);
  
  proptable_VitrovexGlass->AddProperty("FIRSTCOMPONENT",Scnt_PP,Scnt_SLOW,32);
  proptable_VitrovexGlass->AddProperty("SECONDCOMPONENT",Scnt_PP,Scnt_SLOW,32);
  proptable_VitrovexGlass->AddProperty("THIRDCOMPONENT",Scnt_PP,Scnt_SLOW,32);
  
  proptable_VitrovexGlass->AddConstProperty("FIRSTTIME",FirstTime[tempIndex]);
  proptable_VitrovexGlass->AddConstProperty("SECONDTIME",SecondTime[tempIndex]);
  proptable_VitrovexGlass->AddConstProperty("THIRDTIME",ThirdTime[tempIndex]);
  
  proptable_VitrovexGlass->AddConstProperty("RESOLUTIONSCALE", 1.0);


  //----------------_Scintillation-----------
  
  proptable_VitrovexGlass->AddProperty("RINDEX", GeneralRIndPhotonEnergy, VitroVexGlassRInd, 49);
  proptable_VitrovexGlass->AddProperty("ABSLENGTH", VitroVexGlassPhotonEnergy, VitroVexGlassAbsLen,35);
  
  // Chiba glass
  G4MaterialPropertiesTable* proptable_ChibaGlass = new G4MaterialPropertiesTable();
  proptable_ChibaGlass->AddProperty("RINDEX", GeneralRIndPhotonEnergy, ChibaGlassRInd, 49);
  proptable_ChibaGlass->AddProperty("ABSLENGTH", ChibaGlassAbsPhotonEnergy, ChibaGlassAbsLen, 37);
  
  // IceCube glass (aka Benthos, aka Kopp 9500) using Chiba RINDEX !!
  G4MaterialPropertiesTable* proptable_KoppGlass = new G4MaterialPropertiesTable();
  proptable_KoppGlass->AddProperty("RINDEX", GeneralRIndPhotonEnergy, ChibaGlassRInd, 49);
  proptable_KoppGlass->AddProperty("ABSLENGTH", IceCubeGlassPhotonEnergy, IceCubeGlassAbsLen,37);
  
  
  // myVitroVex glass (values deduced from my measurement)
  G4MaterialPropertiesTable* proptable_myVitrovexGlass = new G4MaterialPropertiesTable();
  proptable_myVitrovexGlass->AddProperty("RINDEX", GeneralRIndPhotonEnergy, VitroVexGlassRInd, 49);
  proptable_myVitrovexGlass->AddProperty("ABSLENGTH", myVitroVexGlassPhotonEnergy, myVitroVexGlassAbsLen, 82);
  
  // myChiba glass (values deduced from my measurement)
  G4MaterialPropertiesTable* proptable_myChibaGlass = new G4MaterialPropertiesTable();
  proptable_myChibaGlass->AddProperty("RINDEX", GeneralRIndPhotonEnergy, ChibaGlassRInd, 49);
  proptable_myChibaGlass->AddProperty("ABSLENGTH", myChibaGlassAbsPhotonEnergy, myChibaGlassAbsLen, 23);
  
  // WOM glass (aka Quartz)
  G4MaterialPropertiesTable* proptable_QuartzGlass = new G4MaterialPropertiesTable();
  proptable_QuartzGlass->AddProperty("RINDEX", QuartzRIndPhotonEnergy, QuartzGlassRInd, 31);
  proptable_QuartzGlass->AddProperty("ABSLENGTH", WOMGlassAbsPhotonEnergy, WOMGlassAbsLen,36); 
  
  // Fused Silica from Thorlabs 
  G4MaterialPropertiesTable* proptable_FusedSilica = new G4MaterialPropertiesTable();
  proptable_FusedSilica->AddProperty("RINDEX", QuartzRIndPhotonEnergy, QuartzGlassRInd, 31);
  proptable_FusedSilica->AddProperty("ABSLENGTH", FusedSilicaAbsPhotonEnergy, FusedSilicaGlassAbsLen, 26); 
  
  // choose here: YYY
  Mat_Tube_Glass->SetMaterialPropertiesTable(proptable_TubeGlass); // no choice here so far...
  
  if (Glass_type == "VitroVex"){
    Mat_Vessel_Glass->SetMaterialPropertiesTable(proptable_VitrovexGlass);
  }
  if (Glass_type == "Chiba"){
    Mat_Vessel_Glass->SetMaterialPropertiesTable(proptable_ChibaGlass);
  }
  if (Glass_type == "Kopp"){
    Mat_Vessel_Glass->SetMaterialPropertiesTable(proptable_KoppGlass);
  }
  if (Glass_type == "myVitroVex"){
    Mat_Vessel_Glass->SetMaterialPropertiesTable(proptable_myVitrovexGlass);
  }
  if (Glass_type == "myChiba"){
    Mat_Vessel_Glass->SetMaterialPropertiesTable(proptable_myChibaGlass);
  }
  if (Glass_type == "WOMQuartz"){
    Mat_Vessel_Glass->SetMaterialPropertiesTable(proptable_QuartzGlass);
  }
  if (Glass_type == "fusedSilica"){
    Mat_Vessel_Glass->SetMaterialPropertiesTable(proptable_FusedSilica);
  }
  
  
  // ---------------------- Wacker gel -----------------------------------------------------
  // from old KM3NeT measuremts 
  // used by Björn and Claudio, Christophe uses same values, scaled down by 2
  
  G4double WackerGelRInd[49];
  for (i = 0; i < 49; i++) { WackerGelRInd[i] = 1.404; }
  
  G4double WackerGelPhotonEnergy[34] = {
    
    hc_eVnm / 250.0*eV,
    hc_eVnm / 300.0*eV,
    hc_eVnm / 310.0*eV,
    hc_eVnm / 320.0*eV,
    hc_eVnm / 330.0*eV,
    hc_eVnm / 340.0*eV,
    hc_eVnm / 350.0*eV,
    hc_eVnm / 360.0*eV,
    hc_eVnm / 370.0*eV,
    hc_eVnm / 380.0*eV,
    hc_eVnm / 390.0*eV,
    hc_eVnm / 400.0*eV,
    hc_eVnm / 410.0*eV,
    hc_eVnm / 420.0*eV,
    hc_eVnm / 430.0*eV,
    hc_eVnm / 440.0*eV,
    hc_eVnm / 450.0*eV,
    hc_eVnm / 460.0*eV,
    hc_eVnm / 470.0*eV,
    hc_eVnm / 480.0*eV,
    hc_eVnm / 490.0*eV,
    hc_eVnm / 500.0*eV,
    hc_eVnm / 510.0*eV,
    hc_eVnm / 520.0*eV,
    hc_eVnm / 530.0*eV,
    hc_eVnm / 540.0*eV,
    hc_eVnm / 550.0*eV,
    hc_eVnm / 560.0*eV,
    hc_eVnm / 570.0*eV,
    hc_eVnm / 580.0*eV,
    hc_eVnm / 590.0*eV,
    hc_eVnm / 600.0*eV,
    hc_eVnm / 610.0*eV,
    hc_eVnm / 730.0*eV
  };
  
  G4double WackerGelAbsLen[34] = {
    
    0.0*mm,
    0.0*mm,
    80.0*mm,
    156.0*mm,
    230.8*mm,
    304.9*mm,
    371.4*mm,
    418.8*mm,
    457.1*mm,
    489.6*mm,
    532.9*mm,
    566.4*mm,
    593.8*mm,
    625.3*mm,
    644.8*mm,
    669.1*mm,
    680.5*mm,
    723.1*mm,
    745.5*mm,
    764.8*mm,
    781.8*mm,
    810.8*mm,
    844.9*mm,
    858.8*mm,
    869.5*mm,
    901.0*mm,
    890.9*mm,
    943.6*mm,
    964.2*mm,
    969.0*mm,
    998.9*mm,
    999.4*mm,
    1008.1*mm,
    1008.1*mm
  };
  
  // ---------------------- Chiba gel ----------------------------------------------------------
  G4double ChibaGelRInd[49] = {
    
    1.4659,
    1.4586,
    1.4524,
    1.4469,
    1.4421,
    1.4379,
    1.4341,
    1.4308,
    1.4277,
    1.4250,
    1.4226,
    1.4203,
    1.4183,
    1.4164,
    1.4147,
    1.4132,
    1.4117,
    1.4104,
    1.4092,
    1.4081,
    1.4070,
    1.4060,
    1.4051,
    1.4042,
    1.4034,
    1.4027,
    1.4020,
    1.4013,
    1.4007,
    1.4001,
    1.3995,
    1.3990,
    1.3985,
    1.3980,
    1.3976,
    1.3972,
    1.3968,
    1.3964,
    1.3960,
    1.3957,
    1.3954,
    1.3951,
    1.3948,
    1.3945,
    1.3942,
    1.3939,
    1.3937,
    1.3935,
    1.3932
  };
  
  G4double ChibaGelPhotonEnergy[17] = {
    
    hc_eVnm / 250.0*eV,
    hc_eVnm / 260.0*eV,
    hc_eVnm / 274.0*eV,
    hc_eVnm / 284.0*eV,
    hc_eVnm / 294.0*eV,
    hc_eVnm / 304.0*eV,
    hc_eVnm / 314.0*eV,
    hc_eVnm / 324.0*eV,
    hc_eVnm / 334.0*eV,
    hc_eVnm / 344.0*eV,
    hc_eVnm / 354.0*eV,
    hc_eVnm / 364.0*eV,
    hc_eVnm / 374.0*eV,
    hc_eVnm / 384.0*eV,
    hc_eVnm / 394.0*eV,
    hc_eVnm / 400.0*eV,
    hc_eVnm / 730.0*eV
  };
  
  G4double ChibaGelAbsLen[17] = {
    
    7.8*mm,		//really??
    8.2*mm,
    18.6*mm,
    65.3*mm,
    281.4*mm,
    583.2*mm,
    854.8*mm,
    1192.8*mm,
    1506.8*mm,
    1927.2*mm,
    2376.6*mm,
    2811.9*mm,
    3293.1*mm,
    4262.8*mm,
    4943.5*mm,
    5000.0*mm,
    5000.0*mm
  };
  
  // ----------------------- IceCube gel ---------------------------------------------------
  G4double IceCubeGelRInd[49];
  for (i = 0; i < 49; i++) { IceCubeGelRInd[i] = 1.404; }
  
  G4double IceCubeGelPhotonEnergy[36] = {
    
    hc_eVnm / 250.0*eV,
    hc_eVnm / 270.0*eV,
    hc_eVnm / 280.0*eV,
    hc_eVnm / 290.0*eV,
    hc_eVnm / 300.0*eV,
    hc_eVnm / 310.0*eV,
    hc_eVnm / 320.0*eV,
    hc_eVnm / 330.0*eV,
    hc_eVnm / 340.0*eV,
    hc_eVnm / 350.0*eV,
    hc_eVnm / 360.0*eV,
    hc_eVnm / 370.0*eV,
    hc_eVnm / 380.0*eV,
    hc_eVnm / 390.0*eV,
    hc_eVnm / 400.0*eV,
    hc_eVnm / 410.0*eV,
    hc_eVnm / 420.0*eV,
    hc_eVnm / 430.0*eV,
    hc_eVnm / 440.0*eV,
    hc_eVnm / 450.0*eV,
    hc_eVnm / 460.0*eV,
    hc_eVnm / 470.0*eV,
    hc_eVnm / 480.0*eV,
    hc_eVnm / 490.0*eV,
    hc_eVnm / 500.0*eV,
    hc_eVnm / 510.0*eV,
    hc_eVnm / 520.0*eV,
    hc_eVnm / 530.0*eV,
    hc_eVnm / 540.0*eV,
    hc_eVnm / 550.0*eV,
    hc_eVnm / 560.0*eV,
    hc_eVnm / 570.0*eV,
    hc_eVnm / 580.0*eV,
    hc_eVnm / 590.0*eV,
    hc_eVnm / 600.0*eV,
    hc_eVnm / 730.0*eV
  };
  
  G4double IceCubeGelAbsLen[36] = {
    
    0.0*mm,
    0.0*mm,
    1.8*mm,
    12.3*mm,
    23.4*mm,
    39.0*mm,
    58.4*mm,
    80.1*mm,
    105.4*mm,
    138.8*mm,
    168.0*mm,
    203.3*mm,
    239.5*mm,
    254.0*mm,
    283.2*mm,
    308.5*mm,
    339.8*mm,
    354.7*mm,
    378.1*mm,
    406.5*mm,
    431.7*mm,
    462.3*mm,
    485.2*mm,
    510.4*mm,
    538.5*mm,
    579.8*mm,
    608.5*mm,
    631.9*mm,
    648.6*mm,
    719.6*mm,
    741.3*mm,
    776.2*mm,
    808.0*mm,
    849.7*mm,
    849.7*mm,
    849.7*mm
  };
  
  // ------------------------ Wacker gel with company data at UV ---------------------
  // data below 300 nm taken from transmission measurement of 5mm gel, provided by company 
  
  G4double WackerCompanyGelRIndPhotonEnergy[59] = {
    
    hc_eVnm / 150*eV,
    hc_eVnm / 160*eV,
    hc_eVnm / 170*eV,
    hc_eVnm / 180*eV,
    hc_eVnm / 190*eV,
    hc_eVnm / 200*eV,
    hc_eVnm / 210*eV,
    hc_eVnm / 220*eV,
    hc_eVnm / 230*eV,
    hc_eVnm / 240*eV,
    hc_eVnm / 250*eV,
    hc_eVnm / 260*eV,
    hc_eVnm / 270*eV,
    hc_eVnm / 280*eV,
    hc_eVnm / 290*eV,
    hc_eVnm / 300*eV,
    hc_eVnm / 310*eV,
    hc_eVnm / 320*eV,
    hc_eVnm / 330*eV,
    hc_eVnm / 340*eV,
    hc_eVnm / 350*eV,
    hc_eVnm / 360*eV,
    hc_eVnm / 370*eV,
    hc_eVnm / 380*eV,
    hc_eVnm / 390*eV,
    hc_eVnm / 400*eV,
    hc_eVnm / 410*eV,
    hc_eVnm / 420*eV,
    hc_eVnm / 430*eV,
    hc_eVnm / 440*eV,
    hc_eVnm / 450*eV,
    hc_eVnm / 460*eV,
    hc_eVnm / 470*eV,
    hc_eVnm / 480*eV,
    hc_eVnm / 490*eV,
    hc_eVnm / 500*eV,
    hc_eVnm / 510*eV,
    hc_eVnm / 520*eV,
    hc_eVnm / 530*eV,
    hc_eVnm / 540*eV,
    hc_eVnm / 550*eV,
    hc_eVnm / 560*eV,
    hc_eVnm / 570*eV,
    hc_eVnm / 580*eV,
    hc_eVnm / 590*eV,
    hc_eVnm / 600*eV,
    hc_eVnm / 610*eV,
    hc_eVnm / 620*eV,
    hc_eVnm / 630*eV,
    hc_eVnm / 640*eV,
    hc_eVnm / 650*eV,
    hc_eVnm / 660*eV,
    hc_eVnm / 670*eV,
    hc_eVnm / 680*eV,
    hc_eVnm / 690*eV,
    hc_eVnm / 700*eV,
    hc_eVnm / 710*eV,
    hc_eVnm / 720*eV,
    hc_eVnm / 730*eV
  };
  
  G4double WackerCompanyGelRInd[59] = {
    
    1.7412,
    1.6667,
    1.6155,
    1.5782,
    1.5500,
    1.5280,
    1.5104,
    1.4961,
    1.4843,
    1.4743,
    1.4659,
    1.4586,
    1.4524,
    1.4469,
    1.4421,
    1.4379,
    1.4341,
    1.4308,
    1.4277,
    1.4250,
    1.4226,
    1.4203,
    1.4183,
    1.4164,
    1.4147,
    1.4132,
    1.4117,
    1.4104,
    1.4092,
    1.4081,
    1.4070,
    1.4060,
    1.4051,
    1.4042,
    1.4034,
    1.4027,
    1.4020,
    1.4013,
    1.4007,
    1.4001,
    1.3995,
    1.3990,
    1.3985,
    1.3980,
    1.3976,
    1.3972,
    1.3968,
    1.3964,
    1.3960,
    1.3957,
    1.3954,
    1.3951,
    1.3948,
    1.3945,
    1.3942,
    1.3939,
    1.3937,
    1.3935,
    1.3932
  };
  
  G4double WackerCompanyGelPhotonEnergy[41] = {
    
    hc_eVnm / 150*eV,
    hc_eVnm / 210.2*eV,
    hc_eVnm / 214.0*eV,
    hc_eVnm / 217.8*eV,
    hc_eVnm / 219.1*eV,
    hc_eVnm / 222.9*eV,
    hc_eVnm / 228.0*eV,
    hc_eVnm / 238.1*eV,
    hc_eVnm / 250.8*eV,
    hc_eVnm / 268.6*eV,
    hc_eVnm / 299.1*eV,
    hc_eVnm / 330.0*eV,
    hc_eVnm / 340.0*eV,
    hc_eVnm / 350.0*eV,
    hc_eVnm / 360.0*eV,
    hc_eVnm / 370.0*eV,
    hc_eVnm / 380.0*eV,
    hc_eVnm / 390.0*eV,
    hc_eVnm / 400.0*eV,
    hc_eVnm / 410.0*eV,
    hc_eVnm / 420.0*eV,
    hc_eVnm / 430.0*eV,
    hc_eVnm / 440.0*eV,
    hc_eVnm / 450.0*eV,
    hc_eVnm / 460.0*eV,
    hc_eVnm / 470.0*eV,
    hc_eVnm / 480.0*eV,
    hc_eVnm / 490.0*eV,
    hc_eVnm / 500.0*eV,
    hc_eVnm / 510.0*eV,
    hc_eVnm / 520.0*eV,
    hc_eVnm / 530.0*eV,
    hc_eVnm / 540.0*eV,
    hc_eVnm / 550.0*eV,
    hc_eVnm / 560.0*eV,
    hc_eVnm / 570.0*eV,
    hc_eVnm / 580.0*eV,
    hc_eVnm / 590.0*eV,
    hc_eVnm / 600.0*eV,
    hc_eVnm / 610.0*eV,
    hc_eVnm / 730.0*eV
  };
  
  G4double WackerCompanyGelAbsLen[41] = {
    
    0*mm,
    0.9*mm,
    1.6*mm,
    2.3*mm,
    3.2*mm,
    4.5*mm,
    6.5*mm,
    13.0*mm,
    21.5*mm,
    33.2*mm,
    54.5*mm,
    230.8*mm,
    304.9*mm,
    371.4*mm,
    418.8*mm,
    457.1*mm,
    489.6*mm,
    532.9*mm,
    566.4*mm,
    593.8*mm,
    625.3*mm,
    644.8*mm,
    669.1*mm,
    680.5*mm,
    723.1*mm,
    745.5*mm,
    764.8*mm,
    781.8*mm,
    810.8*mm,
    844.9*mm,
    858.8*mm,
    869.5*mm,
    901.0*mm,
    890.9*mm,
    943.6*mm,
    964.2*mm,
    969.0*mm,
    998.9*mm,
    999.4*mm,
    1008.1*mm,
    1008.1*mm
  };
  
  
  // ----------------------- choosing gel for simulation -------------------------------------
  // Wacker SilGel 612 A/B (from KM3NeT data)
  G4MaterialPropertiesTable* proptable_WackerGel = new G4MaterialPropertiesTable();
  proptable_WackerGel->AddProperty("RINDEX", GeneralRIndPhotonEnergy, WackerGelRInd, 49);
  proptable_WackerGel->AddProperty("ABSLENGTH", WackerGelPhotonEnergy, WackerGelAbsLen, 34);
  
  // awesome mysterious Chiba gel
  G4MaterialPropertiesTable* proptable_ChibaGel = new G4MaterialPropertiesTable();
  proptable_ChibaGel->AddProperty("RINDEX", GeneralRIndPhotonEnergy, ChibaGelRInd,49);
  proptable_ChibaGel->AddProperty("ABSLENGTH", ChibaGelPhotonEnergy, ChibaGelAbsLen,36);
  
  // old IceCube gel
  G4MaterialPropertiesTable* proptable_IceCubeGel = new G4MaterialPropertiesTable();
  proptable_IceCubeGel->AddProperty("RINDEX", GeneralRIndPhotonEnergy, ChibaGelRInd, 49);		// using ChibaGelRInd !!
  proptable_IceCubeGel->AddProperty("ABSLENGTH", IceCubeGelPhotonEnergy, IceCubeGelAbsLen,17);
  
  // Wacker SilGel 612 A/B (from company transmission plot)
  G4MaterialPropertiesTable* proptable_WackerCompanyGel = new G4MaterialPropertiesTable();
  proptable_WackerCompanyGel->AddProperty("RINDEX", WackerCompanyGelRIndPhotonEnergy, WackerCompanyGelRInd, 59);
  proptable_WackerCompanyGel->AddProperty("ABSLENGTH", WackerCompanyGelPhotonEnergy, WackerCompanyGelAbsLen, 41);
  
  // chosen here:
  if (Gel_type == "Wacker"){
    Mat_Gel->SetMaterialPropertiesTable(proptable_WackerGel);
  }
  if (Gel_type == "Chiba"){
    Mat_Gel->SetMaterialPropertiesTable(proptable_ChibaGel);
  }
  if (Gel_type == "IceCube"){
    Mat_Gel->SetMaterialPropertiesTable(proptable_IceCubeGel);
  }
  if (Gel_type == "Wacker_company"){
    Mat_Gel->SetMaterialPropertiesTable(proptable_WackerCompanyGel);
  }
  
  // ------------------------ air & vacuum ---------------------------------------------------
  
  

  
  G4String AirSpectrumFile = "../Detector_construction_files/Scintillation_spectra/air.txt";
  std::vector<double> fileFirstColumn = readColumnDouble(AirSpectrumFile, 1);
  std::vector<double> fileSecondColumn = readColumnDouble(AirSpectrumFile, 2);  
  //G4double* airScintWL = &
  G4int ArraySize = fileFirstColumn.size();
  
  G4double airScintWLD[ArraySize];
  G4double airScintInD[ArraySize];
  
  G4cout << ArraySize << G4endl;
  
  for (unsigned int u = 0; u <fileFirstColumn.size(); u++) {
    airScintWLD[u] = hc_eVnm / fileFirstColumn.at(u)*eV;  
    airScintInD[u] = fileSecondColumn.at(u);
  }

  G4MaterialPropertiesTable* proptable_air = new G4MaterialPropertiesTable();

  proptable_air->AddConstProperty("SCINTILLATIONYIELD",scintYield);
  proptable_air->AddConstProperty("FIRSTAMPLITUDE",1);
  proptable_air->AddConstProperty("SECONDAMPLITUDE",0);
  proptable_air->AddConstProperty("THIRDAMPLITUDE",0);
  
  proptable_air->AddProperty("FIRSTCOMPONENT",airScintWLD,airScintInD, ArraySize);
  proptable_air->AddProperty("SECONDCOMPONENT",airScintWLD,airScintInD, ArraySize);
  proptable_air->AddProperty("THIRDCOMPONENT",airScintWLD,airScintInD, ArraySize);
  
  proptable_air->AddConstProperty("FIRSTTIME",2*ns);
  proptable_air->AddConstProperty("SECONDTIME",0);
  proptable_air->AddConstProperty("THIRDTIME",0);

  proptable_air->AddConstProperty("RESOLUTIONSCALE", 1.0);

  
  G4double VacuumPhotonEnergy[2] = {PHOTON_NRG_MIN, PHOTON_NRG_MAX};
  G4double VacuumRidx[2] = {1.0003, 1.0003};
  G4MaterialPropertiesTable* proptable_vacuum = new G4MaterialPropertiesTable();
  proptable_vacuum->AddProperty("RINDEX", VacuumPhotonEnergy, VacuumRidx, 2);
  Mat_Vacuum->SetMaterialPropertiesTable(proptable_vacuum);
  
  G4double AirRidx[2] = {1.0003, 1.0003};
  
  proptable_air->AddProperty("RINDEX", VacuumPhotonEnergy, AirRidx, 2);
  Mat_LabAir->SetMaterialPropertiesTable(proptable_air);
  
  // Print all the materials defined.
  //	G4cout << *(G4Material::GetMaterialTable());
  
  ////////////////
  //  Geometry  //
  ////////////////
  
  //	The "World": -----------------------------------------------------------------------------------------------------
  //Spherical world
  //World_solid = new G4Orb("World",gworldsize*m);
   //Cylinder world
  
  G4double innerRadiusOfTheTube = 0.*cm;
  G4double startAngleOfTheTube = 0.*deg;
  G4double spanningAngleOfTheTube = 360.*deg;
  gRadius = gworldsize*m;
  gHeight = gworldsize*m;
  
  World_solid
    = new G4Tubs("tracker_tube",
                 innerRadiusOfTheTube, 
                 gRadius,
                 gHeight,
                 startAngleOfTheTube, 
                 spanningAngleOfTheTube); 
  
  if (World_type == "air"){
    World_logical = new G4LogicalVolume(World_solid, Mat_LabAir, "World logical", 0, 0, 0);
  }
  if (World_type == "ice"){
    World_logical = new G4LogicalVolume(World_solid, Mat_Ice, "World logical", 0, 0, 0);	// xxx
  }
  if (World_type == "spice"){
    World_logical = new G4LogicalVolume(World_solid, Mat_Spice, "World logical", 0, 0, 0);
  }
  // 	if (World_type == "water"){
  // 		World_logical = new G4LogicalVolume(World_solid, Mat_Water, "World logical", 0, 0, 0);
  // 	}
  
  G4RotationMatrix* rot;
  
  
  //	Constructing of Hamamatsu R12199 and ETEL 9320KFL "three inch" PMTs: --------------------------------------------------
  
  G4double RefCone_SheetThickness = 0.5*mm;	// aluminum sheet thickness true for all reflective cones
  G4double RefCone_ConeToHolder = 1.55*mm;	// horizontal distance from Kärcher's construction 
  
  //	Hamamatsu R12199: spherical window --------------------------------------------------------------------------------------------------------------
  //	sphere and ellipsoid for front half, cylinder for back
  //	z-semiaxis: 23.3mm, x,y: 40mm
  //	radius: 40 mm
  G4double PMT_12199_OutRad = 52.0*mm;
  G4double PMT_12199_WallThick = 2*mm;
  G4double PC_12199_Rad = PMT_12199_OutRad - PMT_12199_WallThick;
  // 	G4double PMT_12199_SphereAngle = 50*degree; //explicit declaration
  G4double PMT_12199_SphereRad  = 33.5*mm;	// 34*mm found by superimposing circle and ellipse, 33.5*mm due to visualizer restrictions
  G4double PMT_12199_SphereAngle  = asin(PMT_12199_SphereRad / PMT_12199_OutRad);	// implicit using diameter
  
  PMT_12199_ellips_solid = new G4Ellipsoid("PMT_12199 solid bulb ellipsoid", 40*mm, 40*mm, 21.5*mm);
  PMT_12199_sphere_solid = new G4Sphere("PMT_12199 solid bulb sphere", 0.0, PMT_12199_OutRad, 0, 2*pi, 0, PMT_12199_SphereAngle);
  PMT_12199_bulk_solid = new G4Tubs("PMT_12199 solid bulk", 0.0, 0.5*51.9*mm, (97-23.3)*0.5*mm, 0, 2*pi); // Glass Cylinder
  
  G4ThreeVector PMT_12199_ort_1;
  
  PMT_12199_tube_solid = new G4UnionSolid("PMT_12199 tube solid", PMT_12199_ellips_solid, PMT_12199_sphere_solid, 0, G4ThreeVector(0,0,-28*mm));
  PMT_12199_ort_1 = G4ThreeVector(0,0,-(97-23.3)*0.5*mm);	
  PMT_12199_tube_solid = new G4UnionSolid("PMT_12199 tube solid", PMT_12199_tube_solid, PMT_12199_bulk_solid, 0, PMT_12199_ort_1);    
  PMT_12199_tube_logical = new G4LogicalVolume(PMT_12199_tube_solid, Mat_Tube_Glass, "PMT_12199 tube logical");
  /*
   *	G4VSolid* PMT_12199_bulkA_solid = new G4Tubs("PMT_12199 solid bulkA", 0.0, 0.5*47.9*mm, (97-27.3)*0.5*mm, 0, 2*pi); // Semi-Vacuum Cylinder solid !
   *	G4VSolid* PMT_12199_sphereA_solid = new G4Sphere("PMT_12199 solid bulb sphereA", 0.0, PMT_12199_OutRad-2*mm, 0, 2*pi, 0, PMT_12199_SphereAngle); // Semi-Vacuum Sphere !
   *	G4VSolid* PMT_12199_ellipsA_solid = new G4Ellipsoid("PMT_12199 solid bulb ellipsoidA", 38*mm, 38*mm, 19.5*mm); // Semi-Vacuum Ellipsoid solid !	
   *	G4VSolid* PMT_12199_tubeA_solid = new G4UnionSolid("PMT_12199 tubeA solid", PMT_12199_ellipsA_solid, PMT_12199_sphereA_solid, 0, G4ThreeVector(0,0,-28*mm));
   */
  G4double PC_12199_LimiterDZ = 25*mm;
  /*
   *	G4VSolid* PC_12199_aux2_solid = new G4Box("PC_12199 aux2 solid", PMT_12199_OutRad, PMT_12199_OutRad, PC_12199_LimiterDZ); 
   *	
   *	PMT_12199_tubeA_solid = new G4SubtractionSolid("PC_12199 tubeA solid", PMT_12199_tubeA_solid, PC_12199_aux2_solid, 0, G4ThreeVector(0,0,(PC_12199_LimiterDZ + 2.9*mm)));	
   *	PMT_12199_tubeA_solid = new G4UnionSolid("PMT_12199 tubeA solid", PMT_12199_tubeA_solid, PMT_12199_bulkA_solid, 0, PMT_12199_ort_1);
   *	
   * 
   *	G4LogicalVolume* PMT_12199_tubeA_logical = new G4LogicalVolume(PMT_12199_tubeA_solid, Mat_HighVacuum, "PMT_12199 tubeA logical");	// Semi-Vacuum Both logical !
   *	G4VPhysicalVolume* PMT_12199_tubeA_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), PMT_12199_tubeA_logical, "PMT_12199 tubeA physical", PMT_12199_tube_logical, false, 0);
   *	
   *	G4VSolid* PMT_12199_metal_solid = new G4Box("PMT_12199 Metal Box Solid", 15*mm, 12.5*mm, 12.5*mm);
   *	G4LogicalVolume* PMT_12199_metal_logical= new G4LogicalVolume(PMT_12199_metal_solid, Mat_Stahl, "PMT 12199 Metal Box Logical");
   *	G4VPhysicalVolume* PMT_12199_metal_physical = new G4PVPlacement(0, G4ThreeVector(0,0,-30*mm), PMT_12199_metal_logical, "PMT_12199 Metal Box physical", PMT_12199_tubeA_logical, false, 0);
   *	
   *	G4VSolid* PMT_12199_metalshield_solid = new G4Tubs("PC_12199 metal shield solid", 0, 0.5*43.9*mm, 1*mm, 0, 2*pi);
   *	G4LogicalVolume* PMT_12199_metalshield_logical = new G4LogicalVolume(PMT_12199_metalshield_solid, Mat_Stahl, "PC_12199 metal shield logical");
   * 	G4VPhysicalVolume* PMT_12199_metalshield_physical = new G4PVPlacement (0, G4ThreeVector(0,0,(-29+12.5)*mm), PMT_12199_metalshield_logical, "PC_12199 metal shield physical", PMT_12199_tubeA_logical, false, 0);
   *	
   *	
   *	G4VisAttributes* metalAttribute = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5,0.2));
   *        PMT_12199_metal_logical->SetVisAttributes(metalAttribute);
   *	PMT_12199_metalshield_logical->SetVisAttributes(metalAttribute);
   *	
   *	G4VisAttributes* tubeA_Attribute = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.2));
   *        PMT_12199_tubeA_logical->SetVisAttributes(tubeA_Attribute);*/
  
  
  //	12199 photocathode: combination of sphere and ellipsoid, cropped to account for real cathode coverage:
  G4double PC_12199_OpenAngle = 45*degree; // sphere opening angle explicit declaration, alternatively definde using diameter below
  // 	G4double PC_12199_Diameter  = 72*mm;	// guaranteed minimum by Hamamatsu
  // 	G4double PC_12199_OpenAngle  = asin(PC_12199_Diameter/(2*PC_12199_Rad));	// implicit, using effctive diameter
  G4double PC_12199_CathodeLimit = 3*mm;
  G4double PC_12199_ShieldToCathode = 0.15*mm;
  PC_12199_sphere_solid = new G4Sphere("PC_12199 sphere solid", 0.0*mm, PC_12199_Rad, 0, 2*pi, 0, PC_12199_OpenAngle);
  PC_12199_ellips_solid = new G4Ellipsoid("PC_12199 ellips solid", 38*mm, 38*mm, 19.5*mm);
  PC_12199_intermed_solid = new G4UnionSolid("PC_12199 intermed solid", PC_12199_ellips_solid, PC_12199_sphere_solid, 0, G4ThreeVector(0,0,-28*mm));
  PC_12199_aux_solid = new G4Box("PC_12199 aux solid", PMT_12199_OutRad, PMT_12199_OutRad, PC_12199_LimiterDZ); 
  PC_12199_solid = new G4SubtractionSolid("PC_12199 solid", PC_12199_intermed_solid, PC_12199_aux_solid, 0, G4ThreeVector(0,0,-(PC_12199_LimiterDZ - PC_12199_CathodeLimit)));
  PC_12199_logical = new G4LogicalVolume(PC_12199_solid, Mat_BiAlkali, "PC_12199 logical");
  // 	PC80_physical = new G4PVPlacement (0, G4ThreeVector(0,0,(-(50-23.3)-2)*mm), PC80_logical, "80mm photocathode", PMT80_tube_logical, false, 0);
  PC_12199_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), PC_12199_logical, "PC_12199 physical", PMT_12199_tube_logical, false, 0);
  //	G4double PC_12199_ShieldRad = (40.- 0.1)*mm * std::sqrt(1.-std::pow((PC_12199_CathodeLimit - PC_12199_ShieldToCathode)/(23.5*mm - 0.1*mm), 2.));	// automatic shield diameter based on ellipsoid formula, with radial offset
  G4double PC_12199_ShieldRad = 38.*mm * std::sqrt(1.-std::pow((PC_12199_CathodeLimit - PC_12199_ShieldToCathode)/(21.5*mm), 2.));	// // automatic shield diameter based on cathode ellipsoid formula
  PC_12199_shield_solid = new G4Tubs("PC_12199 shield solid", 0, PC_12199_ShieldRad, 0.1*mm, 0, 2*pi);
  PC_12199_shield_logical = new G4LogicalVolume(PC_12199_shield_solid, Mat_Absorber, "PC_12199 shield logical");
  PC_12199_shield_physical = new G4PVPlacement (0, G4ThreeVector(0,0,PC_12199_CathodeLimit - PC_12199_ShieldToCathode), PC_12199_shield_logical, "PC_12199 shield physical", PMT_12199_tube_logical, false, 0);
  
  //	//	12199 RefCone:
  // 	using Kärcher's geometry and real metal sheet thickness
  G4double RefCone_12199_InDiam = 82*mm;
  G4double RefCone_12199_OutDiam = 107.5*mm;
  G4double RefCone_12199_dZ = (RefCone_12199_OutDiam - RefCone_12199_InDiam)*0.25;
  RefCone_12199_solid = new G4Cons("RefCone_12199 solid",
				   0.5*RefCone_12199_InDiam, 0.5*RefCone_12199_InDiam + RefCone_SheetThickness*std::sqrt(2.), 
				   0.5*RefCone_12199_OutDiam, 0.5*RefCone_12199_OutDiam + RefCone_SheetThickness*std::sqrt(2.), 
				   RefCone_12199_dZ, 0, 2*pi);
  RefCone_12199_logical = new G4LogicalVolume(RefCone_12199_solid, Mat_Reflector, "RefCone_12199 logical");
  
  //	12199 "Nests" for real reflective cones:
  G4Cons* RefConeNest_12199_cone_solid = new G4Cons("RefConeNest_12199_cone solid", 0, 0.5*RefCone_12199_InDiam + RefCone_ConeToHolder, 0, 0.5*RefCone_12199_InDiam + RefCone_ConeToHolder + 2*30*mm, 30*mm, 0, 2*pi);
  G4UnionSolid* RefConeNest_12199_solid = new G4UnionSolid("RefConeNest_12199 solid", PMT_12199_tube_solid, RefConeNest_12199_cone_solid, 0, G4ThreeVector(0,0,30*mm));
  
  
  
  //	ETEL 9320KFL: ellipsoidal window -------------------------------------------------------------------------------------------------------------------
  //	ellipsoid for front half, cylinder for back
  //	z-semiaxis: 26.2 mm, x,y: 0.5*86.5 mm
  G4ThreeVector PMT_ETEL_ort_1;
  G4RotationMatrix* PMT_ETEL_rot_1 = new G4RotationMatrix();
  PMT_ETEL_ort_1 = G4ThreeVector(0,0,-(95.-26.2)*0.5*mm);
  PMT_ETEL_bulb_solid = new G4Ellipsoid("PMT_ETEL solid bulb", 0.5*86.5*mm, 0.5*86.5*mm, 26.2*mm);
  PMT_ETEL_bulk_solid = new G4Tubs("PMT_ETEL solid bulk", 0.0, 0.5*52.4*mm, (95-26.2)*0.5*mm, 0, 2*pi);
  PMT_ETEL_tube_solid = new G4UnionSolid("PMT_ETEL tube", PMT_ETEL_bulb_solid, PMT_ETEL_bulk_solid, PMT_ETEL_rot_1, PMT_ETEL_ort_1);
  PMT_ETEL_tube_logical = new G4LogicalVolume(PMT_ETEL_tube_solid, Mat_Tube_Glass, "PMT ETEL tube");
  
  //	ETEL photocathode: ellipsoid, assuming glass thickness of 2mm
  G4double PC_ETEL_CathodeLimit = 4*mm; 		// distance between ideal and real cathode coverage, not always at the same level for indicidual PMTs
  G4double PC_ETEL_ShieldToCathode = 0.15*mm;	// distance between back side of cathode and anti photon shield
  PC_ETEL_solid = new G4Ellipsoid("PC_ETEL solid", (0.5*86.5-2.)*mm, (0.5*86.5-2.)*mm, 24.2*mm, PC_ETEL_CathodeLimit, 24.2*mm);
  PC_ETEL_logical = new G4LogicalVolume(PC_ETEL_solid, Mat_BiAlkali, "PC_ETEL logical");
  PC_ETEL_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), PC_ETEL_logical, "PC_ETEL physical", PMT_ETEL_tube_logical, false, 0);
  G4double PC_ETEL_ShieldRad = (0.5*86.5-0.1)*mm * std::sqrt(1.-std::pow((PC_ETEL_CathodeLimit - PC_ETEL_ShieldToCathode)/(26.2*mm -0.1*mm), 2.));	// automatic shield diameter based on bulb ellipsoid formula with radial offset 
  // 	G4double PC_ETEL_ShieldRad = (0.5*86.5-2.)*mm * std::sqrt(1.-std::pow((PC_ETEL_CathodeLimit - PC_ETEL_ShieldToCathode)/(26.2*mm -2.*mm), 2.));	// automatic shield diameter based on cathode ellipsoid formula
  PC_ETEL_shield_solid = new G4Tubs("PC_ETEL shield solid", 0, PC_ETEL_ShieldRad, 0.1*mm, 0, 2*pi);
  PC_ETEL_shield_logical = new G4LogicalVolume(PC_ETEL_shield_solid, Mat_Absorber, "PC_ETEL shield logical");
  PC_ETEL_shield_physical = new G4PVPlacement (0, G4ThreeVector(0,0, PC_ETEL_CathodeLimit - PC_ETEL_ShieldToCathode), PC_ETEL_shield_logical, "PC_ETEL shield physical", PMT_ETEL_tube_logical, false, 0);
  
  //	ETEL RefCone:
  // 	using Kärcher's geometry and real metal sheet thickness
  G4double RefCone_ETEL_InDiam = 90*mm;
  G4double RefCone_ETEL_OutDiam = 107.5*mm;
  G4double RefCone_ETEL_dZ = (RefCone_ETEL_OutDiam - RefCone_ETEL_InDiam)*0.25;
  RefCone_ETEL_solid = new G4Cons("RefCone_ETEL solid",
				  0.5*RefCone_ETEL_InDiam, 0.5*RefCone_ETEL_InDiam + RefCone_SheetThickness*std::sqrt(2.), 
				  0.5*RefCone_ETEL_OutDiam, 0.5*RefCone_ETEL_OutDiam + RefCone_SheetThickness*std::sqrt(2.), 
				  RefCone_ETEL_dZ, 0, 2*pi);
  RefCone_ETEL_logical = new G4LogicalVolume(RefCone_ETEL_solid, Mat_Reflector, "RefCone_ETEL logical");
  
  //	ETEL "Nests" for real reflective cones:
  G4Cons* RefConeNest_ETEL_cone_solid = new G4Cons("RefConeNest_ETEL_cone solid", 0, 0.5*90*mm + RefCone_ConeToHolder, 0, (0.5*90+2*30)*mm + RefCone_ConeToHolder, 30*mm, 0, 2*pi);
  G4UnionSolid* RefConeNest_ETEL_solid = new G4UnionSolid("RefConeNest_ETEL solid", PMT_ETEL_tube_solid, RefConeNest_ETEL_cone_solid, 0, G4ThreeVector(0,0,30*mm));
  
  
  //	Hamamatsu R12199: ellipsoidal window, deviation from real surface ca. 1 mm (Björn style) -------------------------------------------------------------
  //	ellipsoid for front half, cylinder for back, 
  //	z-semiaxis: 23.3mm, x,y: 40mm
  G4ThreeVector PMT80_ort_1;
  G4RotationMatrix* PMT80_rot_1 = new G4RotationMatrix();
  PMT80_ort_1 = G4ThreeVector(0,0,-(95-23.3)*0.5*mm);
  
  PMT80_solid_1 = new G4Ellipsoid("PMT80 solid bulb", 40*mm, 40*mm, 23.3*mm);
  PMT80_solid_2 = new G4Tubs("PMT80 solid bulk", 0.0, 0.5*51.9*mm, (95-23.3)*0.5*mm, 0, 2*pi);
  PMT80_tube_solid = new G4UnionSolid("PMT80 tube", PMT80_solid_1, PMT80_solid_2, PMT80_rot_1, PMT80_ort_1);
  PMT80_tube_logical = new G4LogicalVolume(PMT80_tube_solid, Mat_Tube_Glass, "PMT 80 tube");
  
  // 	photocathode:
  G4Ellipsoid* PC80_solid;
  PC80_solid = new G4Ellipsoid("PC80 solid", 38*mm, 38*mm, 21.3*mm,0*mm,21.3*mm);
  PC80_logical = new G4LogicalVolume(PC80_solid, Mat_BiAlkali, "Photocathode 80mm logical");
  PC80_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), PC80_logical, "80mm photocathode physical", PMT80_tube_logical, false, 0);
  PC80_shield_solid = new G4Tubs("PC80 shield solid", 0, 39.9*mm, 0.1*mm, 0, 2*pi);
  PC80_shield_logical = new G4LogicalVolume(PC80_shield_solid, Mat_Absorber, "PC80 shield logical");
  PC80_shield_physical = new G4PVPlacement (0, G4ThreeVector(0,0,-0.15*mm), PC80_shield_logical, "PC80 shield physical", PMT80_tube_logical, false, 0);
  
  //	PMT80 "Nests" for reflective cones:
  G4Cons* RefConeNest_80_cone_solid = new G4Cons("RefConeNest_80_cone solid", 0, 0.5*RefCone_12199_InDiam + RefCone_ConeToHolder, 0, 0.5*RefCone_12199_InDiam + RefCone_ConeToHolder + 2*30*mm, 30*mm, 0, 2*pi);
  G4UnionSolid* RefConeNest_80_solid = new G4UnionSolid("RefConeNest_80 solid", PMT_12199_tube_solid, RefConeNest_12199_cone_solid, 0, G4ThreeVector(0,0,30*mm));
  
  
  //	Hamamatsu R8071 10" PMT ---------------------------------------------------------------------------------------------------------------------
  // 	PMT center = center of sphere (PMT10_solid_1)
  
  G4double Tube10_geo_curvature_radius = 136.7*mm;
  G4double Tube10_geo_curvature_angle = asin(222*mm/(2*Tube10_geo_curvature_radius));
  G4double PC10_geo_curvature_angle = asin(220*mm/(2*(Tube10_geo_curvature_radius-2*mm)));
  G4double Tube10_back_radii[] = {110*mm, 42*mm, 42*mm};
  G4double Tube10_back_zplanes[] = {-20,-80*mm,-160*mm};
  G4double Tube10_back_inner[] = {0, 0, 0};
  G4ThreeVector PMT10_ort_1, PMT10_ort_2, PMT10_ort_3;
  G4RotationMatrix* PMT10_rot_1 = new G4RotationMatrix();
  PMT10_ort_1 = G4ThreeVector(0,0,50*mm);
  G4RotationMatrix* PMT10_rot_2 = new G4RotationMatrix();
  PMT10_ort_2 = G4ThreeVector(0,0,50*mm);
  G4RotationMatrix* PMT10_rot_3 = new G4RotationMatrix();
  PMT10_ort_3 = G4ThreeVector(0,0,60*mm);
  
  PMT10_solid_1 = new G4Sphere("PMT10 solid front half", 90*mm, Tube10_geo_curvature_radius, 0, 2*pi, 0, Tube10_geo_curvature_angle);
  PMT10_solid_3 = new G4Ellipsoid("PMT10 solid bulk",127*mm,127*mm,70*mm,0*mm,70*mm);
  PMT10_solid_4 = new G4Ellipsoid("PMT10 solid bulk 2",127*mm,127*mm,60*mm,-60*mm,0*mm);
  PMT10_solid_2 = new G4Polycone("10in tube back", 0 , 2*pi, 3, Tube10_back_zplanes, Tube10_back_inner, Tube10_back_radii);
  PMT10_front_solid = new G4UnionSolid("PMT10 front", PMT10_solid_1, PMT10_solid_3,PMT10_rot_1,PMT10_ort_1);
  PMT10_back_solid = new G4UnionSolid("PMT10 back", PMT10_front_solid, PMT10_solid_4,PMT10_rot_2,PMT10_ort_2);
  PMT10_tube_solid = new G4UnionSolid("PMT10 tube", PMT10_back_solid, PMT10_solid_2, PMT10_rot_3, PMT10_ort_3);
  // PMT damit platzieren:
  PMT10_tube_logical = new G4LogicalVolume(PMT10_tube_solid, Mat_Tube_Glass, "PMT 10 tube");
  
  PC10_solid = new G4Sphere("PC10 solid", Tube10_geo_curvature_radius-2.1*mm, Tube10_geo_curvature_radius-2*mm, 0, 2*pi, 0, PC10_geo_curvature_angle);
  PC10_logical = new G4LogicalVolume(PC10_solid, Mat_BiAlkali, "Photocathode 10in logical");
  PC10_physical = new G4PVPlacement (0, G4ThreeVector(0,0,-2*mm), PC10_logical, "Photocathode", PMT10_tube_logical, false, 0);
  
  PC10_shield_solid = new G4Tubs("PC10 shield solid", 0, 110*mm, 0.1*mm, 0, 2*pi);
  PC10_shield_logical = new G4LogicalVolume(PC10_shield_solid, Mat_Absorber, "PC10 shield logical");
  PC10_shield_physical = new G4PVPlacement(0, G4ThreeVector(0, 0, sqrt((Tube10_geo_curvature_radius-2)*(Tube10_geo_curvature_radius-2)-110.0*110.0)*mm-2.5*mm), PC10_shield_logical, "PC10 shield physical", PMT10_tube_logical, false, 0);
  
  
  //	DOM construction ------------------------------------------------------------------------------------------------------------------------------
  
  //	Glass
  // 	GlassSphereTop_solid = new G4Sphere("GlassSphereTop solid", 0, GlasOutRad, 0, 2*pi, 0, 0.51*pi);
  // 	GlassSphereBottom_solid = new G4Sphere("GlassSphereBottom solid", 0, GlasOutRad, 0, 2*pi, 0.49*pi, pi);
  GlassSphereTop_solid = new G4Ellipsoid("GlassSphereTop solid", GlasOutRad, GlasOutRad, GlasOutRad, -5*mm, GlasOutRad+5*mm);
  GlassSphereBottom_solid = new G4Ellipsoid("GlassSphereBottom solid",GlasOutRad, GlasOutRad, GlasOutRad, -(GlasOutRad+5*mm), 5*mm);
  GlassCylinder_solid = new G4Tubs("GlassCylinder solid", 0, GlasOutRad, CylHigh, 0, 2*pi);
  transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
  G4UnionSolid* temp_union = new G4UnionSolid("temp", GlassCylinder_solid, GlassSphereTop_solid, transformers);
  transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-CylHigh));
  Glass_solid = new G4UnionSolid("OM glass body", temp_union, GlassSphereBottom_solid, transformers);
  
  //  Gel
  // 	GelSphereTop_solid = new G4Sphere("GelSphereTop solid", 0, GlasOutRad - GlasThick, 0, 2*pi, 0, 0.51*pi);
  // 	GelSphereBottom_solid = new G4Sphere("GelSphereBottom solid", 0, GlasOutRad - GlasThick, 0, 2*pi, 0.49*pi, pi);
  GelSphereTop_solid = new G4Ellipsoid("GelSphereTop solid", GlasInRad, GlasInRad, GlasInRad, -5*mm, GlasInRad+5*mm);
  GelSphereBottom_solid = new G4Ellipsoid("GelSphereBottom solid", GlasInRad, GlasInRad, GlasInRad, -(GlasInRad+5*mm), 5*mm);
  GelCylinder_solid = new G4Tubs("GelCylinder solid", 0, GlasOutRad - GlasThick, CylHigh , 0, 2*pi);
  transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
  G4UnionSolid* temp_union2 = new G4UnionSolid("temp2", GelCylinder_solid, GelSphereTop_solid, transformers);
  transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-CylHigh));
  Gel_solid = new G4UnionSolid("gel body", temp_union2, GelSphereBottom_solid, transformers);
  
  //  PMT TubeHolder from "foam" primitives & cutting "nests" for PMTs later
  // 	G4Sphere* FoamSphereTop_solid = new G4Sphere("FoamSphereTop solid", 0, GlasOutRad - GlasThick - GelThick, 0, 2*pi, 0, 0.51*pi);
  // 	G4Sphere* FoamSphereBottom_solid = new G4Sphere("FoamSphereBottom solid", 0, GlasOutRad - GlasThick - GelThick, 0, 2*pi, 0.49*pi, pi);
  G4double FoamRad = GlasOutRad - GlasThick - GelThick;
  G4Ellipsoid* FoamSphereTop_solid = new G4Ellipsoid("FoamSphereTop solid", FoamRad, FoamRad, FoamRad, -5*mm, FoamRad+5*mm);
  G4Ellipsoid* FoamSphereBottom_solid = new G4Ellipsoid("FoamSphereBottom solid", FoamRad, FoamRad, FoamRad, -(FoamRad+5*mm), 5*mm);		
  G4Tubs* FoamCylinder_solid = new G4Tubs("FoamCylinder solid", 0, GlasOutRad - GlasThick - GelThick, CylHigh , 0, 2*pi);
  transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(CylHigh )));
  G4UnionSolid* Foam_TempUnion_solid = new G4UnionSolid("Foam TempUnion solid", FoamCylinder_solid, FoamSphereTop_solid, transformers);
  transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-(CylHigh )));
  G4UnionSolid* Foam_solid = new G4UnionSolid("Foam solid", Foam_TempUnion_solid, FoamSphereBottom_solid, transformers);
  
  //	Reflective cones for better "beaming" of angular acceptance (RefCone)	
  G4double RefCone_IdealInRad;
  G4double PMToffset;
  
  if (PMT_type == "12199s"){ 
    PMToffset = 23.3*mm;	
    // 		RefCone_offset = 1*mm;
    RefCone_IdealInRad = 0.5*RefCone_12199_InDiam;
    RefConeNestCone_solid = new G4Cons("RefConeNestCone", 0, RefCone_IdealInRad + RefCone_ConeToHolder / cos(RefCone_angle), 0, RefCone_IdealInRad + RefCone_ConeToHolder / cos(RefCone_angle) + 2*1.5*RefConeDZ*tan(RefCone_angle), 1.5*RefConeDZ, 0, 2*pi);
    RefConeNest_solid = new G4UnionSolid("RefConeNest", PMT_12199_tube_solid, RefConeNestCone_solid, 0, G4ThreeVector(0,0,1.5*RefConeDZ));
    // RefConeNestCone_solid = new G4Cons("RefConeNest_12199_cone solid", 0, 0.5*RefCone_12199_InDiam + RefCone_ConeToHolder, 0, 0.5*RefCone_12199_InDiam + RefCone_ConeToHolder + 2*30*mm, 30*mm, 0, 2*pi);
    // RefConeNest_solid = new G4UnionSolid("RefConeNest_12199 solid", PMT_12199_tube_solid, RefConeNestCone_solid, 0, G4ThreeVector(0,0,30*mm));
  }
  if (PMT_type == "etel"){ 
    PMToffset = 26.2*mm;
    // 		RefCone_offset = 8*mm;
    RefCone_IdealInRad = 0.5*RefCone_ETEL_InDiam;
    RefConeNestCone_solid = new G4Cons("RefConeNestCone", 0, RefCone_IdealInRad + RefCone_ConeToHolder / cos(RefCone_angle), 0, RefCone_IdealInRad + RefCone_ConeToHolder / cos(RefCone_angle) + 2*1.5*RefConeDZ*tan(RefCone_angle), 1.5*RefConeDZ, 0, 2*pi);
    RefConeNest_solid = new G4UnionSolid("RefConeNest", PMT_ETEL_tube_solid, RefConeNestCone_solid, 0, G4ThreeVector(0,0,1.5*RefConeDZ));		
    // RefConeNestCone_solid = new G4Cons("RefConeNest_ETEL_cone solid", 0, 0.5*90*mm + RefCone_ConeToHolder, 0, (0.5*90+2*30)*mm + RefCone_ConeToHolder, 30*mm, 0, 2*pi);
    // RefConeNest_solid = new G4UnionSolid("RefConeNest_ETEL solid", PMT_ETEL_tube_solid, RefConeNestCone_solid, 0, G4ThreeVector(0,0,30*mm));
  }
  if (PMT_type == "12199e"){ 
    PMToffset = 23.3*mm;
    // 		RefCone_offset = 1*mm;
    RefCone_IdealInRad = 0.5*RefCone_12199_InDiam;
    RefConeNestCone_solid = new G4Cons("RefConeNestCone", 0, RefCone_IdealInRad + RefCone_ConeToHolder / cos(RefCone_angle), 0, RefCone_IdealInRad + RefCone_ConeToHolder / cos(RefCone_angle) + 2*1.5*RefConeDZ*tan(RefCone_angle), 1.5*RefConeDZ, 0, 2*pi);
    RefConeNest_solid = new G4UnionSolid("RefConeNest", PMT80_tube_solid, RefConeNestCone_solid, 0, G4ThreeVector(0,0,1.5*RefConeDZ));		
    // RefConeNestCone_solid = new G4Cons("RefConeNest_80_cone solid", 0, 0.5*RefCone_12199_InDiam + RefCone_ConeToHolder, 0, 0.5*RefCone_12199_InDiam + RefCone_ConeToHolder + 2*30*mm, 30*mm, 0, 2*pi);
    // RefConeNest_solid = new G4UnionSolid("RefConeNest_80 solid", PMT80_tube_solid, RefConeNestCone_solid, 0, G4ThreeVector(0,0,30*mm));
  }
  
  
  if (OM_type == "mDOM") {
    //Air slice between halfspheres
//     G4double AirSliceHigh = 0.5*um;
//     G4Tubs* AirCylinder_solid = new G4Tubs("AirCylinder solid", GlasInRad, GlasOutRad, AirSliceHigh, 0, 2*pi);
    
    // Producing PMT & RefCone coordinates
    G4double PMT_theta[99], PMT_phi[99], PMT_x[99], PMT_y[99], PMT_z[99], RefCone_x[99], RefCone_y[99], RefCone_z[99]; 
    G4double PMT_rho;
    G4double RefCone_rho;
    G4double PMT_z_offset;
    G4double PMT_r;
    G4double RefCone_r;
    
    for (i = 0; i <= 23; i++) {
      PMT_r = GlasInRad - GelPMT - PMToffset;        // radius for PMT positioning
      RefCone_r = GlasInRad - GelPMT - PMToffset + RefConeDZ; // radius for RefCone positioning
      if (i>=0 && i<=3){
        PMT_theta[i]=33.0*deg;
        PMT_phi[i]=i*90.0*deg+45*deg;
        PMT_z_offset = CylHigh;
      }
      if (i>=4 && i<=11){
        PMT_theta[i]=72.0*deg;
        PMT_phi[i]=(22.5+(i-4)*45.0)*deg+45*deg;
//         PMT_phi[i]=(0+(i-4)*45.0)*deg;
        PMT_z_offset = CylHigh - MPMTzoffset;
        PMT_r += MPMTroffset;
           RefCone_r += MPMTroffset;
      }
      if (i>=12 && i<=19){
    PMT_theta[i]=108.0*deg;
    //PMT_phi[i]=(i-12)*45.0*deg; //Verdreht
    PMT_phi[i]=(22.5+(i-12)*45.0)*deg+45*deg;
    PMT_z_offset = - CylHigh + MPMTzoffset;
    PMT_r += MPMTroffset;
    RefCone_r += MPMTroffset;
      }
      if (i>=20 && i<=23){
        PMT_theta[i]=147.0*deg;
        //PMT_phi[i]=(22.5+(i-20)*90.0)*deg;//Verdreht
         PMT_phi[i]=(0+(i-20)*90.0)*deg+45*deg;
        PMT_z_offset = - CylHigh;
      } 
      
      //	  G4cout << i << " " << PMT_theta[i] << " " << PMT_phi[i] << G4endl;
      
      PMT_rho = PMT_r * sin(PMT_theta[i]);
      PMT_x[i] = PMT_rho * cos(PMT_phi[i]);
      PMT_y[i] = PMT_rho * sin(PMT_phi[i]);
      PMT_z[i] = PMT_r * cos(PMT_theta[i]) + PMT_z_offset;
      
      RefCone_rho = RefCone_r * sin(PMT_theta[i]);
      RefCone_x[i] = RefCone_rho * cos(PMT_phi[i]);
      RefCone_y[i] = RefCone_rho * sin(PMT_phi[i]);
      RefCone_z[i] = RefCone_r * cos(PMT_theta[i]) + PMT_z_offset;
      
    }	  
    
    //Cutting the PMT nests from foam...
    for (k = 0; k <=23; k++) {
      rot = new G4RotationMatrix();
      rot->rotateY(PMT_theta[k]);
      rot->rotateZ(PMT_phi[k]);
      transformers = G4Transform3D(*rot, G4ThreeVector(PMT_x[k],PMT_y[k],PMT_z[k]));
      if (k==0){
	TubeHolder_solid = new G4SubtractionSolid("TubeHolder solid", Foam_solid, RefConeNest_solid, transformers);
      } 		
      else { 
	TubeHolder_solid = new G4SubtractionSolid("TubeHolder solid", TubeHolder_solid, RefConeNest_solid, transformers);
      }
    }
    
    // building reflective cones:
    // old style:
    // RefConeBasic_solid = new G4Cons("RefConeBasic", 
    // RefCone_IdealInRad, 
    // RefCone_IdealInRad + RefCone_SheetThickness*std::sqrt(2.), 
    // RefCone_IdealInRad + 2*RefConeDZ, 
    // RefCone_IdealInRad + RefCone_SheetThickness*std::sqrt(2.) + 2*RefConeDZ, 
    // RefConeDZ, 0, 2*pi);
    
    RefConeBasic_solid = new G4Cons("RefConeBasic", RefCone_IdealInRad, RefCone_IdealInRad + RefCone_SheetThickness / cos(RefCone_angle), 
				    RefCone_IdealInRad + 2*RefConeDZ*tan(RefCone_angle), RefCone_IdealInRad + RefCone_SheetThickness / cos(RefCone_angle) + 2*RefConeDZ*tan(RefCone_angle), RefConeDZ, 0, 2*pi);									
    
    // reflective cone type1 for spherical part (upper and lowermost rings):
    // old
    // rot = new G4RotationMatrix();
    // transformers = G4Transform3D(*rot, G4ThreeVector(0,0,-(GlasInRad-GelPMT-PMToffset+RefConeDZ)));
    // RefConeType1_solid = new G4IntersectionSolid("RefConeType1", RefConeBasic_solid, FoamSphereTop_solid, transformers);
    // new
    RefConeType1_solid = new G4IntersectionSolid("RefConeType1", RefConeBasic_solid, FoamSphereTop_solid, 0, G4ThreeVector(0,0,-(GlasInRad-GelPMT-PMToffset+RefConeDZ)));
    
    // reflective cone type 2 for upper central ring:	
    G4double RefConeType2_r = GlasInRad - GelPMT - PMToffset + RefConeDZ + MPMTroffset;
    G4double RefConeType2_rho = RefConeType2_r * sin(72*deg);
    G4double RefConeType2_x = RefConeType2_rho * cos(0*deg);
    G4double RefConeType2_y = RefConeType2_rho * sin(0*deg);
    G4double RefConeType2_z = RefConeType2_r * cos(72*deg) + CylHigh - MPMTzoffset;
    
    rot = new G4RotationMatrix();
    rot->rotateY(72*deg);
    transformers = G4Transform3D(*rot, G4ThreeVector(RefConeType2_x, RefConeType2_y, RefConeType2_z));
    
    RefConeType2_solid = new G4IntersectionSolid("RefConeType2", Foam_solid, RefConeBasic_solid, transformers);
    
    // reflective cone type 3 for lower central ring:	
    G4double RefConeType3_r = GlasInRad - GelPMT - PMToffset + RefConeDZ + MPMTroffset;
    G4double RefConeType3_rho = RefConeType3_r * sin(108*deg);
    G4double RefConeType3_x = RefConeType3_rho * cos(0*deg);
    G4double RefConeType3_y = RefConeType3_rho * sin(0*deg);
    G4double RefConeType3_z = RefConeType3_r * cos(108*deg) - CylHigh + MPMTzoffset;
    
    rot = new G4RotationMatrix();
    rot->rotateY(108*deg);
    transformers = G4Transform3D(*rot, G4ThreeVector(RefConeType3_x, RefConeType3_y, RefConeType3_z));
    
    RefConeType3_solid = new G4IntersectionSolid("RefConeType3", Foam_solid, RefConeBasic_solid, transformers);
    
    // etel cone needs additional cropping to prevent collisions:
    // 	if (PMT_type == "etel"){
    G4Cons* RefConeCutter_solid = new G4Cons("RefConeCutter", 
					     0, RefCone_IdealInRad + RefCone_SheetThickness*std::sqrt(2.) + 2*mm, 
					     0, RefCone_IdealInRad + RefCone_SheetThickness*std::sqrt(2.) + 2*RefConeDZ + 2*mm, 
					     RefConeDZ, 0, 2*pi);
    
    // cropping type 2 cone:
    G4double RefConeCut_r = GlasInRad - GelPMT - PMToffset + RefConeDZ + MPMTroffset;
    G4double RefConeCut_rho = RefConeCut_r * sin(72*deg);
    G4double RefConeCut_x = RefConeCut_rho * cos(45*deg);
    G4double RefConeCut_y = RefConeCut_rho * sin(45*deg);
    G4double RefConeCut_z = RefConeCut_r * cos(72*deg) + CylHigh - MPMTzoffset;
    
    rot = new G4RotationMatrix();
    rot->rotateY(72*deg);
    rot->rotateZ(45*deg);
    transformers = G4Transform3D(*rot, G4ThreeVector(RefConeCut_x, RefConeCut_y, RefConeCut_z));
    
    G4SubtractionSolid* RefConeType2_ETEL_solid = new G4SubtractionSolid("RefConeType2 ETEL solid", 
									 RefConeType2_solid,
									 RefConeCutter_solid,
									 transformers
    );
    
    RefConeCut_x = RefConeCut_rho * cos(-45*deg);
    RefConeCut_y = RefConeCut_rho * sin(-45*deg);
    
    rot = new G4RotationMatrix();
    rot->rotateY(72*deg);
    rot->rotateZ(-45*deg);
    transformers = G4Transform3D(*rot, G4ThreeVector(RefConeCut_x, RefConeCut_y, RefConeCut_z));
    
    RefConeType2_ETEL_solid = new G4SubtractionSolid("RefConeType2 ETEL solid", RefConeType2_ETEL_solid, RefConeCutter_solid, transformers);
    
    // cropping type 3 cone:
    RefConeCut_rho = RefConeCut_r * sin(108*deg);
    RefConeCut_x = RefConeCut_rho * cos(45*deg);
    RefConeCut_y = RefConeCut_rho * sin(45*deg);
    RefConeCut_z = RefConeCut_r * cos(108*deg) - CylHigh + MPMTzoffset;
    
    rot = new G4RotationMatrix();
    rot->rotateY(108*deg);
    rot->rotateZ(45*deg);
    transformers = G4Transform3D(*rot, G4ThreeVector(RefConeCut_x, RefConeCut_y, RefConeCut_z));
    
    G4SubtractionSolid* RefConeType3_ETEL_solid = new G4SubtractionSolid("RefConeType3 ETEL solid", 
									 RefConeType3_solid,
									 RefConeCutter_solid,
									 transformers
    );
    
    RefConeCut_x = RefConeCut_rho * cos(-45*deg);
    RefConeCut_y = RefConeCut_rho * sin(-45*deg);
    
    rot = new G4RotationMatrix();
    rot->rotateY(108*deg);
    rot->rotateZ(-45*deg);
    transformers = G4Transform3D(*rot, G4ThreeVector(RefConeCut_x, RefConeCut_y, RefConeCut_z));
    
    RefConeType3_ETEL_solid = new G4SubtractionSolid("RefConeType3 ETEL solid", RefConeType3_ETEL_solid, RefConeCutter_solid, transformers);
    
    G4LogicalVolume* RefConeType2_ETEL_logical = new G4LogicalVolume(RefConeType2_ETEL_solid, Mat_Reflector, "RefConeType2 ETEL logical");
    G4LogicalVolume* RefConeType3_ETEL_logical = new G4LogicalVolume(RefConeType3_ETEL_solid, Mat_Reflector, "RefConeType3 ETEL logical");
    // 	}
    
    
    //mDOM Harness
		// all called pDOM but we are still in mDOM
        G4double PDOM_Harness_inner[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; //Inner radii of the harness
		//G4double GlasOutRad_abitbigger = GlasOutRad + 1*mm;
		//G4double PDOM_Harness_inner[] = {GlasOutRad_abitbigger, GlasOutRad_abitbigger, GlasOutRad_abitbigger, GlasOutRad_abitbigger, GlasOutRad_abitbigger, GlasOutRad_abitbigger, GlasOutRad_abitbigger, GlasOutRad_abitbigger, GlasOutRad_abitbigger, GlasOutRad_abitbigger}; //Inner radii of the harness
        G4double PDOM_Harness_radii[] = {179 * mm, (196.4 - 8) * mm, (196.4 - 6.35) * mm, 199 * mm, 196.4 * mm, 196.4 * mm, 199 * mm, (196.4 - 6.35) * mm, (196.4 - 8) * mm, 179 * mm}; //Outer radii of the harness
        G4double PDOM_Harness_zplanes[] = {-23 * mm, -22.5 * mm, -22.5 * mm, -10 * mm, -10 * mm, 10 * mm, 10 * mm, 22.5 * mm, 22.5 * mm, 23 * mm};  //Corresponding z values for the radii
        //PDOM_GlassSphere_solid = new G4Orb("PDOM_GlassSphere solid", GlasOutRad); //Radius of the sphere to be cut from the harness

        PDOM_HarnessAux_solid = new G4Polycone("PDOM_HarnessAux solid", 0, 2 * pi, 10, PDOM_Harness_zplanes, PDOM_Harness_inner, PDOM_Harness_radii);   //Creating harness via polycone geometry
        PDOM_Harness_solid = new G4SubtractionSolid("PDOM_Harness solid", PDOM_HarnessAux_solid, Glass_solid);   //Cutting PDOM_GlassSphere_solid from PDOM_HarnessAux_solid
		//PDOM_Harness_solid = new G4SubtractionSolid("PDOM_Harness solid", PDOM_HarnessAux_solid, CylinderToSubstract_solid);   //Cutting PDOM_GlassSphere_solid from PDOM_HarnessAux_solid
        PDOM_Harness_logical = new G4LogicalVolume (PDOM_Harness_solid, Mat_Stahl, "PDOM_Harness logical");
        PDOM_HarnessSurface = new G4LogicalSkinSurface("PDOM_Harness_skin", PDOM_Harness_logical, PDOM_Harness_optical);
		
	//mDOM Ropes!
		G4double ropeThickness = 6.35 / 2 * mm;                                                   //Radius of the rope
		G4double ropeLength = 1000. * mm / 2;         //Half the length of the rope
	
		Rope_solid = new G4Tubs("Rope_solid", 0, ropeThickness, ropeLength, 0, 2 * pi);         //Creating a rope
		Rope_logical = new G4LogicalVolume(Rope_solid, Mat_Absorber, "Rope_logical");           //Logical for the solid
		Rope_surface = new G4LogicalSkinSurface("Rope_surface", Rope_logical, PDOM_Harness_optical); //same optical properties than the harness
		
    
    //  Placing da stuff ----------------------------------------------------------------------------------------------------------------------------------
    rot = new G4RotationMatrix();
    
    G4RotationMatrix* flipOM = new G4RotationMatrix();
    // 			flipOM->rotateY(180*deg); // for hunting the bug...
    
    Glass_logical = new G4LogicalVolume (Glass_solid, Mat_Vessel_Glass, "Glass_log");
	
	//Placing mDOMs
	if (gn_mDOMs <= 1) {
		Glass_physical = new G4PVPlacement (flipOM, G4ThreeVector(0,0,0), Glass_logical, "Glass_phys_0", World_logical, false, 0);
		PlacingHarnessAndRopes(0*m, rot, ropeThickness, ropeLength);
	} else {
		std::stringstream moduleconverter;
		Glass_physical_vector.resize(gn_mDOMs);
		for (k = 0; k < gn_mDOMs; k++){
			moduleconverter.str("");
			moduleconverter << "Glass_phys_" << k ;
			G4double zpos;
			if (gn_mDOMs % 2 == 0) {
				zpos = gmdomseparation*(gn_mDOMs/2-k-1./2.);
			} else {
				zpos = gmdomseparation*(gn_mDOMs/2-k);
			}
			Glass_physical_vector[k] = new G4PVPlacement (flipOM, G4ThreeVector(0,0,zpos), Glass_logical, moduleconverter.str(), World_logical, false, 0);
			PlacingHarnessAndRopes(zpos, rot, ropeThickness, ropeLength);
		}
	}
    
//     G4LogicalVolume* AirCylinder_logical = new G4LogicalVolume (AirCylinder_solid, Mat_LabAir, "Airslice logical");
//     G4PVPlacement* AirCylinder_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), AirCylinder_logical, "Airslice physical", Glass_logical, false, 0);
    
    Gel_logical = new G4LogicalVolume (Gel_solid, Mat_Gel, "Gelcorpus logical");
    Gel_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), Gel_logical, "Gelcorpus physical", Glass_logical, false, 0);
    
    TubeHolder_logical = new G4LogicalVolume (TubeHolder_solid, Mat_Absorber, "TubeHolder logical");
    TubeHolder_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), TubeHolder_logical, "TubeHolder physical", Gel_logical, false, 0);
    
    RefConeType1_logical = new G4LogicalVolume(RefConeType1_solid, Mat_Reflector, "RefConeType1 logical");    
    RefConeType2_logical = new G4LogicalVolume(RefConeType2_solid, Mat_Reflector, "RefConeType2 logical");
    RefConeType3_logical = new G4LogicalVolume(RefConeType3_solid, Mat_Reflector, "RefConeType3 logical");
    
    
    // placing PMTs & RefCones (into gel)
    std::stringstream converter;
    for (k = 0; k <= 23; k++){
      converter.str("");
      converter << "PMT_" << k << "_physical";
      rot = new G4RotationMatrix();
      rot->rotateY(PMT_theta[k]);
      rot->rotateZ(PMT_phi[k]);
      transformers = G4Transform3D(*rot, G4ThreeVector(PMT_x[k],PMT_y[k],PMT_z[k]));
      
      // placing PMTs depending on type:
      if (PMT_type == "12199s"){
	PMT_physical[k] = new G4PVPlacement (transformers, PMT_12199_tube_logical, converter.str(), Gel_logical, false, 0);
      }
      if (PMT_type == "etel"){
	PMT_physical[k] = new G4PVPlacement (transformers, PMT_ETEL_tube_logical, converter.str(), Gel_logical, false, 0);
      }
      if (PMT_type == "12199e"){
	PMT_physical[k] = new G4PVPlacement (transformers, PMT80_tube_logical, converter.str(), Gel_logical, false, 0);
      }
      
      // placing reflective cones:
      // type 1 for all "polar" conees, type 2 for upper central, type 3 for lower central
      converter.str("");
      converter << "RefCone_" << k << "_physical";
      transformers = G4Transform3D(*rot, G4ThreeVector(RefCone_x[k],RefCone_y[k],RefCone_z[k]));
      if (k >= 4 && k <= 11) {
	rot = new G4RotationMatrix();
	rot->rotateZ(PMT_phi[k]);
	transformers = G4Transform3D(*rot, G4ThreeVector(0,0,0));
	if (PMT_type == "etel"){
	  RefCone_physical[k] = new G4PVPlacement (transformers, RefConeType2_ETEL_logical, converter.str(), Gel_logical, false, 0);
	}
	else {
	  RefCone_physical[k] = new G4PVPlacement (transformers, RefConeType2_logical, converter.str(), Gel_logical, false, 0);
	}
      }
      else if (k >= 12 && k <= 19) {
	rot = new G4RotationMatrix();
	rot->rotateZ(PMT_phi[k]);
	transformers = G4Transform3D(*rot, G4ThreeVector(0,0,0));
	if (PMT_type == "etel"){
	  RefCone_physical[k] = new G4PVPlacement (transformers, RefConeType3_ETEL_logical, converter.str(), Gel_logical, false, 0);
	}
	else {
	  RefCone_physical[k] = new G4PVPlacement (transformers, RefConeType3_logical, converter.str(), Gel_logical, false, 0);
	}
      }
      else {
	RefCone_physical[k] = new G4PVPlacement (transformers, RefConeType1_logical, converter.str(), Gel_logical, false, 0);
      }
    }
    
    
    // ------------------- optical border surfaces --------------------------------------------------------------------------------
    G4LogicalSkinSurface* RefConeSurface = new G4LogicalSkinSurface("RefCone_skin", RefConeType1_logical, RefCone_optical);
    RefConeSurface = new G4LogicalSkinSurface("RefCone_skin", RefConeType2_logical, RefCone_optical);
    RefConeSurface = new G4LogicalSkinSurface("RefCone_skin", RefConeType3_logical, RefCone_optical);
    RefConeSurface = new G4LogicalSkinSurface("RefCone_skin", RefConeType2_ETEL_logical, RefCone_optical);
    RefConeSurface = new G4LogicalSkinSurface("RefCone_skin", RefConeType3_ETEL_logical, RefCone_optical);
    RefConeSurface = new G4LogicalSkinSurface("RefCone_skin", RefCone_ETEL_logical, RefCone_optical);
    RefConeSurface = new G4LogicalSkinSurface("RefCone_skin", RefCone_12199_logical, RefCone_optical);
    
    // White Holder
    if (Holder_color == "White") { G4LogicalSkinSurface* HolderSurface = new G4LogicalSkinSurface("Holder_skin", TubeHolder_logical, Holder_optical); }
    
    // ---------------- visualisation attributes --------------------------------------------------------------------------------
    Glass_logical->SetVisAttributes(Glass_vis);
    Gel_logical->SetVisAttributes(Gel_vis);
    TubeHolder_logical->SetVisAttributes(Absorber_vis);
    
    RefConeType1_logical->SetVisAttributes(Alu_vis);
    RefConeType2_logical->SetVisAttributes(Alu_vis);
    RefConeType3_logical->SetVisAttributes(Alu_vis);
    RefConeType2_ETEL_logical->SetVisAttributes(Alu_vis);
    RefConeType3_ETEL_logical->SetVisAttributes(Alu_vis);
    RefCone_ETEL_logical->SetVisAttributes(PhotoCathode_vis);
    RefCone_12199_logical->SetVisAttributes(PhotoCathode_vis);
    
  } // closing mDOM area
  
  if (OM_type == "PDOM") {
    
    G4double PDOM_GelThick = 10*mm;
    G4double PDOM_PMT_z = 0.5*12*25.4*mm-Tube10_geo_curvature_radius-PDOM_GelThick;
    
    PDOM_GlassSphere_solid = new G4Orb("PDOM_GlassSphere solid",0.5*13*25.4*mm);
    PDOM_GelSphere_solid = new G4Orb("PDOM_GelSphere solid",0.5*12*25.4*mm);
    
    PDOM_AirAux_solid = new G4Ellipsoid("PDOM_AirAux solid", 0.5*12*25.4*mm, 0.5*12*25.4*mm, 0.5*12*25.4*mm, -0.5*13*25.4*mm, 50*mm);
    
    PDOM_Air_solid = new G4SubtractionSolid("PDOM_Air solid", PDOM_AirAux_solid, PMT10_tube_solid, 0, G4ThreeVector(0,0, PDOM_PMT_z));
    PDOM_Board_solid = new G4Tubs("PDOM_Board solid", 52*mm, 0.5*11*25.4*mm, 2*mm, 0, 2*pi);
    G4Tubs* PDOM_Base_solid = new G4Tubs("PDOM_Board solid", 0*mm, 6*cm, 2*mm, 0, 2*pi);
    
    // steel harness of PDOM; dimensions by Perry form "gel_and_band.pdf"
    G4double PDOM_Harness_inner[] = {0, 0, 0, 0};
    G4double PDOM_Harness_radii[] = {(0.5*365.76-8.3)*mm, 0.5*365.76*mm, 0.5*365.76*mm, (0.5*365.76-8.3)*mm};
    G4double PDOM_Harness_zplanes[] = {-31.75*mm, -10*mm, 10*mm, 31.75*mm};	
    PDOM_HarnessAux_solid =  new G4Polycone("PDOM_HarnessAux solid", 0, 2*pi, 4, PDOM_Harness_zplanes, PDOM_Harness_inner, PDOM_Harness_radii);
    PDOM_Harness_solid = new G4SubtractionSolid("PDOM_Harness solid", PDOM_HarnessAux_solid, PDOM_GlassSphere_solid);

    
    //  Platzierungen
    
    PDOM_Glass_logical = new G4LogicalVolume (PDOM_GlassSphere_solid, Mat_Vessel_Glass, "PDOM_Glass logical");
    
    PDOM_Gel_logical = new G4LogicalVolume (PDOM_GelSphere_solid, Mat_Gel, "PDOM_Gel logical");
    PDOM_Gel_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), PDOM_Gel_logical, "PDOM_Gel physical", PDOM_Glass_logical, false, 0);
    PDOM_Board_logical = new G4LogicalVolume (PDOM_Board_solid, Mat_Absorber, "PDOM_Board logical");
    G4LogicalVolume* PDOM_Base_logical = new G4LogicalVolume (PDOM_Base_solid, Mat_Absorber, "PDOM_Base logical");
    
    PDOM_Harness_logical = new G4LogicalVolume (PDOM_Harness_solid, Mat_Stahl, "PDOM_Harness logical");
    
    PDOM_Air_logical = new G4LogicalVolume (PDOM_Air_solid, Mat_Vacuum, "PDOM_Air logical");	
    PDOM_Board_physical = new G4PVPlacement (0, G4ThreeVector(0,0,-40*mm), PDOM_Board_logical, "PDOM_Board physical", PDOM_Air_logical, false, 0);
    G4PVPlacement* PDOM_Base_physical = new G4PVPlacement (0, G4ThreeVector(0,0,-105*mm), PDOM_Base_logical, "PDOM_Base physical", PDOM_Air_logical, false, 0);
    PDOM_Air_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), PDOM_Air_logical, "PDOM_Air physical", PDOM_Gel_logical, false, 0);
    
    PMT_physical[0] = new G4PVPlacement (0, G4ThreeVector(0,0, PDOM_PMT_z), PMT10_tube_logical, "PMT_0_physical", PDOM_Gel_logical, false, 0);
    PDOM_Harness_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), PDOM_Harness_logical, "PDOM_Harness physical", World_logical, false, 0);
    
    rot = new G4RotationMatrix();
    // 		rot->rotateY(180*deg); //flip OM downwards here before positioning, corresponds to actual position in ice
    PDOM_Glass_physical = new G4PVPlacement (rot, G4ThreeVector(0,0,0), PDOM_Glass_logical, "Glass_phys", World_logical, false, 0);
    
    // ------------------- optical border surfaces --------------------------------------------------------------------------------
    G4LogicalSkinSurface* PDOM_HarnessSurface = new G4LogicalSkinSurface("PDOM_Harness_skin", PDOM_Harness_logical, PDOM_Harness_optical);
    
    
    // ------------------- visualisation attributes -------------------------------------------------------------------------------
    
    PDOM_Glass_logical->SetVisAttributes(Glass_vis);
    PDOM_Gel_logical->SetVisAttributes(Gel_vis);
    PDOM_Air_logical->SetVisAttributes(World_vis);
    PDOM_Board_logical->SetVisAttributes(Board_vis);
    PDOM_Base_logical->SetVisAttributes(Board_vis);
    PDOM_Harness_logical->SetVisAttributes(stahl_vis);
    
  }
  
  
  
  if (OM_type == "halfmDOM") {
    G4double singlePMT_z = CylHigh + GlasInRad - GelPMT - PMToffset;
    
    // Glass
    G4Ellipsoid* single_GlassSphereTop_solid = new G4Ellipsoid("single_GlassSphereTop solid", GlasOutRad, GlasOutRad, GlasOutRad, -5*mm, GlasOutRad+5*mm);
    G4Tubs* single_GlassCylinder_solid = new G4Tubs("single_GlassCylinder solid", 0, GlasOutRad, CylHigh, 0, 2*pi);
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
    G4UnionSolid* single_Glass_solid = new G4UnionSolid("single_glass_body", single_GlassCylinder_solid, single_GlassSphereTop_solid, transformers);
    //for full module
    //		G4Sphere* single_GlassSphereTop_solid = new G4Sphere("single_GlassSphereTop solid", 0, GlasOutRad, 0, 2*pi, 0, 0.51*pi);
    //		G4Sphere* single_GlassSphereBottom_solid = new G4Sphere("single_GlassSphereBottom solid", 0, GlasOutRad, 0, 2*pi, 0.49*pi, pi);
    //		G4Tubs* single_GlassCylinder_solid = new G4Tubs("single_GlassCylinder solid", 0, GlasOutRad, CylHigh, 0, 2*pi);
    //		transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
    //		G4UnionSolid* single_temp_union = new G4UnionSolid("single_temp", single_GlassCylinder_solid, single_GlassSphereTop_solid, transformers);
    //		transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-CylHigh));
    //		G4UnionSolid* single_Glass_solid = new G4UnionSolid("single_glass_body", single_temp_union, single_GlassSphereBottom_solid, transformers);
    
    //  Gel
    G4Ellipsoid* single_GelSphereTop_solid = new G4Ellipsoid("GelSphereTop solid", GlasInRad, GlasInRad, GlasInRad, -5*mm, GlasInRad+5*mm);
    G4Tubs* single_GelCylinder_solid = new G4Tubs("single_GelCylinder solid", 0, GlasOutRad - GlasThick, CylHigh , 0, 2*pi);
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
    G4UnionSolid* single_Gel_solid = new G4UnionSolid("single_gel_body", single_GelCylinder_solid, single_GelSphereTop_solid, transformers);
    //for full module
    //		G4Sphere* single_GelSphereTop_solid = new G4Sphere("single_GelSphereTop solid", 0, GlasOutRad - GlasThick, 0, 2*pi, 0, 0.51*pi);
    //		G4Sphere* single_GelSphereBottom_solid = new G4Sphere("single_GelSphereBottom solid", 0, GlasOutRad - GlasThick, 0, 2*pi, 0.49*pi, pi);
    //		G4Tubs* single_GelCylinder_solid = new G4Tubs("single_GelCylinder solid", 0, GlasOutRad - GlasThick, CylHigh , 0, 2*pi);
    //		transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
    //		G4UnionSolid* single_temp_union2 = new G4UnionSolid("single_temp2", single_GelCylinder_solid, single_GelSphereTop_solid, transformers);
    //		transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-CylHigh));
    //		G4UnionSolid* single_Gel_solid = new G4UnionSolid("single_gel_body", single_temp_union2, single_GelSphereBottom_solid, transformers);
    
    //  PMT TubeHolder from "foam" primitives & cutting "nests" for PMTs later
    // 		G4Sphere* single_FoamSphereTop_solid = new G4Sphere("single FoamSphereTop solid", 0, GlasOutRad - GlasThick - GelThick, 0, 2*pi, 0, 0.51*pi);
    // 		G4Tubs* single_FoamCylinder_solid = new G4Tubs("single FoamCylinder solid", 0, GlasOutRad - GlasThick - GelThick, CylHigh , 0, 2*pi);
    // 		transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(CylHigh )));
    // 		G4UnionSolid* single_Union_solid = new G4UnionSolid("single Union solid", single_FoamCylinder_solid, single_FoamSphereTop_solid, transformers);
    // 		G4Tubs* single_Cylinder_solid = new G4Tubs("single Cylinder solid", 0, 60*mm, 30*cm , 0, 2*pi);
    // 		transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,0));
    // 		G4IntersectionSolid* single_Foam_solid = new G4IntersectionSolid("single Foam solid",single_Union_solid,single_Cylinder_solid,transformers);
    G4double FoamRad = GlasOutRad - GlasThick - GelThick;
    G4Ellipsoid* FoamSphereTop_solid = new G4Ellipsoid("FoamSphereTop solid", FoamRad, FoamRad, FoamRad, -5*mm, FoamRad+5*mm);
    G4Tubs* single_FoamCylinder_solid = new G4Tubs("single FoamCylinder solid", 0, GlasOutRad - GlasThick, CylHigh , 0, 2*pi);
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(CylHigh )));
    G4UnionSolid* single_Foam_solid = new G4UnionSolid("single Foam solid", single_FoamCylinder_solid, FoamSphereTop_solid, transformers);
    //for full module
    //		G4Sphere* single_FoamSphereTop_solid = new G4Sphere("single FoamSphereTop solid", 0, GlasOutRad - GlasThick, 0, 2*pi, 0, 0.51*pi);
    //		G4Sphere* single_FoamSphereBottom_solid = new G4Sphere("single FoamSphereBottom solid", 0, GlasOutRad - GlasThick, 0, 2*pi, 0.49*pi, pi);
    //		G4Tubs* single_FoamCylinder_solid = new G4Tubs("single FoamCylinder solid", 0, GlasOutRad - GlasThick, CylHigh , 0, 2*pi);
    //		transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(CylHigh )));
    //		G4UnionSolid* single_Foam_TempUnion_solid = new G4UnionSolid("single Foam TempUnion solid", single_FoamCylinder_solid, single_FoamSphereTop_solid, transformers);
    //		transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-(CylHigh )));
    //		G4UnionSolid* single_Foam_solid = new G4UnionSolid("single Foam solid", single_Foam_TempUnion_solid, single_FoamSphereBottom_solid, transformers);
    
    RefConeBasic_solid = new G4Cons("RefConeBasic", 
				    RefCone_IdealInRad, 
				    RefCone_IdealInRad + RefCone_SheetThickness*std::sqrt(2.), 
				    RefCone_IdealInRad + 2*RefConeDZ, 
				    RefCone_IdealInRad + RefCone_SheetThickness*std::sqrt(2.) + 2*RefConeDZ, 
				    RefConeDZ, 0, 2*pi);
    rot = new G4RotationMatrix();
    transformers = G4Transform3D(*rot, G4ThreeVector(0,0,-(GlasInRad-GelPMT-PMToffset+RefConeDZ)));
    RefConeType1_solid = new G4IntersectionSolid("RefConeType1", RefConeBasic_solid, FoamSphereTop_solid, transformers);
    RefConeType1_logical = new G4LogicalVolume(RefConeType1_solid, Mat_Reflector, "RefConeType1 logical");   
    
    G4LogicalVolume* single_Glass_logical = new G4LogicalVolume (single_Glass_solid, Mat_Vessel_Glass, "single_Glasscorpus logical");
    G4LogicalVolume* single_Gel_logical = new G4LogicalVolume (single_Gel_solid, Mat_Gel, "single_Gelcorpus logical");
    
    rot = new G4RotationMatrix();
    transformers = G4Transform3D(*rot, G4ThreeVector(0,0,singlePMT_z));
    G4SubtractionSolid* single_TubeHolder_solid = new G4SubtractionSolid("single TubeHolder solid", single_Foam_solid, RefConeNest_solid, transformers);
    //		G4LogicalVolume* single_TubeHolder_logical = new G4LogicalVolume (single_TubeHolder_solid, Mat_Absorber, "single TubeHolder logical");
    G4LogicalVolume* single_TubeHolder_logical = new G4LogicalVolume (single_TubeHolder_solid, Mat_LabAir, "single TubeHolder logical");
    
    // 		PMT_physical[0] = new G4PVPlacement (0, G4ThreeVector(0,0,singlePMT_z), PMT_12199_tube_logical, "PMT_0_physical", single_Gel_logical, false, 0);
    PMT_physical[0] = new G4PVPlacement (0, G4ThreeVector(0,0,singlePMT_z), PMT80_tube_logical, "PMT_0_physical", single_Gel_logical, false, 0);
    RefCone_physical[0] = new G4PVPlacement (0, G4ThreeVector(0,0,singlePMT_z + RefConeDZ), RefConeType1_logical, "RefCone_0_physical", single_Gel_logical, false, 0);
    
    TubeHolder_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), single_TubeHolder_logical, "TubeHolder physical", single_Gel_logical, false, 0);
    Gel_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), single_Gel_logical, "Gelcorpus physical", single_Glass_logical, false, 0);
    Glass_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), single_Glass_logical, "Glass_phys", World_logical, false, 0);
    G4LogicalSkinSurface* RefConeSurface = new G4LogicalSkinSurface("RefCone_skin", RefConeType1_logical, RefCone_optical);
    
    // ------------------- visualisation attributes -------------------------------------------------------------------------------
    
    single_Glass_logical->SetVisAttributes(Glass_vis);
    single_Gel_logical->SetVisAttributes(Gel_vis);
    // 		single_TubeHolder_logical->SetVisAttributes(Absorber_vis);
    single_TubeHolder_logical->SetVisAttributes(World_vis);//
    //		single_Air_logical->SetVisAttributes(World_vis);
  }
  
  
  
  
  
  
  
  if (OM_type == "halfmDOM2") {
    G4double singlePMT_z = CylHigh + GlasInRad - GelPMT - PMToffset;
    
    // Glass
    G4Ellipsoid* single_GlassSphereTop_solid = new G4Ellipsoid("single_GlassSphereTop solid", GlasOutRad, GlasOutRad, GlasOutRad, -5*mm, GlasOutRad+5*mm);
    G4Tubs* single_GlassCylinder_solid = new G4Tubs("single_GlassCylinder solid", 0, GlasOutRad, CylHigh, 0, 2*pi);
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
    G4UnionSolid* single_Glass_solid = new G4UnionSolid("single_glass_body", single_GlassCylinder_solid, single_GlassSphereTop_solid, transformers);
    G4LogicalVolume* single_Glass_logical = new G4LogicalVolume (single_Glass_solid, Mat_Vessel_Glass, "single_Glasscorpus logical");
    
    
    // two half-vessels:
    //G4LogicalVolume* single_Glass_logical = new G4LogicalVolume (Glass_solid, Mat_Vessel_Glass, "single_Glasscorpus logical");
    
    
    
    //  Gel
    G4Ellipsoid* single_Gel_solid1 = new G4Ellipsoid("GelSphereTop solid", GlasInRad, GlasInRad, GlasInRad, 120*mm, GlasInRad+5*mm);
    //G4Tubs* single_GelCylinder_solid = new G4Tubs("single_GelCylinder solid", 0, GlasOutRad - GlasThick, CylHigh , 0, 2*pi);
    //transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
    //G4UnionSolid* single_Gel_solid = new G4UnionSolid("single_gel_body", single_GelCylinder_solid, single_GelSphereTop_solid, transformers);
    
    
    
    //  Gel
    G4Ellipsoid* single_GelSphereTop_solid2 = new G4Ellipsoid("GelSphereTop solid", GlasInRad, GlasInRad, GlasInRad, -5*mm, GlasInRad+5*mm);
    GelSphereBottom_solid = new G4Ellipsoid("GelSphereBottom solid", GlasInRad, GlasInRad, GlasInRad, -(GlasInRad+5*mm), 5*mm);
    G4Tubs* single_GelCylinder_solid2 = new G4Tubs("single_GelCylinder solid", 0, GlasOutRad - GlasThick-0.01*mm, CylHigh , 0, 2*pi);
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
    G4UnionSolid* single_Gel_solid2 = new G4UnionSolid("single_gel_body", single_GelCylinder_solid2, single_GelSphereTop_solid2, transformers);
    
    // two half-vessels:
    // transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-CylHigh));
    //single_Gel_solid2  = new G4UnionSolid("gel body", single_Gel_solid2 , GelSphereBottom_solid, transformers);
    
    
    
    
    
    G4double FoamRad = GlasOutRad - GlasThick - GelThick;
    G4Ellipsoid* FoamSphereTop_solid  = new G4Ellipsoid("FoamSphereTop solid", GlasInRad, GlasInRad, GlasInRad, -5*mm, GlasInRad+5*mm);
    G4Tubs* single_FoamCylinder_solid = new G4Tubs("single FoamCylinder solid", 0, GlasOutRad - GlasThick, CylHigh , 0, 2*pi);
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(CylHigh )));
    G4UnionSolid* single_Foam_solid = new G4UnionSolid("single Foam solid", single_FoamCylinder_solid, FoamSphereTop_solid, transformers);
    //----------------two vessels:
    // G4Ellipsoid* FoamSphereBottom_solid = new G4Ellipsoid("FoamSphereBottom solid", GlasInRad, GlasInRad, GlasInRad, -(GlasInRad+5*mm), 5*mm);		
    // transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-(CylHigh )));
    // single_Foam_solid = new G4UnionSolid("Foam solid", single_Foam_solid, FoamSphereBottom_solid, transformers);
    //-----------------------------
    
    // REF CONE------------------------------
    RefConeBasic_solid = new G4Cons("RefConeBasic", 
				    RefCone_IdealInRad, 
				    RefCone_IdealInRad + RefCone_SheetThickness*std::sqrt(2.), 
				    RefCone_IdealInRad + 2*RefConeDZ, 
				    RefCone_IdealInRad + RefCone_SheetThickness*std::sqrt(2.) + 2*RefConeDZ, 
				    RefConeDZ, 0, 2*pi);
    rot = new G4RotationMatrix();
    transformers = G4Transform3D(*rot, G4ThreeVector(0,0,-(GlasInRad-GelPMT-PMToffset+RefConeDZ)));
    RefConeType1_solid = new G4IntersectionSolid("RefConeType1", RefConeBasic_solid, FoamSphereTop_solid,transformers);
    
    // ABS CONE------------------------------
    G4VSolid* RefConeBasic_solid22 = new G4Cons("RefConeBasic2", 
						0.5*51.9*mm, 
						0.5*51.9*mm+RefCone_SheetThickness*std::sqrt(2.), 
						RefCone_IdealInRad + 1.7*RefConeDZ+1*RefCone_SheetThickness*std::sqrt(2.), 
						RefCone_IdealInRad + 2*RefCone_SheetThickness*std::sqrt(2.) + 1.7*RefConeDZ, 
						RefConeDZ*2.28, 0, 2*pi);
    transformers = G4Transform3D(*rot, G4ThreeVector(0,0,-GlasInRad+RefConeDZ*2));
    G4VSolid* RefConeBasic_solid2 = new G4IntersectionSolid("RefConeType2", RefConeBasic_solid22, single_Gel_solid1, transformers);
    
    
    
    
    
    //TubeHolder booleans substraction	
    //	rot = new G4RotationMatrix();
    //	transformers = G4Transform3D(*rot, G4ThreeVector(0,0,singlePMT_z));
    //	G4SubtractionSolid* single_TubeHolder_solid2 = new G4SubtractionSolid("single TubeHolder2 solid", single_Foam_solid, RefConeNest_solid, transformers);
    
    
    
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh+0.5*mm));
    G4SubtractionSolid* single_TubeHolder_solid3 = new G4SubtractionSolid("single TubeHolder3 solid",single_Foam_solid, single_Gel_solid1, transformers);
    
    
    //	transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(GlasInRad-GelPMT-PMToffset+RefConeDZ)));
    //	G4SubtractionSolid* single_TubeHolder_solid4 = new G4SubtractionSolid("single TubeHolder4 solid",single_TubeHolder_solid3, RefConeBasic_solid2, transformers);
    
    
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(singlePMT_z)));
    G4SubtractionSolid* single_TubeHolder_solid = new G4SubtractionSolid("single TubeHolder solid",single_TubeHolder_solid3, PMT_12199_tube_solid, transformers);
    
    //Gel booleans substraction
    
    //	transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(-CylHigh+GlasInRad-GelPMT-PMToffset+RefConeDZ)));
    //	G4SubtractionSolid* single_Gel_solid3 = new G4SubtractionSolid("Gel substraction 1",single_Gel_solid1, RefConeBasic_solid2, transformers);
    
    //	transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(singlePMT_z-CylHigh)));
    //	G4SubtractionSolid* single_Gel_solid4 = new G4SubtractionSolid("Gel substraction 2",single_Gel_solid1, PMT80_tube_solid, transformers);
    
    //	transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(singlePMT_z-CylHigh)));
    //	G4SubtractionSolid* single_Gel_solid = new G4SubtractionSolid("Gel substraction 3",single_Gel_solid4, RefConeNest_solid, transformers);
    
    
    
    
    
    
    G4LogicalVolume* RefConeType1_logical = new G4LogicalVolume(RefConeType1_solid, Mat_Reflector, "RefConeType1 logical");
    G4LogicalVolume* RefConeType2_logical = new G4LogicalVolume(RefConeBasic_solid2, Mat_Absorber, "RefConeType2 logical"); 
    G4LogicalVolume* single_Gel_logical = new G4LogicalVolume (single_Gel_solid2, Mat_LabAir, "single_Gelcorpus logical");
    G4LogicalVolume* single_TubeHolder_logical = new G4LogicalVolume (single_TubeHolder_solid,  Mat_LabAir, "single TubeHolder logical");
    
    // 		PMT_physical[0] = new G4PVPlacement (0, G4ThreeVector(0,0,singlePMT_z), PMT_12199_tube_logical, "PMT_0_physical", single_Gel_logical, false, 0);
    PMT_physical[0] = new G4PVPlacement (0, G4ThreeVector(0,0,singlePMT_z), PMT_12199_tube_logical, "PMT_0_physical", single_Gel_logical, true, 0);
    //RefCone_physical[0] = new G4PVPlacement (0, G4ThreeVector(0,0,singlePMT_z + RefConeDZ), RefConeType1_logical, "RefCone_0_physical", single_Gel_logical, true, 0);
    G4VPhysicalVolume* RefCone2_physical = new G4PVPlacement (0, G4ThreeVector(0,0,GlasInRad-RefConeDZ*2+CylHigh), RefConeType2_logical, "RefCone_2_physical", single_Gel_logical, true, 0);
    
    TubeHolder_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), single_TubeHolder_logical, "TubeHolder physical", single_Gel_logical, true, 0);
    Gel_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), single_Gel_logical, "Gelcorpus physical", single_Glass_logical, true, 0);
	
    Glass_physical = new G4PVPlacement (0, G4ThreeVector(0,0,0), single_Glass_logical, "Glass_phys", World_logical, true, 0);
	
    G4LogicalSkinSurface* RefConeSurface = new G4LogicalSkinSurface("RefCone_skin", RefConeType1_logical, RefCone_optical);
    
    // ------------------- visualisation attributes -------------------------------------------------------------------------------
    
    single_Glass_logical->SetVisAttributes(Glass_vis);
    single_Gel_logical->SetVisAttributes(World_vis);
    RefConeType2_logical->SetVisAttributes(Absorber_vis);
    single_TubeHolder_logical->SetVisAttributes(World_vis);
    //		single_Air_logical->SetVisAttributes(World_vis);
  }
  
  
  if (OM_type == "PMTwithSample") {
    G4double singlePMT_z = CylHigh + GlasInRad - GelPMT - PMToffset;

    G4Ellipsoid* single_GlassSphereTop_solid = new G4Ellipsoid("single_GlassSphereTop solid", GlasOutRad, GlasOutRad, GlasOutRad, -5*mm, GlasOutRad+5*mm);
    G4Tubs* single_GlassCylinder_solid = new G4Tubs("single_GlassCylinder solid", 0, GlasOutRad, CylHigh, 0, 2*pi);
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
    G4UnionSolid* single_Glass_solid = new G4UnionSolid("single_glass_body", single_GlassCylinder_solid, single_GlassSphereTop_solid, transformers);
    G4LogicalVolume* single_Glass_logical = new G4LogicalVolume (single_Glass_solid, Mat_LabAir, "single_Glasscorpus logical");

    G4Ellipsoid* single_Gel_solid1 = new G4Ellipsoid("GelSphereTop solid", GlasInRad, GlasInRad, GlasInRad, 120*mm, GlasInRad+5*mm);

    G4Ellipsoid* single_GelSphereTop_solid2 = new G4Ellipsoid("GelSphereTop solid", GlasInRad, GlasInRad, GlasInRad, -5*mm, GlasInRad+5*mm);
    GelSphereBottom_solid = new G4Ellipsoid("GelSphereBottom solid", GlasInRad, GlasInRad, GlasInRad, -(GlasInRad+5*mm), 5*mm);
    G4Tubs* single_GelCylinder_solid2 = new G4Tubs("single_GelCylinder solid", 0, GlasOutRad - GlasThick-0.01*mm, CylHigh , 0, 2*pi);
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh ));
    G4UnionSolid* single_Gel_solid2 = new G4UnionSolid("single_gel_body", single_GelCylinder_solid2, single_GelSphereTop_solid2, transformers);

    
    
    
    
    G4double FoamRad = GlasOutRad - GlasThick - GelThick;
    G4Ellipsoid* FoamSphereTop_solid  = new G4Ellipsoid("FoamSphereTop solid", GlasInRad, GlasInRad, GlasInRad, -5*mm, GlasInRad+5*mm);
    G4Tubs* single_FoamCylinder_solid = new G4Tubs("single FoamCylinder solid", 0, GlasOutRad - GlasThick, CylHigh , 0, 2*pi);
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(CylHigh )));
    G4UnionSolid* single_Foam_solid = new G4UnionSolid("single Foam solid", single_FoamCylinder_solid, FoamSphereTop_solid, transformers);

    RefConeBasic_solid = new G4Cons("RefConeBasic", 
				    RefCone_IdealInRad, 
				    RefCone_IdealInRad + RefCone_SheetThickness*std::sqrt(2.), 
				    RefCone_IdealInRad + 2*RefConeDZ, 
				    RefCone_IdealInRad + RefCone_SheetThickness*std::sqrt(2.) + 2*RefConeDZ, 
				    RefConeDZ, 0, 2*pi);
    rot = new G4RotationMatrix();
    transformers = G4Transform3D(*rot, G4ThreeVector(0,0,-(GlasInRad-GelPMT-PMToffset+RefConeDZ)));
    RefConeType1_solid = new G4IntersectionSolid("RefConeType1", RefConeBasic_solid, FoamSphereTop_solid,transformers);

    G4VSolid* RefConeBasic_solid22 = new G4Cons("RefConeBasic2", 
						0.5*51.9*mm, 
						0.5*51.9*mm+RefCone_SheetThickness*std::sqrt(2.), 
						RefCone_IdealInRad + 1.7*RefConeDZ+1*RefCone_SheetThickness*std::sqrt(2.), 
						RefCone_IdealInRad + 2*RefCone_SheetThickness*std::sqrt(2.) + 1.7*RefConeDZ, 
						RefConeDZ*2.28, 0, 2*pi);
    transformers = G4Transform3D(*rot, G4ThreeVector(0,0,-GlasInRad+RefConeDZ*2));
    G4VSolid* RefConeBasic_solid2 = new G4IntersectionSolid("RefConeType2", RefConeBasic_solid22, single_Gel_solid1, transformers);

    
    
    
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,CylHigh+0.5*mm));
    G4SubtractionSolid* single_TubeHolder_solid3 = new G4SubtractionSolid("single TubeHolder3 solid",single_Foam_solid, single_Gel_solid1, transformers);
    
    

    
    transformers = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,(singlePMT_z)));
    G4SubtractionSolid* single_TubeHolder_solid = new G4SubtractionSolid("single TubeHolder solid",single_TubeHolder_solid3, PMT_12199_tube_solid, transformers);
    

    
    
    G4LogicalVolume* RefConeType1_logical = new G4LogicalVolume(RefConeType1_solid, Mat_Reflector, "RefConeType1 logical");
    G4LogicalVolume* RefConeType2_logical = new G4LogicalVolume(RefConeBasic_solid2, Mat_Absorber, "RefConeType2 logical"); 
    G4LogicalVolume* single_Gel_logical = new G4LogicalVolume (single_Gel_solid2, Mat_LabAir, "single_Gelcorpus logical");
    G4LogicalVolume* single_TubeHolder_logical = new G4LogicalVolume (single_TubeHolder_solid,  Mat_LabAir, "single TubeHolder logical");
    
    G4Tubs* quelle_solid = new G4Tubs("Quelle solid", 0, 25*mm/2., 0.5*mm, 0, 2*pi);
    G4LogicalVolume* quelle_logical = new G4LogicalVolume(quelle_solid, MatDatBase->FindOrBuildMaterial("G4_Fe"), "Quelle logical");
    
    
    G4Tubs* ring_solid = new G4Tubs("Ring solid", 22*mm/2., 26*mm/2., 0.5*mm, 0, 2*pi);
    G4LogicalVolume* ring_logical = new G4LogicalVolume(ring_solid, MatDatBase->FindOrBuildMaterial("G4_Fe"), "Ring logical");
    
    
    G4Box* holder1 = new G4Box("Holder 1", 8*cm/2., 16*mm/2., 16*mm/2.);
    G4Box* holder2 = new G4Box("Holder 2", 8*cm/2., 16*mm/2., 16*mm/2.);
    G4LogicalVolume* holder1_logical = new G4LogicalVolume(holder1, Mat_Absorber, "Holder1 logical");
    G4LogicalVolume* holder2_logical = new G4LogicalVolume(holder2, Mat_Absorber, "Holder2 logical");

    
    G4double zdist = 0.968/2.*cm;
    G4double distToPMT = 7.3*cm-8*mm;
    G4double PMTheight = (48.)*0.5*mm;
    G4Box* mySample = new G4Box("Glass_phys", 1.5044*cm,1.5044*cm, zdist);
    G4LogicalVolume* mySampleLog = new G4LogicalVolume (mySample, Mat_Vessel_Glass, "single_Glasscorpus logical");
    
    
    G4Orb* myWorld = new G4Orb("myAirWorldsolid",gworldsize*m);
    G4LogicalVolume* myWorldLog = new G4LogicalVolume (myWorld, Mat_LabAir, "myairworld");
    G4PVPlacement* myworldphy = new G4PVPlacement (0, G4ThreeVector(0,0,0), myWorldLog, "myworld", World_logical, true, 0);
    
    //G4VPhysicalVolume* RefCone2_physical = new G4PVPlacement (0, G4ThreeVector(0,0,RefConeDZ+0*CylHigh), RefConeType2_logical, "RefCone_2_physical", myWorldLog, true, 0);
    PMT_physical[0] = new G4PVPlacement (0, G4ThreeVector(0,0,0), PMT_12199_tube_logical, "PMT_0_physical", myWorldLog, true, 0);
    
    //Glass_physical = new G4PVPlacement (0, G4ThreeVector(0,0,distToPMT+8*mm-zdist+PMTheight), mySampleLog, "Glass_phys", myWorldLog, true, 0);
    G4PVPlacement* quelle_physical = new G4PVPlacement (0, G4ThreeVector(0,0,distToPMT+8*mm+0.5*mm+1*mm+PMTheight), quelle_logical, "Quelle_phys", myWorldLog, true, 0);
    G4PVPlacement* ring_physical = new G4PVPlacement (0, G4ThreeVector(0,0,distToPMT+8*mm+0.5*mm+PMTheight), ring_logical, "ring_phys", myWorldLog, true, 0);
    
    G4cout << distToPMT+8*mm+0.5*mm+1*mm+PMTheight << G4endl;
    
    //G4PVPlacement* holder1_physical = new G4PVPlacement (0, G4ThreeVector(0,-1.5044*cm-8*mm,distToPMT+PMTheight), holder1_logical, "holder1_phys", myWorldLog, true, 0);
    //G4PVPlacement* holder2_physical = new G4PVPlacement (0, G4ThreeVector(0,1.5044*cm+8*mm,distToPMT+PMTheight), holder2_logical, "holder2_phys", myWorldLog, true, 0);
    
    
    // ------------------- visualisation attributes -------------------------------------------------------------------------------
    
    mySampleLog->SetVisAttributes(Glass_vis);
    myWorldLog->SetVisAttributes(World_vis);
    single_Gel_logical->SetVisAttributes(World_vis);
    holder1_logical->SetVisAttributes(Absorber_vis);
    holder2_logical->SetVisAttributes(Absorber_vis);
    single_TubeHolder_logical->SetVisAttributes(World_vis);
    //		single_Air_logical->SetVisAttributes(World_vis);
  }
  
  
  
  
  
  
  
  // ------------ visualisation attributes ---------------------------
  World_logical->SetVisAttributes(Gel_vis);
  
  PMT80_tube_logical->SetVisAttributes(Glass_vis);
  PMT_12199_tube_logical->SetVisAttributes(Glass_vis);
  PMT_ETEL_tube_logical->SetVisAttributes(Glass_vis);
  PMT10_tube_logical->SetVisAttributes(Glass_vis);
  
  PC80_logical->SetVisAttributes(PhotoCathode_vis);
  PC_12199_logical->SetVisAttributes(PhotoCathode_vis);
  PC_ETEL_logical->SetVisAttributes(PhotoCathode_vis);
  PC10_logical->SetVisAttributes(PhotoCathode_vis);
  
  PC80_shield_logical->SetVisAttributes(Alu_vis);
  PC_12199_shield_logical->SetVisAttributes(Alu_vis);
  PC_ETEL_shield_logical->SetVisAttributes(Alu_vis);
  PC10_shield_logical->SetVisAttributes(Alu_vis);

  
  World_physical = new G4PVPlacement (0, G4ThreeVector(0.,0.,0.), World_logical, "World_phys", 0, false, 0);

  aNavigator->SetWorldVolume(World_physical);
  
  return World_physical;
}


void mdomDetectorConstruction::PlacingHarnessAndRopes(G4double zpos, G4RotationMatrix* rot, G4double ropeThickness, G4double ropeLength){
	if (gmdomharness) {
        PDOM_Harness_physical = new G4PVPlacement (0, G4ThreeVector(0*m,0*m,zpos), PDOM_Harness_logical, "PDOM_Harness physical", World_logical, false, 0);
	}
		
	if (gropes) {
		G4int user_NbOfRopes = 4;                                                                 //Number of ropes
		G4double ropeDistance = 196.4 * mm - ropeThickness / 2.;                                   //Distance of the center of the rope from the center of the mDOM
		G4double ropeAngle = atan2(178, 1000);                                                    //Tilting angle of the rope with respect to the mDOM axis
		G4double radialShift = ropeLength * sin(ropeAngle);                                       //Translation of the center of the rope part due to tilting
		G4double zShift = ropeLength * cos(ropeAngle) - ropeThickness * sin(ropeAngle);           //Z-position of the center of the rope incl. tilting correction and rope thickness correction
		G4double gHicetube_pos_x = 0;
		G4double totalDOMsize = 0;
		G4double gHicetube_pos_y= 0;
		//Now the rope gets placed uniformly distributed around the harness (if 3 ropes then with an angular distance of 120°)
		//Therefore the rope solid gets rotated around the z and y axes and a transformation vector is filled with the stuff just defined above
		//After that the rope solid gets tucked onto the harness -> RopeHarnessUnion
		//Then the same is done with the other part of the rope cuz only one half of the rope was dealt with yet
		G4VPhysicalVolume* Rope_physical[8];
		//Rope_physical = new G4VPhysicalVolume*[2 * user_NbOfRopes]; //Physical for the logical
		G4ThreeVector TransformationVector;
		std::stringstream converter;
		for (int it = 0; it < user_NbOfRopes; it++){
			rot = new G4RotationMatrix();       //Upper rope
			rot->rotateZ(360. / user_NbOfRopes * (it + 1) * deg);
			rot->rotateY(ropeAngle);
			TransformationVector = G4ThreeVector((ropeDistance - radialShift) * cos(360. / user_NbOfRopes * (it + 1) * deg) + gHicetube_pos_x * totalDOMsize, -(ropeDistance - radialShift) * sin(360. / user_NbOfRopes * (it + 1) * deg) + gHicetube_pos_y * totalDOMsize, zShift) + G4ThreeVector(0,0,zpos);
			converter.str("");
			converter << "Rope_solid_" << it;
			Rope_physical[it] = new G4PVPlacement(rot, TransformationVector, Rope_logical, converter.str(), World_logical, false, 0);

			rot = new G4RotationMatrix();       //Lower rope
			rot->rotateZ(360. / user_NbOfRopes * (it + 1) * deg);
			rot->rotateY(-ropeAngle);
			TransformationVector = G4ThreeVector((ropeDistance - radialShift) * cos(360. / user_NbOfRopes * (it + 1) * deg) + gHicetube_pos_x * totalDOMsize, -(ropeDistance - radialShift) * sin(360. / user_NbOfRopes * (it + 1) * deg) + gHicetube_pos_y * totalDOMsize, -zShift) + G4ThreeVector(0,0,zpos);
			converter.str("");
			converter << "Rope_solid_" << user_NbOfRopes + it;
			
			Rope_physical[user_NbOfRopes + it] = new G4PVPlacement(rot, TransformationVector, Rope_logical, converter.str(), World_logical, false, 0);
		} 
	}
}
