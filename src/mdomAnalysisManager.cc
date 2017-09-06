#include "mdomAnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "mdomMinimization.hh"
#include "G4ios.hh"

extern MdomMinimization Minimization;
extern G4double gposX;
extern G4double gposY;
extern G4double gposZ;
extern G4double gtheta;
extern G4double gphi;
extern G4double 	gtMin, 	gfMin;
extern G4double gRadius;
extern G4double gHeight;
extern G4bool gQE;
extern G4bool gQEweigh;
extern G4int 	gSNGun;


MdomAnalysisManager::MdomAnalysisManager(){
  }

MdomAnalysisManager::~MdomAnalysisManager(){
}

void MdomAnalysisManager::ResetEvent()
{
	foundPhoton = false;
	photonTheta = -99.;
	photonPhi = -99.;
	photonR = -99.;
	evtStat.nrHitPMTs = 0;
	evtStat.nrHitTot = 0;
	evtStat.hitsPMTs.clear();
	hitStats.clear();
}

void MdomAnalysisManager::AnalyzeEvent()
{
// 	EvtStat evtStat;
// 	evtStat.nrHitPMTs = 0;
// 	evtStat.nrHitTot = (G4int)hitStats.size();
	evtStat.nrHitTot = (G4int)hitStats.size();

	for  ( int i=0; i<24; i++ ) {
		std::pair<G4int,G4int> hitsPMT;
		hitsPMT.first = i;
		hitsPMT.second = 0;
		G4bool hit = false;
		for ( int j=0; j<(G4int)hitStats.size(); j++ ) {
			if ( hitStats[j].pmtNr == i ) {
				hitsPMT.second += 1;
				if ( ! hit ) {
//					evtStat.nrHitPMTs += 1;
					evtStat.nrHitPMTs += 1;
					hit = true;
				}
			}
		}
//		if ( hitsPMT.second != 0 ) evtStat.hitsPMTs.push_back(hitsPMT);
		if ( hitsPMT.second != 0 ) evtStat.hitsPMTs.push_back(hitsPMT);
	}

//	evtStats.push_back(evtStat);
	
	// Write data to Minimization function
	// Minimization.BuildAcumulateHits(evtStat);
	//  Write data to file
	Write();
}

void MdomAnalysisManager::WriteHeader()
{
G4String part;
G4String nu;
  datafile << "#Generated in a cylinder of r = " << gRadius/m<< "m and total height of "<< 2*gHeight/m << " m"<< G4endl;
  datafile << "#"<< G4endl;
  if (gSNGun == 1){
	datafile << "# Neutrino electron elastic scattering from SN at 10kpc" << G4endl;
	part = "e-";
	nu = "nu";
  } else if (gSNGun == 2){
	datafile << "# Inverse beta decay from SN at 10kpc" << G4endl;
	part = "e+";
	nu = "nubar";
  } else if (gSNGun == 3){
	datafile << "# Solar neutrinos! Weigh take into account the total cross section but not the flux!" << G4endl;
	part = "e-";
	nu = "nu";
  }
  if (gSNGun == 1 || gSNGun == 2  || gSNGun == 3) {
  datafile << "#"<< G4endl;
  datafile << "# Time of Flux | Mean energy of "<<nu<<" | "<<nu<<" energy | costheta of "<<part<<" from z dir | "<<part<<" energy | event weigh | Vertex Position (X, Y, Z) [m] | Primary direction (Px,Py,Pz)  |Total hits | PMTs hit |";
  					if (gQEweigh) {
						datafile << " Total QE prob |";
					}
					datafile <<"...for PMT hit...| PMT number | Hits in that PMT |";
					datafile << "...for Hit...| photon wavelengh |";
					if (gQEweigh) {
						datafile << " QE prob |";
					}
					datafile << " hit time |";
					if (gQEweigh) {
						datafile << "...end of hit loop...| Total QE prob in PMT |" << G4endl;
					}
  datafile << "#"<< G4endl;
  }

}

void MdomAnalysisManager::Write()
{
  
  //To test the generator: 
  /*
  	  if (gSNGun == 3) {	      
		datafile << nuEnergy/MeV << "\t";
		datafile << primaryEnergy/MeV << "\t";
		datafile << primaryTheta << "\t";
		datafile << weigh << "\t\t";
		datafile << G4endl;
	      
	  } */
 // To get real results
  
	if (evtStat.nrHitTot != 0 ) {
	  if (gSNGun == 1 || gSNGun == 2 || gSNGun == 3)  {
		datafile << nuTime << "\t";
		datafile << nuMeanEnergy/MeV<< "\t";
		datafile << nuEnergy/MeV<< "\t";
		datafile << cosTheta<< "\t";
		datafile << primaryEnergy/MeV << "\t";
		datafile << weigh << "\t\t";
	  }
		datafile << primaryX/m <<"\t";
		datafile << primaryY/m << "\t";
		datafile << primaryZ/m << "\t\t";
		datafile << primaryDirX <<"\t";
		datafile << primaryDirY << "\t";
		datafile << primaryDirZ << "\t\t";
		datafile << evtStat.nrHitTot << "\t";
		datafile << evtStat.nrHitPMTs << "\t";
		if (gQEweigh) {
			G4double sum = 0;
			for (unsigned int i = 0 ; i<hitStats.size(); i++) {
				sum = sum + hitStats[i].QEprob;
			}
			datafile << sum << "\t";
		}
		datafile << "\t";
		for ( int j=0; j<(G4int)evtStat.hitsPMTs.size(); j++ ) {
			datafile << evtStat.hitsPMTs[j].first << "\t";
			datafile << evtStat.hitsPMTs[j].second << "\t\t";
			G4double PMTprob= 0;
			for (int i=0; i<(G4int)hitStats.size(); i++) {
				if ( evtStat.hitsPMTs[j].first == hitStats[i].pmtNr ) {
					datafile << hitStats[i].wavelen/nm << "\t";
					if (gQEweigh) {
						PMTprob = PMTprob + hitStats[i].QEprob;
						datafile << hitStats[i].QEprob << "\t";
						}
					datafile << hitStats[i].hit_time/ns << "\t";
					}
				}
			if (gQEweigh) {
				datafile << PMTprob << "\t\t";
				}
			}
		datafile << G4endl;
		} 
		
}
