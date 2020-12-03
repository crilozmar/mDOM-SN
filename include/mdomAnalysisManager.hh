#ifndef mdomAnalysisManager_h
#define mdomAnalysisManager_h 1

#include "G4Types.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"

#include <vector>
#include <fstream>
#include <tuple>

struct HitStat {
    HitStat() {};
    ~HitStat() {};
    G4int mothercode; // 0=GPS, 1=SNGun result, 2=gamma of 2 MeV, 3=gamma of 8 MeV
	G4int moduleNr;
	G4int pmtNr;
	G4double theta;
	G4double phi;
	G4double hit_time;
	G4double wavelen;
	G4ThreeVector position;
	G4double QEprob;
        G4double reltheta; //between LED dir and photon dir
};

struct EvtStat {
    EvtStat() {};
    ~EvtStat() {};
  std::vector<G4int> modulescounter; // helper parameter to get how many modules where hit
  G4int mothercode; // 0=GPS, 1=SNGun result, 2=gamma of 2 MeV, 3=gamma of 8 MeV
  G4int nrHitTot;// Total number of hits
  G4int nrHitPMTs;// Number of PMTs hit
  G4int nrHitMod; //Number of modules hit
  std::vector<std::tuple<G4int,G4int,G4int> > hitsPMTs;// Module Nr | PMT number | Hits in this PMT
};

struct CollectiveStat {
    std::vector<G4int> pmtHit; //pmt that detected this photon
    std::vector<G4int> moduleHit; //same for module
};

struct FamilyTrack {
    FamilyTrack() {};
    ~FamilyTrack() {};
    G4String motherparticle;
    G4double energy;
    G4int grandparentID; // parent ID of all tracks here
    std::vector<G4int> parentstracks;
    std::vector<G4int> photonstracks;
    // separar tracks de posibles padres y de optical photons
    
};

struct Test_HitStats {
    Test_HitStats() {};
    ~Test_HitStats() {};
    G4double globtheta; //global theta angle
    G4double reltheta; //between LED dir and photon dir
};


class MdomAnalysisManager
{
	public:
		MdomAnalysisManager();
		~MdomAnalysisManager();
		void ResetEvent();
		void AnalyzeEvent();
        void Helper_AnalyzeEvent(EvtStat& this_evtStat);
		void WriteHeader();
        void HelpTheHeader(std::fstream& thisfile);
        void Writer_InfoFile();
        void Writer_data(std::fstream& thisfile, EvtStat& this_evtStat);
        void Helper_ResetEvent(EvtStat& this_evtStat);
        void ClasifyTracks_New(G4String particle, G4double Energy, G4int firstID);
        void ClasifyTracks_AddTrack(G4String particle, G4int trackID, G4int parentID);
        void testwritter(std::fstream& thisfile);
        void WriteAccept();


        G4double nuTime;
        G4double nuMeanEnergy;
        G4double nuEnergy;
        G4double cosTheta;
        G4double weigh;
        G4double primaryEnergy;
        // event quantities
        G4bool foundPhoton;
        G4double primaryX;
        G4double primaryY;
        G4double primaryZ;
        G4double primaryR;
        G4double primaryDirX;
        G4double primaryDirY;
        G4double primaryDirZ;
        G4double photonTheta;
        G4double photonPhi;
        G4double photonR;
        std::vector<HitStat> hitStats;
        std::vector<Test_HitStats> Test_hitStats;
        std::vector<FamilyTrack> AllFamilyTracks;
        G4double TotHits[24];
        // run quantities
        G4String outputFilename;
        std::fstream datafile;
        std::fstream maininfofile;
        std::fstream datafile_gammas2MeV;
        std::fstream datafile_gammas8MeV;
	G4double realdistance;
        
        EvtStat evtStat0; //GPS
        EvtStat evtStat1; //SNgun2, positron
        EvtStat evtStat2; //SNgun2, gamma 2 MeV
        EvtStat evtStat3; //SNgun2, gamma 8 MeV
        CollectiveStat collectiveStat;
	private:
	 //there is no privacy here
};

#endif
