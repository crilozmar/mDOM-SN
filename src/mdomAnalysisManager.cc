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
extern G4double	gDistance;
extern G4int	gsimevents;


extern G4int	gn_mDOMs;
extern G4int gneutroncapture;
extern const char* gActivateLED;
extern G4double gInnercolumnradius;
extern G4double gInnercolumn_pos_x;// hicetube x pos relative to dom radius
extern G4double gInnercolumn_pos_y;// hicetube y pos relative to dom radius
//extern G4double gInnercolumn_av_costheta; //Our standard ice has 0.9
extern G4double gInnercolumn_b_inv;
extern G4int gEmitter_mdom;


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
    hitStats.clear();
    Helper_ResetEvent(evtStat0);
    if (gSNGun == 2) {
        Helper_ResetEvent(evtStat1);
        Helper_ResetEvent(evtStat2);
        Helper_ResetEvent(evtStat3);
    }
    AllFamilyTracks.clear();
}

void MdomAnalysisManager::Helper_ResetEvent(EvtStat& this_evtStat)
{
    this_evtStat.nrHitPMTs = 0;
    this_evtStat.nrHitTot = 0;
    this_evtStat.nrHitMod = 0;
    this_evtStat.mothercode = -99;
    this_evtStat.hitsPMTs.clear();
}

void MdomAnalysisManager::ClasifyTracks_New(G4String particle, G4double Energy, G4int firstID) 
{   
    bool boolparam = false;
    for (int i=0; i<(G4int)AllFamilyTracks.size(); i++ ) {
        if (AllFamilyTracks.at(i).grandparentID == firstID) {
            boolparam = true;
            break;
        }
    }
    if (boolparam == false) {
        FamilyTrack thisfamilytrack;
        thisfamilytrack.motherparticle = particle;
        thisfamilytrack.energy = Energy;
        thisfamilytrack.grandparentID = firstID;
        AllFamilyTracks.push_back(thisfamilytrack);
    }
}
void MdomAnalysisManager::ClasifyTracks_AddTrack(G4String particle, G4int trackID, G4int parentID) 
{   //First look only for the main parentID, for speed reasons, then for all IDs
    int this_i;
    if (particle == "opticalphoton") {
        for (int i=0; i<(G4int)AllFamilyTracks.size(); i++ ) {
            if (parentID == AllFamilyTracks.at(i).grandparentID) {
                this_i = i;
                goto opticalphotonend;
            }
        }
        for (int i=0; i<(G4int)AllFamilyTracks.size(); i++ ) {
            for (int j=0; j<(G4int)AllFamilyTracks.at(i).parentstracks.size(); j++ ) {
                if (AllFamilyTracks.at(i).parentstracks.at(j) == parentID) {
                    this_i = i;
                    goto opticalphotonend;
                }
            }
        }
        for (int i=0; i<(G4int)AllFamilyTracks.size(); i++ ) {
            for (int j=0; j<(G4int)AllFamilyTracks.at(i).photonstracks.size(); j++ ) {
                if (AllFamilyTracks.at(i).photonstracks.at(j) == parentID) {
                    this_i = i;
                    goto opticalphotonend;
                }
            }
        }
        opticalphotonend:
            if ((G4int)AllFamilyTracks.at(this_i).photonstracks.size() > 0) {
                if ((std::find(std::begin(AllFamilyTracks.at(this_i).photonstracks), std::end(AllFamilyTracks.at(this_i).photonstracks), trackID) == std::end(AllFamilyTracks.at(this_i).photonstracks)) && (AllFamilyTracks.at(this_i).photonstracks.back() != trackID)) { //to avoid duplicates
                    //find returns the last element if it did not find it
                    AllFamilyTracks.at(this_i).photonstracks.push_back(trackID);
                    }
            } else {
                AllFamilyTracks.at(this_i).photonstracks.push_back(trackID);
            }
    } else {
        for (int i=0; i<(G4int)AllFamilyTracks.size(); i++ ) {
            if (parentID == AllFamilyTracks.at(i).grandparentID) {
                this_i = i;
                goto otherparticlesend;
            }
        }
        for (int i=0; i<(G4int)AllFamilyTracks.size(); i++ ) {
            for (int j=0; j<(G4int)AllFamilyTracks.at(i).parentstracks.size(); j++ ) {
                if (AllFamilyTracks.at(i).parentstracks.at(j) == parentID) {
                    this_i = i;
                    goto otherparticlesend;
                }
            }
        }
        otherparticlesend:
            if ((G4int)AllFamilyTracks.at(this_i).parentstracks.size() > 0) {
                if ((std::find(std::begin(AllFamilyTracks.at(this_i).parentstracks), std::end(AllFamilyTracks.at(this_i).parentstracks), trackID) == std::end(AllFamilyTracks.at(this_i).parentstracks)) && (AllFamilyTracks.at(this_i).parentstracks.back() != trackID)) { //to avoid duplicates
                    //find returns the last element if it did not find it
                    AllFamilyTracks.at(this_i).parentstracks.push_back(trackID);
                    }
            } else {
                AllFamilyTracks.at(this_i).parentstracks.push_back(trackID);
            }
    }
}


void MdomAnalysisManager::AnalyzeEvent() {
    //G4cout << (G4int)hitStats.size() << G4endl;
    /*
    G4int counter = 0;
    for (int i=0; i<(G4int)AllFamilyTracks.size(); i++ ) {
        counter = counter + 1;
        FamilyTrack thisfamilytrack = AllFamilyTracks.at(i);
        G4cout << "** Family track number "<< counter << G4endl;
        G4cout << "Mother Particle -> " << thisfamilytrack.motherparticle << " with size " << thisfamilytrack.tracks.size() <<G4endl;
        /*
        for (int j=0; j<(G4int)thisfamilytrack.tracks.size(); j++ ) {
            G4cout << "track check -- " <<thisfamilytrack.tracks.at(j) << G4endl; 
        } 
    } */
    
    
    //if testing:
    /*
    if ((G4int)Test_hitStats.size() > 0) {
        testwritter(datafile);
    } 
    Test_hitStats.clear();
     */
    // if not:
    if ((G4int)hitStats.size() > 0) {
        Writer_InfoFile();
        if ((gSNGun == 2) && (gneutroncapture > 0)) {
            evtStat1.mothercode = 1;
            Helper_AnalyzeEvent(evtStat1);
            Writer_data(datafile,evtStat1);
            
            evtStat2.mothercode = 2;
            Helper_AnalyzeEvent(evtStat2);
            Writer_data(datafile_gammas2MeV,evtStat2);
            
            evtStat3.mothercode = 3;
            Helper_AnalyzeEvent(evtStat3);
            Writer_data(datafile_gammas8MeV,evtStat3);
        } else {
            evtStat0.mothercode = 0;
            Helper_AnalyzeEvent(evtStat0);
            Writer_data(datafile,evtStat0);
        }
    }
    
}


void MdomAnalysisManager::Helper_AnalyzeEvent(EvtStat& this_evtStat)
{
    std::vector<G4int> modulescounter;
    modulescounter.resize(gn_mDOMs);
    for (int k=0; k<gn_mDOMs; k++) {
        for (int i=0; i<24; i++) {
            std::tuple<G4int,G4int,G4int> hitsPMT;
            std::get<0>(hitsPMT) = k;   
            std::get<1>(hitsPMT) = i;   
            std::get<2>(hitsPMT) = 0; 
            G4bool hit = false;
            for ( int j=0; j<(G4int)hitStats.size(); j++ ) {
                if (hitStats[j].mothercode == this_evtStat.mothercode) {
                    if ( (hitStats[j].moduleNr == k) && (hitStats[j].pmtNr == i) ) {
                        this_evtStat.nrHitTot += 1;
                        std::get<2>(hitsPMT) += 1; 
                        if ( ! hit ) {
                            this_evtStat.nrHitPMTs += 1;
                            if (modulescounter.at(k) == 0) {
                                modulescounter.at(k) += 1;
                                this_evtStat.nrHitMod += 1;
                            }
                            hit = true;
                        } 
                    } 
                }
            } 
            if (std::get<2>(hitsPMT) != 0 ) this_evtStat.hitsPMTs.push_back(hitsPMT);
        }
    }

}

void MdomAnalysisManager::WriteHeader()
{
G4String part;
G4String nu;
  //datafile << "#Generated in a cylinder of r = " << gRadius/m<< "m and total height of "<< 2*gHeight/m << " m"<< G4endl;
  if (gSNGun == 1){
    maininfofile<< "# World with "<< gn_mDOMs << " mDOMs" << G4endl;
	maininfofile << "# Neutrino electron elastic scattering from SN at 10kpc" << G4endl;
	part = "e-";
	nu = "nu";
  } else if (gSNGun == 2){
    maininfofile<< "# World with "<< gn_mDOMs << " mDOMs" << G4endl;
	maininfofile << "# Inverse beta decay from SN at 10kpc" << G4endl;
	part = "e+";
	nu = "nubar";
  } else if (gSNGun == 3){
    maininfofile<< "# World with "<< gn_mDOMs << " mDOMs" << G4endl;
	maininfofile << "# Solar neutrinos! Weigh take into account the total cross section but not the flux!" << G4endl;
	part = "e-";
	nu = "nu";
  }
  
  if (gSNGun == 0) {
    realdistance = gDistance*m;// - ((0.5 * 356.0)/1000.0)*m;
    datafile << "# World with "<< gn_mDOMs << " mDOMs" << G4endl;
    datafile << "# Total photons -> " << gsimevents << G4endl;
    datafile << "# Activated LEDs: "<< gActivateLED << " in mDOM " << gEmitter_mdom << G4endl;
    datafile << "# Column radius [cm] " << gInnercolumnradius/cm << " | Column shift x,y [cm] " << gInnercolumn_pos_x/cm << "," << gInnercolumn_pos_y/cm << " | Effective sca. lenght in column [cm] " << gInnercolumn_b_inv;
    if (gQEweigh) {
        datafile << " | Total QE prob ";
    }
    datafile << G4endl;
    datafile <<"# Module number | PMT number | hit time | theta (wrt LED) [rad]";
    if (gQEweigh) {
            datafile << "| QE prob";
    }
    datafile << G4endl;
  }
  
  else if (gSNGun == 1 || gSNGun == 2  || gSNGun == 3) {
    maininfofile << "#"<< G4endl;
    maininfofile << "# Time of Flux [s] | Mean energy of "<<nu<<" | "<<nu<<" energy | costheta of "<<part<<" from z dir | "<<part<<" energy | event weigh | Vertex Position (X, Y, Z) [m] | Primary direction (Px,Py,Pz)" << G4endl;
    maininfofile << "#" << G4endl;
    HelpTheHeader(datafile);
    if (gSNGun == 2) {
        if (gneutroncapture == 3) {
            HelpTheHeader(datafile_gammas2MeV);
            HelpTheHeader(datafile_gammas8MeV);
        }else if (gneutroncapture == 1) {
            HelpTheHeader(datafile_gammas2MeV);
        }else if (gneutroncapture == 2) {
            HelpTheHeader(datafile_gammas8MeV);
        }
        }
    }
}

void MdomAnalysisManager::HelpTheHeader(std::fstream& thisfile)
{
  thisfile << "# Total hits | PMTs hit |";
                    if (gQEweigh) {
                        thisfile << " Total QE prob |";
                    }
                    thisfile <<"...for PMT hit...| Module number | PMT number | Hits in that PMT |";
                    thisfile << "...for Hit...";
                    if (gQEweigh) {
                        thisfile << " QE prob |";
                    }
                    thisfile << " hit time |";
                    if (gQEweigh) {
                        thisfile << "...end of hit loop...| Total QE prob in PMT |" << G4endl;
                    }
  thisfile << "#"<< G4endl;
}


void MdomAnalysisManager::Writer_InfoFile() {
      if (gSNGun == 1 || gSNGun == 2 || gSNGun == 3)  {
        maininfofile << nuTime/s << "\t";
        maininfofile << nuMeanEnergy/MeV<< "\t";
        maininfofile << nuEnergy/MeV<< "\t";
        maininfofile << cosTheta<< "\t";
        maininfofile << primaryEnergy/MeV << "\t";
        maininfofile << weigh << "\t\t";
      }
     if (gSNGun == 0)  {
        maininfofile << realdistance/m << "\t";
      }
        maininfofile << primaryX/m <<"\t";
        maininfofile << primaryY/m << "\t";
        maininfofile << primaryZ/m << "\t\t";
        maininfofile << primaryDirX <<"\t";
        maininfofile << primaryDirY << "\t";
        maininfofile << primaryDirZ << "\t\t";
    maininfofile << G4endl;
}

void MdomAnalysisManager::Writer_data(std::fstream& thisfile, EvtStat& this_evtStat)
{
    if (gSNGun != 0) {
    thisfile << this_evtStat.nrHitTot << "\t";
    thisfile << this_evtStat.nrHitMod << "\t";
    thisfile << this_evtStat.nrHitPMTs << "\t";
    if (gQEweigh) {
        G4double sum = 0;
        for (unsigned int i = 0 ; i<hitStats.size(); i++) {
            sum = sum + hitStats[i].QEprob;
            }
        thisfile << sum << "\t";
        }
    thisfile << "\t";
    }
    for ( int j=0; j<(G4int)this_evtStat.hitsPMTs.size(); j++ ) {
        thisfile << std::get<0>(this_evtStat.hitsPMTs[j]) << "\t";
        thisfile << std::get<1>(this_evtStat.hitsPMTs[j]) << "\t";
        if (gSNGun != 0) {
            thisfile << std::get<2>(this_evtStat.hitsPMTs[j]) << "\t";
        }
        for (int i=0; i<(G4int)hitStats.size(); i++) {
            if (this_evtStat.mothercode == hitStats[i].mothercode) {
                if ( (std::get<0>(this_evtStat.hitsPMTs[j]) == hitStats[i].moduleNr) && (std::get<1>(this_evtStat.hitsPMTs[j]) == hitStats[i].pmtNr) ) {
                    //thisfile << hitStats[i].wavelen/nm << "\t";
                    thisfile << hitStats[i].hit_time/ns << "\t";
                    thisfile << hitStats[i].reltheta/rad << "\t";
                    if (gQEweigh) {
                        thisfile << hitStats[i].QEprob << "\t";
                        }
                    }
            }
        }
    }
    thisfile << G4endl;
}

void MdomAnalysisManager::testwritter(std::fstream& thisfile)
{
    for (int i=0; i<(G4int)Test_hitStats.size(); i++) {
        thisfile << Test_hitStats[i].reltheta/rad << "\t";
        thisfile << Test_hitStats[i].globtheta/rad << "\t";
    }
    thisfile << G4endl;
}

