#include "mdomSteppingAction.hh"

#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include "mdomAnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4Fragment.hh"
#include "G4LorentzVector.hh"
#include "G4PhotonEvaporation.hh"
#include "Randomize.hh"


extern MdomAnalysisManager gAnalysisManager;
extern G4String gQEfile;
extern G4bool gQE;
extern G4bool gQEweigh;
extern std::vector<double> readColumnDouble (G4String fn, int col);
extern G4int gSNGun;
extern G4int gneutroncapture;

mdomSteppingAction::mdomSteppingAction()
{ 
	if ((gQE) || gQEweigh) {
		QEwavelenght= readColumnDouble(gQEfile, 1);
		QEprobability= readColumnDouble(gQEfile, 2);  
		for (unsigned int u = 0; u <QEwavelenght.size(); u++) {
			QEwavelenght[u] = QEwavelenght.at(u)*nm;
			QEprobability[u] = QEprobability.at(u)/100.;
		}
	}
}


G4double mdomSteppingAction::get_excitation (G4String s)
	{
	  std::stringstream str;
	  G4double d;
	  int i = s.find_first_of("[");
	  int j = s.find_first_of("]");
	  str << s.substr(i+1, j-i-1);
	  str >> d;
	  return d;
	}

void mdomSteppingAction::UserSteppingAction(const G4Step* aStep)
{
	G4Track* aTrack = aStep->GetTrack();
	extern std::vector<G4int>	stats_PMT_hit;
	extern std::vector<G4int>	stats_OM_hit;
//	extern G4long current_event_id;
	std::vector<G4String> n_module; //module number
	std::vector<G4String> n; //PMT number
	extern std::vector<G4String> explode (G4String s, char d);
	G4ThreeVector deltapos;
	G4double Ekin;
	G4double h = 4.136E-15*eV*s;
	G4double c = 2.99792458E17*nm/s;
	G4double lambda;
	
	// Find position of decay
	if ( ! gAnalysisManager.foundPhoton ) {
		if ( aTrack->GetCreatorProcess() ) {
			if ( aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov" ) {
				gAnalysisManager.foundPhoton = true;
				gAnalysisManager.photonTheta = aTrack->GetVertexPosition().getTheta();
				gAnalysisManager.photonPhi = aTrack->GetVertexPosition().getPhi();
				gAnalysisManager.photonR = aTrack->GetVertexPosition().getR();
			}
		}
	}

	// Deexcitation of nuclei after radioactive decay
	if ( aTrack->GetParticleDefinition()->GetParticleType() == "nucleus" ) {
	  if ( get_excitation(aTrack->GetParticleDefinition()->GetParticleName()) != 0.0 ) {
	    if ( aTrack->GetTrackID() != 1 ) {G4LorentzVector aMomentum ( aTrack->GetMomentumDirection(), (aTrack->GetTotalEnergy()+get_excitation(aTrack->GetParticleDefinition()->GetParticleName())) );
	      G4Fragment aNucleus ( aTrack->GetParticleDefinition()->GetAtomicMass(), aTrack->GetParticleDefinition()->GetAtomicNumber(),aMomentum );
	      G4PhotonEvaporation* thePhotonEvap = new G4PhotonEvaporation;
	      thePhotonEvap->BreakItUp(aNucleus);
	    }
	  }
	}
	
	// Check if optical photon is about to hit a photocathode, if so, destroy it and save the hit
	if ( aTrack->GetDefinition()->GetParticleName() == "opticalphoton" ) {
		if ( aTrack->GetTrackStatus() != fStopAndKill ) {
			if ( aStep->GetPostStepPoint()->GetMaterial()->GetName() == "Photocathode" ) {
				Ekin = (aTrack->GetKineticEnergy());
				lambda = h*c/Ekin;
				if ( (!gQE) || ( (QEcheck(lambda)) && (gQE)) || (gQEweigh)) {
					HitStat hitStat;
                    G4int mothercode;
                    G4int trackID = aTrack->GetTrackID();
                    bool debugparameter = false;
                    if ((gSNGun == 2) && (gneutroncapture > 0)) { 
                        for (int i=0; i<(G4int)gAnalysisManager.AllFamilyTracks.size(); i++ ) {
                            if ((G4int)gAnalysisManager.AllFamilyTracks.at(i).photonstracks.size() == 0)
                                {continue;} 
                            if ((std::find(std::begin(gAnalysisManager.AllFamilyTracks.at(i).photonstracks),      std::end(gAnalysisManager.AllFamilyTracks.at(i).photonstracks), trackID) != std::end(gAnalysisManager.AllFamilyTracks.at(i).photonstracks)) || (gAnalysisManager.AllFamilyTracks.at(i).photonstracks.back() == trackID)) {
                                if (gAnalysisManager.AllFamilyTracks.at(i).motherparticle == "e+") {
                                    mothercode = 1;
                                } else if (gAnalysisManager.AllFamilyTracks.at(i).motherparticle == "gamma") {
                                    G4double energy2 = 2*MeV;
                                    G4double energy3 = 8*MeV;
                                    if (gAnalysisManager.AllFamilyTracks.at(i).energy == energy2) {
                                        mothercode = 2;
                                    } else if (gAnalysisManager.AllFamilyTracks.at(i).energy == energy3) {
                                        mothercode = 3;
                                    } else {
                                        G4cout << "ERROR!!! Not energy match for gamma mother!!" << G4endl;
                                        G4cout  << gAnalysisManager.AllFamilyTracks.at(i).energy/MeV << G4endl;
                                    }
                                } else {
                                    G4cout << "ERROR!!! Wrong mother particle type!!" << G4endl;
                                    G4cout << gAnalysisManager.AllFamilyTracks.at(i).motherparticle << G4endl;
                                }
                                debugparameter = true;
                                break;
                            }
                        }
                         if (debugparameter != true) {
                             G4cout << "ERROR!!! I could not found the mothercode for this event! " << G4endl;
                         }
                    } else {
                        mothercode = 0;
                    }
					//G4cout << "DEBUG: " << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
					//G4cout << "DEBUG: " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetName() << G4endl;
					n_module = explode(aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetName(),'_');
					//G4cout <<"DEBUG: " << n_module.at(2) << G4endl;
					n = explode(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName(),'_');
					hitStat.moduleNr = atoi(n_module.at(2));
					hitStat.pmtNr = atoi(n.at(1));
                    hitStat.mothercode = mothercode;
	//				G4cout << hitStat.pmtNr << "\t" << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectTranslation().getTheta() << "\t" << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectTranslation().getPhi() << "\t" << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectTranslation().getR() << G4endl;
					hitStat.position = aTrack->GetPosition();
					hitStat.wavelen = lambda;
					hitStat.hit_time = aTrack->GetGlobalTime();
					if ((gQE) || (gQEweigh)) {
                        QEcheck(lambda);
                        hitStat.QEprob = probability;
					} else { 
                        hitStat.QEprob = 1;
					}
					gAnalysisManager.hitStats.push_back(hitStat);
				}
				aTrack->SetTrackStatus(fStopAndKill);		// kills counted photon to prevent scattering and double-counting 
			}
		}
	}
	//kill everything which hits the photocathode that is not an optical photon...
	if ( aTrack->GetDefinition()->GetParticleName() != "opticalphoton" ) {
		if ( aTrack->GetTrackStatus() != fStopAndKill ) {
			if ( aStep->GetPostStepPoint()->GetMaterial()->GetName() == "Photocathode" ) {
				aTrack->SetTrackStatus(fStopAndKill);
			}
		}
	}
}


bool mdomSteppingAction::QEcheck(G4double lambda) {
	if ((gQE) || (gQEweigh)) {
	if (!QEwavelenght.empty()) {
	  	probability = 0;
		if (lambda < QEwavelenght.at(0) || lambda > QEwavelenght.back()) {
			return false;
		} else {
			G4bool boolparameter = false;
			for (unsigned int u = 0; u <= (QEwavelenght.size()-1); u++) {
				if (QEwavelenght.at(u) == lambda) {
					probability = QEprobability.at(u);
					boolparameter = true;
				}
				else if (QEwavelenght.at(u) > lambda) {
					G4double slope = (QEprobability.at(u)-QEprobability.at(u-1))/(QEwavelenght.at(u)-QEwavelenght.at(u-1));
					probability = (slope*(lambda-QEwavelenght.at(u-1))+QEprobability.at(u-1));
					boolparameter = true;
				}
				if (boolparameter){
					G4double rand = G4UniformRand();
					//G4cout << "wavelenght -> " << lambda/nm << "  points -> " << QEwavelenght.at(u-1)/nm << " & " << QEwavelenght.at(u)/nm << "  points Prob-> " << QEprobability.at(u-1) << " & " << QEprobability.at(u) << "  probabiliy -> " <<probability << "  random "<< rand <<G4endl;  
					if (rand< probability) {
						//G4cout << "PARTY" << G4endl;
						return true;
					} else {
						return false;
					}
				}
			}
		}
	} else {
		G4cout << "ERROR!!! -> Check Quantum efficiency function or data" << G4endl;
		return 0;
	}
	}
}


