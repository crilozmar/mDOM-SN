#include "mdomStackingAction.hh"
#include "mdomSteppingAction.hh"
#include "mdomAnalysisManager.hh"

#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SystemOfUnits.hh"

extern MdomAnalysisManager gAnalysisManager;

mdomStackingAction::mdomStackingAction() {}


mdomStackingAction::~mdomStackingAction() {}


G4ClassificationOfNewTrack mdomStackingAction::ClassifyNewTrack(const G4Track * aTrack){
  /*
	//Sphericity test -> we take the new created photons, get their direction and kill them
	if ( aTrack->GetDefinition()->GetParticleName() == "opticalphoton" ) {
		if ( aTrack->GetTrackStatus() != fStopAndKill ) {
			G4double dirX = aTrack->GetMomentumDirection().x();
			G4double dirY = aTrack->GetMomentumDirection().y();
			G4double dirZ = aTrack->GetMomentumDirection().z();
			gAnalysisManager.photon_dirX.push_back(dirX);
			gAnalysisManager.photon_dirY.push_back(dirY);
			gAnalysisManager.photon_dirZ.push_back(dirZ);
			return fKill;
		}
	}
	*/
    return fUrgent;
	
}


void mdomStackingAction::NewStage() {}



void mdomStackingAction::PrepareNewEvent() {}