/** @file mdomStackingAction.cc
 *  @brief Geant4 tracking action
 * 
 *  @author Lew Classen (lclassen@wwu.de)
 * 
 *  @version Geant4 10.7
 */

#include "mdomEventAction.hh"
#include "mdomTrackingAction.hh"
#include "mdomRunAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//#include "TH1D.h"

mdomTrackingAction::mdomTrackingAction()
:G4UserTrackingAction()
{
}

mdomTrackingAction::~mdomTrackingAction()
{
}

void mdomTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{	
   
   //G4cout << "end of tracking action "<<G4endl;
}

void mdomTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
}
