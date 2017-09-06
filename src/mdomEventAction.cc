#include "mdomEventAction.hh"
#include "mdomRunAction.hh"
#include "mdomTrackingAction.hh"

#include "mdomAnalysisManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"

//#include "TH1.h"

extern MdomAnalysisManager gAnalysisManager;

mdomEventAction::mdomEventAction()
{}

mdomEventAction::~mdomEventAction()
{}

void mdomEventAction::BeginOfEventAction(const G4Event* evt)
{
//	Reset analysis manager
	gAnalysisManager.ResetEvent();
}

void mdomEventAction::EndOfEventAction(const G4Event* evt)
{
//	Analyze event
	G4double Xev = evt->GetPrimaryVertex()->GetX0();
	G4double Yev = evt->GetPrimaryVertex()->GetY0();
	G4double Zev = evt->GetPrimaryVertex()->GetZ0();
	G4double Rev = pow(pow(Xev,2)+pow(Yev,2)+pow(Zev,2),1./2.); //Spherical R
	
	G4double Pxev = evt->GetPrimaryVertex()->GetPrimary()->GetPx();
	G4double Pyev = evt->GetPrimaryVertex()->GetPrimary()->GetPy();
	G4double Pzev = evt->GetPrimaryVertex()->GetPrimary()->GetPz();
	G4double Pmod = pow(pow(Pxev,2)+pow(Pyev,2)+pow(Pzev,2),1./2.);
	
	G4double Pxdir = Pxev/Pmod;
	G4double Pydir = Pyev/Pmod;
	G4double Pzdir = Pzev/Pmod;
	
	gAnalysisManager.primaryX = Xev;
	gAnalysisManager.primaryY = Yev;
	gAnalysisManager.primaryZ = Zev;
	gAnalysisManager.primaryR= Rev;
	gAnalysisManager.primaryDirX = Pxdir;
	gAnalysisManager.primaryDirY = Pydir;
	gAnalysisManager.primaryDirZ = Pzdir;
	
	gAnalysisManager.AnalyzeEvent();
}
