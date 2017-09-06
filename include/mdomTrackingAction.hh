#ifndef mdomTrackingAction_h
#define mdomTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class mdomTrackingAction : public G4UserTrackingAction
{
	public:
		mdomTrackingAction();
		~mdomTrackingAction();
	
		void PreUserTrackingAction(const G4Track*);
		void PostUserTrackingAction(const G4Track*);
		
	private:
};

#endif
