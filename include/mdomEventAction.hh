#ifndef mdomEventAction_h
#define mdomEventAction_h 1

#include "G4UserEventAction.hh"

#include <string>

class G4Event;

class mdomEventAction : public G4UserEventAction
{
	public:
		mdomEventAction();
		~mdomEventAction();

	public:
		void BeginOfEventAction(const G4Event*);
		void EndOfEventAction(const G4Event*);

	private:
};

#endif
