#ifndef mdomPhysicsList_h
#define mdomPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class G4ProductionCuts;

class mdomPhysicsList: public G4VUserPhysicsList
{
	public:
		mdomPhysicsList();
		~mdomPhysicsList();

	protected:
		void ConstructParticle();
		void ConstructProcess();
		void SetCuts();

	private:

};
#endif
