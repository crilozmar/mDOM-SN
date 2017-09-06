#ifndef mdomPrimaryGeneratorAction0_h
#define mdomPrimaryGeneratorAction0_h 1
 
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
 
class G4GeneralParticleSource;
class G4Event;
class mdomPrimaryGeneratorAction0 : public G4VUserPrimaryGeneratorAction
{
public:
	mdomPrimaryGeneratorAction0();
	~mdomPrimaryGeneratorAction0();

public:
	void GeneratePrimaries(G4Event* anEvent);

private:
	G4GeneralParticleSource* particleSource;
};

#endif
