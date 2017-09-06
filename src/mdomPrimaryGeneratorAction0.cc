#include "mdomPrimaryGeneratorAction.hh"
#include "mdomPrimaryGeneratorAction0.hh"
#include "mdomDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTypes.hh"

extern double gworldsize;


mdomPrimaryGeneratorAction0::mdomPrimaryGeneratorAction0()
{
	particleSource = new G4GeneralParticleSource ();
	particleSource->SetParticleDefinition(G4GenericIon::GenericIonDefinition());
}

mdomPrimaryGeneratorAction0::~mdomPrimaryGeneratorAction0()
{
  delete particleSource;
}

void mdomPrimaryGeneratorAction0::GeneratePrimaries(G4Event* anEvent)
{
	particleSource->GeneratePrimaryVertex(anEvent);
}
