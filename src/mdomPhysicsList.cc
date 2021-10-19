/** @file mdomPhysicsList.cc
 *  @brief Geant4 physic list
 * 
 *  @author Lew Classen (lclassen@wwu.de)
 * 
 *  @version Geant4 10.7
 */

#include "mdomPhysicsList.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4StepLimiter.hh"
#include "G4SystemOfUnits.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4RadioactiveDecay.hh"
#include "G4DecayPhysics.hh"
#include "G4Cerenkov.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4eplusAnnihilation.hh"

#include "G4LivermoreIonisationModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"

#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"


mdomPhysicsList::mdomPhysicsList():  G4VUserPhysicsList()
{
	defaultCutValue = 0.1*um;
	SetVerboseLevel(0);
}

mdomPhysicsList::~mdomPhysicsList()
{
}

void mdomPhysicsList::ConstructParticle()
{
	G4Gamma::GammaDefinition();
	G4OpticalPhoton::OpticalPhotonDefinition();
    G4GenericIon::GenericIonDefinition();
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4MuonMinus::MuonMinusDefinition();
    G4MuonPlus::MuonPlusDefinition();
	G4NeutrinoE::NeutrinoEDefinition();
	G4AntiNeutrinoE::AntiNeutrinoEDefinition();
	G4NeutrinoMu::NeutrinoMuDefinition();
	G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

void mdomPhysicsList::ConstructProcess()
{
	AddTransportation();

//	The Cherenkov process
	G4Cerenkov* theCerenkovProcess = new G4Cerenkov("Cerenkov");
    theCerenkovProcess->SetTrackSecondariesFirst(true);
    theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    theCerenkovProcess->SetMaxNumPhotonsPerStep(300); 

//	The Livermore models
			G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
			G4LivermoreComptonModel* theLivermoreComptonModel = new G4LivermoreComptonModel();
			theComptonScattering->SetEmModel(theLivermoreComptonModel);
    
	G4eIonisation* theIonizationModel = new G4eIonisation();
	theIonizationModel->SetEmModel(new G4LivermoreIonisationModel());

	G4PhotoElectricEffect* thePhotoElectricEffectModel = new G4PhotoElectricEffect();
	thePhotoElectricEffectModel->SetEmModel(new G4LivermorePhotoElectricModel());

	G4GammaConversion* theGammaConversionModel = new G4GammaConversion();
	theGammaConversionModel->SetEmModel(new G4LivermoreGammaConversionModel());

//	The Decay process
	G4Decay* theDecayProcess = new G4Decay();
	
//	Now assign processes to generated particles
        auto theParticleIterator=GetParticleIterator(); //new geant4
	theParticleIterator->reset();
	while( (*theParticleIterator)() ){
		
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		G4double particleMass = particle->GetPDGMass();
		G4String particleType = particle->GetParticleType();
		G4double particleCharge = particle->GetPDGCharge();
		
		if (particleName == "opticalphoton") {
			pmanager->AddDiscreteProcess(new G4OpAbsorption());
			pmanager->AddDiscreteProcess(new G4OpBoundaryProcess());
	// 		pmanager->AddDiscreteProcess(new G4OpRayleigh);
			pmanager->AddDiscreteProcess(new G4OpMieHG);
		}	
		else if (particleName == "gamma") {
			pmanager->AddDiscreteProcess(theGammaConversionModel);
			pmanager->AddDiscreteProcess(theComptonScattering);
			//pmanager->AddDiscreteProcess(new G4ComptonScattering());
			pmanager->AddDiscreteProcess(thePhotoElectricEffectModel);
		}
		else if (particleName == "GenericIon") {
			G4RadioactiveDecay* theRadioactiveDecay = new G4RadioactiveDecay();
			pmanager->AddProcess(theRadioactiveDecay);
			pmanager->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
			pmanager->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
		}
		else if (particleName == "e-") {
			pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
			pmanager->AddProcess(theIonizationModel,-1,2,2);
			pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
		}
		else if (particleName == "e+") {
			pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
			pmanager->AddProcess(new G4eIonisation,-1,2,2); // The livermore ionization model is only aplicable to electrons
			pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
			pmanager->AddProcess(new G4eplusAnnihilation, 0,-1, 4);
		}

		if (theCerenkovProcess->IsApplicable(*particle)) {
  			pmanager->AddProcess(theCerenkovProcess);
			pmanager->SetProcessOrdering(theCerenkovProcess, idxPostStep);
        }

	}
	
}

void mdomPhysicsList::SetCuts()
{
	SetCutsWithDefault();
	if (verboseLevel>0) DumpCutValuesTable();
}


