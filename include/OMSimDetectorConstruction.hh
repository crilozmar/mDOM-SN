#ifndef OMSimDetectorConstruction_h
#define OMSimDetectorConstruction_h 1
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "OMSimInputData.hh"
#include "OMSimMDOM.hh"
#include "OMSimPDOM.hh"
#include "OMSimLOM16.hh"
#include "abcDetectorComponent.hh"

extern G4int gDOM;

class OMSimDetectorConstruction : public G4VUserDetectorConstruction
{
    public: OMSimDetectorConstruction();
    ~OMSimDetectorConstruction();
    G4VPhysicalVolume* Construct();
    abcDetectorComponent* mOpticalModule;
    /*
    if (gDOM == 0) {
        mDOM* mOpticalModule;
    } else if (gDOM == 1) {
        pDOM* mOpticalModule; //public so we can call it from the main to get the LED positions
    } else if (gDOM == 2) {
        LOM16* mOpticalModule; //public so we can call it from the main to get the LED positions
    } else {
        G4String mssg = "--dom parameter only defined for mDOM (0), pDOM (1) or LOM (2), "<< gDOM<<" given";
        critical(mssg);
    }
    */
    
private:
    G4Tubs*			mWorldSolid;
    G4LogicalVolume*		mWorldLogical;
    G4VPhysicalVolume*		mWorldPhysical;
    void ConstructWorld();
    OMSimInputData *mData;
};


#endif
//
