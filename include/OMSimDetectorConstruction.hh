#ifndef OMSimDetectorConstruction_h
#define OMSimDetectorConstruction_h 1
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "OMSimInputData.hh"
#include "OMSimMDOM.hh"


class OMSimDetectorConstruction : public G4VUserDetectorConstruction
{
    public:         OMSimDetectorConstruction();
    ~OMSimDetectorConstruction();
    G4VPhysicalVolume* Construct();
    mDOM* mOpticalModule; //public so we can call it from the main to get the LED positions

private:
    G4Tubs*				mWorldSolid;
    G4LogicalVolume*		mWorldLogical;
    G4VPhysicalVolume*		mWorldPhysical;
    void ConstructWorld();
    OMSimInputData *mData;
};


#endif
//
