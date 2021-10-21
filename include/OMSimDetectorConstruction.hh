#ifndef OMSimDetectorConstruction_h
#define OMSimDetectorConstruction_h 1
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "OMSimInputData.hh"


class OMSimDetectorConstruction : public G4VUserDetectorConstruction
{
    public:         OMSimDetectorConstruction();
    ~OMSimDetectorConstruction();
    G4VPhysicalVolume* Construct();
    
private:
    G4Tubs*				mWorldSolid;
    G4LogicalVolume*		mWorldLogical;
    G4VPhysicalVolume*		mWorldPhysical;
    void ConstructWorld(OMSimInputData* pData);
    
};


#endif
//
