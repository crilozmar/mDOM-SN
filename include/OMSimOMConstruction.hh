#ifndef OMSimOMConstruction_h
#define OMSimOMConstruction_h 1


#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"
#include "OMSimInputData.hh"
#include "OMSimPMTConstruction.hh"


class OMSimOMConstruction
{

public:
    OMSimOMConstruction(OMSimInputData* pData);
    void SelectModule(G4String pSelectedModule);
    void PlaceIt(G4ThreeVector pPosition, G4RotationMatrix* pRotation, G4LogicalVolume*& pMother, bool pIncludeHarness, G4String pNameExtension = "");
    G4UnionSolid* GetOMSolid();
    //OMProperties* mProps = new OMProperties();
    
private:
    void pDOMConstruction();
    void MultiPMTopticalModule();
    OMSimInputData* mData;
    OMSimPMTConstruction* mPMTManager;
    G4String mSelectedModule;
    
    G4LogicalVolume* mGlassLogical;
    G4LogicalVolume* mHarnessLogical;
    bool mCheckOverlaps = true;
    
    G4VisAttributes* mGlassVis= new G4VisAttributes(G4Colour(0.7,0.7,0.8,0.2));
    G4VisAttributes* mGelVis= new G4VisAttributes(G4Colour(0.45,0.5,0.35,0.2));
    G4VisAttributes* mAluVis= new G4VisAttributes(G4Colour(0.8,0.8,0.9,1.0));
    G4VisAttributes* mSteelVis= new G4VisAttributes(G4Colour(0.75,0.75,0.85,1.0));
    G4VisAttributes* mAbsorberVis= new G4VisAttributes(G4Colour(0.2,0.2,0.2,1.0));
    G4VisAttributes* mBoardVis= new G4VisAttributes(G4Colour(0,1,0,1));
    G4VisAttributes* mAirVis= new G4VisAttributes(G4Colour(1.0,1.0,1.0,1.0));
    G4VisAttributes mInvisibleVis = G4VisAttributes::GetInvisible();
};


#endif
//