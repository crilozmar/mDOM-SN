/** @file OMSimDetectorConstruction.cc
 *  @brief User defined detector.
 *
 * You should define and construct the detector here...this template is an example for a single mDOM. 
 *
 *  @author 
 *  @date March 2021
 * 
 *  @version Geant4 10.7
 *  
 *  @todo 
 */


#include "OMSimDetectorConstruction.hh"

#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"

#include "OMSimInputData.hh"
#include "OMSimPMTConstruction.hh"
#include "OMSimOMConstruction.hh"
#include "G4Navigator.hh"



extern G4double gworldsize;
extern G4double	gmdomseparation;
extern G4int	gn_mDOMs;
extern G4double gRadius; 
extern G4double gHeight; 
extern G4Navigator* aNavigator;


OMSimDetectorConstruction::OMSimDetectorConstruction()
:mWorldSolid(0), mWorldLogical(0), mWorldPhysical(0)
{}

OMSimDetectorConstruction::~OMSimDetectorConstruction()
{}

/**
 * Construct the world volume
 * @param pData Instance of OMSimInputData for getting materials
 */
void OMSimDetectorConstruction::ConstructWorld(OMSimInputData* pData){
    G4double innerRadiusOfTheTube = 0.*cm;
    G4double startAngleOfTheTube = 0.*deg;
    G4double spanningAngleOfTheTube = 360.*deg;
    gRadius = gworldsize*m;
    gHeight = gworldsize*m;
   
    mWorldSolid = new G4Tubs("World",
                 innerRadiusOfTheTube, 
                 gRadius,
                 gHeight,
                 startAngleOfTheTube, 
                 spanningAngleOfTheTube); 
  
    
    mWorldLogical = new G4LogicalVolume(mWorldSolid, pData->GetMaterial("argWorld"), "World_log", 0, 0, 0);
    //mWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
    G4VisAttributes* World_vis= new G4VisAttributes(G4Colour(0.45,0.5,0.35,0.2));
    mWorldLogical->SetVisAttributes(World_vis);
    mWorldPhysical = new G4PVPlacement (0, G4ThreeVector(0.,0.,0.), mWorldLogical, "World_phys", 0, false, 0);
    aNavigator->SetWorldVolume(mWorldPhysical);

}

/**
 * Construct all solids needed for your study. Call OMSimOMConstruction for simulations with optical modules
 * and OMSimPMTConstruction for simulations with only PMTs.
 * @return World physical for the main
 */
G4VPhysicalVolume* OMSimDetectorConstruction::Construct() {
    
    OMSimInputData* lData = new OMSimInputData();
    lData->SearchFolders();
    
    ConstructWorld(lData);
    
    OMSimOMConstruction* lOpticalModule = new OMSimOMConstruction(lData);
    lOpticalModule->SelectModule("argOM");
    
    if (gn_mDOMs <= 1) {
        lOpticalModule->PlaceIt(G4ThreeVector(0,0,0), new G4RotationMatrix(), mWorldLogical, false, "module_phys_0");
    } else {
        std::stringstream moduleconverter;
        for (unsigned int k = 0; k < gn_mDOMs; k++){
            moduleconverter.str("");
            moduleconverter << "module_phys_" << k ;
            G4double zpos;
            if (gn_mDOMs % 2 == 0) {
                    zpos = gmdomseparation*(gn_mDOMs/2-k-1./2.);
            } else {
                    zpos = gmdomseparation*(gn_mDOMs/2-k);
            }
            lOpticalModule->PlaceIt(G4ThreeVector(0,0,zpos), new G4RotationMatrix(), mWorldLogical, false, moduleconverter.str());
        }
    }
    
    return mWorldPhysical;
}
