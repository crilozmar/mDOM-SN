#ifndef OMSimPMTConstruction_h
#define OMSimPMTConstruction_h 1

#include <boost/property_tree/ptree.hpp>                                        
#include <boost/property_tree/json_parser.hpp> 
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "OMSimInputData.hh"
#include <tuple>
#include <map>
namespace pt = boost::property_tree;


class OMSimPMTConstruction
{

public:
    OMSimPMTConstruction(OMSimInputData* pData);
    //PMTproperties* mProps = new PMTproperties();
    
    void SimulateInternalReflections();
    void SelectPMT(G4String pPMTtoSelect);
    G4UnionSolid* GetPMTSolid();
    void PlaceIt(G4ThreeVector pPosition, G4RotationMatrix* pRotation, G4LogicalVolume*& pMother, G4String pNameExtension = "");
    void PlaceIt(G4Transform3D pTransform, G4LogicalVolume*& pMother, G4String pNameExtension = "");

    G4double mPMTCenterToTip;
    G4double mPMTMaxRad;
    
private:
    OMSimInputData* mData = new OMSimInputData();
    void ConstructIt();
    G4String mSelectedPMT;
    
    
    void DynodeSystemConstruction(G4LogicalVolume* pMother);
    void CathodeBackShield(G4LogicalVolume* pPMTIinner);
    void ReadParameters(G4String pSide,  G4bool pDoubleEllipse, G4bool pFullFit);
    std::tuple<G4UnionSolid*, G4SubtractionSolid*>  BulbConstructionFullFit(G4String pSide, G4bool pDoubleEllipse, G4bool pFullFit);
    std::tuple<G4UnionSolid*, G4SubtractionSolid*>  BulbConstructionOnlyFront(G4String pSide, G4bool pDoubleEllipse);
    std::tuple<G4UnionSolid*, G4SubtractionSolid*>  BulbConstruction(G4String pSide);
    
    G4PVPlacement* mVacuumPhotocathodePlacement;
    G4PVPlacement* mVacuumTubePlacement;
    G4bool mInternalReflections = false;
    G4bool mDynodeSystem = false;
    G4UnionSolid* mPMTSolid;
    G4LogicalVolume* mPMTlogical;
    bool mCheckOverlaps = false;
    
    //Variables from json files are saved in the following members
    G4double mTotalLenght; 
    G4double mTubeWidth;
    G4double mOutRad;
    G4double mEllipseXYaxis;
    G4double mEllipseZaxis;
    G4double mSphereEllipseTransition_r;
    G4double mSpherePos_y;
    G4double mEllipsePos_y;
    G4double mLineFitSlope;
    G4double mEllipseConeTransition_x;
    G4double mEllipseConeTransition_y;
    G4double mConeTorusTransition_x;
    G4double mTorusCircleR;
    G4double mTorusCirclePos_x;
    G4double mTorusCirclePos_y;
    G4double mTorusTubeTransition_y;
    G4double mEllipseXYaxis_2;               
    G4double mEllipsePos_y_2;        
    G4double mEllipseZaxis_2;
    
    //Visual attributes
    G4VisAttributes* mAluVis= new G4VisAttributes(G4Colour(0.8,0.8,0.9,1.0));
    G4VisAttributes* mGlassVis= new G4VisAttributes(G4Colour(0.7,0.7,0.8,0.2));
    G4VisAttributes* mPhotocathodeVis= new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.3));
    G4VisAttributes* mSteelVis= new G4VisAttributes(G4Colour(0.75,0.75,0.85,1.0));
    G4VisAttributes* mAirVis= new G4VisAttributes(G4Colour(1.0,1.0,1.0,1.0));
    G4VisAttributes mInvisibleVis = G4VisAttributes::GetInvisible();
    
    
};


#endif
//