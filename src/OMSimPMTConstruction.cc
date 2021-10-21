/** @file OMSimPMTConstruction.cc
 *  @brief Construction of the PMTs.
 *
 *  This class creates the solids of the PMTs and place them in the detector/OMs.
 * 
 *  @author Lew Classen, Martin Unland
 *  @date March 2021
 * 
 *  @version Geant4 10.7
 * 
 *  @todo 
 * 
 */

#include <dirent.h>

#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4Ellipsoid.hh>
#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include <G4PVPlacement.hh>
#include <G4Sphere.hh>
#include <G4SubtractionSolid.hh>
#include "G4SystemOfUnits.hh"
#include <G4Torus.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4UnitsTable.hh>
#include "G4VisAttributes.hh"

#include "OMSimPMTConstruction.hh"


extern G4bool gVisual;
extern G4int gPMT;

// Hash and mix are needed for switch statemens with strings
uint64_t constexpr mix(char m, uint64_t s)
{
    return ((s<<7) + ~(s>>3)) + ~m;
}

uint64_t constexpr hash(const char * m)
{
    return (*m) ? mix(*m,hash(m+1)) : 0;
}


/**
 * Constructor of the class. The InputData instance has to be passed here in order to avoid loading the input data twice and redifining the same materials.
 * @param pData OMSimInputData instance 
 */
OMSimPMTConstruction::OMSimPMTConstruction(OMSimInputData* pData){
    mData = pData;
}

/**
 * Returns the Solid of the constructed PMT.
 * @return G4UnionSolid of the PMT  
 */
G4UnionSolid* OMSimPMTConstruction::GetPMTSolid(){
    return mPMTSolid;
}

/**
 * The physical volume of the constructed PMT is produced here
 * @param pPosition G4ThreeVector with position of the module (as in G4PVPlacement())
 * @param pRotation G4RotationMatrix with rotation of the module (as in G4PVPlacement())
 * @param pMother G4LogicalVolume where the module is going to be placed (as in G4PVPlacement())
 * @param pNameExtension G4String name of the physical volume. You should not have two physicals with the same name
 */
void OMSimPMTConstruction::PlaceIt(G4ThreeVector pPosition, G4RotationMatrix* pRotation, G4LogicalVolume*& pMother, G4String pNameExtension ){
    G4PVPlacement* lPMTPhysical = new G4PVPlacement(pRotation, pPosition, mPMTlogical, "PMT_"+pNameExtension, pMother, false, 0, mCheckOverlaps);
        if (mDynodeSystem){
        G4cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << G4endl;
        G4OpticalSurface* Photocathode_opsurf = new G4OpticalSurface("Photocathode_opsurf");
        new G4LogicalBorderSurface("Photocathode_out", mVacuumPhotocathodePlacement, lPMTPhysical, Photocathode_opsurf);
        new G4LogicalBorderSurface("Photocathode_in", lPMTPhysical, mVacuumPhotocathodePlacement, Photocathode_opsurf);//*/
        new G4LogicalBorderSurface("PMT_mirrorglass", mVacuumTubePlacement, lPMTPhysical, mData->GetOpticalSurface("Refl_100polished"));//*/
        new G4LogicalBorderSurface("PMT_mirrorglass", lPMTPhysical, mVacuumTubePlacement, mData->GetOpticalSurface("Refl_100polished"));//*/
        
    }
}
/**
 * The physical volume of the constructed PMT is produced here
 * @param pTransform G4Transform3D with position & rotation of PMT
 * @param pMother G4LogicalVolume where the module is going to be placed (as in G4PVPlacement())
 * @param pNameExtension G4String name of the physical volume. You should not have two physicals with the same name
 */
void OMSimPMTConstruction::PlaceIt(G4Transform3D pTransform, G4LogicalVolume*& pMother, G4String pNameExtension ){
    G4PVPlacement* lPMTPhysical = new G4PVPlacement(pTransform, mPMTlogical, "PMT_"+pNameExtension, pMother, false, 0, mCheckOverlaps);
    
    if (mDynodeSystem){
        G4cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << G4endl;
        G4OpticalSurface* Photocathode_opsurf = new G4OpticalSurface("Photocathode_opsurf");
        new G4LogicalBorderSurface("Photocathode_out", mVacuumPhotocathodePlacement, lPMTPhysical, Photocathode_opsurf);
        new G4LogicalBorderSurface("Photocathode_in", lPMTPhysical, mVacuumPhotocathodePlacement, Photocathode_opsurf);//*/
        new G4LogicalBorderSurface("PMT_mirrorglass", mVacuumTubePlacement, lPMTPhysical, mData->GetOpticalSurface("Refl_100polished"));//*/
        new G4LogicalBorderSurface("PMT_mirrorglass", lPMTPhysical, mVacuumTubePlacement, mData->GetOpticalSurface("Refl_100polished"));//*/
        
    }
}

/**
 * The PMT is constructed allowing internal reflections (in the mirrored part opposite to the photocathode). Do not use it if you really don't need it!!
 */
void OMSimPMTConstruction::SimulateInternalReflections(){
    
    
    if (mSelectedPMT=="pmt_Hamamatsu_R15458"){
        mInternalReflections = true;
        mDynodeSystem = true;
    }
    else{
        G4cerr << "Internal reflections inside PMT is on, but this PMT has no Dynode system defined, only vacuum! I will continue with internal reflections anyway..." << G4endl;
        mInternalReflections = true;
        mDynodeSystem = false;
    }
    ConstructIt();
}

/**
 * Constructs the PMT.
 */
void OMSimPMTConstruction::ConstructIt(){
    
    G4UnionSolid* lGlassInside;
    G4UnionSolid* lVacuumMirrorSolid;
    G4SubtractionSolid* lVacuumPhotocathodeSolid;
    G4SubtractionSolid* lVacuumTubeSolid;
    
    std::tie(mPMTSolid, lVacuumPhotocathodeSolid) = BulbConstruction("jOuterShape");
    std::tie(lGlassInside, lVacuumPhotocathodeSolid) = BulbConstruction("jInnerShape");
    
    mPMTlogical = new G4LogicalVolume(mPMTSolid, mData->GetMaterial("RiAbs_Glass_Tube"), "PMT tube logical");
    mPMTlogical->SetVisAttributes(mGlassVis);
    
    if (mInternalReflections){
        lVacuumTubeSolid = new G4SubtractionSolid("Vacuum Tube solid", lGlassInside, lVacuumPhotocathodeSolid, 0, G4ThreeVector(0,0,0));
        
        G4LogicalVolume* lVacuumPhotocathodeLogical      = new G4LogicalVolume(lVacuumPhotocathodeSolid, mData->GetMaterial("Ri_Vacuum"), "Photocathode area vacuum");
        G4LogicalVolume* lVacuumTubeLogical              = new G4LogicalVolume(lVacuumTubeSolid, mData->GetMaterial("Ri_Vacuum"), "Tube vacuum");
        
        mVacuumPhotocathodePlacement    = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lVacuumPhotocathodeLogical, "Vacuum_1", mPMTlogical, false, 0, mCheckOverlaps);
        mVacuumTubePlacement            = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lVacuumTubeLogical, "Vacuum_3", mPMTlogical, false, 0, mCheckOverlaps);
        
        if (mDynodeSystem){
            DynodeSystemConstruction(lVacuumTubeLogical);
        }
        lVacuumPhotocathodeLogical->SetVisAttributes(mPhotocathodeVis);
        
        //lVacuumTubeLogical->SetVisAttributes(mAirVis);
        //lVacuumMirrorLogical->SetVisAttributes(mInvisibleVis);
    }
    
    else{
        G4LogicalVolume* lTubeVacuum = new G4LogicalVolume(lGlassInside, mData->GetMaterial("Ri_Vacuum"), "PMTvacuum");
        G4LogicalVolume* lPhotocathode = new G4LogicalVolume(lVacuumPhotocathodeSolid, mData->GetMaterial("RiAbs_Photocathode"), "Photocathode");
        G4PVPlacement* lVacuumPlacement    = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lTubeVacuum, "VacuumTube", mPMTlogical, false, 0, mCheckOverlaps);
        G4PVPlacement* lPhotocathodePlacement    = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lPhotocathode, "VacuumPhoto", lTubeVacuum, false, 0, mCheckOverlaps);
        CathodeBackShield(lTubeVacuum);
        lPhotocathode->SetVisAttributes(mPhotocathodeVis);
        lTubeVacuum->SetVisAttributes(mInvisibleVis);
        
    }
    
    
}

/**
 * Construction & placement of the dynode system entrance for internal reflections. Currently only geometry for Hamamatsu R15458.
 * @param pMother GloogicalVolume of the mother, where the dynode system entrance is placed (vacuum volume)
 */
void OMSimPMTConstruction::DynodeSystemConstruction(G4LogicalVolume* pMother){
    
    G4double lEntranceH                 = 19.5*mm;
    G4double lEntranceW                 = 21.9*mm;
    G4double lPlateWidth                = 0.7*mm; 
    G4double lPlateRad                  = 46.*0.5*mm; 
    G4double lFirstDynodeH              = 18.06*mm;
    G4double lFirstDynodeCurvatureR     = 14.6439*mm;
    
    G4double lFurthestX = 13.8593179*mm;
    G4double lFurthestY = 2.0620*mm;
    
    G4Tubs* lFirstDynode = new G4Tubs("FirstDynode", 0, lFirstDynodeCurvatureR, lEntranceW*0.5, 8.46*degree, 106.64*degree-8.46*degree);
    G4Tubs* lSubstTube = new G4Tubs("FirstDynode", 0, lFirstDynodeCurvatureR, lEntranceW*1, 0, 114.39*degree);
    G4Box*  lSubBox = new G4Box("box", 30*mm, lFurthestY, lEntranceW);
    
    G4double lleftBoxH = 14.102*mm;
    G4double lleftBoxW = 4.216*mm;
    G4Box*  lleftBox = new G4Box("box", 0.5*lleftBoxW, 0.5*lleftBoxH, lEntranceW*0.5);
    
    G4double lThirdDynodeHeight = 6.2*mm;
    
    G4double lLargeBoxH =  lleftBoxH-lThirdDynodeHeight-lFurthestY;
    G4double lLargeBoxW =  7*mm;
    G4Box*  lLargeBox = new G4Box("box", 0.5*lLargeBoxW, 0.5*lLargeBoxH, lEntranceW*0.5);
    
    G4SubtractionSolid* newSubst = new G4SubtractionSolid("Subst", lFirstDynode, lSubstTube, 0,G4ThreeVector(lFurthestX,lFurthestY, 0)); 
    newSubst = new G4SubtractionSolid("Subst", newSubst, lSubBox, 0,G4ThreeVector(0,0, 0)); 
    
    G4UnionSolid* lUnion = new G4UnionSolid("curvature leftbox", newSubst, lleftBox, 0, G4ThreeVector(-lleftBoxW*0.5, lleftBoxH*0.5,0));
    lUnion = new G4UnionSolid("curvature leftbox", lUnion, lLargeBox, 0, G4ThreeVector(-lLargeBoxW, lLargeBoxH*0.5+lFurthestY,0));
    
    
    
    G4Tubs* lDynodeEntrancePlate = new G4Tubs("DynodeEntrancePlate", 0, lPlateRad, lPlateWidth, 0, 360.*degree);
    G4Box* lEntranceWindowSubst = new G4Box("Entrance window", lEntranceW*0.5,lEntranceH*0.5,30*mm);
    G4SubtractionSolid* lDynodeEntrance = new G4SubtractionSolid("Substracted Plate", lDynodeEntrancePlate, lEntranceWindowSubst);
    
    G4Tubs* lDynodeSystem = new G4Tubs("DynodeEntrancePlate", 0, lPlateRad-2*mm, 0.5*15*mm, 0, 360.*degree);
    G4Tubs* lAbsorber = new G4Tubs("DynodeEntrancePlate", lPlateRad-1.99*mm, (52.2*mm-2*1.351*mm)*0.5, 0.5*10*mm, 0, 360.*degree);
    G4Tubs* lBackGlass = new G4Tubs("DynodeEntrancePlate", lPlateRad-1.99*mm, (52.2*mm-2*1.351*mm)*0.5, 0.5*2*mm, 0, 360.*degree);
    
    G4RotationMatrix* myRotation1 = new G4RotationMatrix();
    myRotation1->rotateX(90.*deg);
    myRotation1->rotateY(90.*deg);
    G4SubtractionSolid* lSubstractedDynodeSystem = new G4SubtractionSolid("Substracteddynode", lDynodeSystem, lUnion, myRotation1, G4ThreeVector(0,-(lFurthestX-lEntranceH*0.5), +lFurthestY+15*0.5*mm+0.1*mm));
    
    G4LogicalVolume* lDynodeSystemLog = new G4LogicalVolume(lSubstractedDynodeSystem, mData->GetMaterial("NoOptic_Reflector"), "PMT_12199 tube logical");
    G4LogicalVolume* lAbsorberLog = new G4LogicalVolume(lAbsorber,mData->GetMaterial("NoOptic_Absorber"), "PMT_12199_absorber");
    G4LogicalVolume* lBackGlassLog = new G4LogicalVolume(lBackGlass, mData->GetMaterial("RiAbs_Glass_Tube"), "PMT_12199_absorber");
    G4LogicalVolume* lDynodeEntranceLog = new G4LogicalVolume(lDynodeEntrance, mData->GetMaterial("NoOptic_Reflector"), "PMT_12199 tube logical");
    
    G4double lDistanceDynodeFrontalEllipsoid = -(29.0-7.9-lPlateWidth)*mm;
    G4PVPlacement* lDynodePlatePhysical = new G4PVPlacement (0,G4ThreeVector(0,0,lDistanceDynodeFrontalEllipsoid), lDynodeEntranceLog, "Frontal_plate", pMother, false, 0, mCheckOverlaps);
    G4PVPlacement* lDynodeSystemPhysical = new G4PVPlacement (0,G4ThreeVector(0, 0,-15*0.5*mm-lPlateWidth+lDistanceDynodeFrontalEllipsoid), lDynodeSystemLog, "FirstDynode", pMother, false, 0, mCheckOverlaps);
    G4PVPlacement* lAbsorberPhys = new G4PVPlacement (0,G4ThreeVector(0, 0,-40*0.5*mm-lPlateWidth+lDistanceDynodeFrontalEllipsoid), lAbsorberLog, "Absorber", pMother, false, 0, mCheckOverlaps);
    G4PVPlacement* lBackGlassPhys = new G4PVPlacement (0,G4ThreeVector(0, 0,-18*0.5*mm-lPlateWidth+lDistanceDynodeFrontalEllipsoid), lBackGlassLog, "BackGlass", pMother, false, 0, mCheckOverlaps);
    new G4LogicalBorderSurface("PMT_platemirror", mVacuumTubePlacement, lDynodePlatePhysical, mData->GetOpticalSurface("Refl_100polished"));
    new G4LogicalBorderSurface("PMT_platemirror", lDynodePlatePhysical, mVacuumTubePlacement, mData->GetOpticalSurface("Refl_100polished"));
    new G4LogicalBorderSurface("PMT_1dynmirror", mVacuumTubePlacement, lDynodeSystemPhysical, mData->GetOpticalSurface("Refl_100polished"));
    new G4LogicalBorderSurface("PMT_1dynmirror", lDynodeSystemPhysical, mVacuumTubePlacement, mData->GetOpticalSurface("Refl_100polished"));
    
}

/**
 * 
 */
std::tuple<G4UnionSolid*, G4SubtractionSolid*>  OMSimPMTConstruction::BulbConstruction(G4String pSide){
    G4String lParameterSet = mData->mPMTdata->mTable.at(mSelectedPMT).get<G4String>("jShape");
    G4UnionSolid* U1;
    G4SubtractionSolid* S1;
    if (mInternalReflections){
        switch(hash(lParameterSet)){
            case hash("jOnlyFront_SphereEllipse")   : std::tie(U1, S1) = BulbConstructionFullFit(pSide, false, false); break;
            case hash("jOnlyFront_2Ellipses")       : std::tie(U1, S1) = BulbConstructionFullFit(pSide, true, false); break;
            case hash("jFullFit_SphereEllipse")     : std::tie(U1, S1) = BulbConstructionFullFit(pSide, false, true); break;
            case hash("jFullFit_2Ellipses")         : std::tie(U1, S1) = BulbConstructionFullFit(pSide, true, true); break;
            default                                 : std::tie(U1, S1) = BulbConstructionFullFit(pSide, false, false); break;
        }
    }
    
    else{
        switch(hash(lParameterSet)){
            case hash("jOnlyFront_SphereEllipse")   : std::tie(U1, S1) = BulbConstructionFullFit(pSide, false, false); break;
            case hash("jOnlyFront_2Ellipses")       : std::tie(U1, S1) = BulbConstructionFullFit(pSide, true, false); break;
            case hash("jFullFit_SphereEllipse")     : std::tie(U1, S1) = BulbConstructionFullFit(pSide, false, false); break;
            case hash("jFullFit_2Ellipses")         : std::tie(U1, S1) = BulbConstructionFullFit(pSide, true, false); break;
            default                                 : std::tie(U1, S1) = BulbConstructionFullFit(pSide, false, false); break;
        }
    }
    return std::make_tuple(U1, S1);
}

/**
 * 
 */
void OMSimPMTConstruction::ReadParameters(G4String pSide,  G4bool pDoubleEllipse, G4bool pFullFit){
    mTotalLenght               = mData->mPMTdata->Get(pSide+".jTotalLenght");
    mTubeWidth                 = mData->mPMTdata->Get(pSide+".jTubeWidth");
    mOutRad                    = mData->mPMTdata->Get(pSide+".jOutRad");
    mEllipseXYaxis             = mData->mPMTdata->Get(pSide+".jEllipseXYaxis");
    mEllipseZaxis              = mData->mPMTdata->Get(pSide+".jEllipseZaxis");
    mSphereEllipseTransition_r = mData->mPMTdata->Get(pSide+".jSphereEllipseTransition_r");
    mSpherePos_y               = mData->mPMTdata->Get(pSide+".jSpherePos_y");
    mEllipsePos_y              = mData->mPMTdata->Get(pSide+".jEllipsePos_y");
    
    if (pFullFit){
        mLineFitSlope              = mData->mPMTdata->Get(pSide+".jLineFitSlope");
        mEllipseConeTransition_x   = mData->mPMTdata->Get(pSide+".jEllipseConeTransition_x");
        mEllipseConeTransition_y   = mData->mPMTdata->Get(pSide+".jEllipseConeTransition_y");
        mConeTorusTransition_x     = mData->mPMTdata->Get(pSide+".jConeTorusTransition_x");
        mTorusCircleR              = mData->mPMTdata->Get(pSide+".jTorusCircleR");
        mTorusCirclePos_x          = mData->mPMTdata->Get(pSide+".jTorusCirclePos_x");
        mTorusCirclePos_y          = mData->mPMTdata->Get(pSide+".jTorusCirclePos_y");
        mTorusTubeTransition_y     = mData->mPMTdata->Get(pSide+".jTorusTubeTransition_y");
    }
    
    if (pDoubleEllipse) {
        mEllipseXYaxis_2             = mData->mPMTdata->Get(pSide+".jEllipseXYaxis_2");
        mEllipseZaxis_2              = mData->mPMTdata->Get(pSide+".jEllipseZaxis_2");
        mEllipsePos_y_2              = mData->mPMTdata->Get(pSide+".jEllipsePos_y_2");
    }
    
    if (pSide=="jOuterShape"){
        mPMTMaxRad = mEllipseXYaxis;
    }
    
}

/**
 * 
 */
void OMSimPMTConstruction::CathodeBackShield(G4LogicalVolume* pPMTinner){
    
    ReadParameters("jInnerShape",  false , false);
    G4double lPCDiameter  = 2 * mEllipseXYaxis;
    G4double lShieldBottomcut = mEllipseZaxis / mEllipseXYaxis * std::sqrt(std::pow(mEllipseXYaxis, 2.) - std::pow(0.5*lPCDiameter, 2.));
    G4double lShieldRad = -0.01*mm + mEllipseXYaxis / mEllipseZaxis * std::sqrt(std::pow(mEllipseZaxis, 2.) - std::pow(0.5*lShieldBottomcut, 2.)); 	
    G4Tubs* lShieldSolid = new G4Tubs("Shield solid", 0, lShieldRad, 0.1*mm, 0, 2*CLHEP::pi);
    G4LogicalVolume* lShieldLogical = new G4LogicalVolume(lShieldSolid, mData->GetMaterial("NoOptic_Absorber"), "Shield logical");
    new G4PVPlacement (0, G4ThreeVector(0,0, -0.1*mm), lShieldLogical, "Shield physical", pPMTinner, false, 0, mCheckOverlaps);
    lShieldLogical->SetVisAttributes(mSteelVis);
}


/**
 * 
 */
std::tuple<G4UnionSolid*, G4SubtractionSolid*>  OMSimPMTConstruction::BulbConstructionFullFit(G4String pSide,  G4bool pDoubleEllipse, G4bool pFullFit){
    G4String lParameterSet = mData->mPMTdata->mTable.at(mSelectedPMT).get<G4String>("jShape");
    
    ReadParameters(pSide, pDoubleEllipse, pFullFit);
    
    G4double SphereAngle  = asin(mSphereEllipseTransition_r / mOutRad);
    
    // PMT frontal glass envelope as union of sphere and ellipse	
    G4Ellipsoid* lBulbEllipsoid = new G4Ellipsoid("Solid Bulb Ellipsoid", mEllipseXYaxis, mEllipseXYaxis, mEllipseZaxis);
    G4Sphere* lBulbSphere = new G4Sphere("Solid Bulb Ellipsoid", 0.0, mOutRad, 0, 2*CLHEP::pi, 0, SphereAngle);
    G4UnionSolid* lBulbSolid = new G4UnionSolid("Solid Bulb", lBulbEllipsoid, lBulbSphere, 0, G4ThreeVector(0,0, mSpherePos_y- mEllipsePos_y));
    
    if (lParameterSet == "jOnlyFront_2Ellipses" || lParameterSet == "jFullFit_2Ellipses") {
        G4Ellipsoid* lBulbEllipsoid_2 = new G4Ellipsoid("Solid Bulb Ellipsoid 2", mEllipseXYaxis_2, mEllipseXYaxis_2, mEllipseZaxis_2);
        G4double lExcess = mEllipsePos_y - mEllipsePos_y_2;
        G4Tubs* lSubtractionTube =  new G4Tubs("substracion_tube_large_ellipsoid", 0.0, mEllipseXYaxis_2*2, 0.5*mTotalLenght,0,2*CLHEP::pi);
        G4SubtractionSolid* lSubstractedLargeEllipsoid = new G4SubtractionSolid("Substracted Bulb Ellipsoid 2", lBulbEllipsoid_2, lSubtractionTube, 0,  G4ThreeVector(0,0,lExcess-mTotalLenght*0.5));  
        lBulbSolid = new G4UnionSolid("Solid Bulb", lBulbSolid, lSubstractedLargeEllipsoid, 0, G4ThreeVector(0,0,mEllipsePos_y_2-mEllipsePos_y));
    }
    else {
        G4double lSphereAngle  = asin( mSphereEllipseTransition_r / mOutRad);
        G4Sphere* lBulbSphere = new G4Sphere("Solid Bulb Ellipsoid", 0.0, mOutRad, 0, 2*CLHEP::pi, 0, lSphereAngle);
        lBulbSolid = new G4UnionSolid("Solid Bulb", lBulbEllipsoid, lBulbSphere, 0, G4ThreeVector(0,0,mSpherePos_y-mEllipsePos_y));   
    }
    
    
    // Defining volume with boundaries of photocathode volume
    G4Tubs* lLargeTube = new G4Tubs("LargeTube", 0, mEllipseXYaxis, 50*cm,0,2*CLHEP::pi);
    G4SubtractionSolid* lPhotocathodeSide = new G4SubtractionSolid("SubstractionPhotocathodeSide", lBulbSolid, lLargeTube, 0, G4ThreeVector(0,0,-50*cm));
    
    //     G4SubtractionSolid* lReflectiveSide = new G4SubtractionSolid("Reflective side", lBulbSolid, lLargeTube, 0, G4ThreeVector(0,0, mOutRad*2));
    
    // Rest of tube
    G4double lFrontToEllipse_y = mOutRad + mSpherePos_y - mEllipsePos_y;
    
    G4double lMissingTubeLength = (mTotalLenght-lFrontToEllipse_y)*0.5*mm;       
    G4Tubs* lBulkSolid = new G4Tubs("Bulb bulk solid", 0.0, 0.5*mTubeWidth, lMissingTubeLength, 0, 2*CLHEP::pi);
    
    if (pSide == "jOuterShape"){
        mPMTCenterToTip = lFrontToEllipse_y;
    }
    
    if (pFullFit){
        // Creating Cone
        G4double lConeLength_x = mEllipseConeTransition_x - mConeTorusTransition_x;
        G4double lConeHalfHeight = mLineFitSlope * lConeLength_x*0.5;
        G4Cons* lCone = new G4Cons("Solid substraction cone", mConeTorusTransition_x, mConeTorusTransition_x+mTubeWidth, mEllipseConeTransition_x, mEllipseConeTransition_x+mTubeWidth, lConeHalfHeight, 0, 2*CLHEP::pi);
        
        // Cone is substracted from frontal volume 
        G4double lConeEllipse_y = mEllipseConeTransition_y - mEllipsePos_y - lConeHalfHeight;
        G4SubtractionSolid* lBulbSolidSubstractions = new G4SubtractionSolid("Substracted solid bulb", lBulbSolid, 
                                                                             lCone, 0, G4ThreeVector(0,0,lConeEllipse_y));
        
        // Creating Torus
        G4Torus* lTorus = new G4Torus("Solid substraction torus", 0.0, mTorusCircleR, mTorusCirclePos_x, 0, 2*CLHEP::pi);
        G4double lTorusToEllipse = mTorusCirclePos_y - mEllipsePos_y;
        G4Tubs* lTubeEdge = new G4Tubs("Solid edge of torus", mTorusCirclePos_x, mEllipseConeTransition_x+mTubeWidth, mTorusCirclePos_x*0.5, 0, 2*CLHEP::pi);
        
        G4UnionSolid* lTorusTubeEdge = new G4UnionSolid("Solid torus with cylindrical edges", lTorus, lTubeEdge, 0, G4ThreeVector(0,0, 0));
        
        // Create Tube for substracting cone and torus
        G4double lSubstractionTubeLength = mEllipseConeTransition_y - mTorusTubeTransition_y;
        G4Tubs* lSubstractionTube = new G4Tubs("substracion_tube", 0.0, mEllipseConeTransition_x, 0.5 * lSubstractionTubeLength,0,2*CLHEP::pi);
        
        G4double lSTubeEllipse_y = mEllipseConeTransition_y - mEllipsePos_y - lSubstractionTubeLength*0.5;
        
        G4SubtractionSolid* lBulbBack = new G4SubtractionSolid("Solid back of PMT", lSubstractionTube, lCone,0, G4ThreeVector(0,0,lConeEllipse_y-lSTubeEllipse_y)); 
        lBulbBack = new G4SubtractionSolid("Solid back of PMT", lBulbBack, lTorusTubeEdge,0, G4ThreeVector(0,0,lTorusToEllipse-lSTubeEllipse_y)); 
        
        //     G4UnionSolid* lReflectiveVolume = new G4UnionSolid("photocathode", lReflectiveSide, lBulbBack, 0, G4ThreeVector(0,0, lSTubeEllipse_y));
        
        // For some parameters this shape is too complicated for the visualizer, so I don't use the boolean in case of visualization
        // if you wanna see the PMT, negate this condition and use RayTracerX
        if (gVisual){
            lBulbSolid = new G4UnionSolid("Bulb tube solid", lBulbSolidSubstractions, lBulbBack, 0,G4ThreeVector(0,0, lSTubeEllipse_y) );
            lBulbSolid = new G4UnionSolid("Bulb tube solid", lBulbSolid, lBulkSolid, 0, G4ThreeVector(0, 0, -lMissingTubeLength));
        }
        else {
            lBulbSolid = new G4UnionSolid("Bulb tube solid", lBulbSolid, lBulkSolid, 0, G4ThreeVector(0, 0, -lMissingTubeLength));
            
        }
        
        
    }
    else {
        lBulbSolid = new G4UnionSolid("Bulb tube solid", lBulbSolid, lBulkSolid, 0, G4ThreeVector(0,0, -lMissingTubeLength));
        
    }
    
    return std::make_tuple(lBulbSolid, lPhotocathodeSide);
    
}

/**
 * 
 */
void OMSimPMTConstruction::SelectPMT(G4String pPMTtoSelect) {
    if (pPMTtoSelect.substr (0,6) == "argPMT"){
        G4String lPMTTypes[] = {"pmt_Hamamatsu_R15458", "pmt_ETEL_9320KFL-KFB", "pmt_HZC_XP82B2F", "pmt_Hamamatsu_4inch"};
        SelectPMT(lPMTTypes[gPMT]);
    }
    else {
        G4int lFound = mData->mPMTdata->mTable.count(pPMTtoSelect);
        if (lFound > 0) {
            mSelectedPMT = pPMTtoSelect;
            mData->mPMTdata->mSelectedKey = pPMTtoSelect;
            G4cout << pPMTtoSelect << " selected." << G4endl;
            ConstructIt();
        }
        else{
            G4cerr <<  "Selected PMT not in PMT tree, please check that requested PMT exists in data folder." <<  G4endl; 
        }
    }
}

