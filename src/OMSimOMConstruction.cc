/** @file OMSimOMConstruction.cc
 *  @brief Construction of the optical module.
 *
 *  This class creates the solids of the optical modules and place them in the detector.
 * 
 *  @author Lew Classen, Martin Unland
 *  @date March 2021
 * 
 *  @version Geant4 10.7
 * 
 * 
 *   @todo  -Harness of the mDOM<br>
 *          -Correct PMT positioning for the LOM<br> 
 *          -Cutting the reflectors MultiPMTopticalModule() when overlapping<br>
 *          -Check material of lHarnessSurface opticalsurface<br>
 */

#include "OMSimOMConstruction.hh"
#include <dirent.h>

#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "G4Polycone.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include <G4UnitsTable.hh>
#include "G4VisAttributes.hh"

extern G4double gRefCone_angle;
extern G4int gDOM;

// Hash and mix are needed for switch statemens with strings
uint64_t constexpr mix(char m, uint64_t s)
{
    return ((s<<7) + ~(s>>3)) + ~m;
}

uint64_t constexpr hash(const char * m)
{
    return (*m) ? mix(*m,hash(m+1)) : 0;
}



OMSimOMConstruction::OMSimOMConstruction(OMSimInputData* pData){
    mData = pData;
    mPMTManager =  new OMSimPMTConstruction(pData);
}

/**
 * Select the module with string input and call the construction function. Add a case in the switch if you want another module.
 * @param pSelectedModule string name of the module to be selected. Currently only [pDOM, mDOM, LOM].
 */
void OMSimOMConstruction::SelectModule(G4String pSelectedModule){
    
    G4String lOM[] = {"mDOM", "pDOM"};
    if (pSelectedModule=="argOM"){
        pSelectedModule = lOM[gDOM];
    mSelectedModule = pSelectedModule;}
    
    switch(hash(pSelectedModule)){
        case hash("pDOM")       : pDOMConstruction();
        break;
        case hash("mDOM")       : 
            mSelectedModule = "om_mDOM";
            mData->mOMdata->mSelectedKey = mSelectedModule;
            MultiPMTopticalModule();
            break;
        case hash("LOM")       : 
            mSelectedModule = "om_LOM";
            mData->mOMdata->mSelectedKey = mSelectedModule;
            MultiPMTopticalModule();
            break;
        default : G4cerr << "Requested Module does not exist!" << G4endl;
    }
}


/**
 * The physical volume of the constructed module is produced here. Also the harness is constructed and placed if requested. In the case of the mDOM, the harness is still not in the code, but its construction/placement should probably be called from here!
 * @param pPosition G4ThreeVector with position of the module (as in G4PVPlacement())
 * @param pRotation G4RotationMatrix with rotation of the module (as in G4PVPlacement())
 * @param pMother G4LogicalVolume where the module is going to be placed (as in G4PVPlacement())
 * @param pIncludeHarness bool Harness is placed if true
 * @param pNameExtension G4String name of the physical volume. You should not have two physicals with the same name
 */
void OMSimOMConstruction::PlaceIt(G4ThreeVector pPosition, G4RotationMatrix* pRotation, G4LogicalVolume*& pMother, bool pIncludeHarness, G4String pNameExtension ){
    
    //Placement of harness
    if(pIncludeHarness){
        switch(hash(mSelectedModule)){
            case hash("pDOM")       : new G4PVPlacement (pRotation, pPosition, mHarnessLogical, "HarnessPhys"+pNameExtension, pMother, false, 0, mCheckOverlaps);
            break;
            case hash("mDOM")       : G4cout << "mDOM Harness to be implemented" << G4endl;
            break;
            default : G4cerr << "Select module before placing it!" << G4endl;
        }
        
    }
    
    //Placement of module
    new G4PVPlacement (pRotation, pPosition, mGlassLogical, "OMPhys"+pNameExtension, pMother, false, 0, mCheckOverlaps);
}

/**
 * Construction of the pDOM, as in the older code.
 */
void OMSimOMConstruction::pDOMConstruction(){
    
    G4UnionSolid* lPMTsolid;
    
    mPMTManager->SelectPMT("pmt_Hamamatsu_R7081");
    
    lPMTsolid = mPMTManager->GetPMTSolid();
    
    G4double lGelThickness = 10*mm;
    
    G4double lPMTz = 0.5*12*25.4*mm - mPMTManager->mPMTCenterToTip - lGelThickness;
    
    G4Orb* lGlassSphereSolid = new G4Orb("PDOM_GlassSphere solid",0.5*13*25.4*mm);
    G4Orb* lGelSphereSolid = new G4Orb("PDOM_GelSphere solid",0.5*12*25.4*mm);
    
    G4Ellipsoid* lAirAuxSolid = new G4Ellipsoid("PDOM_AirAux solid", 0.5*12*25.4*mm, 0.5*12*25.4*mm, 0.5*12*25.4*mm, -0.5*13*25.4*mm, 50*mm);
    
    G4SubtractionSolid* lAirSolid = new G4SubtractionSolid("PDOM_Air solid", lAirAuxSolid, lPMTsolid, 0, G4ThreeVector(0,0, lPMTz));
    G4Tubs* lBoardSolid = new G4Tubs("PDOM_Board solid", 52*mm, 0.5*11*25.4*mm, 2*mm, 0, 2*CLHEP::pi);
    G4Tubs* lBaseSolid = new G4Tubs("PDOM_Board solid", 0*mm, 6*cm, 2*mm, 0, 2*CLHEP::pi);
    
    //Harness
    G4double lHarnessInner[] = {0, 0, 0, 0};
    G4double lHarnessRadii[] = {(0.5*365.76-8.3)*mm, 0.5*365.76*mm, 0.5*365.76*mm, (0.5*365.76-8.3)*mm};
    G4double lHarnessZplanes[] = {-31.75*mm, -10*mm, 10*mm, 31.75*mm};	
    G4Polycone* lHarnessAuxSolid =  new G4Polycone("PDOM_HarnessAux solid", 0, 2*CLHEP::pi, 4, lHarnessZplanes, lHarnessInner, lHarnessRadii);
    G4SubtractionSolid* lHarnessSolid = new G4SubtractionSolid("PDOM_Harness solid", lHarnessAuxSolid, lGlassSphereSolid);
    
    
    //Logicals mData
    mGlassLogical     = new G4LogicalVolume(lGlassSphereSolid,
                                            mData->GetMaterial("argVesselGlass"), 
                                            "PDOM_Glass logical");
    mHarnessLogical   = new G4LogicalVolume(lHarnessSolid,
                                            mData->GetMaterial("NoOptic_Stahl"),
                                            "PDOM_Harness logical");
    
    G4LogicalVolume* lGelLogical       = new G4LogicalVolume(lGelSphereSolid,
                                                             mData->GetMaterial("argGel"), 
                                                             "PDOM_Gel logical");
    
    G4LogicalVolume* lBoardLogical     = new G4LogicalVolume(lBoardSolid,
                                                             mData->GetMaterial("NoOptic_Absorber"), 
                                                             "PDOM_Board logical");
    
    G4LogicalVolume* lBaseLogical      = new G4LogicalVolume(lBaseSolid,
                                                             mData->GetMaterial("NoOptic_Absorber"), 
                                                             "PDOM_Base logical");
    
    G4LogicalVolume* lAirLogical       = new G4LogicalVolume(lAirSolid, 
                                                             mData->GetMaterial("Ri_Vacuum"), 
                                                             "PDOM_Air logical");	
    
    G4PVPlacement* lBoardPhysical = new G4PVPlacement (0, G4ThreeVector(0,0,-40*mm), lBoardLogical, "pDOMBoardPhys", lAirLogical, false, 0);
    G4PVPlacement* lBasePhysical = new G4PVPlacement (0, G4ThreeVector(0,0,-105*mm), lBaseLogical, "pDOMBasePhys", lAirLogical, false, 0);
    
    G4PVPlacement* lAirPhysical = new G4PVPlacement (0, G4ThreeVector(0,0,0), lAirLogical, "pDOMAirPhys ", lGelLogical, false, 0);
    
    mPMTManager->PlaceIt(G4ThreeVector(0,0, lPMTz), new G4RotationMatrix(), lGelLogical);
    
    G4PVPlacement* lGelPhysical = new G4PVPlacement (0, G4ThreeVector(0,0,0), lGelLogical, "pDOMGelPhys", mGlassLogical, false, 0);
    
    // ------------------- optical border surfaces --------------------------------------------------------------------------------
    G4LogicalSkinSurface* lHarnessSurface = new G4LogicalSkinSurface("PDOM_Harness_skin", mHarnessLogical, mData->GetOpticalSurface("Refl_V95Gel"));  
    
    mGlassLogical->SetVisAttributes(mGlassVis);
    mHarnessLogical->SetVisAttributes(mSteelVis);
    lGelLogical->SetVisAttributes(mGelVis);
    lAirLogical->SetVisAttributes(mAirVis);
    lBoardLogical->SetVisAttributes(mBoardVis);
    lBaseLogical->SetVisAttributes(mBoardVis);
}

/**
 * Construction of a segmented module, as in the older code.
 */
void OMSimOMConstruction::MultiPMTopticalModule(){
    
    G4Transform3D lTransformers;
    G4double lGlassOutRad = mData->mOMdata->Get("jGlassOutRad"); // outer radius of galss cylinder (pressure vessel)
    G4double lGlassThick = mData->mOMdata->Get("jGlassThick"); // maximum Glass thickness
    G4double lCylHigh = mData->mOMdata->Get("jCylHigh"); // height of cylindrical part of glass half-vessel
    G4double lGelThicknessFrontPMT = mData->mOMdata->Get("jGelThicknessFrontPMT"); // distance between inner glass surface and tip of PMTs
    G4double lGelThickness = mData->mOMdata->Get("jGelThickness");// distance between inner glass surface and holding structure, filled with gel
    G4double lCylinderAngle =  mData->mOMdata->Get("lCylinderAngle");// Deviation angle of cylindrical part of the pressure vessel
    
    G4double lEqPMTrOffset = mData->mOMdata->Get("jEqPMTrOffset");// middle PMT circles are slightly further out due to lEqPMTzOffset
    G4double lEqPMTzOffset = mData->mOMdata->Get("jEqPMTzOffset");// z-offset of middle PMT circles w.r.t. center of glass sphere
    G4double lRefConeHalfZ = mData->mOMdata->Get("jRefConeHalfZ"); // half-height of reflector (before cutting to right form)
    
    G4double lRefConeSheetThickness = mData->mOMdata->Get("lRefConeSheetThickness"); // aluminum sheet thickness true for all reflective cones
    G4double lRefConeToHolder = mData->mOMdata->Get("jRefConeToHolder"); // horizontal distance from K??rcher's construction 
    
    G4double lThetaPolar = mData->mOMdata->Get("jThetaPolar");
    G4double lThetaEquatorial = mData->mOMdata->Get("jThetaEquatorial");
    G4int lNrPolarPMTs = mData->mOMdata->mTable.at(mSelectedModule).get<G4int>("jNrPolarPMTs");
    G4int lNrEqPMTs = mData->mOMdata->mTable.at(mSelectedModule).get<G4int>("jNrEqPMTs");
    G4double lPolEqPMTPhiPhase = mData->mOMdata->Get("jPolEqPMTPhiPhase");
    

    G4double lGlassInRad   = lGlassOutRad - lGlassThick;
    G4double lRefConeAngle = gRefCone_angle * deg;
    G4int lTotalNrPMTs = (lNrPolarPMTs+lNrEqPMTs)*2;

    mPMTManager->SelectPMT("argPMT");
    //mPMTManager->SimulateInternalReflections();
    
    G4UnionSolid* lPMTsolid = mPMTManager->GetPMTSolid();
    G4double lPMToffset = mPMTManager->mPMTCenterToTip;
    G4double lPMTBulbRad = mPMTManager->mPMTMaxRad;
    
    //	Glass
    G4Ellipsoid* lGlassTopSolid = new G4Ellipsoid("GlassSphereTop solid", lGlassOutRad, lGlassOutRad, lGlassOutRad, -5*mm, lGlassOutRad+5*mm);
    G4Ellipsoid* lGlassBottomSolid = new G4Ellipsoid("GlassSphereBottom solid", lGlassOutRad, lGlassOutRad, lGlassOutRad, -(lGlassOutRad+5*mm), 5*mm);
    
    G4double zCorners[] = {-lCylHigh*1.001,-lCylHigh, 0, lCylHigh, lCylHigh*1.001};
    G4double rCorners[] = {0, lGlassOutRad, lGlassOutRad+lCylHigh*sin(lCylinderAngle), lGlassOutRad, 0 };
    G4Polycone* lGlassCylinderSolid = new G4Polycone("GlassCylinder solid", 0, 2*CLHEP::pi,  5, rCorners, zCorners);
    
    G4UnionSolid* lTempUnion = new G4UnionSolid("temp", lGlassCylinderSolid, lGlassTopSolid, 0,G4ThreeVector(0, 0, lCylHigh) );
    G4UnionSolid* lGlassSolid = new G4UnionSolid("OM glass body", lTempUnion, lGlassBottomSolid, 0, G4ThreeVector(0,0,-lCylHigh));
    
    //  Gel
    G4Ellipsoid* lGelTopSolid = new G4Ellipsoid("GelSphereTop solid", lGlassInRad, lGlassInRad, lGlassInRad, -5*mm, lGlassInRad+5*mm);
    G4Ellipsoid* lGelBottomSolid = new G4Ellipsoid("GelSphereBottom solid", lGlassInRad, lGlassInRad, lGlassInRad, -(lGlassInRad+5*mm), 5*mm);
    
    G4double rCornersGel[] = {0, lGlassOutRad- lGlassThick, lGlassOutRad-lGlassThick+lCylHigh*sin(lCylinderAngle), lGlassOutRad-lGlassThick, 0 };
    G4Polycone* lGelCylinderSolid = new G4Polycone("GelCylinder solid", 0, 2*CLHEP::pi,  5, rCornersGel, zCorners);
    
    G4UnionSolid* lTempUnion2 = new G4UnionSolid("temp2", lGelCylinderSolid, lGelTopSolid, 0, G4ThreeVector(0,0,lCylHigh ));
    G4UnionSolid* lGelSolid = new G4UnionSolid("gel body", lTempUnion2, lGelBottomSolid, 0, G4ThreeVector(0,0,-lCylHigh));
    
    //  PMT support structure primitives & cutting "nests" for PMTs later
    G4double lSupStructureRad = lGlassOutRad - lGlassThick - lGelThickness;
    G4Ellipsoid* lSupStructureTopSolid = new G4Ellipsoid("FoamSphereTop solid", lSupStructureRad, lSupStructureRad, lSupStructureRad, -5*mm, lSupStructureRad+5*mm);
    G4Ellipsoid* lSupStructureBottomSolid = new G4Ellipsoid("FoamSphereBottom solid", lSupStructureRad, lSupStructureRad, lSupStructureRad, -(lSupStructureRad+5*mm), 5*mm);		
    G4Tubs* lSupStructureCylinderSolid = new G4Tubs("FoamCylinder solid", 0, lGlassOutRad - lGlassThick - lGelThickness, lCylHigh , 0, 2*CLHEP::pi);
    
    G4UnionSolid* lSupStructureTempUnion = new G4UnionSolid("Foam TempUnion solid", lSupStructureCylinderSolid, lSupStructureTopSolid, 0, G4ThreeVector(0,0,(lCylHigh )));
    G4UnionSolid* lSupStructureFirstSolid = new G4UnionSolid("Foam solid", lSupStructureTempUnion, lSupStructureBottomSolid, 0, G4ThreeVector(0,0,-(lCylHigh )));
    
    double lPMTtheta[99], lPMTphi[99], lPMTx[99], lPMTy[99], lPMTz[99], lRefConeX[99], lRefConeY[99], lRefConeZ[99]; 
    G4double lPMTr = lGlassInRad - lGelThicknessFrontPMT - lPMToffset; // radius for PMT positioning
    G4double lRefConeR = lGlassInRad - lGelThicknessFrontPMT - lPMToffset + lRefConeHalfZ; // radius for RefCone positioning
    G4double lPMTrho;
    G4double lRefConeRho;
    G4double lPMTzOffset;
    
    
    for (int i = 0; i <= lTotalNrPMTs-1; i++) {
        lPMTr = lGlassInRad - lGelThicknessFrontPMT - lPMToffset;
        lRefConeR = lGlassInRad - lGelThicknessFrontPMT - lPMToffset + lRefConeHalfZ;
        
        if (i>=0 && i<=lNrPolarPMTs-1){
            lPMTtheta[i]=lThetaPolar;
            lPMTphi[i]=(lPolEqPMTPhiPhase+i*360.*deg/lNrPolarPMTs);
            lPMTzOffset = lCylHigh;
        }
        if (i>=lNrPolarPMTs && i<=lNrPolarPMTs+lNrEqPMTs-1){
            lPMTtheta[i]=lThetaEquatorial;
            lPMTphi[i]=(i-4)*360.*deg/lNrEqPMTs;
            lPMTzOffset = lCylHigh - lEqPMTzOffset;
            lPMTr += lEqPMTrOffset;
            lRefConeR += lEqPMTrOffset;
        }
        if (i>=lNrPolarPMTs+lNrEqPMTs && i<=lNrPolarPMTs+2*lNrEqPMTs-1){
            lPMTtheta[i]=180.*deg-lThetaEquatorial;
            lPMTphi[i]=(i-12)*360.*deg/lNrEqPMTs;
            lPMTzOffset = - lCylHigh + lEqPMTzOffset;
            lPMTr += lEqPMTrOffset;
            lRefConeR += lEqPMTrOffset;
        }
        if (i>=lNrPolarPMTs+2*lNrEqPMTs && i<=lTotalNrPMTs-1){
            lPMTtheta[i]=180.*deg-lThetaPolar;
            lPMTphi[i]=(lPolEqPMTPhiPhase+(i-20)*360.*deg/lNrPolarPMTs);
            lPMTzOffset = - lCylHigh;
        }
        
        lPMTrho = lPMTr * sin(lPMTtheta[i]);
        lPMTx[i] = lPMTrho * cos(lPMTphi[i]);
        lPMTy[i] = lPMTrho * sin(lPMTphi[i]);
        lPMTz[i] = lPMTr * cos(lPMTtheta[i]) + lPMTzOffset;
        lRefConeRho = lRefConeR * sin(lPMTtheta[i]);
        lRefConeX[i] = lRefConeRho * cos(lPMTphi[i]);
        lRefConeY[i] = lRefConeRho * sin(lPMTphi[i]);
        lRefConeZ[i] = lRefConeR * cos(lPMTtheta[i]) + lPMTzOffset;
        
        
    }	
    
    //Reflectors nests for support structure substraction
    G4double lRefConeIdealInRad = lPMTBulbRad + 2*mm;
    G4double lRefConeOuterNegRad = lRefConeIdealInRad + lRefConeToHolder / cos(lRefConeAngle);
    G4double lRefConeOuterPosRad = lRefConeIdealInRad + lRefConeToHolder / cos(lRefConeAngle) + 2*1.5*lRefConeHalfZ*tan(lRefConeAngle);
    G4Cons* lRefConeNestConeSolid = new G4Cons("RefConeNestCone", 0, lRefConeOuterNegRad, 0, lRefConeOuterPosRad, 1.5*lRefConeHalfZ, 0, 2*CLHEP::pi);
    G4UnionSolid* lRefConeNestSolid = new G4UnionSolid("RefConeNest", lPMTsolid, lRefConeNestConeSolid, 0, G4ThreeVector(0,0,1.5*lRefConeHalfZ));
    
    //Support structure substraction
    G4RotationMatrix* lRot = new G4RotationMatrix();
    G4SubtractionSolid* lSupStructureSolid;
    for (int k = 0; k <=lTotalNrPMTs-1; k++) { 
        lRot= new G4RotationMatrix();
        lRot->rotateY(lPMTtheta[k]);
        lRot->rotateZ(lPMTphi[k]);
        lTransformers = G4Transform3D(*lRot, G4ThreeVector(lPMTx[k],lPMTy[k],lPMTz[k]));
        if (k==0){
            lSupStructureSolid = new G4SubtractionSolid("TubeHolder solid", lSupStructureFirstSolid, lRefConeNestSolid, lTransformers );
        } 		
        else { 
            lSupStructureSolid = new G4SubtractionSolid("TubeHolder solid", lSupStructureSolid, lRefConeNestSolid, lTransformers);
        }
    }
    
    //Reflector solid
    G4Cons* lRefConeBasicSolid = new G4Cons("RefConeBasic", lRefConeIdealInRad,
                                            lRefConeIdealInRad + lRefConeSheetThickness / cos(lRefConeAngle), 
                                            lRefConeIdealInRad + 2*lRefConeHalfZ*tan(lRefConeAngle), lRefConeIdealInRad + lRefConeSheetThickness / cos(lRefConeAngle) + 2*lRefConeHalfZ*tan(lRefConeAngle), lRefConeHalfZ, 0, 2*CLHEP::pi);									
    
    //Reflective cone for polar PMTs
    G4IntersectionSolid* lRefConePolarSolid = new G4IntersectionSolid("PolarRefCones", lRefConeBasicSolid, lSupStructureTopSolid, 0, G4ThreeVector(0,0,-(lGlassInRad-lGelThicknessFrontPMT-lPMToffset+lRefConeHalfZ)));
    
    //Reflective upper equatorial PMTs
    G4double lRefConeEqUR = lGlassInRad - lGelThicknessFrontPMT - lPMToffset + lRefConeHalfZ + lEqPMTrOffset;
    G4double lRefConeEqURho = lRefConeEqUR * sin(lThetaEquatorial);
    G4double lRefConeEqUX = lRefConeEqURho * cos(0*deg);
    G4double lRefConeEqUY = lRefConeEqURho * sin(0*deg);
    G4double lRefConeEqUZ = lRefConeEqUR * cos(lThetaEquatorial) + lCylHigh - lEqPMTzOffset;
    
    lRot = new G4RotationMatrix();
    lRot->rotateY(lThetaEquatorial);
    lTransformers = G4Transform3D(*lRot, G4ThreeVector(lRefConeEqUX, lRefConeEqUY, lRefConeEqUZ));
    G4IntersectionSolid* lRefConeEqUpSolid = new G4IntersectionSolid("RefConeType2", lSupStructureFirstSolid, lRefConeBasicSolid, lTransformers );
    
    //Reflective lower equatorial PMTs
    G4double lRefConeEqLR = lGlassInRad - lGelThicknessFrontPMT - lPMToffset + lRefConeHalfZ + lEqPMTrOffset;
    G4double lRefConeEqLRho = lRefConeEqLR * sin(180.*deg-lThetaEquatorial);
    G4double lRefConeEqLX = lRefConeEqLRho * cos(0*deg);
    G4double lRefConeEqLY = lRefConeEqLRho * sin(0*deg);
    G4double lRefConeEqLZ = lRefConeEqLR * cos(180.*deg-lThetaEquatorial) - lCylHigh + lEqPMTzOffset;
    
    lRot = new G4RotationMatrix();
    lRot->rotateY(180.*deg-lThetaEquatorial);
    lTransformers = G4Transform3D(*lRot, G4ThreeVector(lRefConeEqLX, lRefConeEqLY, lRefConeEqLZ));
    G4IntersectionSolid* lRefConeEqLoSolid = new G4IntersectionSolid("RefConeType3", lSupStructureFirstSolid, lRefConeBasicSolid, lTransformers);
    
    
    G4Cons* lRefConCutterSolid = new G4Cons("RefConeCutter", 
                                            0, lRefConeIdealInRad + lRefConeSheetThickness*std::sqrt(2.) + 2*mm, 
                                            0, lRefConeIdealInRad + lRefConeSheetThickness*std::sqrt(2.) + 2*lRefConeHalfZ + 2*mm, 
                                            lRefConeHalfZ, 0, 2*CLHEP::pi);
    
    //Cropping upper equatorial cones
    G4double lRefConeCutR = lGlassInRad - lGelThicknessFrontPMT - lPMToffset + lRefConeHalfZ + lEqPMTrOffset;
    G4double lRefConeCutRho = lRefConeCutR * sin(lThetaEquatorial);
    G4double lRefConeCutX = lRefConeCutRho * cos(360.*deg/lNrEqPMTs);
    G4double lRefConeCutY = lRefConeCutRho * sin(360.*deg/lNrEqPMTs);
    G4double lRefConeCutZ = lRefConeCutR * cos(lThetaEquatorial) + lCylHigh - lEqPMTzOffset;
    
    lRot = new G4RotationMatrix();
    lRot->rotateY(lThetaEquatorial);
    lRot->rotateZ(360.*deg/lNrEqPMTs);
    lTransformers = G4Transform3D(*lRot, G4ThreeVector(lRefConeCutX, lRefConeCutY, lRefConeCutZ));
    
    G4SubtractionSolid* lRefconeEqUpCutSolid = new G4SubtractionSolid("RefConeEqUpCut", 
                                                                      lRefConeEqUpSolid,
                                                                      lRefConCutterSolid, 
                                                                      lTransformers);
    
    lRefConeCutX = lRefConeCutRho * cos(-360.*deg/lNrEqPMTs);
    lRefConeCutY = lRefConeCutRho * sin(-360.*deg/lNrEqPMTs);
    
    lRot = new G4RotationMatrix();
    lRot->rotateY(lThetaEquatorial);
    lRot->rotateZ(-360.*deg/lNrEqPMTs);
    lTransformers = G4Transform3D(*lRot, G4ThreeVector(lRefConeCutX, lRefConeCutY, lRefConeCutZ));
    
    lRefconeEqUpCutSolid = new G4SubtractionSolid("RefConeEqLoCut", lRefconeEqUpCutSolid, lRefConCutterSolid, lTransformers);
    
    //Cropping lower equatorial cones
    lRefConeCutRho = lRefConeCutR * sin(180.*deg-lThetaEquatorial);
    lRefConeCutX = lRefConeCutRho * cos(360.*deg/lNrEqPMTs);
    lRefConeCutY = lRefConeCutRho * sin(360.*deg/lNrEqPMTs);
    lRefConeCutZ = lRefConeCutR * cos(180.*deg-lThetaEquatorial) - lCylHigh + lEqPMTzOffset;
    
    lRot = new G4RotationMatrix();
    lRot->rotateY(180.*deg-lThetaEquatorial);
    lRot->rotateZ(360.*deg/lNrEqPMTs);
    lTransformers = G4Transform3D(*lRot, G4ThreeVector(lRefConeCutX, lRefConeCutY, lRefConeCutZ));
    
    G4SubtractionSolid* lRefconeEqLoCutSolid = new G4SubtractionSolid("RefConeType3 ETEL solid", 
                                                                      lRefConeEqLoSolid,
                                                                      lRefConCutterSolid,
                                                                      lTransformers);
    
    lRefConeCutX = lRefConeCutRho * cos(-360.*deg/lNrEqPMTs);
    lRefConeCutY = lRefConeCutRho * sin(-360.*deg/lNrEqPMTs);
    
    lRot = new G4RotationMatrix();
    lRot->rotateY(180.*deg-lThetaEquatorial);
    lRot->rotateZ(-360.*deg/lNrEqPMTs);
    lTransformers = G4Transform3D(*lRot, G4ThreeVector(lRefConeCutX, lRefConeCutY, lRefConeCutZ));
    
    lRefconeEqLoCutSolid = new G4SubtractionSolid("RefConeType3 ETEL solid", lRefconeEqLoCutSolid, lRefConCutterSolid, lTransformers);
    
    
    //Logicals
    lRot = new G4RotationMatrix();
    mGlassLogical = new G4LogicalVolume (lGlassSolid, 
                                         mData->GetMaterial("argVesselGlass"),
                                         "Glass_log");
    G4LogicalVolume* lGelLogical = new G4LogicalVolume (lGelSolid,
                                                        mData->GetMaterial("argGel"), 
                                                        "Gelcorpus logical");
    G4LogicalVolume* lSupStructureLogical = new G4LogicalVolume (lSupStructureSolid, 
                                                                 mData->GetMaterial("NoOptic_Absorber"), 
                                                                 "TubeHolder logical");
    
    G4LogicalVolume* lRefConePolarLogical = new G4LogicalVolume(lRefConePolarSolid,
                                                                mData->GetMaterial("NoOptic_Reflector"), 
                                                                "RefConeType1 logical");    
    
    G4LogicalVolume* lRefconeEqUpCutLogical = new G4LogicalVolume(lRefconeEqUpCutSolid,
                                                                  mData->GetMaterial("NoOptic_Reflector"),
                                                                  "RefConeType2 ETEL logical");
    
    G4LogicalVolume* lRefconeEqLoCutLogical = new G4LogicalVolume(lRefconeEqLoCutSolid,
                                                                  mData->GetMaterial("NoOptic_Reflector"),
                                                                  "RefConeType3 ETEL logical");
    
    
    //Placements
    new G4PVPlacement (0, G4ThreeVector(0,0,0), lSupStructureLogical, "SupportStructure_physical", lGelLogical, false, 0);        
    
    new G4PVPlacement (0, G4ThreeVector(0,0,0), lGelLogical, "Gel_physical", mGlassLogical, false, 0);
    
    
    std::stringstream lConverter;
    for (int k = 0; k <= lTotalNrPMTs-1; k++){
        lConverter.str("");
        lConverter << k << "_physical";
        
        lRot = new G4RotationMatrix();
        lRot->rotateY(lPMTtheta[k]);
        lRot->rotateZ(lPMTphi[k]);
        lTransformers = G4Transform3D(*lRot, G4ThreeVector(lPMTx[k],lPMTy[k],lPMTz[k]));
        mPMTManager->PlaceIt(lTransformers, lGelLogical, lConverter.str());
        
        //Placing reflective cones:
        lConverter.str("");
        lConverter << "RefCone_" << k << "_physical";
        lTransformers = G4Transform3D(*lRot, G4ThreeVector(lRefConeX[k],lRefConeY[k],lRefConeZ[k]));
        if (k >= lNrPolarPMTs && k <= lNrPolarPMTs+lNrEqPMTs-1) {
            lRot = new G4RotationMatrix();
            lRot->rotateZ(lPMTphi[k]);
            lTransformers = G4Transform3D(*lRot, G4ThreeVector(0,0,0));
            new G4PVPlacement (lTransformers, lRefconeEqUpCutLogical, lConverter.str(), lGelLogical, false, 0);
        }
        else if (k >= lNrPolarPMTs+lNrEqPMTs && k <= lNrPolarPMTs+lNrEqPMTs*2-1) {
            lRot = new G4RotationMatrix();
            lRot->rotateZ(lPMTphi[k]);
            lTransformers = G4Transform3D(*lRot, G4ThreeVector(0,0,0));
            new G4PVPlacement (lTransformers, lRefconeEqLoCutLogical, lConverter.str(), lGelLogical, false, 0);
        }
        else {
            new G4PVPlacement (lTransformers, lRefConePolarLogical, lConverter.str(), lGelLogical, false, 0);
        }
    }
    
    // ------------------- optical border surfaces --------------------------------------------------------------------------------
    G4LogicalSkinSurface* lRefConePolarSurface = new G4LogicalSkinSurface("RefCone_skin", lRefConePolarLogical,
                                                                          mData->GetOpticalSurface("argReflector"));
    G4LogicalSkinSurface* lRefconeEqUpCutSurface = new G4LogicalSkinSurface("RefCone_skin", lRefconeEqUpCutLogical, 
                                                                            mData->GetOpticalSurface("argReflector"));
    G4LogicalSkinSurface* lRefconeEqLoCutSurface = new G4LogicalSkinSurface("RefCone_skin", lRefconeEqLoCutLogical, 
                                                                            mData->GetOpticalSurface("argReflector"));
    
    //     // ---------------- visualisation attributes --------------------------------------------------------------------------------
    mGlassLogical->SetVisAttributes(mGlassVis);
    lGelLogical->SetVisAttributes(mGelVis);
    lSupStructureLogical->SetVisAttributes(mAbsorberVis);
    //lSupStructureLogical->SetVisAttributes(mInvisibleVis);
    lRefConePolarLogical->SetVisAttributes(mAluVis);
    lRefconeEqUpCutLogical->SetVisAttributes(mAluVis);
    lRefconeEqLoCutLogical->SetVisAttributes(mAluVis);
    
    
    
}




