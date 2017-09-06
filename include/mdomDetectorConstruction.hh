#ifndef mdomDetectorConstruction_h
#define mdomDetectorConstruction_h 1


#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

class G4Box;
class G4Sphere;
class G4Orb;
class G4Polycone;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Tubs;
class G4Cons;
class G4Ellipsoid;

class mdomDetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		mdomDetectorConstruction();
		~mdomDetectorConstruction();
	
	public:
		G4VPhysicalVolume* Construct();

	private:
		G4double Mie_Scattering(int u, int depth_pos);
		G4double Spice_Absorption(int u, int depth_pos);
		G4double Spice_Temperature(G4double depth);
		G4double GetWaveLenght(float energy);
		G4double Spice_Refraction(int u);

// 		double z_offset(int i);
		G4Tubs*					World_solid;
		G4LogicalVolume*		World_logical;
		G4VPhysicalVolume*		World_physical;

		G4Sphere*				GelSphereTop_solid;
		G4Sphere*				GelSphereBottom_solid;
		G4Tubs*					GelCylinder_solid;
		G4UnionSolid*			Gel_solid;
		G4LogicalVolume*		Gel_logical;
		G4VPhysicalVolume*		Gel_physical;

		G4Sphere*				GlassSphereTop_solid;
		G4Sphere*				GlassSphereBottom_solid;
		G4Tubs*					GlassCylinder_solid;
		G4UnionSolid*			Glass_solid;
		G4LogicalVolume*		Glass_logical;
		G4VPhysicalVolume*		Glass_physical;
		
		G4SubtractionSolid*		TubeHolder_solid;	
		G4LogicalVolume*		TubeHolder_logical;
		G4VPhysicalVolume*		TubeHolder_physical;		
	
//	Hamamatsu R12199 (spherical window)
		G4Sphere* 				PMT_12199_sphere_solid;
		G4Ellipsoid*			PMT_12199_ellips_solid;
		G4Tubs*					PMT_12199_bulk_solid;
		G4UnionSolid*			PMT_12199_tube_solid;
		G4LogicalVolume*		PMT_12199_tube_logical;
//		G4Ellipsoid* 			PC_12199_sphere_solid; // possible alternative, ellipse directly allows z cuts with no boolean operations
		G4Sphere*	 			PC_12199_sphere_solid;
		G4Ellipsoid*			PC_12199_ellips_solid;
		G4UnionSolid*	 		PC_12199_intermed_solid;
		G4Box*					PC_12199_aux_solid;
		G4SubtractionSolid*		PC_12199_solid;
		G4LogicalVolume*		PC_12199_logical;
		G4VPhysicalVolume*		PC_12199_physical;
		G4Tubs*					PC_12199_shield_solid;
		G4LogicalVolume*		PC_12199_shield_logical;
		G4VPhysicalVolume*		PC_12199_shield_physical;
		G4Cons*					RefCone_12199_solid;
		G4LogicalVolume*		RefCone_12199_logical;
		
//	ETEL 9320KFL (ellipsoidal window OK)		
		G4Ellipsoid*			PMT_ETEL_bulb_solid;
		G4Tubs*					PMT_ETEL_bulk_solid;
		G4UnionSolid*			PMT_ETEL_tube_solid;
		G4LogicalVolume*		PMT_ETEL_tube_logical;
		G4Ellipsoid* 			PC_ETEL_solid;
		G4LogicalVolume*		PC_ETEL_logical;
		G4VPhysicalVolume*		PC_ETEL_physical;
		G4Tubs*					PC_ETEL_shield_solid;
		G4LogicalVolume*		PC_ETEL_shield_logical;
		G4VPhysicalVolume*		PC_ETEL_shield_physical;
		G4Cons*					RefCone_ETEL_solid;
		G4LogicalVolume*		RefCone_ETEL_logical;
		G4LogicalVolume*		RefConeType2_ETEL_logical;
		G4LogicalVolume*		RefConeType3_ETEL_logical;
		
//	Hamamatsu R12199 (ellipsoidal window, deviation from real surface ca. 1 mm)
		G4Ellipsoid*			PMT80_solid_1;
		G4Tubs*					PMT80_solid_2;
		G4UnionSolid*			PMT80_tube_solid;
		G4LogicalVolume*		PMT80_tube_logical;
		G4Ellipsoid* 			PC80_solid;
		G4LogicalVolume*		PC80_logical;
		G4VPhysicalVolume*		PC80_physical;
		G4Tubs*					PC80_shield_solid;
		G4LogicalVolume*		PC80_shield_logical;
		G4VPhysicalVolume*		PC80_shield_physical;
		
// Hamamatsu R7081 PMT 	
		G4Sphere*				PMT10_solid_1;
		G4Polycone*				PMT10_solid_2;
		G4Ellipsoid*			PMT10_solid_3;
		G4Ellipsoid*			PMT10_solid_4;
		G4UnionSolid*			PMT10_front_solid;
		G4UnionSolid*			PMT10_back_solid;
		G4UnionSolid*			PMT10_tube_solid;
		G4LogicalVolume*		PMT10_tube_logical;
		
		G4Sphere*				PC10_solid;
		G4LogicalVolume*		PC10_logical;
		G4VPhysicalVolume*		PC10_physical;
		G4Tubs*					PC10_shield_solid;
		G4LogicalVolume*		PC10_shield_logical;
		G4VPhysicalVolume*		PC10_shield_physical;
		
// PDOM parts
		G4Tubs*					PDOM_Board_solid;
		G4LogicalVolume*		PDOM_Board_logical;
		G4VPhysicalVolume*		PDOM_Board_physical;
		
		G4Orb*					PDOM_GlassSphere_solid;
		G4LogicalVolume*		PDOM_Glass_logical;
		G4VPhysicalVolume*		PDOM_Glass_physical;
		
		G4Orb*					PDOM_GelSphere_solid;
		G4LogicalVolume*		PDOM_Gel_logical;
		G4VPhysicalVolume*		PDOM_Gel_physical;
		
		G4Sphere*				PDOM_AirAux_solid;
		G4SubtractionSolid*		PDOM_Air_solid;
		G4LogicalVolume*		PDOM_Air_logical;
		G4VPhysicalVolume*		PDOM_Air_physical;
		
		G4Polycone*				PDOM_HarnessAux_solid;
		G4SubtractionSolid*		PDOM_Harness_solid;
		G4LogicalVolume*		PDOM_Harness_logical;
		G4VPhysicalVolume*		PDOM_Harness_physical;

// actually placed PMTs
		G4VPhysicalVolume*		PMT_test_physical;
		G4VPhysicalVolume*		PMT_physical[24];

// reflective cones
		G4Cons*					RefConeNestCone_solid;
		G4UnionSolid*			RefConeNest_solid;
		G4Cons*					RefConeBasic_solid;
		G4IntersectionSolid*	RefConeType1_solid;
		G4IntersectionSolid*	RefConeType2_solid;
		G4IntersectionSolid*	RefConeType3_solid;
        G4LogicalVolume*		RefConeNest_logical;
        G4LogicalVolume*		RefConeType1_logical;
		G4LogicalVolume*		RefConeType2_logical;
		G4LogicalVolume*		RefConeType3_logical;
		G4VPhysicalVolume*		RefCone_physical[24];

		
};

#endif
