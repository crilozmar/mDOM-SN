    Glass_logical = new G4LogicalVolume (Glass_solid, Mat_Vessel_Glass, "Glass_log");
	G4int numberofstrings = 5;
	//Placing mDOMs
			G4double xpos[5] = {0*m,10*m,0*m,0*m,-10*m};
		G4double ypos[5] = {0*m,0*m,10*m,-10*m,0*m};
	if (numberofstrings > 1) {

		if (gn_mDOMs <= 1) {
			Glass_physical_vector.resize(numberofstrings);
		}else{
			Glass_physical_vector.resize(gn_mDOMs+numberofstrings);
		}
	}
	for (unsigned int j = 0; j < numberofstrings; j++){
	if (gn_mDOMs <= 1) {
		std::stringstream moduleconverter;
			moduleconverter.str("");
		moduleconverter << "Glass_phys_" << j << "_0";
		Glass_physical[j] = new G4PVPlacement (flipOM, G4ThreeVector(xpos[j],ypos[j],0), Glass_logical, moduleconverter.str(), World_logical, false, 0);
		PlacingHarnessAndRopes(0*m, rot, ropeThickness, ropeLength);
	} else {
		std::stringstream moduleconverter;
		for (unsigned int k = 0; k < gn_mDOMs; k++){
			moduleconverter.str("");
			moduleconverter << "Glass_phys_" << j << "_" <<k;
			G4double zpos;
			if (gn_mDOMs % 2 == 0) {
				zpos = gmdomseparation*(gn_mDOMs/2-k-1./2.);
			} else {
				zpos = gmdomseparation*(gn_mDOMs/2-k);
			}
			Glass_physical_vector[k+j] = new G4PVPlacement (flipOM, G4ThreeVector(xpos[j],ypos[j],zpos), Glass_logical, moduleconverter.str(), World_logical, false, 0);
			PlacingHarnessAndRopes(zpos, rot, ropeThickness, ropeLength);
		}
	}
	}
    