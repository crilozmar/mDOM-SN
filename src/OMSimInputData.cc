/** @file OMSimInputData.cc
 *  @brief Input data from external files are read and saved to memory.
 *
 *  All materials for the modules, pmts and environment, as well as the shape parameters of the PMTs and Modules are in the ../data folder.
 *  This class read all files and parse them to tables.
 * 
 *  @author Martin Unland
 *  @date 22. July 2021
 * 
 *  @version Geant4 10.7
 * 
 *  @todo  -One could make OMSimInputData::MaterialFromFile more flexible and shorter by having an array of properties in the file and its units. This way one does not need the case separation with prefixes.
 */

#include <boost/property_tree/ptree.hpp>                                        
#include <boost/property_tree/json_parser.hpp> 

#include "OMSimInputData.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4SystemOfUnits.hh"
#include <G4UnitsTable.hh>
#include <dirent.h>
#include <cmath>


namespace pt = boost::property_tree;
extern G4int gEnvironment;
extern G4int gGlass;
extern G4int gGel;
extern G4int gConeMat;
extern G4int gDepthpos;

// Hash and mix are needed for switch statemens with strings
uint64_t constexpr mix(char pChar, uint64_t pUint)
{
    return ((pUint<<7) + ~(pUint>>3)) + ~pChar;
}

uint64_t constexpr hash(const char * pChar)
{
    return (*pChar) ? mix(*pChar,hash(pChar+1)) : 0;
}



OMSimInputData::OMSimInputData(){
    mSpiceDepth_pos = gDepthpos;
}


/*
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *                    Functions for retrieving materials and surfaces
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

/**
 * Get a G4Material. In order to get custom built materials, method SearchFolders() should have already been called.
 * Standard materials from Geant4 are also transfered directly (if found through FindOrBuildMaterial()).
 * Materials defined by arguments of the main should start with "arg" and then we get the material name from the 
 * arrays depending on the global integers (gGel, gGlass, gEnvironment..) .
 * @param pName name of the (custom) material or "argVesselGlass", "argGel", "argWorld" for argument materials
 * @return G4Material  
 */	
G4Material* OMSimInputData::GetMaterial(G4String pName){
    
    //Check if requested material is an argument material
    if (pName.substr (0,3) == "arg"){
        
        G4String lGlass[] = {"RiAbs_Glass_Vitrovex", "RiAbs_Glass_Chiba", "RiAbs_Glass_Benthos", "RiAbs_Glass_myVitrovex", "RiAbs_Glass_myChiba", "RiAbs_Glass_WOMQuartz"};
        G4String lGel[] = {"RiAbs_Gel_Wacker612Measured", "RiAbs_Gel_Shin-Etsu", "RiAbs_Gel_QGel900", "RiAbs_Gel_Wacker612Company", "Ri_Vacuum"};
        G4String lWorld[] = {"Ri_Air", "IceCubeICE", "IceCubeICE_SPICE"};
        
        switch(hash(pName)){
            case hash("argVesselGlass")  : return GetMaterial(lGlass[gGlass]); break;
            case hash("argGel")          : return GetMaterial(lGel[gGel]); break;
            case hash("argWorld")        : return GetMaterial(lWorld[gEnvironment]); break;
        }
    }
    
    //If it is not an argument material, the material is looked up.
    G4Material* lReturn = G4Material::GetMaterial(pName);
    if (!lReturn)
    {
        lReturn = G4NistManager::Instance()->FindOrBuildMaterial(pName);
    }
    return lReturn;
}

/**
 * Get a G4OpticalSurface. In order to get custom built materials, method SearchFolders() should have already been called. 
 * @param pName name of the optical surface or argument reflectors "argReflector"
 * @return G4OpticalSurface  
 */
G4OpticalSurface* OMSimInputData::GetOpticalSurface(G4String pName){
    
    //Check if requested material is an argument surface
    if (pName.substr (0,12) == "argReflector"){
        G4String lRefCones[] = {"Refl_V95Gel", "Refl_V98Gel", "Refl_Aluminium", "Refl_Total98"};
        return GetOpticalSurface(lRefCones[gConeMat]);
    }
    //If not, we look in the dictionary 
    else{
        G4int lFound = mOpticalSurfaceMap.count(pName);
        if (lFound > 0) {
            return mOpticalSurfaceMap.at(pName);
        }
        else{
            G4cerr <<  "Requested Optical Surface " << pName << " not found. This will cause a segmentation fault. Please check the name!!" <<  G4endl; 
        }
    }
}


/*
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *                        Functions for ice optical properties
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

/**
 * This gives you temperature of ice depending on the depth.
 * Function needed for the calculation of scattering and absorption length of the ice. 
 * @param pDepth Depth in m from where we need the temperature
 * @return Temperature
 */
G4double OMSimInputData::Spice_Temperature(G4double pDepth){
    G4double spice_temp = 221.5-0.00045319/m*pDepth+5.822e-6/m2*pow(pDepth, 2.);
    return spice_temp;
}

/**
 * Calculation of the absorption length of IC-ice for a specific wavelength
 * @param pLambd Wavelength
 * @return Absorption length
 */
G4double OMSimInputData::Spice_Absorption(G4double pLambd){
    G4double lKappa = 1.08;
    G4double lParamA = 6954./m;
    G4double lParamB = 6618*nm;
    G4double lAdust = 1./(mSpice_a400inv[mSpiceDepth_pos])*pow(pLambd/(400.*nm), - lKappa);
    G4double lDeltaTau = Spice_Temperature(mSpice_Depth[mSpiceDepth_pos])-Spice_Temperature(1730.);
    G4double la_inv = 1./(lAdust+lParamA*exp(-lParamB/pLambd)*(1.+0.01*lDeltaTau));
    return la_inv;
}

/**
 * Calculation of the refraction index of IC-ice for a specific wavelength.
 * @param pLambd Wavelength
 * @return Refraction index
 */
G4double OMSimInputData::Spice_Refraction(G4double pLambd){
    // unknown depth. Parametrization by Thomas Kittler.
    G4double lLambd3 = pLambd*1e-3;
    G4double lNphase = 1.55749 - 1.57988/nm * lLambd3 + 3.99993/(nm*nm) * pow(lLambd3,2) - 4.68271/(nm*nm*nm) * pow(lLambd3,3) + 2.09354/(nm*nm*nm*nm) * pow(lLambd3,4);
    return lNphase;	// using this now after discussion with Timo
}
/**
 * Calculation of the mie scattering length of IC-ice for a specific wavelength
 * @param pLambd Wavelength
 * @return Mie scattering length
 */
G4double OMSimInputData::Mie_Scattering(G4double pLambd){
    // depth_pos is the coordinate for the chosen depth in Depth_spice. For example to choose
    // depth=2278.2 m, we use depth_pos = 88
    G4double lAlpha = 0.90;
    G4double lAv_costheta = 0.9;
    G4double lBe_inv = 1./(1./(mSpice_be400inv[mSpiceDepth_pos])*pow((pLambd/(400.*nm)), -lAlpha));
    G4double lB_inv = lBe_inv*(1.-lAv_costheta);
    return lB_inv;
}


/*
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *                Functions for creating G4Materials and G4OpticalSurface
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
/**
 * Defines new surface from data in json-file. Currently it is more flexible than OMSimInputData::MaterialFromFile
 * since you can add any property, since no units are needed. 
 * @param pFilename name of the file with material properties
 */
void OMSimInputData::ReflectiveSurfaceFromFile(G4String pFilename){ 
    
    pt::ptree lJsonTree;
    pt::read_json(pFilename, lJsonTree); //read json file into lJsonTree
    
    G4String lName = lJsonTree.get<G4String>("jName");
    G4String lModelStr = lJsonTree.get<G4String>("jModel");
    G4String lFinishStr = lJsonTree.get<G4String>("jFinish");
    G4String lTypeStr = lJsonTree.get<G4String>("jType");
    
    G4OpticalSurfaceModel lModel;
    G4OpticalSurfaceFinish lFinish;
    G4SurfaceType lType;
    
    switch(hash(lModelStr)){
        case hash("glisur")         : lModel = glisur; break;
        case hash("unified")        : lModel = unified; break;
        case hash("LUT")            : lModel = LUT; break;
    }    
    
    switch(hash(lTypeStr)){
        case hash("dielectric_metal")       : lType = dielectric_metal; break;
        case hash("dielectric_dielectric")  : lType = dielectric_dielectric; break;
        case hash("dielectric_LUT")         : lType = dielectric_LUT; break;
        case hash("firsov")                 : lType = firsov; break;
        case hash("x_ray")                  : lType = x_ray; break;
    }
    
    switch(hash(lFinishStr)){
        case hash("polished") 	            : lFinish = polished 	; break;
        case hash("polishedfrontpainted") 	: lFinish = polishedfrontpainted 	; break;
        case hash("polishedbackpainted") 	: lFinish = polishedbackpainted 	; break;
        case hash("ground") 	            : lFinish = ground 	; break;
        case hash("groundfrontpainted") 	: lFinish = groundfrontpainted 	; break;
        case hash("groundbackpainted") 	    : lFinish = groundbackpainted 	; break;
        case hash("polishedlumirrorair") 	: lFinish = polishedlumirrorair 	; break;
        case hash("polishedlumirrorglue") 	: lFinish = polishedlumirrorglue 	; break;
        case hash("polishedair") 	        : lFinish = polishedair 	; break;
        case hash("polishedteflonair") 	    : lFinish = polishedteflonair 	; break;
        case hash("polishedtioair") 	    : lFinish = polishedtioair 	; break;
        case hash("polishedtyvekair") 	    : lFinish = polishedtyvekair 	; break;
        case hash("polishedvm2000air") 	    : lFinish = polishedvm2000air 	; break;
        case hash("polishedvm2000glue") 	: lFinish = polishedvm2000glue 	; break;
        case hash("etchedlumirrorair") 	    : lFinish = etchedlumirrorair 	; break;
        case hash("etchedlumirrorglue") 	: lFinish = etchedlumirrorglue 	; break;
        case hash("etchedair") 	            : lFinish = etchedair 	; break;
        case hash("etchedteflonair") 	    : lFinish = etchedteflonair 	; break;
        case hash("etchedtioair") 	        : lFinish = etchedtioair 	; break;
        case hash("etchedtyvekair") 	    : lFinish = etchedtyvekair 	; break;
        case hash("etchedvm2000air") 	    : lFinish = etchedvm2000air 	; break;
        case hash("etchedvm2000glue") 	    : lFinish = etchedvm2000glue 	; break;
	    case hash("groundlumirrorair") 	    : lFinish = groundlumirrorair 	; break;
        case hash("groundlumirrorglue") 	: lFinish = groundlumirrorglue 	; break;
        case hash("groundair") 	            : lFinish = groundair 	; break;
        case hash("groundteflonair") 	    : lFinish = groundteflonair 	; break;
        case hash("groundtioair") 	        : lFinish = groundtioair 	; break;
        case hash("groundtyvekair") 	    : lFinish = groundtyvekair 	; break;
        case hash("groundvm2000air") 	    : lFinish = groundvm2000air 	; break;
        case hash("groundvm2000glue")       : lFinish = groundvm2000glue; break;
    }
    
    G4OpticalSurface* lOpticalSurface = new G4OpticalSurface(lName, lModel,  lFinish, lType);
    G4MaterialPropertiesTable* lMPT = new G4MaterialPropertiesTable();
    
    try // Only few materials have jSigmaAlpha defined
    {
        G4double lSigmaAlpha = lJsonTree.get<G4double>("jSigmaAlpha");
        lOpticalSurface -> SetSigmaAlpha(lSigmaAlpha);
    }
    catch (...) {} // not very elegant, I know...
    
    for (pt::ptree::value_type &key : lJsonTree.get_child("jProperties"))
    {   
        G4String lKey = key.second.get_value<G4String>();
        std::vector<G4double> lPhotonEnergy;
        std::vector<G4double> lValues;
        ParseToVector(lValues, lJsonTree, "jValues_"+lKey, 1., false);
        ParseToVector(lPhotonEnergy, lJsonTree, "jWavelength_"+lKey, mHC_eVnm, true);
        lMPT-> AddProperty(lKey,  &lPhotonEnergy[0], &lValues[0], static_cast<int>(lPhotonEnergy.size()));
    }
    
    lOpticalSurface->SetMaterialPropertiesTable(lMPT);
    mOpticalSurfaceMap[lName] = lOpticalSurface;
    
    G4cout << "New Optical Surface: " << lName << G4endl;
}


/**
 * Defines new material from data in json-file. The prefix of the filename tells this function what kind of
 * material we are handling:
 * "RiAbs"  : Material with refraction index and absorption length.
 * "Ri"     : Material with refraction index only
 * "NoOptic": Material has no optics defined
 *  "ICE"   : This is a special file with data from Icecube. Here we call several function for calculating the optical properties
 * @param pFilename name of the file with material properties
 */
void OMSimInputData::MaterialFromFile(G4String pFilename){    
    
    pt::ptree lJsonTree;
    pt::read_json(pFilename, lJsonTree); //read json file into lJsonTree
    
    mMaterialData->AppendParameterTable(pFilename);
    
    G4MaterialPropertiesTable* lMPT  = new G4MaterialPropertiesTable();
    G4NistManager* lMatDatBase = G4NistManager::Instance();
    
    G4String lName = lJsonTree.get<G4String>("jName");
    G4String lDataType = lJsonTree.get<G4String>("jDataType");
    
    mMaterialData->mSelectedKey = lName;
    G4double lDensity = mMaterialData->Get("jDensity");
    
    G4State lState;
    G4String lState_str = lJsonTree.get<G4String>("jState");
    
    std::vector<G4double> lRefractionIndex;
    std::vector<G4double> lRefractionIndexEnergy;
    std::vector<G4double> lAbsLength;
    std::vector<G4double> lAbsLengthEnergy;
    
    switch(hash(lState_str)){
        case hash("kStateSolid")  : lState = kStateSolid; break;
        case hash("kStateLiquid") : lState = kStateLiquid; break;
        case hash("kStateGas")    : lState = kStateGas; break;
        
        default                   : lState = kStateUndefined;
    }
    
    //Defining the material with its density, number of components, state and name
    G4Material* lMaterial = new G4Material(lName, lDensity, lJsonTree.get_child("jComponents").size() , lState);
    
    //Construct material with fractional components (isotopes or G4-Materials)
    for (pt::ptree::value_type &key : lJsonTree.get_child("jComponents"))
    {
        std::string componentName = key.first;
        double componentFraction = key.second.get_value<double>();
        lMaterial->AddMaterial(lMatDatBase->FindOrBuildMaterial(componentName), componentFraction);
    }
    
    //Check which optical properties are included in the json and parse them to vectors
    if (lDataType == "RiAbs") {
        ParseToVector(lRefractionIndex, lJsonTree, "jRefractiveIdx", 1., false);
        ParseToVector(lRefractionIndexEnergy, lJsonTree, "jRefractiveIdxWavelength", mHC_eVnm, true);
        ParseToVector(lAbsLength, lJsonTree, "jAbsLength", 1*mm, false );
        ParseToVector(lAbsLengthEnergy, lJsonTree, "jAbsLengthWavelength", mHC_eVnm, true);
        lMPT->AddProperty("RINDEX", &lRefractionIndexEnergy[0], &lRefractionIndex[0], static_cast<int>(lRefractionIndex.size()));
        lMPT->AddProperty("ABSLENGTH", &lAbsLengthEnergy[0] , &lAbsLength[0], static_cast<int>(lAbsLength.size()));
        lMaterial->SetMaterialPropertiesTable(lMPT);
    }
    else if (lDataType == "Ri") {
        ParseToVector(lRefractionIndex, lJsonTree, "jRefractiveIdx", 1., false);
        ParseToVector(lRefractionIndexEnergy, lJsonTree, "jRefractiveIdxWavelength", mHC_eVnm, true);
        lMPT->AddProperty("RINDEX", &lRefractionIndexEnergy[0], &lRefractionIndex[0], static_cast<int>(lRefractionIndex.size()));
        lMaterial->SetMaterialPropertiesTable(lMPT);
    }
    else if (lDataType == "ICE"){
        G4Material* lIceMie = new G4Material(lName+"_SPICE", lDensity, lMatDatBase->FindOrBuildMaterial("G4_WATER") , lState);
        G4double MIE_spice_const[3]={0.972,0, 1};
        std::vector<G4double> lMieScatteringLength;
        std::vector<G4double> lWavelength;
        ParseToVector(lRefractionIndexEnergy, lJsonTree, "jWavelength_spice", mHC_eVnm, true);
        ParseToVector(lWavelength, lJsonTree, "jWavelength_spice", 1*nm, false);
        ParseToVector(mSpice_be400inv, lJsonTree, "jbe400inv_spice", 1*m, false);
        ParseToVector(mSpice_a400inv, lJsonTree, "ja400inv_spice", 1*m, false );
        ParseToVector(mSpice_Depth, lJsonTree, "jDepth_spice", 1*m, false);
        
        for (int u = 0; u < static_cast<int>(lRefractionIndexEnergy.size()); u++) {
            lRefractionIndex.push_back(Spice_Refraction(lWavelength.at(u)));
            lAbsLength.push_back(Spice_Absorption(lWavelength.at(u)));
            lMieScatteringLength.push_back(Mie_Scattering(lWavelength.at(u)));
        }
        
        lMPT->AddProperty("RINDEX", &lRefractionIndexEnergy[0], &lRefractionIndex[0], static_cast<int>(lRefractionIndex.size()));
        lMaterial->SetMaterialPropertiesTable(lMPT);
        
        lMPT->AddProperty("ABSLENGTH", &lRefractionIndexEnergy[0] , &lAbsLength[0], static_cast<int>(lAbsLength.size()));
        lMPT->AddProperty("MIEHG",&lRefractionIndexEnergy[0], &lMieScatteringLength[0] , static_cast<int>(lRefractionIndex.size()))->SetSpline(true);
        lMPT->AddConstProperty("MIEHG_FORWARD",MIE_spice_const[0]);
        lMPT->AddConstProperty("MIEHG_BACKWARD",MIE_spice_const[1]);
        lMPT->AddConstProperty("MIEHG_FORWARD_RATIO",MIE_spice_const[2]);
        lIceMie->SetMaterialPropertiesTable(lMPT);
        G4cout << "Ice properties at depth " << mSpice_Depth[mSpiceDepth_pos]/m << " m." << G4endl;
    }
    G4cout << "New Material defined: " << lMaterial->GetName()<<G4endl;
    
    
}


/*
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *                Function for reading json files and handling of json-trees
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

/**
 * Calls ScannDataDirectory in the defined folders.
 * @see ScannDataDirectory
 */
void OMSimInputData::SearchFolders(){
    mDataDirectory = "../data/Materials";
    ScannDataDirectory();
    mDataDirectory = "../data/PMTs";
    ScannDataDirectory();
    mDataDirectory = "../data/SegmentedModules";
    ScannDataDirectory();
}

template <typename T>
/**
 * Transforms the values inside a pTree-array to a vector. The values can be also transformed to a G4double.
 * @param pVector  vector where the (transformed) values are saved
 * @param pTree    pTree containing json data
 * @param pKey     json attribute label where values are found
 * @param pScaling Values of array are multiplied by this factor. You can set the fisical unit with this.
 * @param pInverse In case you need the inverse of a value x, 1/x is appended (e.g. transforming from nm to eV)
 */
void OMSimInputData::ParseToVector(std::vector<T>& pVector, pt::ptree pTree, std::basic_string<char> pKey, G4double pScaling, bool pInverse){
    for (pt::ptree::value_type &ridx : pTree.get_child(pKey)){ //get array from element with key "pKey" of the json
        if (pInverse) { // if we need 1/x
            pVector.push_back(pScaling/ridx.second.get_value<T>());
        }
        else { // otherwise we only by scaling factor
            pVector.push_back(ridx.second.get_value<T>()*pScaling); 
        }
    }
}

/**
 * Scann for data files inside mDataDirectory and creates materials.
 * @param pName name of the material
 */
void OMSimInputData::ScannDataDirectory(){
    
    struct dirent* lFile= NULL;
    DIR* lDirectory = NULL;
    
    lDirectory=opendir( mDataDirectory.data());
    if(lDirectory == NULL)
    {
        G4cerr << "Couldn't open directory" << mDataDirectory << G4endl;
        return;
    }
    
    while( (lFile = readdir(lDirectory)) )
    {
        std::string fileName = lFile->d_name;
        if (lFile->d_type == 8){
            if ((fileName.substr (0,5) == "RiAbs"  || fileName.substr (0,10) == "IceCubeICE" || fileName.substr (0,7) == "NoOptic" || fileName.substr (0,2) == "Ri"))
            {
                MaterialFromFile(mDataDirectory+"/"+fileName);
            }
            else if ((fileName.substr (0,4) == "Refl" )) ReflectiveSurfaceFromFile(mDataDirectory+"/"+fileName);
            else if ((fileName.substr (0,4) == "pmt_")) mPMTdata->AppendParameterTable(mDataDirectory+"/"+fileName);
            else if ((fileName.substr (0,3) == "om_")) mOMdata->AppendParameterTable(mDataDirectory+"/"+fileName);
        }
    }
    closedir(lDirectory);
}



/**
 * This class provides handling of table with json Trees
 */
OMSimInputData::ParameterTable::ParameterTable() {}

/**
 * Appends information of json-file containing PMT/OM parameters to a vector of ptrees
 * @param pFileName Name of file containing json
 */
void OMSimInputData::ParameterTable::AppendParameterTable(G4String pFileName) {
    pt::ptree lJsonTree;
    pt::read_json(pFileName, lJsonTree); 
    G4String lName = lJsonTree.get<G4String>("jName");
    mTable[lName] = lJsonTree;
    G4cout << lName << " added to dictionary..." <<  G4endl;
}
/**
 * Get values from a pTree with its unit, transforming it to a G4double
 * @param pParameter Name of parameter in json
 * @return G4double of the value and its unit
 */
G4double OMSimInputData::ParameterTable::Get(G4String pParameter) {
    G4String lUnit =  mTable.at(mSelectedKey).get<G4String>(pParameter+".jUnit");
    G4double lValue = mTable.at(mSelectedKey).get<G4double>(pParameter+".jValue");
    if (lUnit ==  "NULL") {
        return lValue;
    }
    else {
        return lValue*G4UnitDefinition::GetValueOf(lUnit);
    }
    
}




