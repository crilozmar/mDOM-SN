#ifndef OMSimInputData_h
#define OMSimInputData_h 1

#include <boost/property_tree/ptree.hpp>                                        
#include <boost/property_tree/json_parser.hpp> 
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4SystemOfUnits.hh"


namespace pt = boost::property_tree;



class OMSimInputData
{
    
    private:
    class ParameterTable{
    public:
        ParameterTable();
        G4double Get(G4String pParameter);
        G4String mSelectedKey;
        void AppendParameterTable(G4String pFileName);
        std::map<G4String, pt::ptree> mTable;
    };
    
public:
    
    OMSimInputData();
    G4OpticalSurface* GetOpticalSurface(G4String pName);
    G4Material* GetMaterial(G4String name);
    G4int mSpiceDepth_pos; //depth = 2278.2 m
    void SearchFolders();
    
    ParameterTable* mPMTdata= new ParameterTable();
    ParameterTable* mOMdata = new ParameterTable();
    ParameterTable* mMaterialData = new ParameterTable();
    
private:

    void MaterialFromFile(G4String filename);
    template <typename T>
    void ParseToVector(std::vector<T>& pVector, pt::ptree pTree, std::basic_string<char> pKey, G4double pScaling, bool pInverse);
    void ReflectiveSurfaceFromFile(G4String pFilename);
    void ScannDataDirectory();
    
    G4String mDataDirectory;
    
    std::vector<G4double> mSpice_be400inv;
    std::vector<G4double> mSpice_a400inv;
    std::vector<G4double> mSpice_Depth;
    
    std::map<G4String, pt::ptree> mPMTtable;
    std::map<G4String, pt::ptree> mOMtable;
    
    G4double Spice_Temperature(G4double depth);
    G4double Spice_Absorption(G4double pLambd);
    G4double Spice_Refraction(G4double pLambd);
    G4double Mie_Scattering(G4double pLambd);
    
    G4double mHC_eVnm = 1239.84193*eV; // h*c in eV * nm
    
    std::map<G4String, G4OpticalSurface*> mOpticalSurfaceMap;
};

#endif
//
