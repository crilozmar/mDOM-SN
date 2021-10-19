#ifndef mdomSteppingAction_h
#define mdomSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class mdomSteppingAction : public G4UserSteppingAction
{
  public:
    mdomSteppingAction();
   ~mdomSteppingAction(){};

    void UserSteppingAction(const G4Step*);
    G4double get_excitation (G4String s);
    bool QEcheck(G4double lambda);
    std::vector<double> QEwavelenght;
    std::vector<double> QEprobability;
    G4double probability;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
