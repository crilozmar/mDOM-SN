#ifndef mdomStackingAction_H
#define mdomStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class mdomStackingAction : public G4UserStackingAction
{
  public:

    mdomStackingAction();
    virtual ~mdomStackingAction();
 
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
	virtual void PrepareNewEvent();
  private:
};

#endif