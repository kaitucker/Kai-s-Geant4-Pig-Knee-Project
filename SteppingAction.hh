#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4String.hh"

class RunAction;
class G4Step;

class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction(RunAction* runAction, const G4String& targetLogicalName);
  ~SteppingAction() override;
  void UserSteppingAction(const G4Step* step) override;
private:
  RunAction* fRunAction;
  G4String fTargetLogicalName;
};
#endif
