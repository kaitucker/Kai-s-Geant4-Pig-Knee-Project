#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4String.hh"

class G4Step;
class EventAction;

class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction (EventAction* eventAction, const G4String& targetLogicalName);
  ~SteppingAction() override;
  void UserSteppingAction(const G4Step* step) override;

private:
  EventAction* fEventAction;
  G4String fTargetLogicalName;
};
#endif
