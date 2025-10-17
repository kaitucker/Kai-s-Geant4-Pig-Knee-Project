#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4TouchableHandle.hh"
#include "G4LogicalVolume.hh"
#include "EventAction.hh"

SteppingAction::SteppingAction(EventAction* eventAction, const G4String& targetLogicalName)
 : G4UserSteppingAction(), fEventAction(eventAction), fTargetLogicalName(targetLogicalName) {}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
  if (!step) return;
  auto prePoint = step->GetPreStepPoint();
  if (!prePoint) return;
  auto touch = prePoint->GetTouchableHandle();
  if (!touch) return;
  auto volume = touch->GetVolume();
  if (!volume) return;
  G4LogicalVolume* lv = volume->GetLogicalVolume();
  if (!lv) return;

  if (lv->GetName() == fTargetLogicalName) {
    G4double edep = step->GetTotalEnergyDeposit();
    if (edep > 0.) {
      fEventAction->AddEdep(edep);
    }
  }
}
