#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "EventAction.hh"

ActionInitialization::ActionInitialization()
: G4VUserActionInitialization()
{}
ActionInitialization::~ActionInitialization()
{}
void ActionInitialization::BuildForMaster() const
{
// Only the RunAction is needed on the master for summaries
SetUserAction(new RunAction());
}
//ONLY CARE ABOUT THIS:
void ActionInitialization::Build() const
{
  // create the RunAction first so pointers can be passed to other actions
  RunAction* runAction = new RunAction();
  SetUserAction(runAction); // <-- KEEP THIS

  SetUserAction(new PrimaryGeneratorAction());

  // --- REPLACE THE OLD LOGIC WITH THIS ---
  auto eventAction = new EventAction(runAction);
  SetUserAction(eventAction);

  // Pass the eventAction pointer to SteppingAction, not the runAction pointer
  SetUserAction(new SteppingAction(eventAction, "PatellarTendonLog"));
}
