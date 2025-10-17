#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
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
  SetUserAction(new PrimaryGeneratorAction);

  auto runAction = new RunAction;
  SetUserAction(runAction);

  auto eventAction = new EventAction(runAction); // New line
  SetUserAction(eventAction);                     // New line

  SetUserAction(new SteppingAction(eventAction)); // Pass eventAction, not runAction
}
